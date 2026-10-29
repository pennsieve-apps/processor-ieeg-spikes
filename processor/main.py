"""
processor-ieeg-spikes: interictal epileptiform discharge (spike) detection on iEEG.

Reads the BIDS-iEEG derivative of ieeg-prep (or a plain BIDS-iEEG dataset, or loose
EDF/BrainVision/FIF files) from INPUT_DIR and writes to OUTPUT_DIR (or OUTPUT_DIR/<LABEL>/):

  spikes.csv          recording_id, channel, time_s, amplitude_uv, detector, event_id
  spike_rates.csv     recording_id, channel, n_spikes, minutes, rate_per_min, detector
  spike_events.csv    recording_id, event_id, time_s, leader, n_channels, channels, detector
                      (cross-channel events, when GROUP_EVENTS=true)
  provenance.json

One recording at a time, one time window (WINDOW_S) of the good channels at a time,
channels in parallel (THREADS).
"""
from __future__ import annotations

import csv
import json
import logging
import os
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

import numpy as np

from processor import resources as R
from processor import detect as D
from processor import ieeg as I

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-ieeg-spikes")
log.setLevel(logging.INFO)
for noisy in ("mne", "numba", "joblib"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-ieeg-spikes"
REPO_ROOT = Path(__file__).resolve().parent.parent
PAD_S = 10.0     # context read on each side of a window (never across excluded time)
EDGE_S = 1.0     # Janca ignores the first and last second it is given

SPIKES_HEADER = ["recording_id", "channel", "time_s", "amplitude_uv", "detector", "event_id"]
RATES_HEADER = ["recording_id", "channel", "n_spikes", "minutes", "rate_per_min", "detector"]
EVENTS_HEADER = ["recording_id", "event_id", "time_s", "leader", "n_channels", "channels", "detector"]


def version_of(root: Path) -> str:
    try:
        return (root / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for pkg in ("epycom", "mne", "numpy", "scipy", "joblib", "numba"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


# Memory per parallel worker (a loky process holding one channel of a window) and for the main
# process (the window of every good channel), measured on a 186-channel, 1024 Hz demo clip:
# peak 0.75 GB with THREADS=1, 1.7 GB with 4.
WORKER_MB = 320
MAIN_MB = 850


def _bool(name: str, value: str) -> bool:
    v = value.strip().lower()
    if v in ("true", "1", "yes"):
        return True
    if v in ("false", "0", "no"):
        return False
    raise ValueError(f"{name} must be true or false (got {value!r})")


def env_params() -> dict:
    """Parameters from the environment; raises ValueError on a bad value."""
    e = lambda k, d: (os.environ.get(k) or d).strip()  # noqa: E731
    raw = {
        "DETECTOR": e("DETECTOR", "janca").lower(),
        "JANCA_K1": e("JANCA_K1", "3.65"),
        "SEGMENT": e("SEGMENT", "all").lower(),
        "EXCLUDE_SEIZURE": e("EXCLUDE_SEIZURE", "true"),
        "SEIZURE_MARGIN_S": e("SEIZURE_MARGIN_S", "10"),
        "WINDOW_S": e("WINDOW_S", "300"),
        "GROUP_EVENTS": e("GROUP_EVENTS", "true"),
        "GROUP_WINDOW_MS": e("GROUP_WINDOW_MS", "50"),
        "THREADS": e("THREADS", "0"),
        "LABEL": e("LABEL", ""),
    }
    p = dict(raw)
    if p["DETECTOR"] not in ("janca", "barkmeier"):
        raise ValueError(f"DETECTOR must be janca or barkmeier (got {raw['DETECTOR']!r})")
    if p["SEGMENT"] not in ("all", "baseline"):
        raise ValueError(f"SEGMENT must be all or baseline (got {raw['SEGMENT']!r})")
    p["k1"] = float(raw["JANCA_K1"])
    if p["k1"] <= 0:
        raise ValueError("JANCA_K1 must be > 0")
    p["exclude_seizure"] = _bool("EXCLUDE_SEIZURE", raw["EXCLUDE_SEIZURE"])
    p["margin"] = float(raw["SEIZURE_MARGIN_S"])
    if p["margin"] < 0:
        raise ValueError("SEIZURE_MARGIN_S must be >= 0")
    p["window"] = float(raw["WINDOW_S"])
    if p["window"] < 60:
        raise ValueError("WINDOW_S must be >= 60")
    p["group"] = _bool("GROUP_EVENTS", raw["GROUP_EVENTS"])
    p["group_window"] = float(raw["GROUP_WINDOW_MS"]) / 1000.0
    if p["group_window"] <= 0:
        raise ValueError("GROUP_WINDOW_MS must be > 0")
    threads = int(raw["THREADS"])
    if threads < 0:
        raise ValueError("THREADS must be >= 0 (0 = every vCPU)")
    p["threads"], p["threads_why"] = R.resolve_workers(threads, WORKER_MB, MAIN_MB)
    if p["LABEL"] and not all(c.isalnum() or c in "-_" for c in p["LABEL"]):
        raise ValueError("LABEL may hold letters, digits, - and _ only")
    return p


class CsvOut:
    def __init__(self, path: Path, header: list[str]):
        self.fh = open(path, "w", newline="", encoding="utf-8")
        self.w = csv.DictWriter(self.fh, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        self.w.writeheader()

    def rows(self, rows):
        self.w.writerows(rows)
        self.fh.flush()

    def close(self):
        self.fh.close()


def process_recording(rec: I.Recording, params: dict, parallel) -> dict:
    """Detect spikes in one recording. Returns spikes, rates, events and a status summary."""
    from joblib import delayed

    t0 = time.monotonic()
    raw = I.open_raw(rec.path)
    fs = float(raw.info["sfreq"])
    n_times = int(raw.n_times)
    duration = n_times / fs
    if fs < 100:
        raise ValueError(f"sampling rate {fs:g} Hz is too low for spike detection (needs >= 100 Hz)")
    chans = I.channels(rec, list(raw.ch_names), list(raw.get_channel_types()))
    good = [c for c in chans if c.status == "good"]
    if not good:
        raise ValueError(f"no good iEEG channel among {len(chans)}")
    index = {name: i for i, name in enumerate(raw.ch_names)}
    picks = [index[c.name] for c in good]
    plan = I.analysis_plan(rec, duration, params["SEGMENT"], params["exclude_seizure"], params["margin"])

    minutes = {c.name: 0.0 for c in good}
    found: dict[str, list[tuple[float, float]]] = {c.name: [] for c in good}
    errors: dict[str, str] = {}
    windows = 0
    t_read = t_detect = 0.0
    for c0, c1, r0, r1 in I.chunks(plan.segments, params["window"], PAD_S, D.MIN_SECONDS, EDGE_S):
        if r0 is None:
            rec.notes.append(f"segment {c0:.1f}-{c1:.1f} s is shorter than {D.MIN_SECONDS:g} s; not analysed")
            continue
        s, e = max(0, int(round(r0 * fs))), min(n_times, int(round(r1 * fs)))
        ta = time.monotonic()
        data = I.read_window(raw, picks, s, e)
        tb = time.monotonic()
        ratios = [1.0] * len(good)
        if params["DETECTOR"] == "barkmeier":
            levels = np.array(parallel(delayed(D.broadband_level)(data[i], fs) for i in range(len(good))))
            ref = np.nanmedian(levels) if np.any(np.isfinite(levels)) else np.nan
            ratios = [float(v / ref) if np.isfinite(ref) and ref > 0 else 1.0 for v in levels]
        results = parallel(delayed(D.detect_channel)(data[i], fs, params["DETECTOR"], params["k1"], s / fs, (c0, c1),
                                                     ratios[i]) for i in range(len(good)))
        del data
        t_read += tb - ta
        t_detect += time.monotonic() - tb
        windows += 1
        for c, (dets, err) in zip(good, results):
            if err:
                errors.setdefault(c.name, f"{c0:.0f}-{c1:.0f} s: {err}")
                continue
            minutes[c.name] += (c1 - c0) / 60.0
            found[c.name].extend(dets)

    analysed = [c for c in good if minutes[c.name] > 0]
    total_min = max(minutes.values()) if minutes else 0.0
    det = params["DETECTOR"]
    spikes = []
    for c in good:
        for t, amp in sorted(found[c.name]):
            spikes.append({"recording_id": rec.recording_id, "channel": c.name, "time_s": f"{t:.4f}",
                           "amplitude_uv": f"{amp:.2f}", "detector": det, "event_id": "", "_t": t, "_a": amp})
    events = []
    if params["group"] and spikes:
        grouped = D.group_events([(s["_t"], s["channel"], s["_a"]) for s in spikes], params["group_window"])
        for k, ev in enumerate(grouped, start=1):
            for i in ev["members"]:
                spikes[i]["event_id"] = k
            events.append({"recording_id": rec.recording_id, "event_id": k, "time_s": f"{ev['time_s']:.4f}",
                           "leader": ev["leader"], "n_channels": ev["n_channels"],
                           "channels": "|".join(ev["channels"]), "detector": det})
    spikes.sort(key=lambda s: (s["_t"], s["channel"]))
    rates = []
    for c in analysed:
        n = len(found[c.name])
        m = minutes[c.name]
        rates.append({"recording_id": rec.recording_id, "channel": c.name, "n_spikes": n,
                      "minutes": f"{m:.3f}", "rate_per_min": f"{n / m:.4f}", "detector": det})
    summary = {
        "sfreq": fs, "duration_s": round(duration, 3), "n_channels": len(chans), "n_good": len(good),
        "n_analysed": len(analysed), "analysed": plan.description,
        "segments": [[round(a, 3), round(b, 3)] for a, b in plan.segments], "excluded": plan.excluded,
        "minutes": round(total_min, 3), "windows": windows, "n_spikes": len(spikes), "n_events": len(events),
        "channels_not_analysed": {c.name: c.reason for c in chans if c.status != "good"},
        "channel_errors": errors, "notes": rec.notes,
        "seconds": {"read": round(t_read, 2), "detect": round(t_detect, 2), "total": round(time.monotonic() - t0, 2)},
    }
    return {"spikes": spikes, "rates": rates, "events": events, "summary": summary}


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    try:
        params = env_params()
    except ValueError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    out_dir = output_dir / params["LABEL"] if params["LABEL"] else output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    recs, seen, problems = I.find_recordings(input_dir)
    log.info("input: %d recording(s), %d file(s); detector %s, segment %s, seizure %s, %d thread(s)",
             len(recs), len(seen), params["DETECTOR"], params["SEGMENT"],
             f"excluded (margin {params['margin']:g} s)" if params["exclude_seizure"] else "kept", params["threads"])
    log.info("parallel workers: %d (%s)", params["threads"], params["threads_why"])
    for p in problems:
        log.warning("%s: %s", p["input"], p["reason"])

    spikes_out = CsvOut(out_dir / "spikes.csv", SPIKES_HEADER)
    rates_out = CsvOut(out_dir / "spike_rates.csv", RATES_HEADER)
    events_out = CsvOut(out_dir / "spike_events.csv", EVENTS_HEADER) if params["group"] else None
    prov_recs: list[dict] = [{"recording_id": p["recording_id"], "input": p["input"],
                              "status": p.get("status", "failed"), "error": p["reason"]} for p in problems]
    n_ok = n_failed = 0
    n_failed += sum(1 for p in problems if p.get("status", "failed") == "failed")

    from joblib import Parallel

    with Parallel(n_jobs=params["threads"], backend="loky", max_nbytes=None) as parallel:
        for rec in recs:
            try:
                res = process_recording(rec, params, parallel)
            except Exception as exc:  # noqa: BLE001 - one recording's failure must not sink the run
                log.error("%s: FAILED: %s\n%s", rec.recording_id, exc, traceback.format_exc())
                prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": "failed",
                                  "error": f"{type(exc).__name__}: {exc}"})
                n_failed += 1
                continue
            s = res["summary"]
            if s["n_analysed"] == 0:
                status = "skipped: no analysable time" + (
                    f" ({'; '.join(s['notes'])})" if s["notes"] else f" ({s['analysed']})")
                log.warning("%s: %s", rec.recording_id, status)
                prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": status, **s})
                continue
            spikes_out.rows(res["spikes"])
            rates_out.rows(res["rates"])
            if events_out:
                events_out.rows(res["events"])
            n_ok += 1
            for ch, err in s["channel_errors"].items():
                log.warning("%s: channel %s: %s", rec.recording_id, ch, err)
            for note in s["notes"]:
                log.warning("%s: %s", rec.recording_id, note)
            log.info("%s: %d channels (%d good, %d analysed), %.1f Hz, %.2f min analysed (%s), %d spikes "
                     "(%.2f/min/channel)%s, %.1f s (read %.1f s, detect %.1f s)",
                     rec.recording_id, s["n_channels"], s["n_good"], s["n_analysed"], s["sfreq"], s["minutes"],
                     s["analysed"], s["n_spikes"], s["n_spikes"] / max(s["minutes"], 1e-9) / max(s["n_analysed"], 1),
                     f", {s['n_events']} events" if events_out else "", s["seconds"]["total"], s["seconds"]["read"],
                     s["seconds"]["detect"])
            prov_recs.append({"recording_id": rec.recording_id, "input": rec.rel, "status": "ok", **s})
    for f in (spikes_out, rates_out, events_out):
        if f:
            f.close()

    status = "ok"
    if not recs and not problems:
        log.error("no recording found under INPUT_DIR (looked for recordings.csv and %s files); it holds %d file(s): %s",
                  "/".join(I.SIGNAL_EXTENSIONS), len(seen), ", ".join(seen[:40]) + (" …" if len(seen) > 40 else ""))
        status = "failed: no input"
    elif n_ok == 0 and n_failed > 0:
        status = "failed: every recording failed"
    elif n_ok == 0:
        status = "ok: nothing analysable"
    prov = {
        "app": APP_NAME, "version": version_of(REPO_ROOT), "tools": tool_versions(),
        "parameters": {k: params[k] for k in ("DETECTOR", "JANCA_K1", "SEGMENT", "EXCLUDE_SEIZURE",
                                              "SEIZURE_MARGIN_S", "WINDOW_S", "GROUP_EVENTS", "GROUP_WINDOW_MS",
                                              "THREADS", "LABEL")},
        "threads_used": params["threads"],
        "threads_why": params["threads_why"],
        "inputs": seen, "status": status,
        "started": started.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2),
        "recordings": prov_recs,
    }
    with open(out_dir / "provenance.json", "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=2)
    log.info("done: %d recording(s) ok, %d failed, %d skipped, %.1f s -> %s", n_ok, n_failed,
             len(prov_recs) - n_ok - n_failed, prov["seconds"], out_dir)
    return 1 if status.startswith("failed") else 0


if __name__ == "__main__":
    sys.exit(main())
