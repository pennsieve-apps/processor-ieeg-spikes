"""
A tiny synthetic BIDS-iEEG derivative that follows workflow-epilepsy-ieeg/docs/conventions.md
(the `ieeg` output of ieeg-prep), for the tests and `make run`:

  <out>/ieeg/dataset_description.json
  <out>/ieeg/recordings.csv
  <out>/ieeg/sub-01/ses-01/ieeg/sub-01_ses-01_task-ictal_run-01_desc-preproc_{ieeg.edf,ieeg.json,channels.tsv,events.tsv}
  <out>/ieeg/sub-02/ieeg/sub-02_task-rest_run-01_desc-preproc_{ieeg.edf,ieeg.json,channels.tsv,events.tsv}

sub-01: 512 Hz, 180 s, bipolar SEEG, a seizure at 120 s (30 s), baseline 0-110 s. Channel
LA01-LA02 carries a spike every ~2 s (rate ~30/min) and LA02-LA03 the same spikes 10 ms
later (one cross-channel event, led by LA01-LA02); LB01-LB02 is marked bad, EKG is an ECG
channel; the other channels are background. During the seizure every channel has
rhythmic spiking, which EXCLUDE_SEIZURE must keep out of the counts.
sub-02: 200 Hz (like the HUP atlas clips), 90 s, no labels; RA01-RA02 carries spikes at
~12/min.

Usage: python -m tests.make_fixture <out-dir>
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

SEED = 42


def pink(n: int, rng: np.random.Generator, std: float) -> np.ndarray:
    """1/f background, `std` µV."""
    spec = np.fft.rfft(rng.standard_normal(n))
    f = np.fft.rfftfreq(n)
    f[0] = f[1]
    x = np.fft.irfft(spec / np.sqrt(f), n)
    return (x / x.std() * std).astype(np.float64)


def spike_shape(fs: float, amp: float) -> np.ndarray:
    """A sharp biphasic transient (~60 ms) followed by a slow wave (~200 ms), peak `amp` µV."""
    t = np.arange(int(0.3 * fs)) / fs
    sharp = amp * np.exp(-((t - 0.04) / 0.012) ** 2) - 0.35 * amp * np.exp(-((t - 0.075) / 0.015) ** 2)
    slow = -0.3 * amp * np.exp(-((t - 0.18) / 0.05) ** 2)
    return sharp + slow


def spike_times(duration: float, rate_per_min: float, rng: np.random.Generator, start: float = 2.0,
                stop: float | None = None) -> np.ndarray:
    stop = duration - 2.0 if stop is None else stop
    step = 60.0 / rate_per_min
    t = np.arange(start, stop, step)
    return t + rng.uniform(-0.2 * step, 0.2 * step, size=len(t))


def add_at(x: np.ndarray, fs: float, times, shape: np.ndarray, delay: float = 0.0) -> None:
    for t in times:
        i = int(round((t + delay) * fs))
        j = min(len(x), i + len(shape))
        if 0 <= i < len(x):
            x[i:j] += shape[: j - i]


def write_edf(path: Path, data: np.ndarray, names: list[str], fs: float) -> None:
    """EDF+ with per-channel physical range from the data (as ieeg-prep writes it)."""
    from edfio import Edf, EdfSignal

    sigs = []
    for x, name in zip(data, names):
        r = max(1.0, float(np.max(np.abs(x))) * 1.1)
        sigs.append(EdfSignal(x, sampling_frequency=fs, label=name, physical_dimension="uV",
                              physical_range=(-r, r)))
    Edf(sigs).write(path)


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


CHANNEL_COLS = ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "status",
                "status_description", "contacts", "soz", "resected", "onset", "spread"]


def channel_rows(names, types, fs, bad=(), soz=()):
    rows = []
    for n, t in zip(names, types):
        rows.append({"name": n, "type": t, "units": "µV", "low_cutoff": "0.5", "high_cutoff": f"{fs / 2:g}",
                     "sampling_frequency": f"{fs:g}", "status": "bad" if n in bad else "good",
                     "status_description": "flat" if n in bad else "n/a", "contacts": n.replace("-", "|"),
                     "soz": ("true" if n in soz else "false") if soz else "", "resected": "", "onset": "",
                     "spread": ""})
    return rows


def make(out: Path) -> Path:
    rng = np.random.default_rng(SEED)
    root = out / "ieeg"
    root.mkdir(parents=True, exist_ok=True)
    (root / "dataset_description.json").write_text(json.dumps({
        "Name": "synthetic fixture", "BIDSVersion": "1.9.0", "DatasetType": "derivative",
        "GeneratedBy": [{"Name": "ieeg-prep", "Description": "tests/make_fixture.py of processor-ieeg-spikes"}]},
        indent=2))
    truth: dict = {}
    recs = []

    # --- sub-01: ictal clip, 512 Hz ------------------------------------------------
    fs, dur = 512.0, 180.0
    n = int(fs * dur)
    names = ["LA01-LA02", "LA02-LA03", "LA03-LA04", "LB01-LB02", "LB02-LB03", "LC01-LC02", "EKG"]
    types = ["SEEG"] * 6 + ["ECG"]
    data = np.vstack([pink(n, rng, 20.0) for _ in names])
    data[3] *= 0.01  # the "bad" channel
    data[6] = 300 * np.sin(2 * np.pi * 1.2 * np.arange(n) / fs)
    times = spike_times(dur, 30.0, rng, stop=118.0)
    shape = spike_shape(fs, 400.0)
    add_at(data[0], fs, times, shape)
    add_at(data[1], fs, times, 0.6 * shape, delay=0.010)
    sz0, sz1 = 120.0, 150.0
    sz_times = np.arange(sz0 + 0.5, sz1, 0.25)
    for i in range(6):
        add_at(data[i], fs, sz_times, spike_shape(fs, 250.0))
    stem = "sub-01_ses-01_task-ictal_run-01"
    d = root / "sub-01" / "ses-01" / "ieeg"
    d.mkdir(parents=True, exist_ok=True)
    write_edf(d / f"{stem}_desc-preproc_ieeg.edf", data, names, fs)
    write_tsv(d / f"{stem}_desc-preproc_channels.tsv", CHANNEL_COLS,
              channel_rows(names, types, fs, bad=("LB01-LB02",), soz=("LA01-LA02", "LA02-LA03")))
    write_tsv(d / f"{stem}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"], [
        {"onset": "0.0", "duration": "110.0", "trial_type": "baseline", "source": "fixture"},
        {"onset": f"{sz0}", "duration": f"{sz1 - sz0}", "trial_type": "seizure", "source": "fixture"}])
    (d / f"{stem}_desc-preproc_ieeg.json").write_text(json.dumps(
        {"SamplingFrequency": fs, "PowerLineFrequency": 60, "iEEGReference": "bipolar"}, indent=2))
    recs.append({"recording_id": stem, "subject": "01", "session": "01", "task": "ictal", "run": "01",
                 "path": f"sub-01/ses-01/ieeg/{stem}_desc-preproc_ieeg.edf", "sfreq": f"{fs:g}",
                 "n_channels": len(names), "n_good": 5, "duration_s": f"{dur:g}", "montage": "bipolar",
                 "kind": "ictal", "onset_s": f"{sz0:g}", "baseline_start_s": "0", "baseline_end_s": "110",
                 "source_format": "bids-edf"})
    truth[stem] = {"spikes": times.tolist(), "channel": "LA01-LA02", "follower": "LA02-LA03", "seizure": [sz0, sz1]}

    # --- sub-02: interictal clip, 200 Hz -------------------------------------------
    fs, dur = 200.0, 90.0
    n = int(fs * dur)
    names = ["RA01-RA02", "RA02-RA03", "RB01-RB02"]
    data = np.vstack([pink(n, rng, 15.0) for _ in names])
    times2 = spike_times(dur, 12.0, rng)
    add_at(data[0], fs, times2, spike_shape(fs, 300.0))
    stem = "sub-02_task-rest_run-01"
    d = root / "sub-02" / "ieeg"
    d.mkdir(parents=True, exist_ok=True)
    write_edf(d / f"{stem}_desc-preproc_ieeg.edf", data, names, fs)
    write_tsv(d / f"{stem}_desc-preproc_channels.tsv", CHANNEL_COLS, channel_rows(names, ["SEEG"] * 3, fs))
    write_tsv(d / f"{stem}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"], [])
    (d / f"{stem}_desc-preproc_ieeg.json").write_text(json.dumps({"SamplingFrequency": fs}, indent=2))
    recs.append({"recording_id": stem, "subject": "02", "session": "", "task": "rest", "run": "01",
                 "path": f"sub-02/ieeg/{stem}_desc-preproc_ieeg.edf", "sfreq": f"{fs:g}", "n_channels": 3,
                 "n_good": 3, "duration_s": f"{dur:g}", "montage": "bipolar", "kind": "interictal", "onset_s": "",
                 "baseline_start_s": "", "baseline_end_s": "", "source_format": "bids-edf"})
    truth[stem] = {"spikes": times2.tolist(), "channel": "RA01-RA02"}

    header = ["recording_id", "subject", "session", "task", "run", "path", "sfreq", "n_channels", "n_good",
              "duration_s", "montage", "kind", "onset_s", "baseline_start_s", "baseline_end_s", "source_format"]
    with open(root / "recordings.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
        w.writeheader()
        w.writerows(recs)
    (out / "truth.json").write_text(json.dumps(truth))
    return out


if __name__ == "__main__":
    print(make(Path(sys.argv[1] if len(sys.argv) > 1 else "data/fixture")))
