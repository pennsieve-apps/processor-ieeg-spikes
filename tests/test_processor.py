"""
Unit tests on a synthetic fixture (tests/make_fixture.py); no Docker, no real data.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from processor import detect as D
from processor import ieeg as I
from processor import main as M
from tests import make_fixture

REC1 = "sub-01_ses-01_task-ictal_run-01"
REC2 = "sub-02_task-rest_run-01"
PARAMS = ("DETECTOR", "JANCA_K1", "SEGMENT", "EXCLUDE_SEIZURE", "SEIZURE_MARGIN_S", "WINDOW_S", "GROUP_EVENTS",
          "GROUP_WINDOW_MS", "THREADS", "LABEL")


@pytest.fixture(scope="module")
def fixture_dir(tmp_path_factory) -> Path:
    return make_fixture.make(tmp_path_factory.mktemp("fixture"))


def run(monkeypatch, input_dir: Path, output_dir: Path, **env) -> int:
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("THREADS", "1")
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    return M.main()


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def rates(out: Path) -> dict[tuple[str, str], dict]:
    return {(r["recording_id"], r["channel"]): r for r in read_csv(out / "spike_rates.csv")}


@pytest.fixture(scope="module")
def default_run(fixture_dir, tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    mp = pytest.MonkeyPatch()
    try:
        assert run(mp, fixture_dir, out) == 0
    finally:
        mp.undo()
    return out


# --- end to end ---------------------------------------------------------------------

def test_outputs_and_columns(default_run):
    for name, header in (("spikes.csv", M.SPIKES_HEADER), ("spike_rates.csv", M.RATES_HEADER),
                         ("spike_events.csv", M.EVENTS_HEADER)):
        with open(default_run / name) as fh:
            assert next(csv.reader(fh)) == header
    prov = json.loads((default_run / "provenance.json").read_text())
    assert prov["app"] == "processor-ieeg-spikes" and prov["version"] == "0.1.0"
    assert prov["status"] == "ok"
    assert {r["recording_id"]: r["status"] for r in prov["recordings"]} == {REC1: "ok", REC2: "ok"}
    assert "epycom" in prov["tools"] and "mne" in prov["tools"]
    assert any(p.endswith("recordings.csv") for p in prov["inputs"])


def test_rates_on_spiking_and_background_channels(default_run, fixture_dir):
    r = rates(default_run)
    # LA01-LA02 spikes at ~30/min until 118 s; 0-110 s analysed (seizure 120-150 s, margin 10 s)
    assert 25 <= float(r[(REC1, "LA01-LA02")]["rate_per_min"]) <= 35
    assert float(r[(REC1, "LA01-LA02")]["minutes"]) == pytest.approx(108 / 60, abs=0.01)
    assert 9 <= float(r[(REC2, "RA01-RA02")]["rate_per_min"]) <= 15
    for key in ((REC1, "LA03-LA04"), (REC1, "LC01-LC02"), (REC2, "RB01-RB02")):
        assert float(r[key]["rate_per_min"]) < 2
    # bad and non-iEEG channels are not analysed
    assert (REC1, "LB01-LB02") not in r and (REC1, "EKG") not in r


def test_spike_times_match_truth(default_run, fixture_dir):
    truth = json.loads((fixture_dir / "truth.json").read_text())[REC1]
    det = np.array([float(s["time_s"]) for s in read_csv(default_run / "spikes.csv")
                    if s["recording_id"] == REC1 and s["channel"] == "LA01-LA02"])
    peaks = np.array(truth["spikes"]) + 0.04  # the template peaks 40 ms after its start
    peaks = peaks[(peaks > 1) & (peaks < 109)]
    err = np.array([np.min(np.abs(det - p)) for p in peaks])
    assert np.mean(err < 0.03) > 0.9
    amps = [float(s["amplitude_uv"]) for s in read_csv(default_run / "spikes.csv") if s["channel"] == "LA01-LA02"]
    assert 250 < np.median(amps) < 600


def test_seizure_excluded(default_run):
    t = [float(s["time_s"]) for s in read_csv(default_run / "spikes.csv") if s["recording_id"] == REC1]
    assert not [x for x in t if 110 <= x <= 160]


def test_events_led_by_first_channel(default_run):
    ev = [e for e in read_csv(default_run / "spike_events.csv") if e["recording_id"] == REC1]
    multi = [e for e in ev if int(e["n_channels"]) == 2]
    assert len(multi) > 40
    assert all(e["leader"] == "LA01-LA02" for e in multi)
    sp = read_csv(default_run / "spikes.csv")
    assert all(s["event_id"] for s in sp)


def test_keep_seizure(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, EXCLUDE_SEIZURE="false", GROUP_EVENTS="false") == 0
    assert not (tmp_path / "spike_events.csv").exists()
    t = [float(s["time_s"]) for s in read_csv(tmp_path / "spikes.csv")
         if s["recording_id"] == REC1 and s["channel"] == "LA03-LA04"]
    assert len([x for x in t if 120 <= x <= 150]) > 10  # Janca adapts to the dense ictal background
    assert all(s["event_id"] == "" for s in read_csv(tmp_path / "spikes.csv"))


def test_segment_baseline_and_label(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, SEGMENT="baseline", LABEL="baseline") == 0
    out = tmp_path / "baseline"
    r = rates(out)
    assert float(r[(REC1, "LA01-LA02")]["minutes"]) == pytest.approx(108 / 60, abs=0.01)
    prov = json.loads((out / "provenance.json").read_text())
    rec2 = [x for x in prov["recordings"] if x["recording_id"] == REC2][0]
    assert any("no baseline interval" in n for n in rec2["notes"])


def test_barkmeier(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, DETECTOR="barkmeier") == 0
    r = rates(tmp_path)
    assert r[(REC1, "LA01-LA02")]["detector"] == "barkmeier"
    assert 20 <= float(r[(REC1, "LA01-LA02")]["rate_per_min"]) <= 36
    assert float(r[(REC1, "LC01-LC02")]["rate_per_min"]) < 3


def test_plain_bids_without_recordings_csv(fixture_dir, tmp_path, monkeypatch):
    src = fixture_dir / "ieeg" / "sub-02"
    shutil.copytree(src, tmp_path / "in" / "sub-02")
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out") == 0
    r = rates(tmp_path / "out")
    assert (REC2, "RA01-RA02") in r  # the _desc-preproc entity is dropped from the recording id


def test_duplicate_recording_prefers_recordings_csv(fixture_dir, tmp_path, monkeypatch):
    shutil.copytree(fixture_dir / "ieeg", tmp_path / "in" / "ieeg")
    shutil.copytree(fixture_dir / "ieeg" / "sub-02", tmp_path / "in" / "raw" / "sub-02")
    recs, _, problems = I.find_recordings(tmp_path / "in")
    assert [r.recording_id for r in recs] == [REC1, REC2]
    assert all(r.source == "recordings.csv" for r in recs)
    assert len(problems) == 1 and "duplicate" in problems[0]["reason"]


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "notes.txt").write_text("x")
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out") == 1
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    assert prov["status"] == "failed: no input"


def test_bad_parameter(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir, tmp_path, DETECTOR="magic") == 2
    assert run(monkeypatch, fixture_dir, tmp_path, WINDOW_S="10") == 2


def test_failed_recording_is_skipped(fixture_dir, tmp_path, monkeypatch):
    shutil.copytree(fixture_dir / "ieeg", tmp_path / "in")
    edf = next((tmp_path / "in" / "sub-02").rglob("*.edf"))
    edf.write_bytes(b"not an edf")
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out") == 0
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    st = {r["recording_id"]: r["status"] for r in prov["recordings"]}
    assert st[REC1] == "ok" and st[REC2] == "failed"


# --- units ------------------------------------------------------------------------

def test_timing_after_resampling_from_1024_hz():
    """A spike at a known time in a 1024 Hz signal is reported within 10 ms (epycom's own decimation is 2.4 % off)."""
    rng = np.random.default_rng(1)
    fs = 1024.0
    x = make_fixture.pink(int(fs * 120), rng, 20.0)
    times = np.arange(5.0, 115.0, 3.7)
    make_fixture.add_at(x, fs, times, make_fixture.spike_shape(fs, 400.0))
    dets, err = D.detect_channel(x.astype(np.float32), fs, "janca", 3.65, 0.0, (0.0, 120.0))
    assert err is None
    t = np.array([d[0] for d in dets])
    late = times[times > 100] + 0.04
    assert all(np.min(np.abs(t - p)) < 0.01 for p in late)


def test_group_events():
    spikes = [(1.000, "A", 100), (1.020, "B", 50), (1.040, "C", 80), (1.060, "D", 90), (2.0, "B", 10),
              (2.0, "A", 20)]
    ev = D.group_events(spikes, 0.05)
    assert [(e["leader"], e["n_channels"]) for e in ev] == [("A", 3), ("D", 1), ("A", 2)]


def test_subtract_and_chunks():
    assert I.subtract([(0, 100)], [(10, 20), (50, 200)]) == [(0, 10), (20, 50)]
    ch = list(I.chunks([(0, 700), (800, 810)], 300, 10, 30, 1))
    assert ch[-1][2] is None  # the 10 s segment is too short
    cores = [(a, b) for a, b, r, _ in ch if r is not None]
    assert len(cores) == 3 and cores[0][0] == 1 and cores[-1][1] == 699
    assert all(b - a <= 300 for a, b in cores)
    assert ch[0][2] == 0 and ch[1][2] == pytest.approx(cores[1][0] - 10)


def test_recording_id_from_prefix():
    p = Path("sub-HUP1_ses-x_task-ictal1_run-01_desc-preproc_ieeg.edf")
    assert I.recording_id_from_prefix(I.signal_prefix(p)) == "sub-HUP1_ses-x_task-ictal1_run-01"
    assert I.signal_prefix(Path("rec.edf")) == "rec"


def test_plan_with_unknown_seizure_end(tmp_path):
    rec = I.Recording("r", tmp_path / "r.edf", "r.edf", "recordings.csv", onset_s=100.0)
    plan = I.analysis_plan(rec, 300.0, "all", True, 10.0)
    assert plan.segments == [(0.0, 90.0)]
    assert "end unknown" in plan.excluded[0]["reason"]
