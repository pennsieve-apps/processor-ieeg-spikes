"""
Interictal epileptiform discharge (spike) detectors, from epycom 0.3 (Apache-2.0,
Cimbálník et al., FNUSA-ICRC Brno).

Both detectors run on the channel resampled to 200 Hz (the rate Janca's detector is designed
for and resamples to internally). We resample with a polyphase anti-aliasing filter
ourselves because epycom's internal decimation uses an integer factor while labelling the
result 200 Hz (1024 Hz -> factor 5 -> 204.8 Hz treated as 200 Hz, a 2.4 % timing error).
"""
from __future__ import annotations

import warnings
from fractions import Fraction

import numpy as np
from scipy.signal import butter, filtfilt, resample_poly

DETECTOR_FS = 200
# Janca needs >= 20 s per call (5 s windows, 1 s steps, a 5-point smoothing filtfilt of the
# window statistics); we ask for a little more.
MIN_SECONDS = 30.0


def to_detector_rate(x: np.ndarray, fs: float) -> np.ndarray:
    """Resample one channel to 200 Hz (polyphase FIR); unchanged when already at 200 Hz."""
    x = np.asarray(x, dtype=np.float64)
    fs_int = int(round(fs))
    if fs_int == DETECTOR_FS:
        return x
    frac = Fraction(DETECTOR_FS, fs_int)
    return resample_poly(x, frac.numerator, frac.denominator)


def janca(x200: np.ndarray, k1: float = 3.65, bandwidth=(10, 60)) -> list[tuple[float, float]]:
    """
    Janca et al. 2015 detector (epycom.detect_spikes_janca) on a 200 Hz signal in µV.
    Returns (time_s from the start of x200, signed amplitude in µV of the 1 Hz high-passed
    signal at the discharge).
    """
    from epycom.event_detection.spike.janca_detector import detect_spikes_janca

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # log(0) on flat stretches
        out = detect_spikes_janca(x200, fs=DETECTOR_FS, k1=k1, k2=k1, bandwidth=tuple(bandwidth))
    return [(float(peak) / DETECTOR_FS, float(amp)) for peak, _env, amp in out]


def _broad(x200: np.ndarray) -> np.ndarray:
    """The 1-80 Hz signal Barkmeier's detector works on (its own filters)."""
    fs = DETECTOR_FS
    bh, ah = butter(2, 1 / (fs / 2), "highpass")
    bl, al = butter(4, 80 / (fs / 2), "lowpass")
    return filtfilt(bl, al, filtfilt(bh, ah, x200))


def broadband_level(x: np.ndarray, fs: float) -> float:
    """Mean |1-80 Hz signal| of one channel (µV), for Barkmeier's recording-level scaling."""
    try:
        x200 = to_detector_rate(x, fs)
        return float(np.mean(np.abs(_broad(x200)))) if np.all(np.isfinite(x200)) else float("nan")
    except Exception:  # noqa: BLE001
        return float("nan")


def barkmeier(x200: np.ndarray, level_ratio: float = 1.0, union_s: float = 0.12) -> list[tuple[float, float]]:
    """
    Barkmeier et al. 2012 detector (epycom.detect_spikes_barkmeier, default thresholds) on a
    200 Hz signal in µV.

    The detector rescales the data to a 70 µV reference (epycom: scale /
    np.median(np.mean(np.abs(x))), which reads as a median over channels) and then applies fixed
    amplitude/slope thresholds, so a channel with larger discharges than its neighbours keeps
    them. Given one channel, as epycom's function is, the median is of a single value and
    each channel is normalised to itself, which removes exactly that difference; we pass
    scale = 70 × (this channel's level / the median level of the recording's channels in
    the window), which turns epycom's per-channel scaling back into the paper's
    recording-level one. epycom also compares the gap to the previous detection in samples
    with 0.005 (seconds), so it can report one discharge at neighbouring samples:
    detections closer than union_s (Janca's polyspike union time) are merged here, keeping
    the largest. Returns (time_s, amplitude in µV of the 1-80 Hz signal at the peak).
    """
    from epycom.event_detection.spike.barkmeier_detector import detect_spikes_barkmeier

    fs = DETECTOR_FS
    scale = 70.0 * level_ratio if np.isfinite(level_ratio) and level_ratio > 0 else 70.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out = detect_spikes_barkmeier(x200, fs=fs, scale=scale)
    if not out:
        return []
    broad = _broad(x200)
    peaks = sorted(int(o[0]) for o in out)
    merged: list[int] = []
    for p in peaks:
        if merged and p - merged[-1] <= union_s * fs:
            if broad[p] > broad[merged[-1]]:
                merged[-1] = p
            continue
        merged.append(p)
    return [(p / fs, float(broad[p])) for p in merged]


def detect_channel(x: np.ndarray, fs: float, detector: str, k1: float, offset_s: float,
                   core: tuple[float, float], level_ratio: float = 1.0) -> tuple[list[tuple[float, float]], str | None]:
    """
    Run one detector over one channel window. x starts at offset_s (seconds from the
    recording start); only detections inside core = (start_s, end_s) are returned, as
    (time_s from the recording start, amplitude µV). level_ratio is Barkmeier's
    channel / recording broadband level. Returns (detections, error).
    """
    try:
        x200 = to_detector_rate(x, fs)
        if not np.all(np.isfinite(x200)):
            return [], "non-finite samples"
        if np.ptp(x200) == 0:
            return [], None  # flat channel: nothing to detect
        found = janca(x200, k1=k1) if detector == "janca" else barkmeier(x200, level_ratio)
    except Exception as exc:  # noqa: BLE001 - one channel's failure must not sink the recording
        return [], f"{type(exc).__name__}: {exc}"
    out = []
    for t, amp in found:
        t_abs = offset_s + t
        if core[0] <= t_abs < core[1]:
            out.append((round(t_abs, 4), round(amp, 2)))
    return out, None


def group_events(spikes: list[tuple[float, str, float]], window_s: float) -> list[dict]:
    """
    Cross-channel grouping: spikes sorted by time; an event starts at a spike and takes
    every spike within window_s of that first spike. The leader is the channel that
    fires first (ties: the larger |amplitude|). Returns one dict per event with
    time_s, leader, n_channels, channels, and the index list of its spikes.
    """
    order = sorted(range(len(spikes)), key=lambda i: (spikes[i][0], -abs(spikes[i][2])))
    events: list[dict] = []
    cur: list[int] = []
    for i in order:
        if cur and spikes[i][0] - spikes[cur[0]][0] > window_s:
            events.append(_event(spikes, cur))
            cur = []
        cur.append(i)
    if cur:
        events.append(_event(spikes, cur))
    return events


def _event(spikes, idx: list[int]) -> dict:
    lead = idx[0]
    chans: list[str] = []
    for i in idx:
        if spikes[i][1] not in chans:
            chans.append(spikes[i][1])
    return {"time_s": spikes[lead][0], "leader": spikes[lead][1], "n_channels": len(chans),
            "channels": chans, "members": idx}
