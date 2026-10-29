# processor-ieeg-spikes

A step of the **Interictal Biomarker Mapping** workflow
(`ieeg-prep → {ieeg-spikes, ieeg-bandpower → ieeg-normative, ieeg-hfo} → ieeg-report`). It
finds interictal epileptiform discharges (IEDs, "spikes") on every good channel of an
intracranial EEG recording (SEEG or ECoG) and counts them per channel per minute.

**For the epileptologist.** Interictal spikes are most frequent in and around the
epileptogenic zone, and contacts with high spike rates are a long-standing part of the
evaluation of a surgical candidate. Marking them by hand is slow and reviewers disagree
(Barkmeier et al. 2012 report this inter-reviewer variability). This step marks them
automatically with a validated detector, so every channel gets a rate (spikes/min) over the
same interictal time, seizures left out, and cross-channel discharges are grouped into one
event with the channel that fires first (the "leader"), which is where a propagating
discharge most likely starts.

**For the researcher.** The default detector is Janca et al. 2015 (envelope distribution
modelling of the 10–60 Hz band, adaptive to the background in 5 s windows), as implemented in
epycom 0.3 (FNUSA-ICRC Brno), on the signal resampled to 200 Hz. Barkmeier et al. 2012 (slope
and amplitude criteria) is available for comparison. Outputs are long-format CSV keyed by
`recording_id` and `channel` (workflow-epilepsy-ieeg `docs/conventions.md`), ready to join
with band power, HFO rates, normative abnormality and labels (`soz`, `resected`) from
`_channels.tsv`.

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `ieeg` | input | root | The BIDS-iEEG derivative of ieeg-prep: `recordings.csv`, `*_desc-preproc_ieeg.edf` with `_channels.tsv`, `_events.tsv`. Also any BIDS-iEEG dataset, or plain EDF/BDF/BrainVision/FIF/EEGLAB files, found anywhere under `INPUT_DIR` (folders and symlinks followed). |
| `spikes` | output | `OUTPUT_DIR/` (or `OUTPUT_DIR/<LABEL>/`) | `spikes.csv`, `spike_rates.csv`, `spike_events.csv`, `provenance.json` |

## Parameters

| Name | Default | What it does |
|---|---|---|
| `DETECTOR` | `janca` | `janca` (Janca et al. 2015) or `barkmeier` (Barkmeier et al. 2012). |
| `JANCA_K1` | `3.65` | Janca's threshold factor k1 (the published default). Higher: fewer, more obvious spikes. |
| `SEGMENT` | `all` | `all`: the whole recording. `baseline`: only the baseline interval of `recordings.csv` (`baseline_start_s`–`baseline_end_s`) or a `baseline` event, e.g. the pre-ictal baseline of a seizure clip. Recordings without one are analysed whole (noted in `provenance.json`). |
| `EXCLUDE_SEIZURE` | `true` | Leave known seizures out: `onset_s` of `recordings.csv` and `seizure` rows of `_events.tsv` (also `ictal`, for plain BIDS datasets), with `SEIZURE_MARGIN_S` on both sides. A seizure without a known end is cut to the end of the recording. `artifact` events are always left out. |
| `SEIZURE_MARGIN_S` | `10` | Seconds left out before the onset and after the end of each seizure. |
| `WINDOW_S` | `300` | Windows the recording is read and analysed in: all good channels at once, plus 10 s of context on each side (never across excluded time). Bounds memory: channels × (`WINDOW_S` + 20 s) × sampling rate × 4 bytes. |
| `GROUP_EVENTS` | `true` | Group spikes of different channels into events (`spike_events.csv`, `event_id` in `spikes.csv`). |
| `GROUP_WINDOW_MS` | `50` | An event takes every spike within this many ms of its first spike. |
| `THREADS` | `0` | Channels analysed in parallel (processes); `0` = one per vCPU. |
| `LABEL` | (empty) | When set, write to `OUTPUT_DIR/<LABEL>/` (a workflow running this step twice, e.g. once per detector). |

## Outputs

**`spikes.csv`**: one row per detected spike.

| Column | |
|---|---|
| `recording_id` | the BIDS entity chain of the source recording (`sub-HUP257_ses-postimplant_task-ictal387681_run-01`) |
| `channel` | the `_channels.tsv` name (`LA01-LA02`) |
| `time_s` | the spike's peak, seconds from the recording start (4 decimals; 5 ms resolution at 200 Hz) |
| `amplitude_uv` | signed peak amplitude in µV: Janca, of the 1 Hz high-passed signal at the discharge; Barkmeier, of its 1–80 Hz signal at the peak |
| `detector` | `janca` or `barkmeier` |
| `event_id` | the cross-channel event it belongs to (per recording, from 1; empty when `GROUP_EVENTS=false`) |

**`spike_rates.csv`**: one row per analysed channel (good channels of an iEEG type), zeros included.

| Column | |
|---|---|
| `recording_id`, `channel` | as above |
| `n_spikes` | spikes counted |
| `minutes` | minutes of that channel actually analysed (excluded seizures and artifacts, and the first and last second of each analysed stretch, which Janca cannot see, are not counted) |
| `rate_per_min` | `n_spikes / minutes` |
| `detector` | |

Bad channels (`status=bad`) and non-iEEG channels (ECG, EMG, EOG, MISC, …) get no row; why
each was left out is in `provenance.json`.

**`spike_events.csv`**: one row per cross-channel event: `recording_id`, `event_id`, `time_s`
(the leader's spike), `leader` (the channel that fires first; ties go to the larger
amplitude), `n_channels`, `channels` (`|`-separated, in firing order), `detector`. An event
starts at a spike and takes every spike within `GROUP_WINDOW_MS` (50 ms) of it; isolated
spikes are events with `n_channels = 1`.

**`provenance.json`**: app and version, tool versions (epycom, MNE, NumPy, SciPy, joblib,
numba, Python), parameters, every input file seen, timings, and per recording: status
(`ok`, `failed` with the error, `skipped: …`), sampling rate, channel counts, the analysed
segments and what was excluded and why, minutes, spike and event counts, channels not
analysed with the reason, notes.

## How it works

1. **Finds the recordings**: every `recordings.csv` (its `path` column, relative to the CSV),
   plus any signal file no CSV lists, with the sidecars next to it. Two recordings with the
   same id: the one in a `recordings.csv` wins (a derivative next to its source dataset).
2. **Per recording** (one at a time): opens it lazily with MNE; picks good channels of an
   iEEG type from `_channels.tsv` (without one: MNE's EEG/SEEG/ECoG channels whose names do
   not look like ECG/EMG/EOG/…); works out the analysed time (`SEGMENT`, minus seizures and
   artifacts) and splits it into `WINDOW_S` windows.
3. **Per window**: reads all good channels (float32 µV), and analyses each channel in a
   worker process: resample to 200 Hz with a polyphase anti-aliasing filter, run the
   detector, keep detections inside the window's core.
4. **Groups** spikes across channels and writes the rows.

Why we resample ourselves: epycom's Janca implementation decimates by an integer factor but
labels the result 200 Hz (1024 Hz ÷ 5 = 204.8 Hz), which shifts detections by 2.4 % (seconds
over a few minutes; `tests/test_processor.py::test_timing_after_resampling_from_1024_hz`).
Given a 200 Hz signal it does not decimate. We also keep each call under 600 s (epycom's own
buffering of longer signals refers to a field that does not exist) and use `k2 = k1` (its
ambiguous-spike branch is unfinished).

Barkmeier's detector rescales the data to a reference amplitude of 70 µV and then applies
fixed slope and amplitude thresholds. epycom's scaling, `scale / np.median(np.mean(np.abs(x)))`,
reads as a median over channels of each channel's mean |amplitude| (the paper describes a
multi-channel algorithm), but epycom's function takes one channel, so every channel is
normalised to itself, which removes the amplitude differences between channels that the fixed
thresholds would otherwise see. We
pass a per-channel `scale` that turns this back into a recording-level scaling (median over
the good channels, per window), and merge detections closer than 120 ms (epycom compares a
gap in samples with 0.005 s and can report one discharge at neighbouring samples). On the
atlas this raised Barkmeier's agreement with `spike_24h` from ρ = 0.01 to 0.07.

## Runtime

`computeTypes: [standard]`, **2 vCPU / 4 GiB**, 4 h timeout. Channels run in parallel
(`THREADS`, default one per vCPU); memory is bounded by the window, not the recording length.

Measured in Docker (linux/amd64 **emulated** on an M-series Mac, `--cpus 2 --memory 4g`) on 4
seizure clips of the demo dataset (Multi-Expert Seizure Annotation, 1024 Hz, as-recorded
referential montage: HUP224 230 channels × 338 s, HUP230 106 × 384 s, HUP250 170 × 296 s,
HUP257 118 × 463 s; 439 MB of EDF):

| Run | Analysed | CPU time | Wall time | Peak memory (cgroup `memory.peak`) |
|---|---|---|---|---|
| whole recordings (`EXCLUDE_SEIZURE=false`) | 24.6 min, 3,710 channel-minutes | 113 s | 78 s | 1.71 GB |
| 120 s pre-ictal baselines (`SEGMENT=baseline`) | 7.9 min, 1,232 channel-minutes | 53 s | 34 s | 1.43 GB |

`memory.peak` includes the page cache of the EDF files read; the largest window (230 channels
× 320 s × 1024 Hz) is 300 MB as float32. About 0.03 CPU-seconds per channel-minute
(emulated), so a 1 h, 128-channel recording at 1024 Hz takes ~2 min on 2 vCPU and ~0.6 GB for
its window. Recordings of 256+ channels at 2 kHz fit in 4 GiB with the default window; lower
`WINDOW_S` if they do not. The image is 199 MB compressed (829 MB unpacked).

## Results on real data

**HUP interictal atlas** (Bernabei et al. 2022; Discover 179): 60 patients, 3,431 bipolar
channels, one 60 s wake clip each at 200 Hz (0.5–80 Hz), converted to a derivative with the
atlas' `spike_24h` (spikes counted over 24 h) in `_electrodes.tsv`. Spike rate per channel
from the 60 s clip vs `spike_24h`:

| Detector | Spearman ρ (all channels) | per patient, median ρ (IQR) | SOZ vs other channels, median rate |
|---|---|---|---|
| `janca` | **0.35** (p = 3e-98, n = 3,431) | 0.35 (0.22–0.51); ρ > 0 in 56 of 58 | 2.07 vs 0.00 /min (Mann-Whitney p = 3e-30) |
| `barkmeier` | 0.07 (p = 5e-5) | 0.14 (−0.04–0.26); ρ > 0 in 42 of 59 | 2.07 vs 1.03 /min (p = 5e-9) |

A 60 s clip is a small sample of a 24 h count, so a moderate correlation is what to expect;
Janca's ranks agree, Barkmeier's barely do (it flags a similar number of events on every
channel; see Limitations).

**Seizure clips, 120 s pre-ictal baseline** (`SEGMENT=baseline SEIZURE_MARGIN_S=0`, the 4
clips above, referential montage): baseline spike rate in the channels of the experts'
consensus unequivocal electrographic onset (`ueo_consensus`) vs all other channels:

| Recording | UEO channels | Janca: UEO vs others (median /min), AUC | Barkmeier AUC |
|---|---|---|---|
| HUP224 ictal339135 | 6 | 1.02 vs 0.00, **0.83** | 0.81 |
| HUP230 ictal296207 | 3 | 41.7 vs 0.00, **0.99** | 0.88 |
| HUP250 ictal209153 | 2 | 5.08 vs 0.00, **0.99** | 0.64 |
| HUP257 ictal387681 | 5 | 6.61 vs 0.00, **0.98** | 0.92 |

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python        # unit tests on a synthetic fixture
make run PYTHON=.venv/bin/python         # Docker over the fixture -> data/output/
make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg   # an ieeg-prep output
```

Or by hand:

```sh
docker build --platform linux/amd64 -t processor-ieeg-spikes:dev .
docker run --rm --platform linux/amd64 --cpus 2 --memory 4g --env-file dev.env \
  -v /path/to/ieeg:/data/input:ro -v "$PWD/data/output":/data/output processor-ieeg-spikes:dev
```

## Limitations

- A detector, not a diagnosis. Spike rates rank channels; they do not by themselves delimit
  the epileptogenic zone, and high rates also occur in irritative cortex outside it.
- Janca's thresholds adapt to each channel's background in 5 s windows: a channel that is
  spiking continuously (e.g. during a seizure, or periodic discharges) raises its own
  background and reports fewer spikes than it has. Leave seizures out (the default).
- Montage matters: the counts depend on how the channel was referenced. ieeg-prep's
  bipolar montage is the intended input; on a common reference, a spike at the reference
  shows on every channel.
- The detectors look for spikes and sharp waves; they do not classify them (no
  spike-and-wave vs polyspike, no morphology), and they cannot tell an epileptic discharge
  from a sharp artifact that looks like one. Mark artifacts as `artifact` events.
- epycom's Janca notches 50 Hz (and harmonics up to 66 Hz) internally, whatever the power
  line frequency; 60 Hz line noise must be removed upstream (ieeg-prep does this).
- Barkmeier: epycom also accepts peaks whose slopes and amplitude are all *below* the
  thresholds (a second condition the paper does not have); it matters little in practice
  (≤ 0.1 % of detections on the atlas), but together with per-channel narrow-band
  thresholds it makes Barkmeier rates nearly flat across channels (ρ = 0.07 with `spike_24h`).
  Use it for comparison, not as the primary rate.
- The first and last second of each analysed stretch, and stretches shorter than 30 s, are
  not analysed (Janca needs ≥ 20 s of signal); the minutes in `spike_rates.csv` say exactly
  how much time each rate is based on.
- Cross-channel grouping is purely temporal (50 ms from the first spike); it does not check
  that the channels are neighbours.

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles epycom (Apache-2.0), MNE-Python
(BSD-3-Clause) and their dependencies; see `THIRD_PARTY.md`. Nothing GPL is linked in-process.

## Citations

Spike detectors:

> Janca R, Jezdik P, Cmejla R, Tomasek M, Worrell GA, Stead M, Wagenaar J, Jefferys JGR,
> Krsek P, Komarek V, Jiruska P, Marusic P. **Detection of Interictal Epileptiform Discharges
> Using Signal Envelope Distribution Modelling: Application to Epileptic and Non-Epileptic
> Intracranial Recordings.** Brain Topogr 28(1):172–183 (2015). doi:10.1007/s10548-014-0379-1,
> PMID 24970691.

> Barkmeier DT, Shah AK, Flanagan D, Atkinson MD, Agarwal R, Fuerst DR, Jafari-Khouzani K,
> Loeb JA. **High inter-reviewer variability of spike detection on intracranial EEG addressed
> by an automated multi-channel algorithm.** Clin Neurophysiol 123(6):1088–1095 (2012).
> doi:10.1016/j.clinph.2011.09.023, PMID 22033028.

Implementation: epycom 0.3, ElectroPhYsiology COmputational Module, FNUSA-ICRC Brno,
https://github.com/ICRC-BME/epycom (Apache-2.0).

Reading the recordings:

> Gramfort A, Luessi M, Larson E, Engemann DA, Strohmeier D, Brodbeck C, Goj R, Jas M,
> Brooks T, Parkkonen L, Hämäläinen M. **MEG and EEG data analysis with MNE-Python.** Front
> Neurosci 7:267 (2013). doi:10.3389/fnins.2013.00267, PMID 24431986.

> Holdgraf C, Appelhoff S, Bickel S, et al. (39 authors; all in `CITATION.cff`).
> **iEEG-BIDS, extending the Brain Imaging Data Structure specification to human
> intracranial electrophysiology.** Sci Data 6:102 (2019). doi:10.1038/s41597-019-0105-7.

Validation data:

> Bernabei JM, Sinha N, Arnold TC, Conrad E, Ong I, Pattnaik AR, Stein JM, Shinohara RT,
> Lucas TH, Bassett DS, Davis KA, Litt B. **Normative intracranial EEG maps epileptogenic
> tissues in focal epilepsy.** Brain 145(6):1949–1961 (2022). doi:10.1093/brain/awab480,
> PMID 35640886. Data: Bernabei J, Litt B. HUP Interictal iEEG Atlas. Pennsieve Discover
> (2021), doi:10.26275/lzhy-ox20 (ODC-By).

> Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. **Multi-Expert Seizure
> Annotation.** Pennsieve Discover (2026), doi:10.26275/dely-fzwx (CC BY 4.0).

Every DOI was checked against Crossref (articles) or DataCite (datasets): title, author list
and year; PMIDs against PubMed. Crossref lists only the first author of the MNE-Python paper;
its full author list here and in `CITATION.cff` is PubMed's.
