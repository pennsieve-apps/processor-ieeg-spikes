# sample-input

Real input for this app is the `ieeg` output of **ieeg-prep**: a BIDS-iEEG derivative with
`recordings.csv`, and per recording `*_desc-preproc_ieeg.edf`, `_channels.tsv`,
`_events.tsv` (see `docs/conventions.md` in workflow-epilepsy-ieeg). Any BIDS-iEEG dataset
or loose EDF/BrainVision/FIF files work too.

The workflows' demo data are public on Pennsieve Discover:

> Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. Multi-Expert Seizure
> Annotation. Pennsieve Discover 530 (2026), doi:10.26275/dely-fzwx (CC BY 4.0).
> BIDS-iEEG SEEG seizure clips at 1024 Hz with 120 s before the clinically marked onset.

> Bernabei J, Litt B. HUP Interictal iEEG Atlas. Pennsieve Discover 179 (2021),
> doi:10.26275/lzhy-ox20 (ODC-By). 60 s wake clips at 200 Hz of 3,431 bipolar channels of
> 60 patients, with a 24 h spike count per channel (`spike_24h`).

```sh
cd ../workflow-epilepsy-ieeg && scripts/fetch-sample-data.sh        # -> data/seizures, data/atlas
# run ieeg-prep on data/seizures (-> data/work/prep/ieeg), then:
cd ../processor-ieeg-spikes && make run-demo DEMO=../workflow-epilepsy-ieeg/data/work/prep/ieeg
```

The seizure clips also run directly as a plain BIDS dataset (`ictal` events are treated as
seizures): mount `data/seizures/DATA` as the input. For the 120 s pre-ictal baseline, set
`SEGMENT=baseline` on an ieeg-prep output (its `recordings.csv` carries the baseline), or
`SEIZURE_MARGIN_S=0` on the raw clips (everything before the onset is analysed; the
post-ictal part after the annotated end is analysed too unless the seizure has no end).

For tests and `make run`, `tests/make_fixture.py` writes a small synthetic derivative
instead (`make fixture` -> `data/fixture/ieeg/`): a 512 Hz bipolar SEEG clip with planted
spikes on two channels (one leading the other by 10 ms), a seizure, a bad channel and an
ECG channel, and a 200 Hz clip without labels (`data/fixture/truth.json` holds the planted
times).
