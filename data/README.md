# Data

All files are MATLAB v5 `.mat` files. They load with `load` in base MATLAB
and with `scipy.io.loadmat` in Python. They were built by
[`tools/make_data.py`](../tools/make_data.py).

| File | Content | Size |
|---|---|---|
| `synthetic_noise.mat` | `noise`: 1501 × 100 unit-variance Gaussian noise (single precision, seed 20210926). The noisy gathers are `clean + σ · noise`. | 0.6 MB |
| `chirp_noise.mat` | `noise`: 2001 × 1 unit-variance Gaussian noise for the chirp test (same generator, drawn after the field above). | 8 kB |
| `earthquake_IU_DAV_BHZ.mat` | `counts`: 26,400 samples (660 s at 40 Hz) of a vertical broadband record, raw counts (no instrument correction); `fs_hz`, `station`, `start_time`, `units`, `source`. | 64 kB |

The clean synthetic gather isn't stored. It is computed exactly from the
config (`hyperbolic_gather`), so it is identical in both languages.

## Provenance of the earthquake record

The record is variable `sig1` of the original project file `Sigs.mat`. The
original script `Kalman_Earth.m` read the text file
`IU.DAV.00.BHZ.M.2021-09-05T152300.019538.txt`. The station code
(IU.DAV, location 00, channel BHZ: Davao, Philippines) and the start time
(2021-09-05 15:23:00.0195 UTC) are taken from that file name. The sampling
rate (40 Hz) matches the BH channel and the 660 s record length used in
`kalmancode.m`.

The data come from the Global Seismograph Network, distributed openly by the
EarthScope Consortium (formerly IRIS) Data Management Center. If you use them,
please cite the network:

> Albuquerque Seismological Laboratory/USGS (1988). *Global Seismograph Network
> (GSN – IRIS/USGS)* [Data set]. International Federation of Digital
> Seismograph Networks. https://doi.org/10.7914/SN/IU

The variable `sig2` of `Sigs.mat` (int16) is almost a copy of `sig1`
(correlation 0.986) and could not be identified, so it is not included.

## Not included

* The real shot gathers `shot1.sgy` and `shot20.sgy` used by
  `legacy/Kalman_Real.m` and `legacy/kalman_Real_22_sep.m` were not available.
* `data.rar` from the original project (nine SEISAN event files from
  2010-12-24) isn't used by any script.
