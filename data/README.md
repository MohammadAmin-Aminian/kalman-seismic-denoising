# Data

All files are MATLAB v5 `.mat` files. They load with `load` in base MATLAB
and with `scipy.io.loadmat` in Python. They were built by
[`tools/make_data.py`](../tools/make_data.py).

| File | Content | Size |
|---|---|---|
| `synthetic_noise.mat` | `noise`: 1501 × 100 unit-variance Gaussian noise (single precision, seed 20210926). The noisy gathers are `clean + σ · noise`. | 0.6 MB |
| `chirp_noise.mat` | `noise`: 2001 × 1 unit-variance Gaussian noise for the chirp test (same generator, drawn after the field above). | 8 kB |
| `earthquake_IU_DAV_BHZ.mat` | `counts`: 26,400 samples (660 s at 40 Hz) of a vertical broadband record, raw counts (no instrument correction); `fs_hz`, `station`, `start_time`, `units`, `source`. | 64 kB |
| `modelled_shot20.mat` | `d`: 1000 × 120 noise-free modelled shot gather with ground roll (shot 20 of the original `shot20.sgy`, 2 ms); `offset_m` (−1140 … 1260 m, 20 m spacing), `dt_s`, `source_x_m`. | 0.45 MB |
| `modelled_noise.mat` | `noise`: 1000 × 120 unit-variance Gaussian noise for the modelled shot (seed 20210927). | 0.44 MB |
| `field_shot_ffid3043.mat` | `d`: 751 × 447 field shot record FFID 3043 from the original `shot1.sgy` (4 ms, IBM floats as recorded); `offset_m` (−4762 … 6383 m, 25 m spacing), `dt_s`, `ffid`. | 0.89 MB |

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

## Shot gathers

Both come from the original project and were converted with
`tools/make_data.py` / `tools/segy.py`. Only the record used is stored.

* `shot20.sgy` holds 20 modelled shots (120 traces × 1000 samples). The file
  is little-endian, and its binary header declares IBM floats although the
  samples are IEEE floats; the converter handles this. Shot 20 (source at
  x = 2240 m) is the one stored.
* `shot1.sgy` holds field records FFID 3043 (447 traces) and 3044 (28
  traces). Its text header gives area Tang-e-Homam, line ns02 (2012). FFID
  3043 is stored as recorded, without the text header.

## Not included

* The second field record (FFID 3044, 28 traces) and the other 19 modelled shots.
* `data.rar` from the original project (nine SEISAN event files from
  2010-12-24) isn't used by any script.
