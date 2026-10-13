"""Minimal SEG-Y rev 0/1 reader used once by tools/make_data.py.

Supports fixed-length traces, 4-byte IBM or IEEE samples, either byte order.
The byte order and the sample format can be forced, because the original
modelled file (shot20.sgy) is little-endian and declares format code 1 (IBM)
while its samples are in fact IEEE floats.
"""
from __future__ import annotations

import struct

import numpy as np


def ibm_to_float(u):
    """Convert 32-bit IBM floats (as uint32) to float64."""
    u = np.asarray(u, dtype=np.uint32)
    sign = 1.0 - 2.0 * ((u >> 31) & 1)
    expo = ((u >> 24) & 0x7F).astype(np.int32) - 64
    mant = (u & 0x00FFFFFF).astype(np.float64) / float(1 << 24)
    return np.where(u == 0, 0.0, sign * mant * np.power(16.0, expo))


def read_segy(path, endian=">", sample_format=None):
    """Return (data[ns, ntraces], headers dict, binary-header dict, text header)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    text = raw[:3200]
    try:
        text = text.decode("cp500") if text[:1] == b"\xc3" else text.decode("ascii", "replace")
    except Exception:  # pragma: no cover
        text = ""
    bh = raw[3200:3600]
    dt_us = struct.unpack(endian + "H", bh[16:18])[0]
    ns = struct.unpack(endian + "H", bh[20:22])[0]
    fmt = struct.unpack(endian + "H", bh[24:26])[0] if sample_format is None else sample_format
    nb = 240 + 4 * ns
    ntr = (len(raw) - 3600) // nb
    block = np.frombuffer(raw[3600:3600 + ntr * nb], dtype=np.uint8).reshape(ntr, nb)
    samples = block[:, 240:].copy()
    if fmt == 1:
        data = ibm_to_float(samples.view(endian + "u4"))
    elif fmt == 5:
        data = samples.view(endian + "f4").astype(np.float64)
    else:
        raise ValueError(f"unsupported sample format {fmt}")
    data = data.reshape(ntr, ns).T

    def field(off, code="i", size=4):
        return np.array([struct.unpack(endian + code, bytes(block[i, off:off + size]))[0] for i in range(ntr)])

    headers = {
        "fldr": field(8), "tracf": field(12), "cdp": field(20), "offset": field(36),
        "scalco": field(70, "h", 2), "sx": field(72), "sy": field(76), "gx": field(80), "gy": field(84),
    }
    return data, headers, {"dt_us": dt_us, "ns": ns, "format": fmt}, text
