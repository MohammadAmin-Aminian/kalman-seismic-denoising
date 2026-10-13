"""Quality metrics."""
from __future__ import annotations

import numpy as np

__all__ = ["snr_db", "snr_db_original"]


def snr_db(clean, estimate):
    """SNR of an estimate against the known clean signal (dB).

    SNR = 10 log10( sum(clean^2) / sum((estimate - clean)^2) ).
    """
    c = np.asarray(clean, dtype=float)
    e = np.asarray(estimate, dtype=float)
    return 10.0 * np.log10(np.sum(c ** 2) / np.sum((e - c) ** 2))


def snr_db_original(clean, estimate):
    """SNR as computed in the original scripts (kept only for comparison).

    10 log10( sum(estimate^2) / sum((estimate - clean)^2) ) -- the numerator
    uses the estimate instead of the clean signal, which rewards outputs
    with extra energy and is therefore not a valid quality measure.
    """
    c = np.asarray(clean, dtype=float)
    e = np.asarray(estimate, dtype=float)
    return 10.0 * np.log10(np.sum(e ** 2) / np.sum((e - c) ** 2))
