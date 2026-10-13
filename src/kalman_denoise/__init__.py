"""Kalman-filter random-noise attenuation of seismic data.

Main entry points: :func:`kalman_denoise.kalman.kalman1d`,
:func:`kalman_denoise.kalman.kalman_rts` and
:func:`kalman_denoise.kalman.steady_state`.
"""
from .kalman import cutoff_frequency, ewma_response, kalman1d, kalman_rts, steady_state  # noqa: F401

__version__ = "1.0.0"
