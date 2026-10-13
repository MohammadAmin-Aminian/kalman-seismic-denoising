"""Publication figures (matplotlib, PNG at 200 dpi)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .kalman import ewma_response, steady_state  # noqa: E402
from .spectra import amplitude_spectrum, fk_spectrum  # noqa: E402
from .synthetic import ricker  # noqa: E402

INK = "#0b0b0b"
MUTED = "#8a8986"
COLORS = {  # fixed identity per method (validated categorical order)
    "Kalman": "#2a78d6",
    "Savitzky-Golay": "#eb6834",
    "f-x decon": "#1baf7a",
    "SVD": "#eda100",
    "Kalman RTS": "#6250d6",
    "Clean": INK,
    "Noisy": MUTED,
}

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 200, "font.size": 9, "axes.titlesize": 9.5,
    "axes.labelsize": 9, "axes.edgecolor": "#52514e", "axes.linewidth": 0.6,
    "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.grid": False,
    "legend.frameon": False, "lines.linewidth": 1.2, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


def _gather(ax, d, t, h, clip, title):
    ax.imshow(d, cmap="gray_r", vmin=-clip, vmax=clip, aspect="auto",
              extent=[h[0] / 1000, h[-1] / 1000, t[-1], t[0]], interpolation="bilinear")
    ax.set_title(title)
    ax.set_xlabel("Offset (km)")


def fig_gathers(syn_lo, syn_hi, path):
    t, h = syn_lo["t"], syn_lo["h"]
    clip = 0.5 * np.abs(syn_lo["clean"]).max()
    fig, axs = plt.subplots(1, 3, figsize=(9, 4.2), sharey=True)
    _gather(axs[0], syn_lo["clean"], t, h, clip, "(a) Clean")
    _gather(axs[1], syn_lo["noisy"], t, h, clip, f"(b) Noisy, σ = {syn_lo['sigma']:g}  "
            f"(SNR {syn_lo['snr']['input']:.1f} dB)")
    _gather(axs[2], syn_hi["noisy"], t, h, clip, f"(c) Noisy, σ = {syn_hi['sigma']:g}  "
            f"(SNR {syn_hi['snr']['input']:.1f} dB)")
    axs[0].set_ylabel("Time (s)")
    fig.savefig(path)
    plt.close(fig)


def fig_methods(syn, path):
    t, h = syn["t"], syn["h"]
    clip = 0.5 * np.abs(syn["clean"]).max()
    names = list(syn["out"])
    fig, axs = plt.subplots(2, len(names), figsize=(2.2 * len(names), 6.4), sharex=True, sharey=True)
    for j, m in enumerate(names):
        e = syn["out"][m]
        _gather(axs[0, j], e, t, h, clip, f"{m}\nSNR {syn['snr'][m]:.1f} dB")
        _gather(axs[1, j], syn["noisy"] - e, t, h, clip, "removed (noisy − output)")
    for i in range(2):
        axs[i, 0].set_ylabel("Time (s)")
    for ax in axs[0]:
        ax.set_xlabel("")
    fig.suptitle(f"Denoising of the synthetic gather, σ = {syn['sigma']:g} "
                 f"(input SNR {syn['snr']['input']:.1f} dB)", y=1.0)
    fig.savefig(path)
    plt.close(fig)


def fig_trace(syn, j, path, window=(0.35, 1.25)):
    t = syn["t"]
    sel = (t >= window[0]) & (t <= window[1])
    c = syn["clean"][:, j]
    rows = [("Noisy", syn["noisy"][:, j]), ("Kalman", syn["out"]["Kalman"][:, j]),
            ("Kalman RTS", syn["out"]["Kalman RTS"][:, j]),
            ("Savitzky-Golay", syn["out"]["Savitzky-Golay"][:, j])]
    fig, axs = plt.subplots(len(rows), 1, figsize=(7.5, 6.2), sharex=True, sharey=True)
    for ax, (m, e) in zip(axs, rows):
        ax.plot(t[sel], c[sel], color=INK, lw=1.0, label="Clean")
        ax.plot(t[sel], e[sel], color=COLORS[m], lw=1.2, label=m)
        ax.axhline(0, color="#d8d7d2", lw=0.5, zorder=0)
        ax.legend(loc="upper right", ncol=2, fontsize=8)
        ax.set_ylabel("Amplitude")
    axs[-1].set_xlabel("Time (s)")
    fig.suptitle(f"Trace {j + 1} (offset {syn['h'][j]:.0f} m), σ = {syn['sigma']:g}")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_spectra(syn, dt, fmax, path):
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    f, A = amplitude_spectrum(syn["clean"], dt)
    ax.plot(f, A, color=INK, lw=1.6, label="Clean")
    f, A = amplitude_spectrum(syn["noisy"], dt)
    ax.plot(f, A, color=MUTED, lw=1.0, label="Noisy")
    for m in ("Kalman", "Kalman RTS", "Savitzky-Golay", "f-x decon"):
        f, A = amplitude_spectrum(syn["out"][m], dt)
        ax.plot(f, A, color=COLORS[m], lw=1.1, label=m)
    ax.set_yscale("log")
    ax.set_xlim(0, fmax)
    ax.set_ylim(bottom=1e-2)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Mean amplitude spectrum")
    ax.set_title(f"Amplitude spectra, σ = {syn['sigma']:g}")
    ax.legend(ncol=3, fontsize=8, loc="upper right")
    fig.savefig(path)
    plt.close(fig)


def fig_fk(syn, dt, dx, fmax, path):
    panels = [("(a) Clean", syn["clean"]), ("(b) Noisy", syn["noisy"]),
              ("(c) Kalman", syn["out"]["Kalman"])]
    fig, axs = plt.subplots(1, 3, figsize=(9, 3.6), sharey=True)
    for ax, (ttl, d) in zip(axs, panels):
        k, f, S = fk_spectrum(d, dt, dx)
        im = ax.imshow(S, origin="lower", aspect="auto", cmap="Blues", vmin=-40, vmax=0,
                       extent=[k[0], k[-1], f[0], f[-1]])
        ax.set_ylim(0, fmax)
        ax.set_title(ttl)
        ax.set_xlabel("Wavenumber (1/km)")
    axs[0].set_ylabel("Frequency (Hz)")
    cb = fig.colorbar(im, ax=axs, shrink=0.9, pad=0.015)
    cb.set_label("Amplitude (dB)")
    fig.savefig(path)
    plt.close(fig)


def fig_sweep(sw, sigmas, Q0, path):
    fig, axs = plt.subplots(1, len(sigmas), figsize=(9, 3.5))
    for ax, s in zip(np.atleast_1d(axs), sigmas):
        ax.semilogx(sw["Q"], sw[f"forward_{s:g}"], "-o", ms=3, color=COLORS["Kalman"], label="Kalman (forward)")
        ax.semilogx(sw["Q"], sw[f"rts_{s:g}"], "-s", ms=3, color=COLORS["Kalman RTS"], label="Kalman RTS")
        ax.axhline(sw[f"input_{s:g}"], color=MUTED, ls="--", lw=1, label="Input (no filtering)")
        ax.axvline(Q0, color="#c3c2b7", lw=1)
        ax.annotate(f"Q = {Q0:g} (original)", xy=(Q0, 0.0), xycoords=("data", "axes fraction"),
                    xytext=(3, 4), textcoords="offset points", fontsize=7.5, color="#52514e")
        ax.set_xlabel("Process-noise variance Q (R = 0.05)")
        ax.set_ylabel("Output SNR (dB)")
        ax.set_title(f"σ = {s:g}")
    np.atleast_1d(axs)[0].legend(fontsize=8, loc="center left")
    fig.savefig(path)
    plt.close(fig)


def fig_response(Qs, R, dt, f0, path):
    f = np.linspace(0, 0.5 / dt, 1000)
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    nt = 4096
    w = ricker((np.arange(nt) - nt // 2) * dt, f0)
    W = np.abs(np.fft.rfft(w))
    ax.fill_between(np.fft.rfftfreq(nt, dt), 0, W / W.max(), color="#e8e7e2", label=f"Ricker {f0:g} Hz (normalised)")
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"]
    for c, Q in zip(pal, Qs):
        K = steady_state(Q, R)["K"]
        ax.plot(f, np.abs(ewma_response(K, f, dt)), color=c, label=f"Q = {Q:g}  (K = {K:.3f})")
    ax.set_xlim(0, 0.5 / dt)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("|H(f)| of steady-state filter")
    ax.set_title("The random-walk Kalman filter is a first-order low-pass filter")
    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.savefig(path)
    plt.close(fig)


def fig_chirp(ch, path, windows=((0.0, 0.6), (3.4, 3.6))):
    t = ch["t"]
    fig = plt.figure(figsize=(7.5, 5.8))
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.25)
    for col, (w0, w1) in enumerate(windows):
        ax = fig.add_subplot(gs[0, col])
        sel = (t >= w0) & (t <= w1)
        ax.plot(t[sel], ch["noisy"][sel], color=MUTED, lw=0.6, label="Chirp + noise")
        ax.plot(t[sel], ch["clean"][sel], color=INK, lw=1.0, label="Chirp")
        ax.plot(t[sel], ch["filtered_noisy"][sel], color=COLORS["Kalman"], lw=1.1, label="Kalman output")
        fa, fb = np.interp([w0, w1], t, ch["f_inst"])
        ax.set_title(f"({'ab'[col]}) {w0:g}–{w1:g} s  ({fa:.0f}–{fb:.0f} Hz)")
        ax.set_xlabel("Time (s)")
        ax.set_ylim(-2.4, 2.4)
        if col == 0:
            ax.set_ylabel("Amplitude")
            ax.legend(fontsize=7.5, loc="lower left", ncol=1)
    ax = fig.add_subplot(gs[1, :])
    f = ch["f_inst"]
    sel = (f >= 5) & (f <= 115)
    ax.plot(f[sel], ch["gain_measured"][sel], color=COLORS["Kalman"], lw=1.4, label="Measured (envelope ratio)")
    ax.plot(f[sel], ch["gain_theory"][sel], "--", color=INK, lw=1.0,
            label=r"Theory $|K / (1-(1-K)e^{-i\omega\Delta t})|$")
    ax.axvline(ch["f3db_hz"], color="#c3c2b7", lw=1)
    ax.text(ch["f3db_hz"] + 1.5, 0.9, f"−3 dB at {ch['f3db_hz']:.1f} Hz", fontsize=8, color="#52514e")
    ax.set_xlabel("Instantaneous frequency (Hz)")
    ax.set_ylabel("Amplitude gain")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"(c) Gain of the filter on the noise-free chirp (Q = 0.004, R = 0.05, K = {ch['K']:.3f})")
    ax.legend(fontsize=8)
    fig.savefig(path)
    plt.close(fig)


def fig_earthquake(eq, fmax, path):
    t = eq["t"]
    x = eq["x"]
    fig = plt.figure(figsize=(8, 6.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.15], hspace=0.55)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0, sharey=ax0)
    ax0.plot(t, x, color=MUTED, lw=0.5, label="Raw (demeaned)")
    ax0.plot(t, eq["alt"]["y"], color="#eb6834", lw=0.6, label=f"Kalman Q = {eq['alt']['Q']:g}")
    ax1.plot(t, x, color=MUTED, lw=0.5, label="Raw (demeaned)")
    ax1.plot(t, eq["main"]["y"], color=COLORS["Kalman"], lw=0.7, label=f"Kalman Q = {eq['main']['Q']:g}")
    for ax in (ax0, ax1):
        ax.set_ylabel("Counts")
        ax.legend(ncol=2, fontsize=8, loc="upper right")
        ax.set_xlim(t[0], t[-1])
    ax1.set_xlabel("Time since " + eq["meta"]["start_time"][:19] + " UTC (s)")
    ax0.set_title(f"{eq['meta']['station']} vertical record (40 Hz)")
    ax2 = fig.add_subplot(gs[2])
    dt = eq["dt"]
    for lab, y, c in (("Raw", x, MUTED), (f"Q = {eq['alt']['Q']:g}", eq["alt"]["y"], "#eb6834"),
                      (f"Q = {eq['main']['Q']:g}", eq["main"]["y"], COLORS["Kalman"])):
        f, A = amplitude_spectrum(y, dt)
        ax2.loglog(f[1:], A[1:], color=c, lw=0.8, label=lab)
    for k, c in (("alt", "#eb6834"), ("main", COLORS["Kalman"])):
        ax2.axvline(eq[k]["f3db_hz"], color=c, lw=0.8, ls="--")
    ax2.set_xlim(1e-2, 0.5 / dt)
    ax2.set_xlabel("Frequency (Hz)  (dashed: −3 dB cut-off of each filter)")
    ax2.set_ylabel("Amplitude")
    ax2.legend(fontsize=8, loc="lower left")
    fig.savefig(path)
    plt.close(fig)
