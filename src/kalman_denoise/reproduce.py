"""Regenerate every figure and number of the study.

    python scripts/reproduce_all.py            # or: kalman-denoise-reproduce
    python scripts/reproduce_all.py --no-figures
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from . import pipeline as P
from .io import load_config, repo_root


def _table(rows, sigmas):
    head = "| Method | " + " | ".join(f"σ = {s:g}" for s in sigmas) + " |"
    sep = "|---|" + "---:|" * len(sigmas)
    lines = [head, sep]
    for name, vals in rows:
        lines.append(f"| {name} | " + " | ".join(f"{v:.2f}" for v in vals) + " |")
    return "\n".join(lines)


def run(argv=None):
    """Run everything and return the metrics dictionary."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--out", default=None, help="output root (default: repository root)")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    root = Path(a.out) if a.out else repo_root()
    (root / "results").mkdir(parents=True, exist_ok=True)
    (root / "figures").mkdir(parents=True, exist_ok=True)

    sigmas = cfg["noise"]["sigma_levels"]
    syn = {s: P.run_synthetic(cfg, s) for s in sigmas}
    sweep = P.run_sweep(cfg)
    chirp = P.run_chirp(cfg)
    eq = P.run_earthquake(cfg)
    mo = P.run_modelled(cfg)
    fd = P.run_field(cfg)

    metrics = {
        "synthetic_snr_db": {f"{s:g}": syn[s]["snr"] for s in sigmas},
        "synthetic_snr_db_original_definition": {f"{s:g}": syn[s]["snr_original_definition"] for s in sigmas},
        "kalman_Q_sweep": {
            "Q": sweep["Q"].tolist(),
            **{f"best_forward_{s:g}": {"Q": float(sweep["Q"][np.argmax(sweep[f"forward_{s:g}"])]),
                                       "snr_db": float(np.max(sweep[f"forward_{s:g}"]))} for s in sigmas},
            **{f"best_rts_{s:g}": {"Q": float(sweep["Q"][np.argmax(sweep[f"rts_{s:g}"])]),
                                   "snr_db": float(np.max(sweep[f"rts_{s:g}"]))} for s in sigmas},
        },
        "chirp": {"K": chirp["K"], "f3db_hz": chirp["f3db_hz"], "max_gain_error": chirp["max_gain_error"]},
        "earthquake": {k: {kk: float(vv) for kk, vv in eq[k].items() if kk != "y"} for k in ("main", "alt")},
        "modelled_shot": {
            "sigma_x_rms": mo["sigma"], "Q": mo["Q"],
            "snr_db": {k: float(v) for k, v in mo["snr"].items()},
            "best_forward": {"Q": float(mo["sweep_Q"][np.argmax(mo["sweep_forward"])]),
                             "snr_db": float(np.max(mo["sweep_forward"]))},
            "best_rts": {"Q": float(mo["sweep_Q"][np.argmax(mo["sweep_rts"])]),
                         "snr_db": float(np.max(mo["sweep_rts"]))},
        },
        "field_record": {"Q": fd["Q"], "K": float(fd["K"]), "f3db_hz": float(fd["f3db_hz"]),
                         "removed_energy_fraction": fd["removed_fraction"]},
    }
    with open(root / "results" / "metrics.json", "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, indent=2, default=float)

    rows = [("Input (noisy)", [syn[s]["snr"]["input"] for s in sigmas])]
    rows += [(m, [syn[s]["snr"][m] for s in sigmas]) for m in P.METHODS]
    rows += [(f"Kalman, best Q (sweep)", [metrics["kalman_Q_sweep"][f"best_forward_{s:g}"]["snr_db"] for s in sigmas])]
    rows += [(f"Kalman RTS, best Q (sweep)", [metrics["kalman_Q_sweep"][f"best_rts_{s:g}"]["snr_db"] for s in sigmas])]
    table = _table(rows, sigmas)
    mrows = [("Input (noisy)", mo["snr"]["input"])] + [(m, mo["snr"][m]) for m in P.METHODS]
    mrows += [("Kalman, best Q (sweep)", metrics["modelled_shot"]["best_forward"]["snr_db"]),
              ("Kalman RTS, best Q (sweep)", metrics["modelled_shot"]["best_rts"]["snr_db"])]
    mtable = "| Method | Output SNR (dB) |\n|---|---:|\n" + "\n".join(f"| {n} | {v:.2f} |" for n, v in mrows)
    ftable = "| Method | Energy removed |\n|---|---:|\n" + "\n".join(
        f"| {m} | {100 * v:.0f} % |" for m, v in fd["removed_fraction"].items())
    (root / "results" / "snr_table.md").write_text(
        "Output SNR (dB) = 10 log10(Σ clean² / Σ (output − clean)²).\n\n## Synthetic gather\n\n" + table
        + f"\n\n## Modelled shot gather (noise σ = {mo['sigma']:g} × RMS, Kalman Q = {mo['Q']:g})\n\n" + mtable
        + f"\n\n## Field record FFID 3043 (no ground truth; Kalman Q = {fd['Q']:g})\n\n" + ftable + "\n",
        encoding="utf-8")
    print(table)
    print(mtable)
    print(ftable)

    if not a.no_figures:
        from . import plotting as G
        fdir = root / "figures"
        s_lo, s_hi = sigmas[0], sigmas[-1]
        dt = cfg["synthetic"]["dt_s"]
        G.fig_gathers(syn[s_lo], syn[s_hi], fdir / "fig01_synthetic_gathers.png")
        G.fig_methods(syn[s_hi], fdir / "fig02_methods_sigma_high.png")
        G.fig_methods(syn[s_lo], fdir / "fig03_methods_sigma_low.png")
        G.fig_trace(syn[s_hi], cfg["display"]["trace_index"] - 1, fdir / "fig04_trace_comparison.png")
        G.fig_spectra(syn[s_hi], dt, cfg["display"]["fmax_hz"], fdir / "fig05_amplitude_spectra.png")
        G.fig_fk(syn[s_hi], dt, cfg["synthetic"]["offset_step_m"], cfg["display"]["fmax_hz"],
                 fdir / "fig06_fk_spectra.png")
        G.fig_sweep(sweep, sigmas, cfg["kalman"]["Q"], fdir / "fig07_snr_vs_Q.png")
        G.fig_response([1e-4, cfg["kalman"]["Q"], 0.05, 0.5], cfg["kalman"]["R"], dt,
                       cfg["synthetic"]["f0_hz"], fdir / "fig08_filter_response.png")
        G.fig_chirp(chirp, fdir / "fig09_chirp_validation.png")
        G.fig_earthquake(eq, cfg["earthquake"]["fmax_hz"], fdir / "fig10_earthquake.png")
        G.fig_modelled(mo, fdir / "fig11_modelled_shot.png")
        G.fig_field(fd, fdir / "fig12_field_record.png")
        G.fig_field_spectra(fd, fdir / "fig13_field_spectra.png")
        print("figures written to", fdir)
    return metrics


def main(argv=None):
    """Command-line entry point (exit status 0 on success)."""
    run(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
