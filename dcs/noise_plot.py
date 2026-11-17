"""Paper Figure 1: amplitude spectral density, not characteristic strain."""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dcs.fim import get_detecror_parameters
from dcs.scans import DETECTORS
from dcs.io import run_dir, atomic_csv, provenance


def run_noise(date=None, results_root="results"):
    directory = run_dir(date, results_root)
    fig, ax = plt.subplots(figsize=(8.5, 5))
    rows = []
    for detector in DETECTORS:
        low, high, _, _, _, _, noise, *_ = get_detecror_parameters(detector)
        frequency = np.geomspace(low, high, 1500)
        asd = np.sqrt(noise(frequency))
        ax.loglog(frequency, asd, label=detector)
        rows.extend(dict(detector=detector, frequency_Hz=float(f), ASD_Hz_minus_half=float(s))
                    for f, s in zip(frequency, asd))
    ax.set(xlabel="Frequency [Hz]", ylabel=r"ASD [$\mathrm{Hz}^{-1/2}$]")
    ax.set_ylim(1e-25, 1e-13)
    ax.legend(ncol=2)
    ax.grid(alpha=.25)
    fig.tight_layout()
    for extension in ("pdf", "png"):
        fig.savefig(directory / "fig" / f"detector_sensitivity.{extension}", dpi=250)
    plt.close(fig)
    path = directory / "data" / "detector_sensitivity.csv"
    atomic_csv(path, rows, list(rows[0]))
    provenance(path.with_suffix(".provenance.json"), dict(plot="ASD", detectors=DETECTORS),
               sorted((Path(__file__).resolve().parents[1] / "Sn/data").glob("*.txt")))
