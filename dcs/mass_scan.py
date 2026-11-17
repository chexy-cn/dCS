"""Mass-dependent horizons (paper Figures 2 and 3), evaluated at 1 Mpc."""

from itertools import product
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from tqdm import tqdm

from dcs.fim import covariance, covariance_two, get_detecror_parameters, snr, snr_two
from dcs.io import atomic_csv, provenance, run_dir
from dcs.scans import DETECTORS
from getresults.q_sensitivity import build_parameters, validity_scale_km


def mass_point(task):
    channel, mass, q = task
    detectors = channel.split("+")
    # Historical mass scans use the ground detector's mass support/steps.
    if len(detectors) == 2:
        detectors = [detectors[1], detectors[0]]
    p1, s1, delta, _, n1 = build_parameters(detectors[0], mass, 1., .25, q)
    if p1["f_min"] >= p1["f_max"]:
        return None
    if len(detectors) == 1:
        cov = covariance(p1, delta, n1)
        rho = snr(s1, n1)
    else:
        p2, s2, _, _, n2 = build_parameters(detectors[1], mass, 1., .25, q)
        # Keep the original multiband finite-difference convention explicit.
        delta["delta_eta"] = get_detecror_parameters(detectors[0])[8]
        delta["delta_x1"] = get_detecror_parameters(detectors[0])[9]
        delta["delta_x2"] = get_detecror_parameters(detectors[0])[10]
        cov = covariance_two(p1, p2, delta, n1, n2)
        rho = np.sqrt(snr_two(s1, n1, n2)**2 + snr_two(s2, n1, n2)**2)
    if cov is None or not np.isfinite(cov[7, 7].real) or cov[7, 7].real <= 0:
        return None
    return dict(detector=channel, q=q, M_total_Msun=mass, D_l_Mpc=1.,
                d_sqrt_alpha_dCS_km=float(np.sqrt(cov[7, 7].real)**.25 * 1e-3),
                alpha_dCS_bar_square_max_km=validity_scale_km(mass/(1+q)),
                SNR=float(rho))


def run_mass(date=None, results_root="results", workers=1, points=1001, quick=False, detectors=None):
    directory = run_dir(date, results_root)
    if quick:
        points = 3
    path = directory / "data" / "mass_sensitivity.csv"
    if path.exists():
        raise FileExistsError(f"{path} exists")
    tasks = []
    channels = detectors or DETECTORS + ["TianQin+CE", "TianQin+ET"]
    for channel, q in product(channels, [1.2, 8.9]):
        detector = channel.split("+")[-1]
        low, high = get_detecror_parameters(detector)[2:4]
        tasks.extend((channel, float(m), q) for m in np.geomspace(low, high, points))
    provenance(path.with_suffix(".provenance.json"), dict(scan="mass", points=points, detectors=channels,
               workers=workers, reference_distance_Mpc=1, T_year=.25, x1=.4, x2=.2, iota=.9))
    with Pool(workers) as pool:
        rows = [row for row in tqdm(pool.imap(mass_point, tasks), total=len(tasks)) if row is not None]
    if not rows:
        raise RuntimeError("No finite mass points")
    atomic_csv(path, rows, list(rows[0]))
    plot_mass(rows, directory / "fig")


def plot_mass(rows, output):
    """Exact distance scaling replaces the old million-cell contour search.

    rho scales as 1/D; the reported length scales as D**(1/4).
    """
    groups = [(["aLIGO", "TianQin", "CE", "LISA", "ET", "Taiji", "BBO", "DECIGO"],
               (4, 2), "all_detectors"),
              (["TianQin+CE", "TianQin+ET"], (1, 2), "TianQin_combined")]
    for channels, shape, name in groups:
        fig, axes = plt.subplots(*shape, figsize=(12, 3.7*shape[0]), squeeze=False)
        for ax, channel in zip(axes.flat, channels):
            for q, style in ((1.2, "-"), (8.9, "--")):
                selected = sorted((r for r in rows if r["detector"] == channel and float(r["q"]) == q),
                                  key=lambda r: float(r["M_total_Msun"]))
                if not selected:
                    continue
                m = np.array([float(r["M_total_Msun"]) for r in selected])
                constraint = np.array([float(r["d_sqrt_alpha_dCS_km"]) for r in selected])
                scale = np.array([float(r["alpha_dCS_bar_square_max_km"]) for r in selected])
                rho = np.array([float(r["SNR"]) for r in selected])
                for y, color in ((rho/8, "royalblue"), ((scale/constraint)**4, "black"),
                                 ((8.5/constraint)**4, "seagreen"), ((1/constraint)**4, "firebrick")):
                    ax.plot(m, y, style, color=color, linewidth=1.2)
            ax.set(xscale="log", yscale="log", title=channel, xlabel=r"$M\ [M_\odot]$", ylabel=r"$D_L\ [\mathrm{Mpc}]$")
            ax.set_xlim(10 if channel in {"TianQin", "LISA", "Taiji"} else 1, 1000)
            ax.set_ylim(1, 1e4 if channel in {"TianQin", "LISA", "Taiji"} else 1e6)
            ax.grid(alpha=.25)
        handles = [Line2D([], [], color=c, label=l) for c, l in
                   [("royalblue", "SNR = 8"), ("black", "Validity"), ("seagreen", "8.5 km"), ("firebrick", "1 km")]]
        axes.flat[0].legend(handles=handles, fontsize=9)
        fig.tight_layout()
        for extension in ("pdf", "png"):
            fig.savefig(output / f"{name}.{extension}", dpi=200, bbox_inches="tight")
        plt.close(fig)
