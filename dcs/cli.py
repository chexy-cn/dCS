"""Public CLI; computational commands belong on an allocated compute node."""

import argparse
import os

from dcs.io import run_dir, today


def main(argv=None):
    parser = argparse.ArgumentParser(prog="dcs-gw")
    commands = parser.add_subparsers(dest="command", required=True)
    scan = commands.add_parser("scan", help="Recompute fixed-source sensitivity scans")
    scan.add_argument("scan", choices=["tc", "q", "component", "inclination", "mass", "all"])
    scan.add_argument("--detectors", nargs="+", default=None)
    scan.add_argument("--workers", type=int, default=int(os.environ.get("DCS_NCPU", "1")))
    scan.add_argument("--points", type=int)
    scan.add_argument("--quick", action="store_true")
    scan.add_argument("--resume", action="store_true")
    plot = commands.add_parser("plot", help="Replot saved data without rerunning Fisher")
    plot.add_argument("--data-dir", help="Directory containing newly computed figure CSVs")
    commands.add_parser("noise", help="Plot detector ASD curves")
    for name in ("scan", "plot", "noise"):
        commands.choices[name].add_argument("--date", default=today())
        commands.choices[name].add_argument("--results-root", default="results")
    args, remaining = parser.parse_known_args(argv)
    if remaining:
        parser.error(f"Unrecognized arguments: {' '.join(remaining)}")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    os.environ["OMP_PROC_BIND"] = "false"
    if args.command == "scan":
        if args.workers < 1 or (args.points is not None and args.points < 2):
            parser.error("workers >= 1 and points >= 2 required")
        from dcs.scans import DETECTORS, run_scan
        detectors = args.detectors or DETECTORS
        if set(detectors).difference(DETECTORS):
            parser.error("Unknown detector name")
        scans = ["tc", "q", "component", "inclination", "mass"] if args.scan == "all" else [args.scan]
        for scan_name in scans:
            if scan_name == "mass":
                from dcs.mass_scan import run_mass
                run_mass(args.date, args.results_root, args.workers, args.points or 1001, args.quick, args.detectors)
            else:
                run_scan(scan_name, detectors, args.date, args.results_root, args.workers,
                         args.points, args.quick, args.resume)
    elif args.command == "plot":
        from getresults.plot_sensitivity import plot_all
        from dcs.io import provenance
        from pathlib import Path
        directory = run_dir(args.date, args.results_root)
        data = Path(args.data_dir) if args.data_dir else directory / "data"
        if not data.is_dir() or not list(data.glob("*.csv")):
            parser.error(f"No CSV data in {data}")
        plot_all(args.date, args.results_root, str(data))
        mass_path = data / "mass_sensitivity.csv"
        if mass_path.exists():
            import csv
            from dcs.mass_scan import plot_mass
            with mass_path.open(newline="") as handle:
                plot_mass(list(csv.DictReader(handle)), directory / "fig")
        provenance(directory / "data" / "plot.provenance.json", vars(args), sorted(data.glob("*.csv")))
    elif args.command == "noise":
        from dcs.noise_plot import run_noise
        run_noise(args.date, args.results_root)
