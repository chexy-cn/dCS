"""Deterministic scan grids shared by local and Slurm entry points."""

from itertools import product
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm

from dcs.io import atomic_csv, provenance, run_dir
from getresults import q_sensitivity, spin_sensitivity

DETECTORS = ["aLIGO", "TianQin", "CE", "LISA", "ET", "Taiji", "BBO", "DECIGO"]


def evaluate(task):
    scan, detector, value, value2 = task
    if scan == "component":
        cs, ca = (value + value2) / 2, (value - value2) / 2
        return spin_sensitivity.get_result(detector, 100., 100., .25, 1.2, cs, ca, .9)
    return q_sensitivity.get_result(
        detector, 100., 100., value if scan == "tc" else .25,
        value if scan == "q" else 1.2, .4, .2, value if scan == "inclination" else .9,
    )


def grid(scan, detectors, points=None, quick=False):
    if scan == "tc":
        values = [1/365.25, 7/365.25, 1/12, .25, .5, 1., 2., 3., 4., 5.]
        if quick:
            values = [.25, 1.]
    elif scan == "q":
        values = np.linspace(1.2, 8.9, points or (2 if quick else 78))
    elif scan == "inclination":
        values = np.linspace(0., np.pi/2, points or (3 if quick else 41))
    elif scan == "component":
        values = np.linspace(-.9, .9, points or (3 if quick else 41))
        return [(scan, d, float(x), float(y)) for d, x, y in product(detectors, values, values)]
    else:
        raise ValueError(f"Unknown figure scan: {scan}")
    return [(scan, d, float(x), None) for d, x in product(detectors, values)]


def run_scan(scan, detectors, date=None, results_root="results", workers=1,
             points=None, quick=False, resume=False):
    directory = run_dir(date, results_root)
    filename = "spin_component_sensitivity" if scan == "component" else f"{scan}_sensitivity"
    path = directory / "data" / f"{filename}.csv"
    configuration = dict(scan=scan, detectors=detectors, points=points, quick=quick,
                         workers=workers, resume=resume, source=dict(M_Msun=100, Dl_Mpc=100,
                         q=1.2, T_year=.25, x1=.4, x2=.2, iota=.9))
    tasks = grid(scan, detectors, points, quick)
    existing = []
    if path.exists():
        if not resume:
            raise FileExistsError(f"{path} exists; choose another date or use --resume")
        import json
        previous_path = path.with_suffix('.provenance.json')
        if not previous_path.exists():
            raise ValueError('Resume requires the original run provenance')
        previous = json.loads(previous_path.read_text())['configuration']
        for name in ('scan', 'detectors', 'points', 'quick', 'source'):
            if previous.get(name) != configuration.get(name):
                raise ValueError(f'Resume configuration mismatch: {name}; use a new run directory')
        import csv
        existing = list(csv.DictReader(path.open()))
        def key(row):
            if scan == "component":
                return (row['detector'], round(float(row['x1']), 12), round(float(row['x2']), 12))
            return (row['detector'], round(float(row[{'tc': 'T_tc_year', 'q': 'q', 'inclination': 'iota'}[scan]]), 12), None)
        completed = {key(row) for row in existing}
        tasks = [t for t in tasks if (t[1], round(t[2], 12), None if t[3] is None else round(t[3], 12)) not in completed]
    provenance(path.with_suffix(".provenance.json"), configuration)
    results = []
    if tasks:
        with Pool(max(1, min(workers, len(tasks)))) as pool:
            for result in tqdm(pool.imap(evaluate, tasks), total=len(tasks)):
                if result is not None:
                    results.append(result)
    results = existing + results
    if not results:
        raise RuntimeError("No scan points produced a result")
    fields = list(results[0])
    results.sort(key=lambda r: (r["detector"], float(r.get("x1", 0)),
                               float(r.get("x2", 0)), float(r.get("q", 0)),
                               float(r.get("T_tc_year", 0)), float(r.get("iota", 0))))
    atomic_csv(path, results, fields)
    print(f"Saved {len(results)} points: {path}")
    return path
