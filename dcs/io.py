"""Portable result layout, atomic files and per-run provenance."""

import csv
from datetime import datetime
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import subprocess
import tempfile
from zoneinfo import ZoneInfo


def today():
    return datetime.now(ZoneInfo("Asia/Shanghai")).strftime("%Y%m%d")


def run_dir(date=None, results_root="results"):
    date = date or today()
    if len(date) != 8 or not date.isdigit():
        raise ValueError("date must be YYYYMMDD")
    datetime.strptime(date, "%Y%m%d")
    path = Path(results_root) / date
    (path / "data").mkdir(parents=True, exist_ok=True)
    (path / "fig").mkdir(exist_ok=True)
    return path


def atomic_csv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", newline="", dir=path.parent,
                                         delete=False, encoding="utf-8") as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance(path, configuration, inputs=()):
    """Store exact configuration, dependency versions, git state and input hashes."""
    packages = {}
    for package in ("numpy", "scipy", "matplotlib", "lalsuite", "tqdm"):
        try:
            packages[package] = version(package)
        except PackageNotFoundError:
            packages[package] = None
    root = Path(__file__).resolve().parents[1]
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root,
                                         stderr=subprocess.DEVNULL, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"],
                                            cwd=root, text=True).strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    record = {
        "created_at": datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(),
        "configuration": configuration, "python": platform.python_version(),
        "packages": packages, "git_commit": commit, "git_dirty": dirty,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "inputs": [{"file": str(p), "sha256": sha256(p)} for p in inputs],
        "source_sha256": {
            str(p.relative_to(root)): sha256(p)
            for folder in ("dcs", "getresults", "Sn", "script")
            for p in sorted((root / folder).glob("*.py"))
        },
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(record, indent=2, default=str) + "\n", encoding="utf-8")
