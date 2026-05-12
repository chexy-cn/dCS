"""运行观测时间 T_tc 的单变量敏感性计算。

只改变 T_tc，固定质量、光度距离、质量比、自旋和倾角。结果写入
`data/<date>/results/tc_sensitivity.csv`。
"""

from pathlib import Path
import csv
import itertools
import os
import sys
from multiprocessing import Pool, cpu_count

from tqdm import tqdm

_PROJECT_ROOT = next(
    path for path in Path(__file__).resolve().parents
    if (path / 'getresults').is_dir() and (path / 'math' / 'fim.py').is_file()
)
for _path in (str(_PROJECT_ROOT), str(_PROJECT_ROOT / 'math')):
    if _path not in sys.path:
        sys.path.insert(0, _path)
os.chdir(_PROJECT_ROOT)

os.environ["OMP_NUM_THREADS"] = "1"

from getresults import tc_sensitivity


def process_task(params):
    """解包 T_tc 敏感性任务参数并调用实际计算函数。"""
    return tc_sensitivity.task(params)


def write_results(file_path, results):
    """把 T_tc 敏感性结果写成 CSV 文件。"""
    fieldnames = [
        'detector',
        'T_tc_year',
        'M_total_Msun',
        'D_l_Mpc',
        'q',
        'eta',
        'x1',
        'x2',
        'iota',
        'f_low_Hz',
        'f_high_Hz',
        'f_pn_low_Hz',
        'f_isco_Hz',
        'deltaF_Hz',
        'sigma_alpha_dCS_bar_square',
        'd_sqrt_alpha_dCS_km',
        'alpha_dCS_bar_square_max_km',
        'SNR',
    ]
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


def main():
    """运行所有探测器和观测时间组合。"""
    date = 20260506
    detectors = ['aLIGO', 'TianQin', 'CE', 'LISA', 'ET', 'Taiji', 'BBO', 'DECIGO']
    T_tc_years = [
        1 / 365.25,
        7 / 365.25,
        1 / 12,
        0.25,
        0.5,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
    ]

    M_total_n = 100.0
    D_l_mpc = 100.0
    q = 1.2
    x1 = 0.40
    x2 = 0.20
    iota = 0.9

    param_space = list(itertools.product(
        detectors,
        [M_total_n],
        [D_l_mpc],
        [q],
        T_tc_years,
        [x1],
        [x2],
        [iota],
    ))
    ncpu = max(1, min(30, len(param_space)))

    with Pool(processes=ncpu) as pool:
        results = list(tqdm(pool.imap_unordered(process_task, param_space), total=len(param_space)))

    results.sort(key=lambda item: (item['detector'], item['T_tc_year']))
    write_results(f'./data/{date}/results/tc_sensitivity.csv', results)


if __name__ == '__main__':
    main()
