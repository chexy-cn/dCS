"""运行质量比 q 的单变量敏感性计算。

只改变 q，固定质量、光度距离、观测时间、自旋和倾角。结果写入
`data/<date>/results/q_sensitivity.csv`。
"""

from pathlib import Path
import csv
import itertools
import os
import sys
from multiprocessing import Pool, cpu_count

import numpy as np
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

from getresults import q_sensitivity


def process_task(params):
    """解包 q 敏感性任务参数并调用实际计算函数。"""
    return q_sensitivity.task(params)


def write_results(file_path, results):
    """把 q 敏感性结果写成 CSV 文件。"""
    fieldnames = [
        'detector',
        'q',
        'eta',
        'M_total_Msun',
        'D_l_Mpc',
        'T_tc_year',
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
    """运行所有探测器和质量比组合。"""
    date = 20260506
    detectors = ['aLIGO', 'TianQin', 'CE', 'LISA', 'ET', 'Taiji', 'BBO', 'DECIGO']
    q_values = np.linspace(1.0, 10.0, 91)

    M_total_n = 100.0
    D_l_mpc = 100.0
    T_tc_year = 0.25
    x1 = 0.40
    x2 = 0.20
    iota = 0.9

    param_space = list(itertools.product(
        detectors,
        [M_total_n],
        [D_l_mpc],
        [T_tc_year],
        q_values,
        [x1],
        [x2],
        [iota],
    ))
    ncpu = max(1, min(30, len(param_space)))

    with Pool(processes=ncpu) as pool:
        results = list(tqdm(pool.imap_unordered(process_task, param_space), total=len(param_space)))

    results.sort(key=lambda item: (item['detector'], item['q']))
    write_results(f'./data/{date}/results/q_sensitivity.csv', results)


if __name__ == '__main__':
    main()
