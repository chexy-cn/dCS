"""批量扫描探测器参数并绘制约束等高线。

脚本入口会并行计算单探测器与 TianQin+地面探测器联合观测结果，
写入 `data/<date>/results/`，随后调用绘图函数生成 PDF。
"""

from pathlib import Path
import os
import sys

_PROJECT_ROOT = next(
    path for path in Path(__file__).resolve().parents
    if (path / 'getresults').is_dir() and (path / 'math' / 'fim.py').is_file()
)
for _path in (str(_PROJECT_ROOT), str(_PROJECT_ROOT / 'math')):
    if _path not in sys.path:
        sys.path.insert(0, _path)
os.chdir(_PROJECT_ROOT)

from tqdm import tqdm
import numpy as np
import os
os.environ["OMP_NUM_THREADS"] = "1"      # 控制numba、scikit-learn中的OpenMP
# os.environ["MKL_NUM_THREADS"] = "1"      # 控制Intel MKL (numpy, scipy可能用)
# os.environ["OPENBLAS_NUM_THREADS"] = "1"    # 控制OpenBLAS (如果编译时选择)
from multiprocessing import Pool, cpu_count
import itertools
from getresults import onedetector, twodetectors, plot
from getresults import plot

def process_task1(params):
    """解包单探测器任务参数并调用实际计算函数。"""
    return onedetector.task(*params)

def process_task2(params):
    """解包双探测器任务参数并调用实际计算函数。"""
    return twodetectors.task(*params)

def main1():
    """批量计算多个单探测器的质量网格约束并绘图。"""
    # 计算单个探测器
    detectors = ['aLIGO','TianQin', 'CE', 'LISA', 'ET', 'Taiji','BBO','DECIGO']
    # detectors = ['CE']
    q = [1.2, 8.9]
    param_space = list(itertools.product(detectors, range(1001), q))

    total_tasks = len(param_space)
    with Pool(processes = ncpu) as pool:
        results = list(tqdm(pool.imap_unordered(process_task1, param_space), total=total_tasks))

    for result in results:
        p_detector, p_l, p_q, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = result
        result_str = ','.join(map(str, [p_l, p_q, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR]))
        directory = f'./data/{date}/results/'
        os.makedirs(directory, exist_ok=True)
        file_path = os.path.join(directory, f'results for {p_detector} with q={p_q}.txt')
        with open(file_path, 'a') as f:
            f.write(result_str + '\n')
    
    # for detector in detectors:
    #     plot.plot_one_simplify(detector, q, date)
    plot.plot_all_detectors_grid(detectors, q, date)


def main2():
    """批量计算 TianQin 与地面探测器联合观测的质量网格约束并绘图。"""
    # 计算TianQin与地面的联合探测
    detectors1 = ['CE', 'ET']
    detectors2 = ['TianQin']
    q = [1.2, 8.9]
    param_space = list(itertools.product(detectors1, detectors2, range(1001), q))
    total_tasks = len(param_space)

    with Pool(processes = ncpu) as pool:
        results = list(tqdm(pool.imap_unordered(process_task2, param_space), total=total_tasks))

    for result in results:
        p_detector1, p_detector2, p_l, p_q, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = result
        result_str = ','.join(map(str, [p_l, p_q, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR]))
        directory = f'./data/{date}/results/'
        os.makedirs(directory, exist_ok=True)
        file_path = os.path.join(directory, f'results for {p_detector1}+{p_detector2} with q={p_q}.txt')
        with open(file_path, 'a') as f:
            f.write(result_str + '\n')

    # for detector in detectors1:
    #     plot.plot_two_simplify(detector, date)
    plot.plot_two_detectors_side_by_side(detectors1, date)

if __name__ == '__main__':
    date = 20251217
    ncpu = cpu_count()-10
    
    main1()
    main2()
    
