"""事件样本的 SNR 批处理入口。

本脚本读取事件参数文件，计算单探测器或 TianQin+地面探测器联合
观测 SNR，并筛选达到 `SNR >= 8` 的事件。
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

import numpy as np
import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
import os
os.environ["OMP_NUM_THREADS"] = "1"
import itertools
from multiprocessing import Pool, cpu_count ,TimeoutError
from getresults import snr, readresult
import time as time
import csv
import psutil 
import logging
logging.basicConfig(
    filename=os.path.expanduser('~/events_snr.log'),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def process_task1(params):
    """解包双探测器 SNR 任务参数并调用实际计算函数。"""
    return snr.task_two(*params)

def process_task2(params):
    """解包单探测器 SNR 任务参数并调用实际计算函数。"""
    return snr.task_one(*params)

def generate_parameter_space1(p_detectors1, p_detectors2, p_file_path):
    """根据事件文件生成双探测器 SNR 任务参数列表。"""
    p_parameters = readresult.read_events(p_file_path)
    parameter_space = []
    for p_parameter in p_parameters:
        parameter_space.extend(list(itertools.product(p_detectors1, p_detectors2, [p_parameter])))
    return parameter_space

def generate_parameter_space2(p_detector, p_file_path):
    """根据事件文件生成单探测器 SNR 任务参数列表。"""
    p_parameters = readresult.read_events(p_file_path)
    parameter_space = []
    for p_parameter in p_parameters:
        parameter_space.extend(list(itertools.product(p_detector, [p_parameter])))
    return parameter_space

def main1(p_param_space, i):
    """并行执行双探测器 SNR 任务并写入结果文件。"""
    with Pool(processes=ncpu) as pool:
        results = list(pool.imap_unordered(process_task1, p_param_space))

    for result in results:
        p_detector1, p_detector2, p_parameters, SNR = result
        result_str = ','.join(map(str, [*p_parameters, SNR]))
        directory = f'/home/chexy/code/dCS/data/{date}/results/'
        os.makedirs(directory, exist_ok=True)
        file_path = os.path.join(directory, f'results for {p_detector2}+{p_detector1}_{i}.txt')
        with open(file_path, 'a') as f:
            f.write(result_str + '\n')

def main2(p_param_space, i):
    """并行执行单探测器 SNR 任务并写入结果文件。"""
    with Pool(processes = ncpu) as pool:
        results = list(pool.imap_unordered(process_task2, p_param_space))

    for result in results:
        p_detector, p_parameters, SNR = result
        result_str = ','.join(map(str, [*p_parameters, SNR]))
        directory = f'/home/chexy/code/dCS/data/{date}/results/'
        os.makedirs(directory, exist_ok=True)
        file_path = os.path.join(directory, f'results for {p_detector}_{i}.txt')
        with open(file_path, 'a') as f:
            f.write(result_str + '\n')

def validity_events(p_detectors, i):
        """从结果文件中筛选 SNR 不低于 8 的事件。"""
        for p_detector in p_detectors:
            with open(f'/home/chexy/code/dCS/data/{date}/results/results for {p_detector}_{i}.txt', 'r') as f:
                lines = f.readlines()
            # 将列表转换为NumPy数组
            Dl = np.array([float(line.strip().split(',')[0]) for line in lines])
            m1 = np.array([float(line.strip().split(',')[1]) for line in lines])
            m2 = np.array([float(line.strip().split(',')[2]) for line in lines])
            Tc = np.array([float(line.strip().split(',')[3]) for line in lines])
            theta_s = np.array([float(line.strip().split(',')[4]) for line in lines])
            SNR = np.array([float(line.strip().split(',')[5]) for line in lines])
            mask = (SNR >= 8)
            filter_Dl = Dl[mask]
            filter_m1 = m1[mask]
            filter_m2 = m2[mask]
            filter_Tc = Tc[mask]
            filter_theta_s = theta_s[mask]
            filter_SNR = SNR[mask]
            directory = f'/home/chexy/code/dCS/data/{date}/validity/'
            os.makedirs(directory, exist_ok=True)
            file_path = os.path.join(directory, f'{p_detector}_{i}.txt')
            file_path2 = os.path.join(directory, f'n_{p_detector}.txt')
            with open(file_path2, 'a') as f:
                f.write(f'{i}\t{np.sum(mask)}\n')
            with open(file_path, 'w') as f:
                for i in range(len(filter_Dl)):
                    line = f"{filter_Dl[i]}\t{filter_m1[i]}\t{filter_m2[i]}\t{filter_Tc[i]}\t{filter_theta_s[i]}\t{filter_SNR[i]}\n"
                    f.write(line)


def events_TianQin_CE(i):
    """处理第 i 批 TianQin+CE 事件 SNR。"""
    detectors1 = ['CE']
    detectors2 = ['TianQin']
    file_path = f'/data/chexy/P2P/parameters-final/TianQin+CE/parameters-{i}.txt'
    param_space = generate_parameter_space1(detectors1, detectors2, file_path)
    main1(param_space, i)
    validity_events(['TianQin+CE'], i)


def events_TianQin_ET(i):
    """处理第 i 批 TianQin+ET 事件 SNR。"""
    detectors1 = ['ET']
    detectors2 = ['TianQin']
    file_path = f'/data/chexy/P2P/parameters-final/TianQin+ET/parameters-{i}.txt'
    param_space = generate_parameter_space1(detectors1, detectors2, file_path)
    main1(param_space, i)
    validity_events(['TianQin+ET'], i)


def events_BBO(i):
    """处理第 i 批 BBO 事件 SNR。"""
    detector = ['BBO']
    file_path = f'/data/chexy/P2P/parameters-final/BBO/parameters-{i}.txt'
    param_space = generate_parameter_space2(detector, file_path)
    main2(param_space, i)
    validity_events(['BBO'], i)


def events_DECIGO(i):
    """处理第 i 批 DECIGO 事件 SNR。"""
    detector = ['DECIGO']
    file_path = f'/data/chexy/P2P/parameters-final/DECIGO/parameters-{i}.txt'
    param_space = generate_parameter_space2(detector, file_path)
    main2(param_space, i)
    validity_events(['DECIGO'], i)


def events_CE(i):
    """处理第 i 批 CE 事件 SNR。"""
    detector = ['CE']
    file_path = f'/data/chexy/P2P/parameters-final/CE/parameters-{i}.txt'
    param_space = generate_parameter_space2(detector, file_path)
    main2(param_space, i)
    validity_events(['CE'], i)


def events_ET(i):
    """处理第 i 批 ET 事件 SNR。"""
    detector = ['ET']
    file_path = f'/data/chexy/P2P/parameters-final/ET/parameters-{i}.txt'
    param_space = generate_parameter_space2(detector, file_path)
    main2(param_space, i)
    validity_events(['ET'], i)


if __name__ == '__main__':
    date = 20251226

    ncpu = max(1, cpu_count() - 10)
    # ncpu=25
    p = psutil.Process()
    p.cpu_affinity(list(range(ncpu)))


    with open('log-events_snr.csv', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteration', 'Function', 'Execution_Time(s)'])
        
        for i in range(100):
            functions = [
                events_CE, 
                # events_ET, 
                # events_TianQin_CE, 
                # events_TianQin_ET,
                # events_DECIGO, 
                # events_BBO
            ]
            
            for func in functions:
                print(f'Starting {func.__name__} for iteration {i}')
                start_time = time.time()
                func(i)
                end_time = time.time()
                print(f'Completed {func.__name__} for iteration {i}','time:', end_time - start_time)
                # 计算执行时间并写入文件
                execution_time = end_time - start_time
                writer.writerow([i, func.__name__, execution_time])
