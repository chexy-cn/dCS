"""事件样本的 dCS 约束批处理入口。

本脚本读取服务器上的事件参数文件，分别计算单探测器和
TianQin+地面探测器联合观测的 dCS 约束，并按有效性条件筛选事件。
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
import itertools
from multiprocessing import Pool, cpu_count
import os
os.environ["OMP_NUM_THREADS"] = "1"      # 控制numba、scikit-learn中的OpenMP
# os.environ["MKL_NUM_THREADS"] = "1"      # 控制Intel MKL (numpy, scipy可能用)
# os.environ["OPENBLAS_NUM_THREADS"] = "1"    # 控制OpenBLAS (如果编译时选择)
from getresults import filter, readresult, plot
import time as time
import csv

def process_task1(params):
    """解包双探测器事件任务参数并调用约束计算。"""
    return filter.task_two(*params)

def process_task2(params):
    """解包单探测器事件任务参数并调用约束计算。"""
    return filter.task_one(*params)

def generate_parameter_space1(p_detectors1, p_detectors2, p_file_path):
    """根据事件文件生成双探测器联合观测的任务参数列表。"""
    p_parameters = readresult.read_events(p_file_path)
    parameter_space = []
    for p_parameter in p_parameters:
        parameter_space.extend(list(itertools.product(p_detectors1, p_detectors2, [p_parameter])))
    return parameter_space

def generate_parameter_space2(p_detector, p_file_path):
    """根据事件文件生成单探测器观测的任务参数列表。"""
    p_parameters = readresult.read_events(p_file_path)
    parameter_space = []
    for p_parameter in p_parameters:
        parameter_space.extend(list(itertools.product(p_detector, [p_parameter])))
    return parameter_space

def main1(p_param_space, i):
    """并行执行双探测器事件约束任务并写入结果文件。"""
    # n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', cpu_count()))
    n_cores = min(len(p_param_space), 150)
    with Pool(processes = n_cores) as pool:
        results = list(pool.imap_unordered(process_task1, p_param_space))

    for result in results:
        if result is None:
            print(f"Warning: Skipping a None result")
            continue  # 跳过这个结果，继续处理下一个
        
        try:
            p_detector1, p_detector2, p_parameters, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = result
            
            # 确保有足够的元素解包
            if len(array_delta) < 8:  # 根据你的需求调整这个数字
                print(f"Warning: array_delta has only {len(array_delta)} elements, expected at least 8")
                continue
            
            result_str = ','.join(map(str, [*p_parameters, array_delta[7], alpha_dCS_bar_square_max, SNR]))
            directory = f'/home/{name}/code/dCS/data/{date}/results/'
            os.makedirs(directory, exist_ok=True)
            file_path = os.path.join(directory, f'results for {p_detector2}+{p_detector1}-{i}.txt')
            with open(file_path, 'a') as f:
                f.write(result_str + '\n')
        
        except (TypeError, ValueError, IndexError) as e:
            print(f"Error processing result: {e}")
            print(f"Result data: {result}")
            continue

def main2(p_param_space, i):
    """并行执行单探测器事件约束任务并写入结果文件。"""
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', cpu_count()))
    with Pool(processes = n_cores) as pool:
        results = list(pool.imap_unordered(process_task2, p_param_space))

    for result in results:
        if result is None:
            print(f"Warning: Skipping a None result")
            continue  # 跳过这个结果，继续处理下一个
        try:
            p_detector, p_parameters, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = result
            # 确保有足够的元素解包
            if len(array_delta) < 8:  # 根据你的需求调整这个数字
                print(f"Warning: array_delta has only {len(array_delta)} elements, expected at least 8")
                continue
            result_str = ','.join(map(str, [*p_parameters, array_delta[7], alpha_dCS_bar_square_max, SNR]))
            directory = f'/home/{name}/code/dCS/data/{date}/results/'
            os.makedirs(directory, exist_ok=True)
            file_path = os.path.join(directory, f'results for {p_detector}-{i}.txt')
            with open(file_path, 'a') as f:
                f.write(result_str + '\n')
        except (TypeError, ValueError, IndexError) as e:
            print(f"Error processing result: {e}")
            print(f"Result data: {result}")
            continue

def validity_events(p_detectors, i):
        """筛选同时满足有效性、SNR 和 8.5 km 约束的事件。"""
        for p_detector in p_detectors:
            try:
                with open(f'/home/{name}/code/dCS/data/{date}/results/results for {p_detector}-{i}.txt', 'r') as f:
                    lines = f.readlines()
            except:
                print(f"No result for {p_detector} iteration {i}")
                return
            # 将列表转换为NumPy数组
            Dl = np.array([float(line.strip().split(',')[0]) for line in lines])
            m1 = np.array([float(line.strip().split(',')[1]) for line in lines])
            m2 = np.array([float(line.strip().split(',')[2]) for line in lines])
            Tc = np.array([float(line.strip().split(',')[3]) for line in lines])
            theta_s = np.array([float(line.strip().split(',')[4]) for line in lines])
            d_alpha_dCS_bar_square = np.array([float(line.strip().split(',')[5])**0.25*1e-3 for line in lines])
            # alpha_dCS_bar_square_max = np.array([((float(line.strip().split(',')[2])*MSUN_SI)**4 * G**4 / c**8 / 32 / pi)**0.25*1e-3 for line in lines])
            alpha_dCS_bar_square_max = np.array([float(line.strip().split(',')[6]) for line in lines])
            SNR = np.array([float(line.strip().split(',')[7]) for line in lines])
            mask = (d_alpha_dCS_bar_square < alpha_dCS_bar_square_max) & (SNR > 8) & (d_alpha_dCS_bar_square < 8.5)
            filter_Dl = Dl[mask]
            filter_m1 = m1[mask]
            filter_m2 = m2[mask]
            filter_Tc = Tc[mask]
            filter_theta_s = theta_s[mask]
            filter_d_alpha_dCS_bar_square = d_alpha_dCS_bar_square[mask]
            directory = f'/home/{name}/code/dCS/data/{date}/validity3/'
            os.makedirs(directory, exist_ok=True)
            file_path = os.path.join(directory, f'{p_detector}-{i}.txt')
            file_path2 = os.path.join(directory, f'n_{p_detector}.txt')
            with open(file_path2, 'a') as f:
                f.write(f'{i}\t{np.sum(mask)}\n')
            with open(file_path, 'w') as f:
                for i in range(len(filter_Dl)):
                    line = f"{filter_Dl[i]}\t{filter_m1[i]}\t{filter_m2[i]}\t{filter_Tc[i]}\t" \
                        f"{filter_theta_s[i]}\t{filter_d_alpha_dCS_bar_square[i]}\n"
                    f.write(line)

def validity_events1(p_detectors, i):
        """筛选满足有效性和 SNR 条件的事件。"""
        for p_detector in p_detectors:
            try:
                with open(f'/home/{name}/code/dCS/data/{date}/results/results for {p_detector}-{i}.txt', 'r') as f:
                    lines = f.readlines()
            except:
                print(f"No result for {p_detector} iteration {i}")
                return
            # 将列表转换为NumPy数组
            Dl = np.array([float(line.strip().split(',')[0]) for line in lines])
            m1 = np.array([float(line.strip().split(',')[1]) for line in lines])
            m2 = np.array([float(line.strip().split(',')[2]) for line in lines])
            Tc = np.array([float(line.strip().split(',')[3]) for line in lines])
            theta_s = np.array([float(line.strip().split(',')[4]) for line in lines])
            d_alpha_dCS_bar_square = np.array([float(line.strip().split(',')[5])**0.25*1e-3 for line in lines])
            # alpha_dCS_bar_square_max = np.array([((float(line.strip().split(',')[2])*MSUN_SI)**4 * G**4 / c**8 / 32 / pi)**0.25*1e-3 for line in lines])
            alpha_dCS_bar_square_max = np.array([float(line.strip().split(',')[6]) for line in lines])
            SNR = np.array([float(line.strip().split(',')[7]) for line in lines])
            mask = (d_alpha_dCS_bar_square < alpha_dCS_bar_square_max) & (SNR > 8)
            filter_Dl = Dl[mask]
            filter_m1 = m1[mask]
            filter_m2 = m2[mask]
            filter_Tc = Tc[mask]
            filter_theta_s = theta_s[mask]
            filter_d_alpha_dCS_bar_square = d_alpha_dCS_bar_square[mask]
            directory = f'/home/{name}/code/dCS/data/{date}/validity2/'
            os.makedirs(directory, exist_ok=True)
            file_path = os.path.join(directory, f'{p_detector}-{i}.txt')
            file_path2 = os.path.join(directory, f'm_{p_detector}.txt')
            with open(file_path2, 'a') as f:
                f.write(f'{i}\t{np.sum(mask)}\n')
            with open(file_path, 'w') as f:
                for i in range(len(filter_Dl)):
                    line = f"{filter_Dl[i]}\t{filter_m1[i]}\t{filter_m2[i]}\t{filter_Tc[i]}\t" \
                        f"{filter_theta_s[i]}\t{filter_d_alpha_dCS_bar_square[i]}\n"
                    f.write(line)

def events_TianQin_CE(i):
    """处理第 i 批 TianQin+CE 事件参数。"""
    detectors1 = ['CE']
    detectors2 = ['TianQin']
    try:
        file_path = f'/data/{name}/BP2P/parameters-validity/TianQin+CE/parameters-{i}.txt'
        param_space = generate_parameter_space1(detectors1, detectors2, file_path)
    except:
        print(f"Error reading file for TianQin+CE iteration {i}")
        return
    main1(param_space, i)
    validity_events(['TianQin+CE'], i)
    validity_events1(['TianQin+CE'], i)

def events_TianQin_ET(i):
    """处理第 i 批 TianQin+ET 事件参数。"""
    detectors1 = ['ET']
    detectors2 = ['TianQin']
    try:
        file_path = f'/data/{name}/BP2P/parameters-validity/TianQin+ET/parameters-{i}.txt'
        param_space = generate_parameter_space1(detectors1, detectors2, file_path)
    except:
        print(f"Error reading file for TianQin+ET iteration {i}")
        return
    main1(param_space, i)
    validity_events(['TianQin+ET'], i)
    validity_events1(['TianQin+ET'], i)

def events_BBO(i):
    """处理第 i 批 BBO 事件参数。"""
    detector = ['BBO']
    try:
        file_path = f'/data/{name}/BP2P/parameters-validity/BBO/parameters-{i}.txt'
        param_space = generate_parameter_space2(detector, file_path)
    except:
        print(f"Error reading file for BBO iteration {i}")
        return
    main2(param_space, i)
    validity_events(['BBO'], i)
    validity_events1(['BBO'], i)

def events_DECIGO(i):
    """处理第 i 批 DECIGO 事件参数。"""
    detector = ['DECIGO']
    try:
        file_path = f'/data/{name}/BP2P/parameters-validity/DECIGO/parameters-{i}.txt'
        param_space = generate_parameter_space2(detector, file_path)
    except:
        print(f"Error reading file for DECIGO iteration {i}")
        return
    main2(param_space, i)
    validity_events(['DECIGO'], i)
    validity_events1(['DECIGO'], i)

def events_CE(i):
    """处理第 i 批 CE 事件参数。"""
    detector = ['CE']
    try:
        file_path = f'/data/{name}/BP2P/parameters-validity/CE/parameters-{i}.txt'
        param_space = generate_parameter_space2(detector, file_path)
    except:
        print(f"Error reading file for CE iteration {i}")
        return
    main2(param_space, i)
    validity_events(['CE'], i)
    validity_events1(['CE'], i)

def events_ET(i):
    """处理第 i 批 ET 事件参数。"""
    detector = ['ET']
    try:
        file_path = f'/data/{name}/BP2P/parameters-validity/ET/parameters-{i}.txt'
        param_space = generate_parameter_space2(detector, file_path)
    except:
        print(f"Error reading file for ET iteration {i}")
        return
    main2(param_space, i)
    validity_events(['ET'], i)
    validity_events1(['ET'], i)

if __name__ == '__main__':
    date = 20251224
    name = 'chexy'  # tq5


    with open('log-events.csv', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Iteration', 'Function', 'Execution_Time(s)'])

        for n in range(100):
            functions = [
                events_CE, 
                events_ET, 
                events_TianQin_CE, 
                events_TianQin_ET,
                events_DECIGO, 
                events_BBO
            ]

            for func in functions:
                i = n
                print(f'Starting {func.__name__} for iteration {i}')
                start_time = time.time()
                func(i)
                end_time = time.time()
                print(f'Completed {func.__name__} for iteration {i}','time:', end_time - start_time)
                # 计算执行时间并写入文件
                execution_time = end_time - start_time
                writer.writerow([i, func.__name__, execution_time])
    # plot_ev()
