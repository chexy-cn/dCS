"""双探测器联合观测参数网格计算封装。

本文件把两个探测器的频段、噪声谱和源参数组合起来，计算联合
Fisher 矩阵、协方差和 SNR。
"""

import numpy as np
from fim import G_c_3, q_to_eta, covariance_two, snr_two, get_detecror_parameters
from lal import MSUN_SI, PC_SI
from scipy.constants import c, G, year, pi

def get_result_simplify(p_detector1, p_detector2, p_l, p_q):
    """在固定距离缩放假设下计算双探测器联合约束。

    Args:
        p_detector1: 第一个探测器名称。
        p_detector2: 第二个探测器名称。
        p_l: 质量网格索引。
        p_q: 质量比 q=m1/m2。

    Returns:
        协方差矩阵、总质量、参考光度距离、dCS 有效性上界和联合 SNR。
    """
    f_D_low1, f_D_high1, M_low_n, M_high_n, _, _, Sn_dector1, delta_M_n, delta_eta_n, delta_x1_n, delta_x2_n = get_detecror_parameters(p_detector1)
    f_D_low2, f_D_high2, _, _, _, _, Sn_dector2, _, _, _, _ = get_detecror_parameters(p_detector2)
    
    M_total_n = pow(10, np.log10(M_low_n) + p_l * (np.log10(M_high_n) - np.log10(M_low_n)) / 1000)

    # 参数
    a = 4
    b = -1
    eta = q_to_eta(p_q)
    x1 = 0.40
    x2 = 0.20
    alpha_dCS_bar_square = 0
    iota = 0.9
    T_tc = 0.25 * year

    D_l = 1e6 * PC_SI
    M_total = M_total_n * MSUN_SI

    f_PN_low = pow(5 * pow(eta, 3 / 5) * G_c_3 * M_total / T_tc, 3 / 8) / (8 * pi * pow(eta, 3 / 5) * G_c_3 * M_total)
    f_ISCO = 1 / (pow(6, 1.5) * pi * G_c_3 * M_total)
    f_high1 = min(f_ISCO, f_D_high1)
    f_low1 = max(f_PN_low, f_D_low1)
    f_high2 = min(f_ISCO, f_D_high2)
    f_low2 = max(f_PN_low, f_D_low2)
    deltaF1 = f_low1 / 50
    deltaF2 = f_low2 / 50
    # deltaF = 1/ T_tc
    # print(f_low1, f_high1, f_low2, f_high2, f_PN_low, f_ISCO)
    
    parameters1 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low1, 'f_max': f_high1, 'deltaF': deltaF1}
    parameters11 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low1, 'f_max': f_D_high1, 'deltaF': deltaF1}

    parameters2 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low2, 'f_max': f_high2, 'deltaF': deltaF2}
    parameters21 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low2, 'f_max': f_D_high2, 'deltaF': deltaF2}

    # 差分间隔
    delta_M = M_total * delta_M_n
    delta_eta = delta_eta_n
    delta_x1 = delta_x1_n
    delta_x2 = delta_x2_n
    delta_para = {'delta_M':delta_M, 'delta_eta':delta_eta, 'delta_x1':delta_x1,'delta_x2':delta_x2}

    cov = covariance_two(parameters1, parameters2, delta_para, Sn_dector1, Sn_dector2)
    SNR = (snr_two(parameters11, Sn_dector1, Sn_dector2)**2 + snr_two(parameters21, Sn_dector1, Sn_dector2)**2)**0.5
    alpha_dCS_bar_square_max = pow(pow(M_total / (p_q+1), 4) * G**4 / c**8 / 32 / pi, 0.25) * 1e-3  # km

    return cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR

def get_result(p_detector1, p_detector2, p_l1, p_l2, p_q):
    """在质量和距离二维网格上计算双探测器联合约束。

    Args:
        p_detector1: 第一个探测器名称。
        p_detector2: 第二个探测器名称。
        p_l1: 质量网格索引。
        p_l2: 距离网格索引。
        p_q: 质量比 q=m1/m2。

    Returns:
        协方差矩阵、总质量、距离网格值、dCS 有效性上界和联合 SNR。
    """
    f_D_low1, f_D_high1, M_low_n, M_high_n, D_l_low_n, D_l_len_n, Sn_dector1, delta_M_n, delta_eta_n, delta_x1_n, delta_x2_n = get_detecror_parameters(p_detector1)
    f_D_low2, f_D_high2, _, _, _, _, Sn_dector2, _, _, _, _ = get_detecror_parameters(p_detector2)
    
    M_total_n = pow(10, np.log10(M_low_n) + p_l1 * (np.log10(M_high_n) - np.log10(M_low_n)) / 100)
    # M_total_n = pow(10, np.log10(M_high_n) - (np.log10(M_high_n) - np.log10(M_low_n)) * p_l1 / 100)
    D_l_n = p_l2 / 100 * D_l_len_n

    # 参数
    a = 4
    b = -1
    eta = q_to_eta(p_q)
    x1 = 0.40
    x2 = 0.20
    alpha_dCS_bar_square = 0
    iota = 0.9
    T_tc = 0.25 * year

    # D_l = 1e6 * PC_SI
    D_l = 10**(D_l_n + D_l_low_n) * PC_SI
    M_total = M_total_n * MSUN_SI

    f_PN_low = pow(5 * pow(eta, 3 / 5) * G_c_3 * M_total / T_tc, 3 / 8) / (8 * pi * pow(eta, 3 / 5) * G_c_3 * M_total)
    f_ISCO = 1 / (pow(6, 1.5) * pi * G_c_3 * M_total)
    f_high1 = min(f_ISCO, f_D_high1)
    f_low1 = max(f_PN_low, f_D_low1)
    f_high2 = min(f_ISCO, f_D_high2)
    f_low2 = max(f_PN_low, f_D_low2)
    deltaF1 = f_low1 / 50
    deltaF2 = f_low2 / 50
    # deltaF = 1/ T_tc
    # print(f_low1, f_high1, f_low2, f_high2, f_PN_low, f_ISCO)
    
    parameters1 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low1, 'f_max': f_high1, 'deltaF': deltaF1}

    parameters11 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low1, 'f_max': f_D_high1, 'deltaF': deltaF1}

    parameters2 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low2, 'f_max': f_high2, 'deltaF': deltaF2}

    parameters21 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low2, 'f_max': f_D_high2, 'deltaF': deltaF2}

    # 差分间隔
    delta_M = M_total * delta_M_n
    delta_eta = eta * delta_eta_n
    delta_x1 = x1 * delta_x1_n
    delta_x2 = x2 * delta_x2_n
    delta_para = {'delta_M':delta_M, 'delta_eta':delta_eta, 'delta_x1':delta_x1,'delta_x2':delta_x2}

    cov = covariance_two(parameters1, parameters2, delta_para, Sn_dector1, Sn_dector2)
    SNR = (snr_two(parameters11, Sn_dector1, Sn_dector2)**2 + snr_two(parameters21, Sn_dector1, Sn_dector2)**2)**0.5
    alpha_dCS_bar_square_max = pow(pow(M_total / (p_q+1), 4) * G**4 / c**8 / 32 / pi, 0.25) * 1e-3  # km

    return cov, M_total_n, D_l_n, alpha_dCS_bar_square_max, SNR

# def task(p_detector1, p_detector2, p_l1, p_l2, p_q):
#     cov, M_total_n, D_l_n, alpha_dCS_bar_square_max, SNR = get_result(p_detector1, p_detector2, p_l1, p_l2, p_q)
#     array_delta = []
#     for i in range(8):
#         array_delta.append(np.sqrt(cov[i, i]).real)
#     return p_detector1, p_detector2, p_l1, p_q, *array_delta, M_total_n, D_l_n, alpha_dCS_bar_square_max, SNR

def task(p_detector1, p_detector2, p_l, p_q):
    """供多进程调用的双探测器任务包装函数。"""
    cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = get_result_simplify(p_detector1, p_detector2, p_l, p_q)
    array_delta = []
    for i in range(8):
        array_delta.append(np.sqrt(cov[i, i]).real)
    return p_detector1, p_detector2, p_l, p_q, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR
