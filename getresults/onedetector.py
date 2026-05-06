"""单探测器参数网格计算封装。

本文件把探测器名称、质量网格位置和质量比转换为完整波形参数，
并调用核心 Fisher/SNR 函数生成单探测器约束结果。
"""

import numpy as np
from fim import G_c_3, q_to_eta, covariance, snr, get_detecror_parameters
from lal import MSUN_SI, PC_SI
from scipy.constants import c, G, year, pi

def get_result(p_detector, p_l, p_q):
    """计算单个探测器在给定质量网格和质量比下的约束。

    Args:
        p_detector: 探测器名称。
        p_l: 质量网格索引，范围按当前脚本约定为 0 到 1000。
        p_q: 质量比 q=m1/m2。

    Returns:
        协方差矩阵、总质量、光度距离、dCS 有效性上界和 SNR。
    """
    f_D_low, f_D_high, M_low_n, M_high_n, _, _, Sn_dector, delta_M_n, delta_eta_n, delta_x1_n, delta_x2_n = get_detecror_parameters(p_detector)
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
    f_high = min(f_ISCO, f_D_high)
    f_low = max(f_PN_low, f_D_low)
    deltaF = f_low / 50
    # deltaF = 1/ T_tc

    parameters = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low, 'f_max': f_high, 'deltaF': deltaF}

    parameters1 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low, 'f_max': f_D_high, 'deltaF': deltaF}

    # 差分间隔
    delta_M = M_total * delta_M_n
    delta_eta = eta * delta_eta_n
    delta_x1 = x1 * delta_x1_n
    delta_x2 = x2 * delta_x2_n
    # delta_M = M_total * 1e-8
    # delta_eta = eta * 1e-9
    # delta_x1 = x1 * 1e-5
    # delta_x2 = x2 * 1e-5
    delta_para = {'delta_M':delta_M, 'delta_eta':delta_eta, 'delta_x1':delta_x1,'delta_x2':delta_x2}

    cov = covariance(parameters, delta_para, Sn_dector)
    SNR = snr(parameters1, Sn_dector)
    alpha_dCS_bar_square_max = pow(pow(M_total / (p_q+1), 4) * G**4 / c**8 / 32 / pi, 0.25) * 1e-3  # km

    return cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR

def task(p_detector, p_l, p_q):
    """供多进程调用的单探测器任务包装函数。"""
    cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = get_result(p_detector, p_l, p_q)
    array_delta = []
    for i in range(8):
        array_delta.append(np.sqrt(cov[i, i]).real)
    return p_detector, p_l, p_q, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR
