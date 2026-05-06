"""事件样本的信噪比计算封装。

该模块只计算 SNR，不求 Fisher 协方差，适合先对大量事件做可探测性
筛选。
"""

import numpy as np
from fim import G_c_3, q_to_eta, covariance_two, snr_two, get_detecror_parameters, covariance, snr
from lal import MSUN_SI, PC_SI
from scipy.constants import c, G, year, pi

def get_result_two(p_detector1, p_detector2, p_parameters):
    """计算单个事件在双探测器联合观测下的 SNR。

    Args:
        p_detector1: 第一个探测器名称。
        p_detector2: 第二个探测器名称。
        p_parameters: 事件参数 `[Dl, m1, m2, Tc, theta_s]`。

    Returns:
        联合 SNR；异常时返回大数占位值。
    """
    Dl, m1, m2, Tc, theta_s = p_parameters[0], p_parameters[1], p_parameters[2], p_parameters[3], p_parameters[4]
    f_D_low1, f_D_high1, _, _, _, _, Sn_dector1, delta_M_n, delta_eta_n, delta_x1_n, delta_x2_n = get_detecror_parameters(p_detector1)
    f_D_low2, f_D_high2, _, _, _, _, Sn_dector2, _, _, _, _ = get_detecror_parameters(p_detector2)
    
    M_total_n = m1 + m2
    q = m1 / m2

    # 参数
    a = 4
    b = -1
    eta = q_to_eta(q)
    x1 = 0.40
    x2 = 0.20
    alpha_dCS_bar_square = 0
    iota = theta_s
    T_tc = Tc

    M_total = M_total_n * MSUN_SI
    D_l = Dl * 1e6 * PC_SI

    f_PN_low = pow(5 * pow(eta, 3 / 5) * G_c_3 * M_total / T_tc, 3 / 8) / (8 * pi * pow(eta, 3 / 5) * G_c_3 * M_total)
    f_ISCO = 1 / (pow(6, 1.5) * pi * G_c_3 * M_total)
    f_high1 = min(f_ISCO, f_D_high1)
    f_low1 = max(f_PN_low, f_D_low1)
    f_high2 = min(f_ISCO, f_D_high2)
    f_low2 = max(f_PN_low, f_D_low2)
    deltaF1 = f_low1 / 50
    deltaF2 = f_low2 / 50
    # deltaF = 1/ T_tc
    
    parameters1 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low1, 'f_max': f_high1, 'deltaF': deltaF1}
    parameters2 = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low2, 'f_max': f_high2, 'deltaF': deltaF2}

    try:
        SNR = (snr_two(parameters1, Sn_dector1, Sn_dector2)**2 + snr_two(parameters2, Sn_dector1, Sn_dector2)**2)**0.5
        return SNR
    except:
        print(p_parameters,'异常')
        return 9999999999

def get_result_one(p_detector, p_parameters):
    """计算单个事件在单探测器观测下的 SNR。

    Args:
        p_detector: 探测器名称。
        p_parameters: 事件参数 `[Dl, m1, m2, Tc, theta_s]`。

    Returns:
        单探测器 SNR；异常时返回大数占位值。
    """
    Dl, m1, m2, Tc, theta_s = p_parameters[0], p_parameters[1], p_parameters[2], p_parameters[3], p_parameters[4]
    f_D_low, f_D_high, M_low_n, M_high_n, _, _, Sn_dector, delta_M_n, delta_eta_n, delta_x1_n, delta_x2_n = get_detecror_parameters(p_detector)

    M_total_n = m1 + m2
    q = m1 / m2

    # 参数
    a = 4
    b = -1
    eta = q_to_eta(q)
    x1 = 0.40
    x2 = 0.20
    alpha_dCS_bar_square = 0
    iota = theta_s
    T_tc = Tc

    M_total = M_total_n * MSUN_SI
    D_l = Dl * 1e6 * PC_SI

    f_PN_low = pow(5 * pow(eta, 3 / 5) * G_c_3 * M_total / T_tc, 3 / 8) / (8 * pi * pow(eta, 3 / 5) * G_c_3 * M_total)
    f_ISCO = 1 / (pow(6, 1.5) * pi * G_c_3 * M_total)
    # f_high = min(f_ISCO, f_D_high)
    f_high = f_D_high
    f_low = max(f_PN_low, f_D_low)
    deltaF = f_low / 50
    # deltaF = 1/ T_tc

    parameters = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low, 'f_max': f_high, 'deltaF': deltaF}

    try:
        SNR = snr(parameters, Sn_dector)
        return SNR
    except:
        print(p_parameters,'异常')
        return 9999999999
    
def task_two(p_detector1, p_detector2, p_parameters):
    """供多进程调用的双探测器 SNR 任务包装函数。"""
    SNR = get_result_two(p_detector1, p_detector2, p_parameters)
    return p_detector1, p_detector2, p_parameters, SNR

def task_one(p_detector, p_parameters):
    """供多进程调用的单探测器 SNR 任务包装函数。"""
    SNR = get_result_one(p_detector, p_parameters)
    return p_detector, p_parameters, SNR
