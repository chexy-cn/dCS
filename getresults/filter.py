"""事件样本的 dCS 约束计算封装。

该模块读取单个事件参数后，构造单探测器或双探测器联合观测的波形
参数，并返回 Fisher 协方差、dCS 约束尺度和 SNR。
"""

import numpy as np
from fim import G_c_3, q_to_eta, covariance_two, snr_two, get_detecror_parameters, covariance, snr
from lal import MSUN_SI, PC_SI
from scipy.constants import c, G, year, pi

def get_result_two(p_detector1, p_detector2, p_parameters):
    """计算单个事件在双探测器联合观测下的 dCS 约束。

    Args:
        p_detector1: 第一个探测器名称。
        p_detector2: 第二个探测器名称。
        p_parameters: 事件参数 `[Dl, m1, m2, Tc, theta_s]`。

    Returns:
        协方差矩阵、总质量、光度距离、dCS 有效性上界和联合 SNR。
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

    # 差分间隔
    delta_M = M_total * delta_M_n
    delta_eta = delta_eta_n
    delta_x1 = delta_x1_n
    delta_x2 = delta_x2_n
    delta_para = {'delta_M':delta_M, 'delta_eta':delta_eta, 'delta_x1':delta_x1,'delta_x2':delta_x2}

    try:
        cov = covariance_two(parameters1, parameters2, delta_para, Sn_dector1, Sn_dector2)
        SNR = (snr_two(parameters1, Sn_dector1, Sn_dector2)**2 + snr_two(parameters2, Sn_dector1, Sn_dector2)**2)**0.5
        alpha_dCS_bar_square_max = ((m2*MSUN_SI)**4 * G**4 / c**8 / 32 / pi)**0.25*1e-3  # km
        return cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR
    except:
        print(p_detector1, p_detector2, Dl, m1, m2, '参数错误！')

def get_result_one(p_detector, p_parameters):
    """计算单个事件在单探测器观测下的 dCS 约束。

    Args:
        p_detector: 探测器名称。
        p_parameters: 事件参数 `[Dl, m1, m2, Tc, theta_s]`。

    Returns:
        协方差矩阵、总质量、光度距离、dCS 有效性上界和 SNR。
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
    f_high = min(f_ISCO, f_D_high)
    f_low = max(f_PN_low, f_D_low)
    deltaF = f_low / 50
    # deltaF = 1/ T_tc

    parameters = {'eta': eta, 'x1': x1, 'x2': x2, 'alpha_dCS_bar_square': alpha_dCS_bar_square, 'T_tc': T_tc, 'iota': iota, 'a': a, 'b': b, 'M': M_total, 'D_l': D_l, 'f_min': f_low, 'f_max': f_high, 'deltaF': deltaF}

    # 差分间隔
    delta_M = M_total * delta_M_n
    delta_eta = delta_eta_n
    delta_x1 = delta_x1_n
    delta_x2 = delta_x2_n
    delta_para = {'delta_M':delta_M, 'delta_eta':delta_eta, 'delta_x1':delta_x1,'delta_x2':delta_x2}
    
    try:
        cov = covariance(parameters, delta_para, Sn_dector)
        SNR = snr(parameters, Sn_dector)
        alpha_dCS_bar_square_max = ((m2*MSUN_SI)**4 * G**4 / c**8 / 32 / pi)**0.25*1e-3  # km
        return cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR
    except:
        print(p_detector, Dl, m1, m2, '参数错误！')
    
def task_two(p_detector1, p_detector2, p_parameters):
    """供多进程调用的双探测器事件任务包装函数。"""
    try:
        cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = get_result_two(p_detector1, p_detector2, p_parameters)
        array_delta = []
        for i in range(8):
            array_delta.append(np.sqrt(cov[i, i]).real)
        return p_detector1, p_detector2, p_parameters, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR
    except:
        print(p_detector1, p_detector2, p_parameters, '参数错误！')

def task_one(p_detector, p_parameters):
    """供多进程调用的单探测器事件任务包装函数。"""
    try:
        cov, M_total_n, D_l, alpha_dCS_bar_square_max, SNR = get_result_one(p_detector, p_parameters)
        array_delta = []
        for i in range(8):
            array_delta.append(np.sqrt(cov[i, i]).real)
        return p_detector, p_parameters, *array_delta, M_total_n, D_l, alpha_dCS_bar_square_max, SNR
    except:
        print(p_detector, p_parameters, '参数错误！')
