"""质量比 q 对 dCS 约束影响的单变量计算。

本模块固定探测器、总质量、光度距离、观测时间、自旋和倾角，只改变
质量比 q，并返回 dCS 参数约束、SNR 和频率截止信息。
"""

import numpy as np
from dcs.fim import G_c_3, covariance, get_detecror_parameters, q_to_eta, snr
from lal import MSUN_SI, PC_SI
from scipy.constants import G, c, pi, year


def validity_scale_km(m2_msun):
    """Small-coupling length scale of the lighter component, in kilometres."""
    return ((m2_msun * MSUN_SI)**4 * G**4 / c**8 / (32*pi))**0.25 * 1e-3


def build_parameters(
    p_detector,
    p_M_total_n,
    p_D_l_mpc,
    p_T_tc_year,
    p_q,
    p_x1=0.40,
    p_x2=0.20,
    p_iota=0.9,
):
    """构造固定源参数、指定质量比下的波形和差分参数。

    Args:
        p_detector: 探测器名称。
        p_M_total_n: 总质量，单位太阳质量。
        p_D_l_mpc: 光度距离，单位 Mpc。
        p_T_tc_year: 观测时间，单位年。
        p_q: 质量比 q=m1/m2。
        p_x1: 第一颗致密天体的无量纲自旋。
        p_x2: 第二颗致密天体的无量纲自旋。
        p_iota: 倾角，单位弧度。

    Returns:
        `(parameters, parameters_snr, delta_para, metadata, Sn_dector)`。
    """
    if not all(np.isfinite(v) for v in (p_M_total_n, p_D_l_mpc, p_T_tc_year, p_q, p_x1, p_x2, p_iota)):
        raise ValueError('Source parameters must be finite')
    if min(p_M_total_n, p_D_l_mpc, p_T_tc_year) <= 0 or p_q < 1:
        raise ValueError('Positive mass/distance/time and q >= 1 required')
    if max(abs(p_x1), abs(p_x2)) > 1 or not 0 <= p_iota <= np.pi:
        raise ValueError('Invalid aligned spin or inclination')
    f_D_low, f_D_high, _, _, _, _, Sn_dector, delta_M_n, delta_eta_n, delta_x1_n, delta_x2_n = get_detecror_parameters(p_detector)

    eta = q_to_eta(p_q)
    T_tc = p_T_tc_year * year
    M_total = p_M_total_n * MSUN_SI
    D_l = p_D_l_mpc * 1e6 * PC_SI

    f_PN_low = pow(5 * pow(eta, 3 / 5) * G_c_3 * M_total / T_tc, 3 / 8) / (8 * pi * pow(eta, 3 / 5) * G_c_3 * M_total)
    f_ISCO = 1 / (pow(6, 1.5) * pi * G_c_3 * M_total)
    f_high = min(f_ISCO, f_D_high)
    f_low = max(f_PN_low, f_D_low)
    deltaF = f_low / 50

    parameters = {
        'eta': eta,
        'x1': p_x1,
        'x2': p_x2,
        'alpha_dCS_bar_square': 0,
        'T_tc': T_tc,
        'iota': p_iota,
        'a': 4,
        'b': -1,
        'M': M_total,
        'D_l': D_l,
        'f_min': f_low,
        'f_max': f_high,
        'deltaF': deltaF,
    }
    parameters_snr = parameters.copy()
    parameters_snr['f_max'] = f_D_high

    delta_para = {
        'delta_M': M_total * delta_M_n,
        'delta_eta': eta * delta_eta_n,
        'delta_x1': (p_x1 if abs(p_x1) >= 1e-6 else 1.0) * delta_x1_n,
        'delta_x2': (p_x2 if abs(p_x2) >= 1e-6 else 1.0) * delta_x2_n,
    }
    metadata = {
        'eta': eta,
        'f_low': f_low,
        'f_high': f_high,
        'f_detector_low': f_D_low,
        'f_detector_high': f_D_high,
        'f_isco': f_ISCO,
        'f_pn_low': f_PN_low,
        'deltaF': deltaF,
    }
    return parameters, parameters_snr, delta_para, metadata, Sn_dector


def get_result(
    p_detector,
    p_M_total_n,
    p_D_l_mpc,
    p_T_tc_year,
    p_q,
    p_x1=0.40,
    p_x2=0.20,
    p_iota=0.9,
):
    """计算指定质量比下的 dCS 约束精度。

    Returns:
        包含输入参数、频率截止、SNR 和 dCS 约束精度的字典。
    """
    parameters, parameters_snr, delta_para, metadata, Sn_dector = build_parameters(
        p_detector,
        p_M_total_n,
        p_D_l_mpc,
        p_T_tc_year,
        p_q,
        p_x1,
        p_x2,
        p_iota,
    )
    if parameters['f_min'] >= parameters['f_max']:
        return None
    cov = covariance(parameters, delta_para, Sn_dector)
    if cov is None:
        return None
    variance = float(cov[7, 7].real)
    if not np.isfinite(variance) or variance <= 0:
        return None
    sigma_alpha_dCS_bar_square = np.sqrt(variance)
    d_sqrt_alpha_dCS_km = sigma_alpha_dCS_bar_square**0.25 * 1e-3
    m2 = p_M_total_n / (p_q + 1)
    alpha_dCS_bar_square_max_km = validity_scale_km(m2)
    SNR = snr(parameters_snr, Sn_dector)

    return {
        'detector': p_detector,
        'q': p_q,
        'eta': metadata['eta'],
        'M_total_Msun': p_M_total_n,
        'D_l_Mpc': p_D_l_mpc,
        'T_tc_year': p_T_tc_year,
        'x1': p_x1,
        'x2': p_x2,
        'iota': p_iota,
        'f_low_Hz': metadata['f_low'],
        'f_high_Hz': metadata['f_high'],
        'f_pn_low_Hz': metadata['f_pn_low'],
        'f_isco_Hz': metadata['f_isco'],
        'deltaF_Hz': metadata['deltaF'],
        'sigma_alpha_dCS_bar_square': sigma_alpha_dCS_bar_square,
        'd_sqrt_alpha_dCS_km': d_sqrt_alpha_dCS_km,
        'alpha_dCS_bar_square_max_km': alpha_dCS_bar_square_max_km,
        'SNR': SNR,
    }


def task(params):
    """供多进程调用的 q 敏感性任务包装函数。"""
    return get_result(*params)
