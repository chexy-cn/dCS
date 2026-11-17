"""对齐/反对齐自旋对 dCS 约束的二维敏感性计算。"""

from getresults import q_sensitivity
from dcs.fim import snr


def get_result(p_detector, p_M_total_n, p_D_l_mpc, p_T_tc_year, p_q, p_chi_s, p_chi_a, p_iota=0.9):
    """计算指定对称和反对称自旋下的约束，跳过非物理自旋。"""
    x1 = p_chi_s + p_chi_a
    x2 = p_chi_s - p_chi_a
    # Permit round-off at the stated physical boundary and clamp it exactly.
    if abs(x1) > 0.9 + 1e-12 or abs(x2) > 0.9 + 1e-12:
        return None
    x1 = min(0.9, max(-0.9, x1))
    x2 = min(0.9, max(-0.9, x2))

    nonspinning = abs(x1) < 1e-12 and abs(x2) < 1e-12
    if nonspinning:
        # Keep the coordinate in the output even though no finite constraint exists.
        parameters, parameters_snr, _, meta, noise = q_sensitivity.build_parameters(
            p_detector, p_M_total_n, p_D_l_mpc, p_T_tc_year, p_q, x1, x2, p_iota)
        result = dict(detector=p_detector, q=p_q, eta=meta['eta'],
                      M_total_Msun=p_M_total_n, D_l_Mpc=p_D_l_mpc,
                      T_tc_year=p_T_tc_year, x1=x1, x2=x2, iota=p_iota,
                      f_low_Hz=meta['f_low'], f_high_Hz=meta['f_high'],
                      f_pn_low_Hz=meta['f_pn_low'], f_isco_Hz=meta['f_isco'],
                      deltaF_Hz=meta['deltaF'], sigma_alpha_dCS_bar_square=float('nan'),
                      d_sqrt_alpha_dCS_km=float('nan'),
                      alpha_dCS_bar_square_max_km=q_sensitivity.validity_scale_km(p_M_total_n/(1+p_q)),
                      SNR=float(snr(parameters_snr, noise)))
    else:
        result = q_sensitivity.get_result(
            p_detector, p_M_total_n, p_D_l_mpc, p_T_tc_year, p_q, x1, x2, p_iota
        )
    if result is None:
        return None
    if abs(x1) < 1e-12 and abs(x2) < 1e-12:
        # The dCS waveform derivative vanishes for a non-spinning binary.
        # A finite inverse-Fisher value here is numerical round-off, not a constraint.
        result['sigma_alpha_dCS_bar_square'] = float('nan')
        result['d_sqrt_alpha_dCS_km'] = float('nan')
    result['chi_s'] = p_chi_s
    result['chi_a'] = p_chi_a
    return result


def task(params):
    """供多进程调用的自旋敏感性任务包装函数。"""
    return get_result(*params)
