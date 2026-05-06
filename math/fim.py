"""dCS 引力波约束的核心计算模块。

本文件集中实现质量参数变换、dCS 修正项、LALSuite 频域波形调用、
偏导数、SNR、Fisher 矩阵、协方差矩阵和探测器参数配置。
"""

import numpy as np
from Sn import Sn
import lalsimulation as lalsim
from scipy.constants import c, G, year, pi
from numpy import linalg
from swsh import swsh_lm

# 常数
G_c_3 = G / pow(c, 3)
c_8_G_4 = pow(c, 8) / pow(G, 4)

# 参数函数
def u(p_f, p_eta, p_M):
    """计算后牛顿展开中使用的无量纲速度参数 u(f)。"""
    return pow(pi * pow(p_eta, 3 / 5) * G_c_3 * p_M * p_f, 1 / 3)

def x_s(p_x1, p_x2):
    """计算两个自旋参数的对称组合。"""
    return (p_x1 + p_x2) / 2

def x_a(p_x1, p_x2):
    """计算两个自旋参数的反对称组合。"""
    return (p_x1 - p_x2) / 2

def zeta_dCS(p_M, p_alpha_dCS_bar_square):
    """由总质量和 dCS 参数平方计算无量纲耦合 zeta_dCS。"""
    return 16 * pi * c_8_G_4 * p_alpha_dCS_bar_square / pow(p_M, 4)

def eta_M_to_m1_m2(p_M, p_eta):
    """由总质量和对称质量比反解两个分量质量。"""
    q = -1 + (np.sqrt(1 - 4 * p_eta) + 1) / (2 * p_eta)
    return p_M * q / (1 + q), p_M / (1 + q)

def q_to_eta(p_q):
    """将质量比 q=m1/m2 转换为对称质量比 eta。"""
    return p_q / pow(1 + p_q, 2)

def eta_to_q(p_eta):
    """将对称质量比 eta 转换为质量比 q。"""
    return -1 + (np.sqrt(1 - 4 * p_eta) + 1) / (2 * p_eta)

def delta_m(p_eta):
    """由对称质量比计算归一化质量差。"""
    q = -1 + (np.sqrt(1 - 4 * p_eta) + 1) / (2 * p_eta)
    return (q - 1) / (q + 1)

def alpha_dCS(p_eta, p_zeta_dCS, p_delta_m, p_x_s, p_x_a):
    """计算 dCS 振幅修正系数 alpha。"""
    return p_zeta_dCS * 57713 / 344064 / pow(p_eta, 14 / 5) * (-2 * p_delta_m * p_x_a * p_x_s + (1 - 14976 *
                                                                                                 p_eta / 57713) * pow(
        p_x_a, 2) + (1 - 215876 * p_eta / 57713) * pow(p_x_s, 2))

def beta_dCS(p_eta, p_zeta_dCS, p_delta_m, p_x_s, p_x_a):
    """计算 dCS 相位修正系数 beta。"""
    return p_zeta_dCS * 481525 / 3670016 / pow(p_eta, 14 / 5) * (-2 * p_delta_m * p_x_a * p_x_s + (1 - 4992 *
                                                                                                   p_eta / 19261) * pow(
        p_x_a, 2) + (1 - 72052 * p_eta / 19261) * pow(p_x_s, 2))

def alpha_dCS_ctb(p_f, p_para):
    """计算波形振幅中的 dCS 乘性修正项。"""
    p_M = p_para['M']
    p_eta = p_para['eta']
    p_x1 = p_para['x1']
    p_x2 = p_para['x2']
    p_alpha_dCS_bar_square = p_para['alpha_dCS_bar_square']
    p_a = p_para['a']
    return (1 + alpha_dCS(p_eta, zeta_dCS(p_M, p_alpha_dCS_bar_square), delta_m(p_eta), x_s(p_x1, p_x2), x_a(p_x1, p_x2)) *
            pow(u(p_f, p_eta, p_M), p_a))

def beta_dCS_ctb(p_f, p_para, p_l, p_m):
    """计算指定模式的 dCS 相位指数修正项。"""
    p_b = p_para['b']
    p_M = p_para['M']
    p_eta = p_para['eta']
    p_x1 = p_para['x1']
    p_x2 = p_para['x2']
    p_alpha_dCS_bar_square = p_para['alpha_dCS_bar_square']
    return np.exp(complex(0, 1) * beta_dCS(p_eta, zeta_dCS(p_M, p_alpha_dCS_bar_square), delta_m(p_eta), x_s(p_x1, p_x2),
                                      x_a(p_x1, p_x2)) * pow(u(p_f, p_eta, p_M), p_b) * pow(2 / p_m, p_b / 3 -1))

def beta_dCS_ctb_ast(p_f, p_para, p_l, p_m):
    """计算 dCS 相位指数修正项的复共轭。"""
    return beta_dCS_ctb(p_f, p_para, p_l, p_m).conjugate()

def alpha_dCS_par(f, p_para):
    """计算波形对 alpha_dCS_bar_square 的振幅修正偏导部分。"""
    p_M = p_para['M']
    p_eta = p_para['eta']
    p_x1 = p_para['x1']
    p_x2 = p_para['x2']
    p_a = p_para['a']
    p_alpha_dCS_bar_square = p_para['alpha_dCS_bar_square']
    return (57713 / 344064 * 16 * pi  * c_8_G_4 / pow(p_M, 4) / pow(p_eta, 14 / 5) *
            (-2 * delta_m(p_eta) * x_a(p_x1, p_x2) * x_s(p_x1, p_x2) + (1 - 14976 * p_eta / 57713) *
            pow(x_a(p_x1, p_x2), 2) + (1 - 215876 * p_eta / 57713) * pow(x_s(p_x1, p_x2), 2)) *
            pow(u(f, p_eta, p_M), p_a))
def beta_dCS_par(f, p_para, p_l, p_m):
    """计算波形对 alpha_dCS_bar_square 的相位修正偏导部分。"""
    p_b = p_para['b']
    p_M = p_para['M']
    p_eta = p_para['eta']
    p_x1 = p_para['x1']
    p_x2 = p_para['x2']
    p_alpha_dCS_bar_square = p_para['alpha_dCS_bar_square']
    return (np.exp(complex(0, 1) *
            beta_dCS(p_eta, zeta_dCS(p_M, p_alpha_dCS_bar_square), delta_m(p_eta), x_s(p_x1, p_x2), x_a(p_x1, p_x2))
            * pow(u(f, p_eta, p_M), p_b) * pow(2 / p_m, p_b / 3 -1)) * complex(0, 1) * 481525 / 3670016 * 16 * pi * c_8_G_4 / pow(p_M, 4) /
            pow(p_eta, 14 / 5) * (-2 * delta_m(p_eta) * x_a(p_x1, p_x2) * x_s(p_x1, p_x2) + (1 - 4992 * p_eta / 19261) *
            pow(x_a(p_x1, p_x2), 2) + (1 - 72052 * p_eta / 19261) * pow(x_s(p_x1, p_x2), 2))
            * pow(u(f, p_eta, p_M), p_b) * pow(2 / p_m, p_b / 3 -1))

def beta_dCS_par_ast(f, para, p_l, p_m):
    """计算 beta_dCS_par 的复共轭。"""
    return beta_dCS_par(f, para, p_l, p_m).conjugate()

# 波形获取
def phenomXHM_lm_fd(p_parameters, ell, emm):
    """调用 LALSuite 生成指定 (ell, m) 模式的频域波形。

    Args:
        p_parameters: 包含质量、距离、频段和步长的参数字典。
        ell: 球谐模式 l。
        emm: 球谐模式 m。

    Returns:
        二元组 `(hp_array, f_array)`，分别为非零复波形和对应频率。
    """
    M, eta = p_parameters['M'], p_parameters['eta']
    m1_SI, m2_SI = eta_M_to_m1_m2(M, eta)
    distance = p_parameters['D_l']
    chi1L = p_parameters['x1']
    chi2L = p_parameters['x2']
    f_min = p_parameters['f_min']
    f_max = p_parameters['f_max']
    deltaF = p_parameters['deltaF']
    phi_Ref = 0.0
    fRef_In = 0.0
    lalParams = dict()

    htildelm = lalsim.SimIMRPhenomXHMGenerateFDOneMode(m1_SI, m2_SI, chi1L, chi2L, ell, emm, distance, f_min, f_max, deltaF, phi_Ref, fRef_In, lalParams)
    hp_array = np.asarray(htildelm.data.data, dtype=np.complex128)
    nonzero_mask = hp_array != 0
    return hp_array[nonzero_mask], deltaF * np.arange(htildelm.data.length)[nonzero_mask]

# 波形调用与定义
def get_waveform_lm(p_para):
    """生成包含 dCS 修正和自旋权球谐投影的总波形及分模波形。

    Args:
        p_para: 波形和 dCS 修正所需的参数字典。

    Returns:
        总波形、修正后分模、原始分模和频率数组。
    """
    hp22, f = phenomXHM_lm_fd(p_para, 2, 2)
    hp21, _ = phenomXHM_lm_fd(p_para, 2, 1)
    hp33, _ = phenomXHM_lm_fd(p_para, 3, 3)
    hp32, _ = phenomXHM_lm_fd(p_para, 3, 2)
    hp44, _ = phenomXHM_lm_fd(p_para, 4, 4)

    p_iota = p_para['iota']
    
    alpha_values = np.array([alpha_dCS_ctb(freq, p_para) for freq in f])
    beta22_values = np.array([beta_dCS_ctb(freq, p_para, 2, 2) for freq in f])
    beta21_values = np.array([beta_dCS_ctb(freq, p_para, 2, 1) for freq in f])
    beta33_values = np.array([beta_dCS_ctb(freq, p_para, 3, 3) for freq in f])
    beta32_values = np.array([beta_dCS_ctb(freq, p_para, 3, 2) for freq in f])
    beta44_values = np.array([beta_dCS_ctb(freq, p_para, 4, 4) for freq in f])
    
    hp_22 = hp22 * swsh_lm(2, 2, p_iota) * beta22_values * alpha_values
    hp_21 = hp21 * swsh_lm(2, 1, p_iota) * beta21_values * alpha_values
    hp_33 = hp33 * swsh_lm(3, 3, p_iota) * beta33_values * alpha_values
    hp_32 = hp32 * swsh_lm(3, 2, p_iota) * beta32_values * alpha_values
    hp_44 = hp44 * swsh_lm(4, 4, p_iota) * beta44_values * alpha_values

    hp_total = hp_22 + hp_21 + hp_33 + hp_32 + hp_44
    
    return hp_total, hp_22, hp_21, hp_33, hp_32, hp_44, hp22, hp21, hp33, hp32, hp44, f

# 偏导定义
def par(p_a, p_b, p_delta):
    """用前向差分计算 `(p_a - p_b) / p_delta`。"""
    return (p_a - p_b) / p_delta

def get_h_par(p_para, p_delta):
    """计算 Fisher 矩阵所需的波形参数偏导。

    参数顺序为 M、eta、D_l、t_c、phi_c、x1、x2、
    alpha_dCS_bar_square。

    Args:
        p_para: 基准波形参数。
        p_delta: M、eta、x1、x2 的差分步长字典。

    Returns:
        `(p_h_par, p_h_par_ast, p_f)`，其中前两项为偏导和复共轭列表。
    """
    p_delta_M = p_delta['delta_M']
    p_delta_eta = p_delta['delta_eta']
    p_delta_x1 = p_delta['delta_x1']
    p_delta_x2 = p_delta['delta_x2']


    para_M = p_para.copy()
    para_eta = p_para.copy()
    para_x1 = p_para.copy()
    para_x2 = p_para.copy()

    para_M['M'] += p_delta_M
    para_eta['eta'] += p_delta_eta
    para_x1['x1'] += p_delta_x1
    para_x2['x2'] += p_delta_x2

    # para_M0 = p_para.copy()
    # para_eta0 = p_para.copy()
    # para_x10 = p_para.copy()
    # para_x20 = p_para.copy()

    # para_M0['M'] -= p_delta_M
    # para_eta0['eta'] -= p_delta_eta
    # para_x10['x1'] -= p_delta_x1
    # para_x20['x2'] -= p_delta_x2

    p_D_l = p_para['D_l']
    p_iota = p_para['iota']
    
    p_h, _, _, _, _, _, hp22, hp21, hp33, hp32, hp44, p_f = get_waveform_lm(p_para)
    p_h_M, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_M)
    p_h_eta, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_eta)
    p_h_x1, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_x1)
    p_h_x2, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_x2)
    # p_h_M0, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_M0)
    # p_h_eta0, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_eta0)
    # p_h_x10, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_x10)
    # p_h_x20, _, _, _, _, _, _, _, _, _, _, _ = get_waveform_lm(para_x20)
    
    p_h_par_M = par(p_h_M, p_h, p_delta_M)
    p_h_par_eta = par(p_h_eta, p_h, p_delta_eta)
    # p_h_par_M = par(p_h_M, p_h_M0, p_delta_M)
    # p_h_par_eta = par(p_h_eta, p_h_eta0, p_delta_eta)
    # p_h_par_M = par(p_h_M, p_h_M0, 2*p_delta_M)
    # p_h_par_eta = par(p_h_eta, p_h_eta0, 2*p_delta_eta)
    # p_h_par_M = median_complex(par(p_h_M, p_h, p_delta_M), par(p_h_M, p_h_M0, p_delta_M), par(p_h_M, p_h_M0, 2*p_delta_M))
    # p_h_par_eta = median_complex(par(p_h_eta, p_h, p_delta_eta), par(p_h_eta, p_h_eta0, p_delta_eta), par(p_h_eta, p_h_eta0, 2*p_delta_eta))
    p_h_par_D_l = p_h * (-p_D_l)**(-1)
    p_h_par_t_c = p_h * (1j * 2 * np.pi * p_f)
    p_h_par_phi_c = p_h * (-1j)
    p_h_par_x1 = par(p_h_x1, p_h, p_delta_x1)
    p_h_par_x2 = par(p_h_x2, p_h, p_delta_x2)
    # p_h_par_x1 = par(p_h_x1, p_h_x10, p_delta_x1)
    # p_h_par_x2 = par(p_h_x2, p_h_x20, p_delta_x2)
    # p_h_par_x1 = par(p_h_x1, p_h_x10, 2*p_delta_x1)
    # p_h_par_x2 = par(p_h_x2, p_h_x20, 2*p_delta_x2)
    # p_h_par_x1 = median_complex(par(p_h_x1, p_h, p_delta_x1), par(p_h_x1, p_h_x10, p_delta_x1), par(p_h_x1, p_h_x10, 2*p_delta_x1))
    # p_h_par_x2 = median_complex(par(p_h_x2, p_h, p_delta_x2), par(p_h_x2, p_h_x20, p_delta_x2), par(p_h_x2, p_h_x20, 2*p_delta_x2))
    alpha_dCS_ctb_vec = np.vectorize(alpha_dCS_ctb)(p_f, p_para)  # 矢量化工具装饰
    alpha_dCS_par_vec = np.vectorize(alpha_dCS_par)(p_f, p_para)
    beta_dCS_par22 = np.vectorize(beta_dCS_par)(p_f, p_para, 2, 2)
    beta_dCS_par21 = np.vectorize(beta_dCS_par)(p_f, p_para, 2, 1)
    beta_dCS_par33 = np.vectorize(beta_dCS_par)(p_f, p_para, 3, 3)
    beta_dCS_par32 = np.vectorize(beta_dCS_par)(p_f, p_para, 3, 2)
    beta_dCS_par44 = np.vectorize(beta_dCS_par)(p_f, p_para, 4, 4)
    beta_dCS_ctb22 = np.vectorize(beta_dCS_ctb)(p_f, p_para, 2, 2)
    beta_dCS_ctb21 = np.vectorize(beta_dCS_ctb)(p_f, p_para, 2, 1)
    beta_dCS_ctb33 = np.vectorize(beta_dCS_ctb)(p_f, p_para, 3, 3)
    beta_dCS_ctb32 = np.vectorize(beta_dCS_ctb)(p_f, p_para, 3, 2)
    beta_dCS_ctb44 = np.vectorize(beta_dCS_ctb)(p_f, p_para, 4, 4)
    term22 = hp22 * (alpha_dCS_ctb_vec * beta_dCS_par22 + alpha_dCS_par_vec * beta_dCS_ctb22) * swsh_lm(2, 2, p_iota)
    term21 = hp21 * (alpha_dCS_ctb_vec * beta_dCS_par21 + alpha_dCS_par_vec * beta_dCS_ctb21) * swsh_lm(2, 1, p_iota)
    term33 = hp33 * (alpha_dCS_ctb_vec * beta_dCS_par33 + alpha_dCS_par_vec * beta_dCS_ctb33) * swsh_lm(3, 3, p_iota)
    term32 = hp32 * (alpha_dCS_ctb_vec * beta_dCS_par32 + alpha_dCS_par_vec * beta_dCS_ctb32) * swsh_lm(3, 2, p_iota)
    term44 = hp44 * (alpha_dCS_ctb_vec * beta_dCS_par44 + alpha_dCS_par_vec * beta_dCS_ctb44) * swsh_lm(4, 4, p_iota)
    p_h_par_alpha_dCS_bar_square = (term22 + term21 + term33 + term32 + term44)


    p_h_par_M_ast = [h.conjugate() for h in p_h_par_M]
    p_h_par_eta_ast = [h.conjugate() for h in p_h_par_eta]
    p_h_par_D_l_ast = [h.conjugate() for h in p_h_par_D_l]
    p_h_par_t_c_ast = [h.conjugate() for h in p_h_par_t_c]
    p_h_par_phi_c_ast = [h.conjugate() for h in p_h_par_phi_c]
    p_h_par_x1_ast = [h.conjugate() for h in p_h_par_x1]
    p_h_par_x2_ast = [h.conjugate() for h in p_h_par_x2]
    p_h_par_alpha_dCS_bar_square_ast = [h.conjugate() for h in p_h_par_alpha_dCS_bar_square]

    p_h_par = [p_h_par_M, p_h_par_eta, p_h_par_D_l, p_h_par_t_c, p_h_par_phi_c, p_h_par_x1, p_h_par_x2, p_h_par_alpha_dCS_bar_square]
    p_h_par_ast = [p_h_par_M_ast, p_h_par_eta_ast, p_h_par_D_l_ast, p_h_par_t_c_ast, p_h_par_phi_c_ast, p_h_par_x1_ast, p_h_par_x2_ast, p_h_par_alpha_dCS_bar_square_ast]

    return p_h_par, p_h_par_ast, p_f

# 内积
def quaddiy(func1, func2, p_f, l1, l2, p_para, p_Sn_dector):
    """计算单探测器下两个波形偏导之间的噪声加权内积。"""
    q = 0
    lenth = len(p_f)
    p_deltaF = p_para['deltaF']
    for i in range(lenth):
        q += (func1[l1][i]*func2[l2][i]+func2[l1][i]*func1[l2][i]) * 2 / p_Sn_dector(p_f[i]) * p_deltaF
    return q.real

def quaddiy_two(func1, func2, p_f, l1, l2, p_para, p_Sn_dector1, p_Sn_dector2):
    """计算双探测器联合噪声权重下的波形偏导内积。"""
    q = 0
    lenth = len(p_f)
    p_deltaF = p_para['deltaF']
    for i in range(lenth):
        q += (func1[l1][i]*func2[l2][i]+func2[l1][i]*func1[l2][i]) * 2 * (1 / p_Sn_dector1(p_f[i]) + 1 / p_Sn_dector2(p_f[i])) * p_deltaF
    return q.real

def quaddiypara(func1, func2, p_f, p_para, p_Sn_dector):
    """计算两个一维波形数组之间的单探测器噪声加权内积。"""
    q = 0
    lenth = len(p_f)
    p_deltaF = p_para['deltaF']
    for i in range(lenth):
        q += (func1[i]*func2[i]+func2[i]*func1[i]) * 2 / p_Sn_dector(p_f[i]) * p_deltaF
    return q.real

# 信噪比
def snr(p_para, p_Sn_dector):
    """计算单探测器匹配滤波信噪比。"""
    p_h, _, _, _, _, _, _, _, _, _, _, p_f = get_waveform_lm(p_para)
    q = 0
    lenth = len(p_f)
    p_deltaF = p_para['deltaF']
    for i in range(lenth):
        q += (p_h[i] * p_h[i].conjugate()) * 4 / p_Sn_dector(p_f[i]) * p_deltaF
    return np.sqrt(q.real)

def snr_two(p_para, p_Sn_dector1, p_Sn_dector2):
    """计算两个探测器联合噪声权重下的信噪比。"""
    p_h, _, _, _, _, _, _, _, _, _, _, p_f = get_waveform_lm(p_para)
    q = 0
    lenth = len(p_f)
    p_deltaF = p_para['deltaF']
    for i in range(lenth):
        q += (p_h[i] * p_h[i].conjugate()) * 4 * (1 / p_Sn_dector1(p_f[i]) + 1 / p_Sn_dector2(p_f[i])) * p_deltaF
    return np.sqrt(q.real)

# fisher矩阵
def fisher_matrix(para, para_delta, p_Sn_dector):
    """计算单探测器 Fisher 信息矩阵。"""
    h_par, h_par_ast, f = get_h_par(para, para_delta)
    pp = len(h_par)
    fm = np.empty((pp, pp), dtype=complex)
    for p1 in range(pp):
        for p2 in range(pp):
            if p1 <= p2:
                fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector)
            else:
                fm[p1, p2] = fm[p2, p1]
    return fm
# def fisher_matrix(para, para_delta, p_Sn_dector):
#     h_par, h_par_ast, f = get_h_par(para, para_delta)
#     pp = len(h_par)
#     fm = np.empty((pp, pp), dtype=complex)
#     for p1 in range(pp):
#         for p2 in range(pp):
#             if p1 < p2:
#                 if p1 == 0 or p2 == 0:
#                     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector) * MSUN_SI
#                 elif p1 == 2 or p2 == 2:
#                     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector) * PC_SI
#                 # elif p1 == 7 or p2 == 7:
#                 #     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector) * 1e12
#                 else:
#                     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector)
#             elif p1 == p2:
#                 if p1 == 0:
#                     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector) * MSUN_SI**2
#                 elif p1 == 2:
#                     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector) * PC_SI**2
#                 # elif p1 == 7:
#                 #     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector) * 1e24
#                 else:
#                     fm[p1, p2] = quaddiy(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector)
#             else:
#                 fm[p1, p2] = fm[p2, p1]
#     return fm

def fisher_matrix_two(para, para_delta, p_Sn_dector1, p_Sn_dector2):
    """计算双探测器联合 Fisher 信息矩阵。"""
    h_par, h_par_ast, f = get_h_par(para, para_delta)
    pp = len(h_par)
    fm = np.empty((pp, pp), dtype=complex)
    for p1 in range(pp):
        for p2 in range(pp):
            if p1 <= p2:
                fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2)
            else:
                fm[p1, p2] = fm[p2, p1]
    return fm
# def fisher_matrix_two(para, para_delta, p_Sn_dector1, p_Sn_dector2):
#     h_par, h_par_ast, f = get_h_par(para, para_delta)
#     pp = len(h_par)
#     fm = np.empty((pp, pp), dtype=complex)
#     for p1 in range(pp):
#         for p2 in range(pp):
#             if p1 < p2:
#                 if p1 == 0 or p2 == 0:
#                     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2) * MSUN_SI
#                 elif p1 == 2 or p2 == 2:
#                     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2) * PC_SI
#                 # elif p1 == 7 or p2 == 7:
#                 #     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2) * 1e12
#                 else:
#                     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2)
#             elif p1 == p2:
#                 if p1 == 0:
#                     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2) * MSUN_SI**2
#                 elif p1 == 2:
#                     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2) * PC_SI**2
#                 # elif p1 == 7:
#                 #     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2) * 1e24
#                 else:
#                     fm[p1, p2] = quaddiy_two(h_par, h_par_ast, f, p1, p2, para, p_Sn_dector1, p_Sn_dector2)
#             else:
#                 fm[p1, p2] = fm[p2, p1]
#     return fm

# 协方差矩阵
def covariance(para, para_delta, p_Sn_dector):
    """由单探测器 Fisher 矩阵求协方差矩阵。"""
    fm = fisher_matrix(para, para_delta, p_Sn_dector)
    try:
        covar = linalg.inv(fm)
    except linalg.LinAlgError:
        print(f"Parameters: {para}, {para_delta}, {p_Sn_dector}", "fm奇异")
        covar = None
    return covar

def covariance_two(para1, para2, para_delta, p_Sn_dector1, p_Sn_dector2):
    """由两个频段/探测器的 Fisher 矩阵之和求协方差矩阵。"""
    fm1 = fisher_matrix_two(para1, para_delta, p_Sn_dector1, p_Sn_dector2)
    fm2 = fisher_matrix_two(para2, para_delta, p_Sn_dector1, p_Sn_dector2)
    try:
        covar = linalg.inv(fm1+fm2)
    except linalg.LinAlgError:
        print(f"Parameters: {para1}, {para2}, {para_delta}, {p_Sn_dector1}, {p_Sn_dector2}", "fm奇异")
        covar = None
    return covar

def get_detecror_parameters(p_detector):
    """返回指定探测器的频段、扫描范围、噪声函数和差分步长。

    Args:
        p_detector: 探测器名称，如 TianQin、LISA、Taiji、DECIGO、
            BBO、ET、CE 或 aLIGO。

    Returns:
        包含频率上下限、质量范围、距离范围、噪声谱函数和差分步长的元组。
    """
    if p_detector == 'TianQin':
        p_f_D_low = 1e-4
        p_f_D_high = 1
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 5
        p_D_l_len_n = 5
        p_Sn_dector = Sn.sensitivity_curves_TianQin_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'LISA':
        p_f_D_low = 1e-4
        p_f_D_high = 1
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 5
        p_D_l_len_n = 5
        p_Sn_dector = Sn.sensitivity_curves_LISA_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'Taiji':
        p_f_D_low = 1e-4
        p_f_D_high = 1
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 5
        p_D_l_len_n = 5
        p_Sn_dector = Sn.sensitivity_curves_Taiji_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'DECIGO':
        p_f_D_low = 1e-3
        p_f_D_high = 1e2
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 6
        p_D_l_len_n = 6
        p_Sn_dector = Sn.sensitivity_curves_DECIGO_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'BBO':
        p_f_D_low = 1e-3
        p_f_D_high = 1e2
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 6
        p_D_l_len_n = 6
        p_Sn_dector = Sn.sensitivity_curves_BBO_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'ET':
        p_f_D_low = 4
        p_f_D_high = 4999
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 6
        p_D_l_len_n = 6
        p_Sn_dector = Sn.sensitivity_curves_ET_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'CE':
        p_f_D_low = 4
        p_f_D_high = 4999
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 6
        p_D_l_len_n = 6
        p_Sn_dector = Sn.sensitivity_curves_CE_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-10
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    elif p_detector == 'aLIGO':
        p_f_D_low = 4
        p_f_D_high = 4999
        p_M_low_n = 1
        p_M_high_n = 1000
        p_D_l_low_n = 6
        p_D_l_len_n = 4
        p_Sn_dector = Sn.sensitivity_curves_aLIGO_sn
        p_delta_M_n = 1e-7
        p_delta_eta_n = 1e-9
        p_delta_x1_n = 1e-5
        p_delta_x2_n = 1e-5
    else:
        print('请重新选择探测器')
    return p_f_D_low, p_f_D_high, p_M_low_n, p_M_high_n, p_D_l_low_n, p_D_l_len_n, p_Sn_dector, p_delta_M_n, p_delta_eta_n, p_delta_x1_n, p_delta_x2_n


def median_complex(c1, c2, c3):
    """逐元素取三个复数组实部和虚部的中位数。"""
    comp = []
    for i in range(len(c1)):
        median_real = sorted([c1[i].real, c2[i].real, c3[i].real])[1]
        median_imag = sorted([c1[i].imag, c2[i].imag, c3[i].imag])[1]
        comp.append(complex(median_real, median_imag))
    return comp
