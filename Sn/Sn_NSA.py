"""非全天平均探测器噪声谱模型。

本文件与 Sn.py 的接口基本一致，但空间探测器噪声谱不乘全天平均
因子，适用于非全天平均近似下的灵敏度估计。
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.constants import c, pi
import numpy as np
import os
mpl.use('Agg')

# 非全天平均

def sensitivity_curves_TianQin_sn(p_f):
    """计算非全天平均 TianQin 噪声功率谱密度。

    Args:
        p_f: 频率或频率数组，单位 Hz。

    Returns:
        与输入频率同形状的 S_n(f)。
    """
    L = pow(3, 0.5) * 1e8  # m
    f_ast = c / (2 * pi * L)
    S_a = 1e-30
    S_x = 1e-24
    S_N = 1 / pow(L, 2) * (4 * S_a / pow((2 * pi * p_f), 4) * (1 + 1e-4 / p_f) + S_x)
    S_n = S_N * (1 + 0.6 * pow((p_f / f_ast), 2))
    # a1, a2, a3, a4, f_k = [0.138, -221, 521, 1680, 0.00113]  # 4year
    # # a1, a2, a3, a4, f_k = [0.165,299,611,1340,0.00173]  # 2year
    # # a1, a2, a3, a4, f_k = [0.171,292,1020,1680,0.00215]  # 1year
    # # a1, a2, a3, a4, f_k = [0.133, 243, 482, 917, 0.00258]  # 0.5year
    # S_c = 9e-45 * p_f ** (-7 / 3) * np.exp(-p_f ** a1 + a2 * p_f * np.sin(a3 * p_f)) * (1 + np.tanh(a4 * (f_k - p_f)))
    return S_n



def sensitivity_curves_LISA_sn(f):
    """计算非全天平均 LISA 噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        与输入频率同形状的 S_n(f)，当前默认不加入混淆噪声。
    """
    P_OMS = 1.5e-11 ** 2 * (1 + (0.002 / f) ** 4)
    P_acc = 3e-15 ** 2 * (1 + (0.0004 / f) ** 2) * (
            1 + (f / 0.008) ** 4)
    L = 2.5e9
    f_ast = c / (2 * pi * L)
    # f_ast = 0.01909
    a1, a2, a3, a4, f_k = [0.138, -221, 521, 1680, 0.00113]  # 4year
    # a1, a2, a3, a4, f_k = [0.165,299,611,1340,0.00173]  # 2year
    # a1, a2, a3, a4, f_k = [0.171,292,1020,1680,0.00215]  # 1year
    # a1, a2, a3, a4, f_k = [0.133, 243, 482, 917, 0.00258]  # 0.5year
    S_c = 9e-45 * f ** (-7 / 3) * np.exp(-f ** a1 + a2 * f * np.sin(a3 * f)) * (1 + np.tanh(a4 * (f_k - f)))
    S_n = 1 / L ** 2 * (P_OMS + 2 * (1 + (np.cos(f / f_ast)) ** 2) * P_acc / (2 * pi * f) ** 4) * (
            1 + 6 / 10 * (f / f_ast) ** 2)
    # return S_n + S_c
    return S_n



def sensitivity_curves_Taiji_sn(f):
    """计算非全天平均 Taiji 噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        与输入频率同形状的 S_n(f)，当前默认不加入混淆噪声。
    """
    P_OMS = 8e-12 ** 2 * (1 + (0.002 / f) ** 4)
    P_acc = 3e-15 ** 2 * (1 + (0.0004 / f) ** 2) * (1 + (f / 0.008) ** 4)
    L = 3e9
    f_ast = c / (2 * pi * L)
    a1, a2, a3, a4, f_k = [0.138, -221, 521, 1680, 0.00113]  # 4year
    # a1, a2, a3, a4, f_k = [0.165,299,611,1340,0.00173]  # 2year
    # a1, a2, a3, a4, f_k = [0.171,292,1020,1680,0.00215]  # 1year
    # a1, a2, a3, a4, f_k = [0.133, 243, 482, 917, 0.00258]  # 0.5year
    S_c = 9e-45 * f ** (-7 / 3) * np.exp(-f ** a1 + a2 * f * np.sin(a3 * f)) * (1 + np.tanh(a4 * (f_k - f)))
    S_n = 1 / L ** 2 * (P_OMS + 2 * (1 + (np.cos(f / f_ast)) ** 2) * P_acc / (2 * pi * f) ** 4) * (
            1 + 6 / 10 * (f / f_ast) ** 2)
    # return S_n + S_c
    return S_n
    # S_n = (4 * S_a * (1 + 1e-4 / f) / math.pow(2 * math.pi * f, 4) + S_x) / math.pow(L, 2)


def sensitivity_curves_DECIGO_sn(f):
    """计算非全天平均 DECIGO 噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        与输入频率同形状的 S_n(f)。
    """
    f_xing = 7.36
    S_n = 7.05e-48 * (1 + (f / f_xing) ** 2) + ((4.8e-51) * f ** (-4)) / (1 + (f / f_xing) ** 2) + 5.33E-52 * f ** (-4)
    S_n = S_n * 3 / 4 * 1
    return S_n


def sensitivity_curves_BBO_sn(f):
    """计算非全天平均 BBO 噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        与输入频率同形状的 S_n(f)。
    """
    S_n = 2e-49 * f ** 2 + 4.58e-49 + 1.26e-52 * f ** (-4)
    S_n = S_n * 3 / 4 * 1
    return S_n


# Sn.py扩展部分
# 预加载ET敏感曲线数据
_ET_freq, _ET_Sn = [], []
_aLIGO_freq, _aLIGO_Sn = [], []
_CE_freq, _CE_Sn = [], []
# ▶只在模块导入时读取一次文件◀
with open('/home/chexinyi/code/dCS/Sn/data/ET_D.txt','r') as f:
    for line in f:
        parts = line.split()
        _ET_freq.append(float(parts[0]))
        _ET_Sn.append(float(parts[1]))

# 生成插值器(全局唯一)
_ET_interpolator = interpolate.interp1d(
    _ET_freq, _ET_Sn, 
    bounds_error=False, 
    fill_value="extrapolate"
)

with open('/home/chexinyi/code/dCS/Sn/data/aLIGO_D.txt','r') as f:
    for line in f:
        parts = line.split()
        _aLIGO_freq.append(float(parts[0]))
        _aLIGO_Sn.append(float(parts[1]))

_aLIGO_interpolator = interpolate.interp1d(
    _aLIGO_freq, _aLIGO_Sn, 
    fill_value="extrapolate",
    bounds_error=False
)

with open('/home/chexinyi/code/dCS/Sn/data/CE_D40.txt','r') as f:
    for line in f:
        parts = line.split()
        _CE_freq.append(float(parts[0]))
        _CE_Sn.append(float(parts[1]))

_CE_interpolator = interpolate.interp1d(
    _CE_freq, _CE_Sn, 
    fill_value="extrapolate",
    bounds_error=False

)

def sensitivity_curves_ET_sn(f):
    """通过 ET 灵敏度数据插值得到噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        ET 的 S_n(f)。超出有效频段时返回较大的占位噪声。
    """
    f = np.asarray(f)
    mask = (f < 1) | (f > 1e4)
    Sensitivity = np.where(mask, 100, _ET_interpolator(f))
    return Sensitivity**2


def sensitivity_curves_aLIGO_sn(f):
    """通过 aLIGO 灵敏度数据插值得到噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        aLIGO 的 S_n(f)。超出有效频段时返回较大的占位噪声。
    """
    f = np.asarray(f)  # 确保输入转为NumPy数组
    mask = (f < 6) | (f > 4999)
    Sensitivity = np.where(mask, 100, _aLIGO_interpolator(f))  # 元素级条件选择
    return Sensitivity**2


def sensitivity_curves_CE_sn(f):
    """通过 CE 灵敏度数据插值得到噪声功率谱密度。

    Args:
        f: 频率或频率数组，单位 Hz。

    Returns:
        CE 的 S_n(f)。超出有效频段时返回较大的占位噪声。
    """
    f = np.asarray(f)  # 确保输入转为NumPy数组
    mask = (f < 6) | (f > 4999)
    Sensitivity = np.where(mask, 100, _CE_interpolator(f))  # 元素级条件选择
    return Sensitivity**2


if __name__ == '__main__':
    date = 20250603

    f1 = np.logspace(-4, 0, 1000)
    TianQin = np.sqrt(sensitivity_curves_TianQin_sn(f1))
    LISA = np.sqrt(sensitivity_curves_LISA_sn(f1))
    Taiji = np.sqrt(sensitivity_curves_Taiji_sn(f1))
    DECIGO = np.sqrt(sensitivity_curves_DECIGO_sn(f1))
    BBO = np.sqrt(sensitivity_curves_BBO_sn(f1))
    directory = f'/home/chexinyi/code/dCS/data/{date}/results/'
    os.makedirs(directory, exist_ok=True)
    file_path = os.path.join(directory, f'Sn.txt')
    with open(file_path, 'w') as f:
        for i in range(1000):
            line = f'{f1[i]}\t{TianQin[i]}\t{Taiji[i]}\t{LISA[i]}\t{DECIGO[i]}\t{BBO[i]}\n'
            f.write(line)



