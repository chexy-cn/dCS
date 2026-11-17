"""自旋权球谐函数取值工具。

当前实现显式给出 l=2、3、4 且常用 m 模式的表达式，用于把
IMRPhenomXHM 的分模波形投影到给定倾角 iota。
"""

import numpy as np
from scipy.constants import pi

def swsh_lm(p_l, p_m, p_iota):
    """计算指定 (l, m) 模式的自旋权球谐函数值。

    Args:
        p_l: 球谐模式的 l。
        p_m: 球谐模式的 m。
        p_iota: 倾角，单位弧度。

    Returns:
        复数形式的自旋权球谐函数值；未实现的模式返回 None。
    """
    if p_l==2:
        if p_m==-2:
            return np.sqrt(5.0/(64.0*pi))*(1.0-np.cos(p_iota))*(1.0-np.cos(p_iota))*np.exp(-1j*p_m*pi/2)
        elif p_m==-1:
            return np.sqrt(5.0/(16.0*pi))*np.sin(p_iota)*(1.0-np.cos(p_iota))*np.exp(-1j*p_m*pi/2)
        elif p_m== 0:
            return np.sqrt(15.0/(32.0*pi))*np.sin(p_iota)*np.sin(p_iota)
        elif p_m== 1:
            return np.sqrt(5.0/(16.0*pi))*np.sin(p_iota)*(1.0+np.cos(p_iota))*np.exp(-1j*p_m*pi/2)
        elif p_m== 2:
            return np.sqrt(5.0/(64.0*pi))*(1.0+np.cos(p_iota))*(1.0+np.cos(p_iota))*np.exp(-1j*p_m*pi/2)
    elif p_l==3:
        if p_m==-3:
            return np.sqrt(21.0/(2.0*pi))*np.cos(p_iota/2.0)*pow(np.sin(p_iota/2.0),5.0)*np.exp(-1j*p_m*pi/2)
        elif p_m==-2:
            return np.sqrt(7.0/(4.0*pi))*(2.0+3.0*np.cos(p_iota))*pow(np.sin(p_iota/2.0),4.0)*np.exp(-1j*p_m*pi/2)
        elif p_m==-1:
            return np.sqrt(35.0/(2.0*pi))*(np.sin(p_iota)+4.0*np.sin(2.0*p_iota)-3.0*np.sin(3.0*p_iota))/32.0*np.exp(-1j*p_m*pi/2)
        elif p_m==0:
            return (np.sqrt(105.0/(2.0*pi))*np.cos(p_iota)*pow(np.sin(p_iota),2.0))/4.0
        elif p_m==1:
            return -np.sqrt(35.0/(2.0*pi))*(np.sin(p_iota)-4.0*np.sin(2.0*p_iota)-3.0*np.sin(3.0*p_iota))/32.0*np.exp(-1j*p_m*pi/2)
        elif p_m==2:
            return np.sqrt(7.0/pi)*pow(np.cos(p_iota/2.0),4.0)*(-2.0+3.0*np.cos(p_iota))/2.0*np.exp(-1j*p_m*pi/2)
        elif p_m==3:
            return -np.sqrt(21.0/(2.0*pi))*pow(np.cos(p_iota/2.0),5.0)*np.sin(p_iota/2.0)*np.exp(-1j*p_m*pi/2)
    elif p_l==4:
        if p_m==-4:
            return 3.0*np.sqrt(7.0/pi)*pow(np.cos(p_iota/2.0),2.0)*pow(np.sin(p_iota/2.0),6.0)*np.exp(-1j*p_m*pi/2)
        elif p_m==-3:
            return 3.0*np.sqrt(7.0/(2.0*pi))*np.cos(p_iota/2.0)*(1.0+2.0*np.cos(p_iota))*pow(np.sin(p_iota/2.0),5.0)*np.exp(-1j*p_m*pi/2)
        elif p_m==-2:
            return (3.0*(9.0+14.0*np.cos(p_iota)+7.0*np.cos(2.0*p_iota))*pow(np.sin(p_iota/2.0),4.0))/(4.0*np.sqrt(pi))*np.exp(-1j*p_m*pi/2)
        elif p_m==-1:
            return (3.0*(3.0*np.sin(p_iota)+2.0*np.sin(2.0*p_iota)+7.0*np.sin(3.0*p_iota)-7.0*np.sin(4.0*p_iota)))/(32.0*np.sqrt(2.0*pi))*np.exp(-1j*p_m*pi/2)
        elif p_m==0:
            return (3.0*np.sqrt(5.0/(2.0*pi))*(5.0+7.0*np.cos(2.0*p_iota))*pow(np.sin(p_iota),2.0))/16.0
        elif p_m==1:
            return (3.0*(3.0*np.sin(p_iota)-2.0*np.sin(2.0*p_iota)+7.0*np.sin(3.0*p_iota)+7.0*np.sin(4.0*p_iota)))/(32.0*np.sqrt(2.0*pi))*np.exp(-1j*p_m*pi/2)
        elif p_m==2:
            return (3.0*pow(np.cos(p_iota/2.0),4.0)*(9.0-14.0*np.cos(p_iota)+7.0*np.cos(2.0*p_iota)))/(4.0*np.sqrt(pi))*np.exp(-1j*p_m*pi/2)
        elif p_m==3:
            return -3.0*np.sqrt(7.0/(2.0*pi))*pow(np.cos(p_iota/2.0),5.0)*(-1.0+2.0*np.cos(p_iota))*np.sin(p_iota/2.0)*np.exp(-1j*p_m*pi/2)
        elif p_m==4:
            return 3.0*np.sqrt(7.0/pi)*pow(np.cos(p_iota/2.0),6.0)*pow(np.sin(p_iota/2.0),2.0)*np.exp(-1j*p_m*pi/2)

if __name__ == "__main__":
    iota = pi
    print(f"SWSH Y_{2}^{2} at iota={iota}: {swsh_lm(2, 2, iota)}")
    print(f"SWSH Y_{2}^{1} at iota={iota}: {swsh_lm(2, 1, iota)}")
    print(f"SWSH Y_{3}^{3} at iota={iota}: {swsh_lm(3, 3, iota)}")
    print(f"SWSH Y_{3}^{2} at iota={iota}: {swsh_lm(3, 2, iota)}")
    print(f"SWSH Y_{4}^{4} at iota={iota}: {swsh_lm(4, 4, iota)}")
