"""结果文件读取和数组整理工具。

本文件集中处理批量计算输出的逗号/制表符文本结果，并提供排序、
二维重排、事件参数读取等辅助函数。
"""

import numpy as np

def parse_number(s):
    """把普通浮点字符串或带括号复数字符串解析为实数。"""
    s = s.strip()
    # 处理带括号的复数（例如 "(23478.045329151766+0j)"）
    if s.startswith('(') and s.endswith(')'):
        return complex(s[1:-1]).real  # 去括号后转换复数取实部
    return float(s)  # 常规数值直接转换

def read_results1(file_path):
    """读取主结果文件并返回网格坐标、8 个不确定度和 SNR 等数组。"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"文件 {file_path} 不存在")

    # 初始化存储数组
    p_x_list = []
    p_l_list = []
    array_delta = [[] for _ in range(8)]  # 8个delta参数
    M_list = []
    D_l_list = []
    alpha_list = []
    SNR_list = []

    for line in lines:
        parts = list(map(parse_number, line.strip().split(',')))
        # parts = list(map(parse_number, line.strip().split()))
        if len(parts) != 14:  # p_x, p_l, 8 deltas, M_total, D_l, alpha
            continue

        p_x_list.append(parts[0])
        p_l_list.append(parts[1])
        
        # 8个delta参数
        for i in range(8):
            array_delta[i].append(parts[2+i])
        M_list.append(parts[10])
        D_l_list.append(parts[11])
        alpha_list.append(parts[12])
        SNR_list.append(parts[13])

    return (
        np.array(p_x_list),
        np.array(p_l_list),
        [np.array(d) for d in array_delta],
        np.array(M_list),
        np.array(D_l_list),
        np.array(alpha_list),
        np.array(SNR_list)
    )

def read_results2(file_path):
    """读取三列结果文件，返回 x、D 和二次型/内积数组。"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"文件 {file_path} 不存在")

    # 初始化存储数组
    x_list = []
    D_list = []
    quad_list = []

    for line in lines:
        parts = list(map(parse_number, line.strip().split(',')))

        x_list.append(parts[0])
        D_list.append(parts[1])
        quad_list.append(parts[2])

    return (np.array(x_list), np.array(D_list), np.array(quad_list))

# def read_events(file_path):
#     try:
#         with open(file_path, 'r') as f:
#             lines = f.readlines()
#     except FileNotFoundError:
#         raise ValueError(f"文件 {file_path} 不存在")
#
#     Dl, m1, m2, Tc, theta_s = [[] for _ in range(5)]
#
#     for line in lines:
#         parts = line.strip().split()
#         Dl.append(float(parts[0]))
#         m1.append(float(parts[1]))
#         m2.append(float(parts[2]))
#         Tc.append(float(parts[3]))
#         theta_s.append(float(parts[4]))
#
#     return (np.array(Dl), np.array(m1), np.array(m2), np.array(Tc), np.array(theta_s))


def reshape_array(arr):
    """把长度为 10201 的一维数组重排为 101x101 网格。"""
    if len(arr) != 10201:
        raise ValueError(f"数组长度必须为10201，当前为{len(arr)}")
    dim = int(np.sqrt(len(arr)))
    return np.asarray(arr).reshape((dim, dim))

def sort_multiple_arrays(*arrays, key_array=None):
    """
    按照第一个数组（或指定的key数组）的升序排列多个数组
    
    参数:
        *arrays: 多个numpy数组或列表，所有数组长度必须相同
        key_array (optional): 如果指定，则使用此数组作为排序基准而非第一个数组
    
    返回:
        tuple: 按第一个数组（或key_array）升序排列后的所有数组
        
    示例用法:
        a = [3,1,2]
        b = ['c','a','b']
        c = [30,10,20]
        sa, sb, sc = sort_multiple_arrays(a, b, c)
        # 或指定key_array:
        sb, sc = sort_multiple_arrays(b, c, key_array=a)
    """
    if not arrays:
        return ()
        
    # 如果指定了key_array，使用它作为排序基准
    if key_array is not None:
        sort_key = np.array(key_array)
    else:
        sort_key = np.array(arrays[0])
    
    # 确保所有数组都有相同长度
    lengths = [len(arr) for arr in arrays]
    if len(set(lengths)) > 1:
        raise ValueError("所有输入数组必须具有相同的长度")
    
    # 获取排序索引
    sort_indices = np.argsort(sort_key)
    
    # 对所有数组应用相同的排序
    sorted_arrays = tuple(np.array(arr)[sort_indices] for arr in arrays)
    
    return sorted_arrays[0] if len(sorted_arrays) == 1 else sorted_arrays

def flatten_2d_array(two_d_array):
    """
    将二维数组转换为一维数组
    
    参数:
        two_d_array (list): 输入的二维数组
    
    返回:
        list: 扁平化后的一维数组
    
    示例:
        >>> flatten_2d_array([[1, 2], [3, 4]])
        [1, 2, 3, 4]
    """
    return [item for sublist in two_d_array for item in sublist]

def read_events(p_flie_path):
    """读取事件参数文件并返回形状为 `(N, 5)` 的参数数组。

    每行格式为 `Dl m1 m2 Tc theta_s`。
    """
    with open(p_flie_path, 'r') as f:
        lines = f.readlines()

    Dl, m1, m2, Tc, theta_s = [[] for _ in range(5)]

    for line in lines:
        parts = line.strip().split()
        Dl.append(float(parts[0]))
        m1.append(float(parts[1]))
        m2.append(float(parts[2]))
        Tc.append(float(parts[3]))
        theta_s.append(float(parts[4]))        
        
    Dl, m1, m2, Tc, theta_s = np.array(Dl), np.array(m1), np.array(m2), np.array(Tc), np.array(theta_s)

    parameters = np.empty([len(Dl), 5])
    for i in range(len(Dl)):
        parameters[i][0] = Dl[i]
        parameters[i][1] = m1[i]
        parameters[i][2] = m2[i]
        parameters[i][3] = Tc[i]
        parameters[i][4] = theta_s[i]
    
    return parameters

def read_results3(file_path):
    """读取两列制表符分隔结果，返回偏心率和约束量数组。"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"文件 {file_path} 不存在")

    # 初始化存储数组
    e_list = []
    delta_list = []

    for line in lines:
        parts = list(map(parse_number, line.strip().split('\t')))

        e_list.append(parts[0])
        delta_list.append(parts[1])

    return (np.array(e_list), np.array(delta_list))

def read_results_events(file_path):
    """读取事件约束结果，返回事件参数、约束值、有效性上界和 SNR。"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"文件 {file_path} 不存在")

    parameters = [[] for _ in range(5)]
    alpha_list = []
    alpha_max_list = []
    SNR_list = []

    for line in lines:
        parts = line.strip().split(',')
        # parts = list(map(parse_number, line.strip().split(',')))
        # parts = list(map(parse_number, line.strip().split()))

        for i in range(5):
            parameters[i].append(float(parts[i]))
        alpha_list.append(float(parts[5]))
        alpha_max_list.append(float(parts[6]))
        SNR_list.append(float(parts[7]))

    return (
        [np.array(d) for d in parameters],
        np.array(alpha_list),
        np.array(alpha_max_list),
        np.array(SNR_list)
    )

def read_results_events0(file_path):
    """读取空格分隔的事件结果，返回事件参数和 alpha 数组。"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"文件 {file_path} 不存在")

    Dl_list = []
    m1_list = []
    m2_list = []
    Tc_list = []
    theta_s_list = []
    alpha_list = []

    for line in lines:
        parts = line.strip().split()

        Dl_list.append(float(parts[0]))
        m1_list.append(float(parts[1]))
        m2_list.append(float(parts[2]))
        Tc_list.append(float(parts[3]))
        theta_s_list.append(float(parts[4]))
        alpha_list.append(float(parts[5]))


    return (
        np.array(Dl_list),
        np.array(m1_list),
        np.array(m2_list),
        np.array(Tc_list),
        np.array(theta_s_list),
        np.array(alpha_list)
    )

def read_results_events01(file_path):
    """读取逗号分隔的事件结果，返回事件参数和 alpha 数组。"""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"文件 {file_path} 不存在")

    Dl_list = []
    m1_list = []
    m2_list = []
    Tc_list = []
    theta_s_list = []
    alpha_list = []

    for line in lines:
        parts = line.strip().split(',')

        Dl_list.append(float(parts[0]))
        m1_list.append(float(parts[1]))
        m2_list.append(float(parts[2]))
        Tc_list.append(float(parts[3]))
        theta_s_list.append(float(parts[4]))
        alpha_list.append(float(parts[6]))


    return (
        np.array(Dl_list),
        np.array(m1_list),
        np.array(m2_list),
        np.array(Tc_list),
        np.array(theta_s_list),
        np.array(alpha_list)
    )


if __name__ == '__main__':
    # a = [3,1,2]
    # b = ['b','a','c']
    # c = [30,10,20]
    # sa, sb, sc = sort_multiple_arrays(a, b, c)
    # print(sa, sb, sc)
    # k = np.logspace(0, 3, 101)
    # for i in range(101):
    #     a = np.append(a, k[i])
    # print(a)

    # array1, array2, array3 = [ [[0 for _ in range(3)] for _ in range(3)] 
    #       for _ in range(3) ]
    # print(array1, array2, array3)

    # shape = (3, 5, 5)  # 3个5x5数组
    array1, array2, array3 = np.zeros((3,5,5))
    print(array1, array2, array3)
