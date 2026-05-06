"""批量计算结果的绘图工具。

本文件读取 `data/<date>/results/` 下的文本结果，绘制 dCS 参数约束
等高线、SNR=8 阈值线和有效性边界，并保存到 `fig/` 目录。
"""

import os
import numpy as np
from fim import get_detecror_parameters
from getresults.readresult import read_results1, read_results2, reshape_array, read_results3, sort_multiple_arrays, read_results_events, read_results_events0,read_results_events01
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def plot_one_simplify(p_detector, q, date):
    """绘制单个探测器在多个质量比下的简化约束图。"""
    p_path = f'./data/{date}/fig/'
    os.makedirs(p_path, exist_ok=True)
    p_l = 1001
    # _, _, _, _, D_l_low_n, D_l_len_n, _, _, _, _, _ = get_detecror_parameters(p_detector)
    fig, ax = plt.subplots(figsize=(8, 6))
    for p_q in q:
        if p_q == q[0]:
            linestyle = '-'
            co = 'lime'
        else:
            linestyle = '--'
            co = 'red'
        data_path = f'./data/{date}/results/results for {p_detector} with q={p_q}.txt'
        _, _, deltas, M_n, _, alpha_dCS_bar_square_max, SNR = read_results1(data_path)
        d_alpha_dCS_bar_square = deltas[7]**0.25 * 1e-3

        M_n_sort, d_alpha_dCS_bar_square_sort, alpha_dCS_bar_square_max_sort, SNR_sort = sort_multiple_arrays(M_n, d_alpha_dCS_bar_square, alpha_dCS_bar_square_max, SNR)
        M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, alpha_dCS_bar_square_max_sort_twodim, SNR_sort_twodim=np.zeros((5, p_l, p_l))

        k = np.logspace(0, 6, p_l)
        for i in range(p_l):
            for j in range(p_l):
                M_n_sort_twodim[i][j] = M_n_sort[i]
                D_l_n_sort_twodim[i][j] = k[j]
                d_alpha_dCS_bar_square_sort_twodim[i][j] = d_alpha_dCS_bar_square_sort[i] * k[j]**0.25
                alpha_dCS_bar_square_max_sort_twodim[i][j] = alpha_dCS_bar_square_max_sort[i]
                SNR_sort_twodim[i][j] = SNR_sort[i] / k[j]

        M_n_validity, D_l_n_validity = [[] for _ in range(2)]
        for i in range(p_l):
            for j in range(p_l):
                # if d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1 < 1e-2 and d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1 > 0:
                if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < 1e-2:
                    M_n_validity.append(M_n_sort_twodim[i][j])
                    D_l_n_validity.append(D_l_n_sort_twodim[i][j])
        
        if p_detector == 'aLIGO':
            for lm in range(p_l):
                for ld in range(p_l):
                    if M_n_sort_twodim[lm][ld]>=600:
                        d_alpha_dCS_bar_square_sort_twodim[lm][ld]=999999999
        elif p_detector == 'CE':
            for lm in range(p_l):
                for ld in range(p_l):
                    if M_n_sort_twodim[lm][ld]>=700:
                        d_alpha_dCS_bar_square_sort_twodim[lm][ld]=999999999
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, levels=[1,8.5], colors=['darkgreen', 'deeppink'],linestyles=linestyle) 
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, SNR_sort_twodim, levels=[8], colors='blue', linewidths=1.5, linestyles=linestyle, antialiased=True)
        # ax.scatter(M_n_validity, D_l_n_validity, c='lime', edgecolor='black', s=80, zorder=5, alpha=0.9, label='Raw Data Points')
        ax.plot(M_n_validity, D_l_n_validity, 'o--', color=co, markersize=1, 
        linewidth=1, alpha=0.7, label='Connected Line')
    
    legend_elements = [Line2D([0], [0], 
                        color='blue', 
                        linestyle='-', 
                        label=r'$\mathrm{SNR}=8$') 
                    ]
    legend_elements.extend([
        Line2D([0], [0], 
            color='lime', 
            linestyle='-', 
            label=r'$\mathrm{Validity}~for~q=1.2$')
    ])
    legend_elements.extend([
        Line2D([0], [0], 
            color='red', 
            linestyle='-', 
            label=r'$\mathrm{Validity}~for~q=8.9$')
    ])
    legend_elements.extend([
        Line2D([0], [0], 
            color='deeppink', 
            linestyle='-', 
            label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=8.5~\mathrm{km}$')
    ])
    legend_elements.extend([
        Line2D([0], [0], 
            color='darkgreen', 
            linestyle='-', 
            label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=1.0~\mathrm{km}$')
    ])

    fig.legend(handles=legend_elements, 
                loc='upper right', 
                framealpha=0.8,
                bbox_to_anchor=(1.2, 0.6),
                )

    
    ax.set_title(f'{p_detector}', fontsize=14, pad=20)
    ax.set_xlabel(r'$M~~[\mathrm{M_{\odot}}]$', fontsize=12)
    ax.set_ylabel(r'$D_{L}~~[\mathrm{Mpc}]$', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if p_detector == 'TianQin' or p_detector == 'LISA' or p_detector == 'Taiji':
        ax.set_xlim(10, 1e3)
        ax.set_ylim(1, 1e4)
    else:
        ax.set_xlim(1, 1e3)
        ax.set_ylim(1, 1e6)
    ax.grid(True, which='major', linestyle=':', color='gray', alpha=1, linewidth=1)
    fig.tight_layout()
    fig.savefig(f'./data/{date}/fig/{p_detector}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('finished plot for', p_detector)
    return 0

def plot_two_simplify(p_detector, p_date):
    """绘制 TianQin 与一个地面探测器联合观测的简化约束图。"""
    
    p_path = f'./data/{p_date}/fig/'
    os.makedirs(p_path, exist_ok=True)
    p_l = 1001
    _, _, _, _, D_l_low_n, D_l_len_n, _, _, _, _, _ = get_detecror_parameters(p_detector)
    q = [1.2, 8.9]
    fig, ax = plt.subplots(figsize=(8, 6))
    for p_q in q:
        if p_q == 1.2:
            linestyle = '-'
            co = 'lime'
        else:
            linestyle = '--'
            co = 'red'
        data_path = f'./data/{p_date}/results/results for {p_detector}+TianQin with q={p_q}.txt'
        _, _, deltas, M_n, _, alpha_dCS_bar_square_max, SNR = read_results1(data_path)
        d_alpha_dCS_bar_square = deltas[7]**0.25 * 1e-3

        M_n_sort, d_alpha_dCS_bar_square_sort, alpha_dCS_bar_square_max_sort, SNR_sort = sort_multiple_arrays(M_n, d_alpha_dCS_bar_square, alpha_dCS_bar_square_max, SNR)
        M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, alpha_dCS_bar_square_max_sort_twodim, SNR_sort_twodim=np.zeros((5, p_l, p_l))

        k = np.logspace(0, 6, p_l)
        for i in range(p_l):
            for j in range(p_l):
                M_n_sort_twodim[i][j] = M_n_sort[i]
                D_l_n_sort_twodim[i][j] = k[j]
                d_alpha_dCS_bar_square_sort_twodim[i][j] = d_alpha_dCS_bar_square_sort[i] * k[j]**0.25
                alpha_dCS_bar_square_max_sort_twodim[i][j] = alpha_dCS_bar_square_max_sort[i]
                SNR_sort_twodim[i][j] = SNR_sort[i] / k[j]

        M_n_validity, D_l_n_validity = [[] for _ in range(2)]
        for i in range(p_l):
            for j in range(p_l):
                if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < 1e-2:
                    M_n_validity.append(M_n_sort_twodim[i][j])
                    D_l_n_validity.append(D_l_n_sort_twodim[i][j])

        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, levels=[1,8.5], colors=['darkgreen', 'deeppink'],linestyles=linestyle) 
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, SNR_sort_twodim, levels=[8], colors='blue', linewidths=1.5, linestyles=linestyle, antialiased=True)
        ax.plot(M_n_validity, D_l_n_validity, 'o--', color=co, markersize=1, 
        linewidth=1, alpha=0.7, label='Connected Line')
    
    legend_elements = [Line2D([0], [0], 
                        color='blue', 
                        linestyle='-', 
                        label=r'$\mathrm{SNR}=8$') 
                    ]
    legend_elements.extend([
        Line2D([0], [0], 
            color='lime', 
            linestyle='-', 
            label=r'$\mathrm{Validity}~for~q=1.2$')
    ])
    legend_elements.extend([
        Line2D([0], [0], 
            color='red', 
            linestyle='-', 
            label=r'$\mathrm{Validity}~for~q=8.9$')
    ])
    legend_elements.extend([
        Line2D([0], [0], 
            color='deeppink', 
            linestyle='-', 
            label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=8.5~\mathrm{km}$')
    ])
    legend_elements.extend([
        Line2D([0], [0], 
            color='darkgreen', 
            linestyle='-', 
            label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=1.0~\mathrm{km}$')
    ])


    fig.legend(handles=legend_elements, 
                loc='upper right', 
                framealpha=0.8,
                bbox_to_anchor=(1.2, 0.6),
                )
    # ax.add_artist(first_legend)


    ax.set_title(f'TianQin+{p_detector}', fontsize=14, pad=20)
    ax.set_xlabel(r'$M~~[\mathrm{M_{\odot}}]$', fontsize=12)
    ax.set_ylabel(r'$D_{L}~~[\mathrm{Mpc}]$', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 1e3)
    ax.set_ylim(1, 1e6)
    ax.grid(True, which='major', linestyle=':', color='gray', alpha=1, linewidth=1)
    fig.tight_layout()
    fig.savefig(f'./data/{p_date}/fig/TianQin+{p_detector}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('finished plot for TianQin+', p_detector)
    return 0

def plot_two(p_detector):
    """绘制旧版 TianQin+探测器二维质量-距离网格等高线图。"""
    p_path = './data/20250401'
    os.makedirs(p_path, exist_ok=True)

    _, _, _, _, D_l_low_n, D_l_len_n, _, _, _, _, _ = get_detecror_parameters(p_detector)
    q = [1.2, 8.9]
    fig, ax = plt.subplots(figsize=(8, 6))
    for p_q in q:
        if p_q == 1.2:
            linestyle = '-'
        else:
            linestyle = '--'
        data_path = f'./data/20250401/results for {p_detector}+TianQin with q={p_q}.txt'
        _, _, deltas, M_n, D_l_n, alpha_dCS_bar_square_max, SNR = read_results1(data_path)
        d_alpha_dCS_bar_square = pow(deltas[7], 0.25)

        sort_keys = (M_n.astype('float64'), D_l_n.astype('float64'))
        multi_sort_idx = np.lexsort(sort_keys)
        M_n_sort = reshape_array(M_n[multi_sort_idx])
        D_l_n_sort = reshape_array(D_l_n[multi_sort_idx])
        alpha_dCS_bar_square_max_sort = reshape_array(alpha_dCS_bar_square_max[multi_sort_idx])
        d_alpha_dCS_bar_square_sort = reshape_array(d_alpha_dCS_bar_square[multi_sort_idx])
        SNR_sort = reshape_array(SNR[multi_sort_idx])

        M_n_validity, D_l_n_validity = [[] for _ in range(2)]
        for i in range(101):
            for j in range(101):
                if abs(d_alpha_dCS_bar_square_sort[i][j] / alpha_dCS_bar_square_max_sort[i][j] - 1) < 1e-2:
                    M_n_validity.append(M_n_sort[i][j])
                    D_l_n_validity.append(D_l_n_sort[i][j])

        ax.contour(M_n_sort, D_l_n_sort, d_alpha_dCS_bar_square_sort, levels=[0.1,1,8.5], colors=['darkorange', 'darkgreen', 'deeppink'],linestyles=linestyle) 
        ax.contour(M_n_sort, D_l_n_sort, SNR_sort, levels=[8], colors='blue', linewidths=1.5, linestyles=linestyle, antialiased=True)
        ax.scatter(M_n_validity, D_l_n_validity, c='lime', edgecolor='black', s=80, zorder=5, alpha=0.9, label='Raw Data Points')
    
        legend_elements = [Line2D([0], [0], 
                                color='darkorange',
                                linestyle=linestyle,
                                label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=0.1~\mathrm{km}$') 
                        ]
        legend_elements.extend([
            Line2D([0], [0], 
                color='darkgreen', 
                linestyle=linestyle, 
                label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=1.0~\mathrm{km}$')
        ])
        legend_elements.extend([
            Line2D([0], [0], 
                color='deeppink', 
                linestyle=linestyle, 
                label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=8.5~\mathrm{km}$')
        ])
        legend_elements.extend([
            Line2D([0], [0], 
                color='blue', 
                linestyle=linestyle, 
                label=r'$\mathrm{SNR}=8$')
        ])
        legend_elements.extend([
            Line2D([0], [0], 
                color='red', 
                linestyle=linestyle, 
                label=r'$\mathrm{Validity}$')
        ])

        if p_q == 1.2:
            first_legend = fig.legend(handles=legend_elements, 
                    title=r'q=1.2',
                    loc='upper right', 
                    framealpha=0.8,
                    bbox_to_anchor=(1.1, 0.9),
                    )
            ax.add_artist(first_legend)
        else:
            fig.legend(handles=legend_elements, 
                    title=r'q=8.9',
                    loc='lower right', 
                    framealpha=0.8,
                    bbox_to_anchor=(0.95, 0.35),
                    )

    
    ax.set_title(f'{p_detector}+TianQin', fontsize=14, pad=20)
    ax.set_xlabel(r'$M~~[\mathrm{M_{\odot}}]$', fontsize=12)
    ax.set_ylabel(r'$D_{L}~~[\mathrm{Mpc}]$', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 1e3)
    ax.set_ylim(1, 1e6)
    ax.grid(True, which='major', linestyle=':', color='gray', alpha=1, linewidth=1)
    fig.tight_layout()
    fig.savefig(f'./data/20250401/{p_detector}+TianQin.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return 0

def plot_delta(p_detector, p_delta, p_q):
    """绘制不同差分步长下参数误差和偏导内积的诊断图。"""
    p_path = f'./data/20250401/{p_detector}/{p_q}/'
    os.makedirs(p_path, exist_ok=True)
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    for i in range(10):
        data_path = f'./data/20250401/q-{p_q}/{p_detector}/{i}/results for {p_delta}.txt'
        x, D, quad = read_results2(data_path)
        sort_keys = (x.astype('float64'))
        sort_idx = np.argsort(sort_keys)
        x_sort = x[sort_idx]
        y1_sort = D[sort_idx]
        y2_sort = quad[sort_idx]
        ax1.plot(x_sort, y1_sort, label=rf'M_{i}')
        ax2.plot(x_sort, y2_sort, label=rf'M_{i}')
    ax1.set_xlabel(rf'$d{p_delta}$', fontsize=12)
    ax1.set_ylabel(rf'$\frac{{\Delta {p_delta}}}{{{p_delta}}}$', fontsize=12)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    # ax1.set_xlim(1e-10, 1e-8)
    # ax1.set_ylim(2.8e16, 3e16)
    fig1.legend()
    fig1.tight_layout()
    fig1.savefig(f'./data/20250401/{p_detector}/{p_q}/results_{p_delta}1.png', dpi=300, bbox_inches='tight', facecolor='white')

    ax2.set_xlabel(rf'$d{p_delta}$', fontsize=12)
    ax2.set_ylabel(rf'$\left(\frac{{\partial h}}{{\partial {{{p_delta}}}}}|\frac{{\partial h}}{{\partial {{{p_delta}}}}}\right)$', fontsize=12)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    # ax2.set_xlim(1e-10, 1e-8)
    # ax2.set_ylim(2.8e16, 3e16)
    fig2.legend()
    fig2.tight_layout()
    fig2.savefig(f'./data/20250401/{p_detector}/{p_q}/results_{p_delta}2.png', dpi=300, bbox_inches='tight', facecolor='white')
    return 0

def plot_e(p_detector,p_date):
    """绘制偏心率与 dCS 约束量之间的关系曲线。"""
    dir_path = f'./data/{p_date}/fig/'
    file_path = os.path.join(dir_path,f'TianQin+{p_detector}.png')
    os.makedirs(dir_path, exist_ok=True)
    data_path = f'./data/{p_date}/results/results for {p_detector}+TianQin.txt'
    e, delta = read_results3(data_path)
    e, delta = sort_multiple_arrays(e, delta)
    
    plt.figure(figsize=(10, 6))
    plt.plot(e, delta, color='red')
    plt.xlabel('eccentricity')
    plt.ylabel(r'$\Delta \sqrt{\bar{\alpha}_{dCS}}~[\mathrm{km}]$')
    plt.title(f'TianQin + {p_detector}')
    plt.savefig(file_path,dpi=300,bbox_inches='tight', facecolor='white')


def plot_two_detectors_side_by_side(detectors, date):
    """把两个 TianQin+地面探测器联合观测子图并排保存。"""
    # 创建一个1行2列的大图
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # 为每个探测器调用绘图函数
    for ax, detector in zip(axes, detectors):
        plot_single_detector(ax, detector, date)
    
    # 调整整体布局
    plt.tight_layout()
    # 保存合并后的图像
    output_path = f'./data/{date}/fig/TianQin_combined.pdf'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)
def plot_single_detector(ax, p_detector, p_date):
    """在指定 Matplotlib 坐标轴上绘制一个联合探测器结果。"""
    # 获取探测器参数
    _, _, _, _, D_l_low_n, D_l_len_n, _, _, _, _, _ = get_detecror_parameters(p_detector)
    q = [1.2, 8.9]
    p_l = 1001
    for p_q in q:
        if p_q == 1.2:
            linestyle = '-'
        else:
            linestyle = '--'
            
        data_path = f'./data/{p_date}/results/results for {p_detector}+TianQin with q={p_q}.txt'
        _, _, deltas, M_n, _, alpha_dCS_bar_square_max, SNR = read_results1(data_path)
        d_alpha_dCS_bar_square = deltas[7]**0.25 * 1e-3
        M_n_sort, d_alpha_dCS_bar_square_sort, alpha_dCS_bar_square_max_sort, SNR_sort = sort_multiple_arrays(M_n, d_alpha_dCS_bar_square, alpha_dCS_bar_square_max, SNR)
        M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, alpha_dCS_bar_square_max_sort_twodim, SNR_sort_twodim = np.zeros((5, p_l, p_l))
        k = np.logspace(0, 6, p_l)
        for i in range(p_l):
            for j in range(p_l):
                M_n_sort_twodim[i][j] = M_n_sort[i]
                D_l_n_sort_twodim[i][j] = k[j]
                d_alpha_dCS_bar_square_sort_twodim[i][j] = d_alpha_dCS_bar_square_sort[i] * k[j]**0.25
                alpha_dCS_bar_square_max_sort_twodim[i][j] = alpha_dCS_bar_square_max_sort[i]
                SNR_sort_twodim[i][j] = SNR_sort[i] / k[j]
        M_n_validity, D_l_n_validity = [[] for _ in range(2)]

        for i in range(p_l):
                for j in range(p_l):
                    if p_q==1.2:
                        if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < 1e-2:
                            M_n_validity.append(M_n_sort_twodim[i][j])
                            D_l_n_validity.append(D_l_n_sort_twodim[i][j])
                    else:
                        if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < pow(10, -3.3):
                            M_n_validity.append(M_n_sort_twodim[i][j])
                            D_l_n_validity.append(D_l_n_sort_twodim[i][j])
        
        # if p_detector == 'CE':
        #     for i in range(p_l):
        #         for j in range(p_l):
        #             if p_q==1.2:
        #                 if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < 1e-2:
        #                     M_n_validity.append(M_n_sort_twodim[i][j])
        #                     D_l_n_validity.append(D_l_n_sort_twodim[i][j])
        #             else:
        #                 if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < pow(10, -3.3):
        #                     M_n_validity.append(M_n_sort_twodim[i][j])
        #                     D_l_n_validity.append(D_l_n_sort_twodim[i][j])
        #     # for lm in range(p_l):
        #     #     for ld in range(p_l):
        #     #         if M_n_sort_twodim[lm][ld]>=733:
        #     #             d_alpha_dCS_bar_square_sort_twodim[lm][ld]=999999999
        #     #             SNR_sort_twodim[lm][ld]=0
        
        # else:
        #     for i in range(p_l):
        #         for j in range(p_l):
        #             if p_q==1.2:
        #                 if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < 1e-2:
        #                     M_n_validity.append(M_n_sort_twodim[i][j])
        #                     D_l_n_validity.append(D_l_n_sort_twodim[i][j])
        #             else:
        #                 if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < pow(10, -3.3):
        #                     M_n_validity.append(M_n_sort_twodim[i][j])
        #                     D_l_n_validity.append(D_l_n_sort_twodim[i][j])
    
        
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, 
                  levels=[1,8.5], colors=['seagreen', 'darkorange'], linestyles=linestyle) 
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, SNR_sort_twodim, 
                  levels=[8], colors='blue', linewidths=1.5, linestyles=linestyle, antialiased=True)
        ax.plot(M_n_validity, D_l_n_validity, linestyle, color="black", markersize=1, 
               linewidth=1, alpha=0.7)
    
    size = 16
    size1 = 20
    size2 = 12
    if p_detector == 'ET':
        legend_elements = [
            Line2D([0], [0], color='blue', linestyle='-', label=r'$\mathrm{SNR}=8$'),
            Line2D([0], [0], color='black', linestyle='-', label=r'$\mathrm{Validity}$'),
            Line2D([0], [0], color='darkorange', linestyle='-', label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=8.5~\mathrm{km}$'),
            Line2D([0], [0], color='seagreen', linestyle='-', label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=1.0~\mathrm{km}$')
        ]
        ax.legend(handles=legend_elements, fontsize=size2, loc='upper right', framealpha=0.8,bbox_to_anchor=(1, 0.77))
    
    ax.set_title(f'TianQin+{p_detector}', fontsize=size1, pad=20)
    ax.set_xlabel(r'$M~~[\mathrm{M_{\odot}}]$', fontsize=size)
    ax.set_ylabel(r'$D_{L}~~[\mathrm{Mpc}]$', fontsize=size)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1, 1e3)
    ax.set_ylim(1, 1e6)
    ax.tick_params(axis='x', labelsize=size)
    ax.tick_params(axis='y', labelsize=size)
    ax.grid(True, which='major', linestyle=':', color='gray', alpha=1, linewidth=1)


def plot_all_detectors_grid(detectors, q, date):
    """把多个单探测器结果绘制到 4x2 网格并保存为 PDF。"""
    # 创建2行4列的子图
    fig, axes = plt.subplots(4, 2, figsize=(16, 24))
    
    # 展平axes数组便于遍历
    axes = axes.ravel()
    
    # 为每个探测器绘制图像
    for ax, detector in zip(axes, detectors):
        plot_single_detector_on_ax(ax, detector, q, date)
    
    # 调整整体布局
    plt.tight_layout()
    # 保存合并后的图像
    output_path = f'./data/{date}/fig/all_detectors.pdf'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_single_detector_on_ax(ax, p_detector, q, date):
    """在指定 Matplotlib 坐标轴上绘制一个单探测器结果。"""
    p_l = 1001
    
    for p_q in q:
        if p_q == 1.2:
            linestyle = '-'
        else:
            linestyle = '--'
            
        data_path = f'./data/{date}/results/results for {p_detector} with q={p_q}.txt'
        _, _, deltas, M_n, _, alpha_dCS_bar_square_max, SNR = read_results1(data_path)
        d_alpha_dCS_bar_square = deltas[7]**0.25 * 1e-3

        M_n_sort, d_alpha_dCS_bar_square_sort, alpha_dCS_bar_square_max_sort, SNR_sort = sort_multiple_arrays(
            M_n, d_alpha_dCS_bar_square, alpha_dCS_bar_square_max, SNR)
        
        M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, \
        alpha_dCS_bar_square_max_sort_twodim, SNR_sort_twodim = np.zeros((5, p_l, p_l))

        k = np.logspace(0, 6, p_l)
        for i in range(p_l):
            for j in range(p_l):
                M_n_sort_twodim[i][j] = M_n_sort[i]
                D_l_n_sort_twodim[i][j] = k[j]
                d_alpha_dCS_bar_square_sort_twodim[i][j] = d_alpha_dCS_bar_square_sort[i] * k[j]**0.25
                alpha_dCS_bar_square_max_sort_twodim[i][j] = alpha_dCS_bar_square_max_sort[i]
                SNR_sort_twodim[i][j] = SNR_sort[i] / k[j]

        M_n_validity, D_l_n_validity = [[] for _ in range(2)]
        for i in range(p_l):
            for j in range(p_l):
                if p_q==1.2:
                    if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < 1e-2:
                        M_n_validity.append(M_n_sort_twodim[i][j])
                        D_l_n_validity.append(D_l_n_sort_twodim[i][j])
                else:
                    if abs(d_alpha_dCS_bar_square_sort_twodim[i][j] / alpha_dCS_bar_square_max_sort_twodim[i][j] - 1) < pow(10, -3.3):
                        M_n_validity.append(M_n_sort_twodim[i][j])
                        D_l_n_validity.append(D_l_n_sort_twodim[i][j])
        
        # 特殊处理aLIGO和CE的质量上限
        if p_detector == 'aLIGO':
            for lm in range(p_l):
                for ld in range(p_l):
                    if M_n_sort_twodim[lm][ld] >= 600:
                        d_alpha_dCS_bar_square_sort_twodim[lm][ld] = 999999999
        elif p_detector == 'CE':
            for lm in range(p_l):
                for ld in range(p_l):
                    if M_n_sort_twodim[lm][ld] >= 700:
                        d_alpha_dCS_bar_square_sort_twodim[lm][ld] = 999999999

        # 绘制等高线
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, d_alpha_dCS_bar_square_sort_twodim, 
                  levels=[1,8.5], colors=['seagreen', 'darkorange'], linestyles=linestyle) 
        ax.contour(M_n_sort_twodim, D_l_n_sort_twodim, SNR_sort_twodim, 
                  levels=[8], colors='blue', linewidths=1.5, linestyles=linestyle, antialiased=True)
        ax.plot(M_n_validity, D_l_n_validity, linestyle, color="black", markersize=1, 
               linewidth=1, alpha=0.7)

    
    size = 16
    size1 = 20
    size2 = 12
    if p_detector == 'aLIGO':
        legend_elements = [
            Line2D([0], [0], color='blue', linestyle='-', label=r'$\mathrm{SNR}=8$'),
            Line2D([0], [0], color='black', linestyle='-', label=r'$\mathrm{Validity}$'),
            Line2D([0], [0], color='darkorange', linestyle='-', label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=8.5~\mathrm{km}$'),
            Line2D([0], [0], color='seagreen', linestyle='-', label=r'$\Delta\sqrt{\overline{\alpha}_{dCS}}=1.0~\mathrm{km}$')
        ]
        ax.legend(handles=legend_elements,fontsize=size2, loc='upper right', framealpha=0.8,bbox_to_anchor=(1, 1))

    # 设置子图标题和坐标轴
    ax.set_title(f'{p_detector}', fontsize=size1, pad=10)
    ax.set_xlabel(r'$M~~[\mathrm{M_{\odot}}]$', fontsize=size)
    ax.set_ylabel(r'$D_{L}~~[\mathrm{Mpc}]$', fontsize=size)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='x', labelsize=size)
    ax.tick_params(axis='y', labelsize=size)
    
    # 设置坐标轴范围
    if p_detector in ['TianQin', 'LISA', 'Taiji']:
        ax.set_xlim(10, 1e3)
        ax.set_ylim(1, 1e4)
    else:
        ax.set_xlim(1, 1e3)
        ax.set_ylim(1, 1e6)
    
    ax.grid(True, which='major', linestyle=':', color='gray', alpha=1, linewidth=0.5)
