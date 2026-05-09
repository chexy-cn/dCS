"""绘制 T_tc 和 q 敏感性分析结果。

本模块读取 `tc_sensitivity.csv` 和 `q_sensitivity.csv`，生成适合论文
使用的约束变化图，以及用于检查低频截止和 SNR 的诊断图。
"""

import csv
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


DETECTOR_ORDER = ['TianQin', 'LISA', 'Taiji', 'DECIGO', 'BBO', 'aLIGO', 'ET', 'CE']
DETECTOR_STYLE = {
    'TianQin': {'color': 'royalblue', 'marker': 'o', 'linestyle': '-'},
    'LISA': {'color': 'limegreen', 'marker': 's', 'linestyle': '--'},
    'Taiji': {'color': 'black', 'marker': '^', 'linestyle': ':'},
    'DECIGO': {'color': 'seagreen', 'marker': 'D', 'linestyle': '-.'},
    'BBO': {'color': 'tomato', 'marker': 'v', 'linestyle': '-'},
    'aLIGO': {'color': 'brown', 'marker': 'P', 'linestyle': '--'},
    'ET': {'color': 'purple', 'marker': 'X', 'linestyle': ':'},
    'CE': {'color': 'gray', 'marker': '*', 'linestyle': '-.'},
}


def read_csv_rows(file_path):
    """读取敏感性 CSV，并把数值列转换为 float。"""
    with open(file_path, 'r') as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key, value in list(row.items()):
            if key != 'detector':
                row[key] = float(value)
    return rows


def group_by_detector(rows, sort_key):
    """按探测器名称组织行数据，并按自变量排序。"""
    grouped = {detector: [] for detector in DETECTOR_ORDER}
    for row in rows:
        grouped.setdefault(row['detector'], []).append(row)
    for detector, detector_rows in grouped.items():
        if detector_rows:
            detector_rows.sort(key=lambda item: item[sort_key])
    return grouped


def style_axes(ax):
    """应用与现有图一致的坐标轴样式。"""
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, which='major', linestyle=':', color='gray', alpha=0.8, linewidth=0.8)
    ax.grid(True, which='minor', linestyle=':', color='gray', alpha=0.25, linewidth=0.5)


def make_legend(ax):
    """在图外放置统一探测器图例。"""
    handles = [
        Line2D([0], [0], label=detector, linewidth=2.0, markersize=6, **DETECTOR_STYLE[detector])
        for detector in DETECTOR_ORDER
    ]
    ax.legend(
        handles=handles,
        fontsize=11,
        ncol=4,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.22),
        framealpha=0.9,
    )


def plot_constraint_sensitivity(date):
    """绘制 dCS 约束随 T_tc 和 q 变化的主图。"""
    tc_rows = read_csv_rows(f'./data/{date}/results/tc_sensitivity.csv')
    q_rows = read_csv_rows(f'./data/{date}/results/q_sensitivity.csv')
    tc_grouped = group_by_detector(tc_rows, 'T_tc_year')
    q_grouped = group_by_detector(q_rows, 'q')

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax_tc, ax_q = axes

    for detector in DETECTOR_ORDER:
        style = DETECTOR_STYLE[detector]
        rows = tc_grouped.get(detector, [])
        if rows:
            x = np.array([row['T_tc_year'] for row in rows])
            y = np.array([row['d_sqrt_alpha_dCS_km'] for row in rows])
            ax_tc.plot(x, y, linewidth=2.0, markersize=6, label=detector, **style)

        rows = q_grouped.get(detector, [])
        if rows:
            x = np.array([row['q'] for row in rows])
            y = np.array([row['d_sqrt_alpha_dCS_km'] for row in rows])
            ax_q.plot(x, y, linewidth=2.0, markersize=3.2, markevery=10, label=detector, **style)

    ax_tc.set_xlabel(r'$T_{\mathrm{obs}}~[\mathrm{yr}]$', fontsize=16)
    ax_tc.set_ylabel(r'$\Delta\sqrt{\overline{\alpha}_{\mathrm{dCS}}}~[\mathrm{km}]$', fontsize=16)
    ax_tc.set_xscale('log')
    ax_tc.set_yscale('log')
    ax_tc.set_xticks([0.25, 0.5, 1.0, 5.0])
    ax_tc.set_xticklabels(['0.25', '0.5', '1', '5'])
    ax_tc.set_title(r'Observation time', fontsize=18, pad=12)
    ax_tc.text(0.03, 0.94, r'$(a)$', transform=ax_tc.transAxes, fontsize=16)
    style_axes(ax_tc)

    ax_q.set_xlabel(r'$q=m_1/m_2$', fontsize=16)
    ax_q.set_ylabel(r'$\Delta\sqrt{\overline{\alpha}_{\mathrm{dCS}}}~[\mathrm{km}]$', fontsize=16)
    ax_q.set_yscale('log')
    ax_q.set_xlim(1.2, 8.9)
    ax_q.set_title(r'Mass ratio', fontsize=18, pad=12)
    ax_q.text(0.03, 0.94, r'$(b)$', transform=ax_q.transAxes, fontsize=16)
    style_axes(ax_q)

    make_legend(ax_q)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.26, wspace=0.28)

    output_dir = f'./data/{date}/fig/'
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, 'sensitivity_constraints.pdf'), dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(output_dir, 'sensitivity_constraints.png'), dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_diagnostics(date):
    """绘制低频截止和 SNR 的诊断图，辅助解释敏感性变化来源。"""
    tc_rows = read_csv_rows(f'./data/{date}/results/tc_sensitivity.csv')
    q_rows = read_csv_rows(f'./data/{date}/results/q_sensitivity.csv')
    tc_grouped = group_by_detector(tc_rows, 'T_tc_year')
    q_grouped = group_by_detector(q_rows, 'q')

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    ax_f, ax_snr = axes

    for detector in DETECTOR_ORDER:
        style = DETECTOR_STYLE[detector]
        rows = tc_grouped.get(detector, [])
        if rows:
            x = np.array([row['T_tc_year'] for row in rows])
            y = np.array([row['f_low_Hz'] for row in rows])
            ax_f.plot(x, y, linewidth=2.0, markersize=6, label=detector, **style)

        rows = q_grouped.get(detector, [])
        if rows:
            x = np.array([row['q'] for row in rows])
            y = np.array([row['SNR'] for row in rows])
            ax_snr.plot(x, y, linewidth=2.0, markersize=3.2, markevery=10, label=detector, **style)

    ax_f.set_xlabel(r'$T_{\mathrm{obs}}~[\mathrm{yr}]$', fontsize=16)
    ax_f.set_ylabel(r'$f_{\mathrm{low}}~[\mathrm{Hz}]$', fontsize=16)
    ax_f.set_xscale('log')
    ax_f.set_yscale('log')
    ax_f.set_xticks([0.25, 0.5, 1.0, 5.0])
    ax_f.set_xticklabels(['0.25', '0.5', '1', '5'])
    ax_f.set_title(r'Low-frequency cutoff', fontsize=18, pad=12)
    ax_f.text(0.03, 0.94, r'$(a)$', transform=ax_f.transAxes, fontsize=16)
    style_axes(ax_f)

    ax_snr.set_xlabel(r'$q=m_1/m_2$', fontsize=16)
    ax_snr.set_ylabel(r'$\mathrm{SNR}$', fontsize=16)
    ax_snr.set_yscale('log')
    ax_snr.set_xlim(1.2, 8.9)
    ax_snr.set_title(r'Signal-to-noise ratio', fontsize=18, pad=12)
    ax_snr.text(0.03, 0.94, r'$(b)$', transform=ax_snr.transAxes, fontsize=16)
    style_axes(ax_snr)

    make_legend(ax_snr)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.26, wspace=0.28)

    output_dir = f'./data/{date}/fig/'
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, 'sensitivity_diagnostics.pdf'), dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(output_dir, 'sensitivity_diagnostics.png'), dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_all(date):
    """生成所有敏感性分析图。"""
    plot_constraint_sensitivity(date)
    plot_diagnostics(date)
