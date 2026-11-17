"""绘制 T_tc 和 q 敏感性分析结果。

本模块读取 `tc_sensitivity.csv` 和 `q_sensitivity.csv`，生成适合论文
使用的约束变化图，以及用于检查低频截止和 SNR 的诊断图。
"""

import csv
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
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
T_OBS_TICKS = [1 / 365.25, 7 / 365.25, 1 / 12, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
T_OBS_LABELS = ['1 d', '1 wk', '1 mo', '3 mo', '0.5', '1', '2', '3', '4', '5']
GROUND_DETECTORS = {'aLIGO', 'ET', 'CE'}
SPACE_LIKE_DETECTORS = {'TianQin', 'LISA', 'Taiji', 'DECIGO', 'BBO'}


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


def spin_contours(ax, x, y, grid, detector):
    """Label only levels present in this panel, including 50 km where requested."""
    finite = grid[np.isfinite(grid) & (grid > 0)]
    if len(finite) < 4:
        return
    low, high = finite.min(), finite.max()
    levels = list(10. ** np.arange(np.ceil(np.log10(low)), np.floor(np.log10(high)) + 1))
    if detector in {'TianQin', 'LISA', 'Taiji', 'aLIGO'} and low < 50 < high:
        levels.append(50.)
    levels = sorted({level for level in levels if low < level < high})
    if levels:
        contour = ax.contour(x, y, np.ma.masked_invalid(grid.T), levels=levels,
                             colors='white', linewidths=1.1, alpha=.9)
        ax.clabel(contour, fmt=lambda value: f'{value:g}', fontsize=7,
                  colors='black', inline=True, inline_spacing=2)


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


def make_combined_legends(detector_ax, flow_ax):
    """为合并图添加探测器图例和低频截止分组图例。"""
    detector_handles = [
        Line2D([0], [0], label=detector, linewidth=2.0, markersize=6, **DETECTOR_STYLE[detector])
        for detector in DETECTOR_ORDER
    ]
    detector_ax.legend(
        handles=detector_handles,
        fontsize=9,
        ncol=2,
        loc='upper right',
        framealpha=0.9,
    )

    quantity_handles = [
        Line2D([0], [0], color='dimgray', linestyle=':', linewidth=2.2,
               label=r'$f_{\mathrm{low}}$: ground-based'),
        Line2D([0], [0], color='saddlebrown', linestyle=':', linewidth=2.2,
               label=r'$f_{\mathrm{low}}$: space/BBO/DECIGO'),
    ]
    flow_ax.legend(
        handles=quantity_handles,
        fontsize=9,
        ncol=1,
        loc='upper right',
        framealpha=0.9,
    )


def plot_observation_time_flow_groups(ax, tc_grouped):
    """在观测时间 panel 的右轴上绘制两类低频截止线。"""
    space_rows = next((tc_grouped[detector] for detector in DETECTOR_ORDER
                       if detector in SPACE_LIKE_DETECTORS and tc_grouped.get(detector)), [])
    ground_rows = next((tc_grouped[detector] for detector in DETECTOR_ORDER
                        if detector in GROUND_DETECTORS and tc_grouped.get(detector)), [])

    if space_rows:
        x = np.array([row['T_tc_year'] for row in space_rows])
        y = np.array([row['f_low_Hz'] for row in space_rows])
        ax.plot(x, y, color='saddlebrown', linestyle=':', linewidth=2.8, alpha=0.85)

    if ground_rows:
        x = np.array([row['T_tc_year'] for row in ground_rows])
        y = np.array([row['f_low_Hz'] for row in ground_rows])
        ax.plot(x, y, color='dimgray', linestyle=':', linewidth=2.8, alpha=0.85)


def plot_combined_sensitivity(date, results_root="results", data_dir=None):
    """绘制一张包含约束和诊断量的双纵轴敏感性图。"""
    tc_rows = read_csv_rows(f'{data_dir or str(results_root) + "/" + str(date) + "/data"}/tc_sensitivity.csv')
    q_rows = read_csv_rows(f'{data_dir or str(results_root) + "/" + str(date) + "/data"}/q_sensitivity.csv')
    tc_grouped = group_by_detector(tc_rows, 'T_tc_year')
    q_grouped = group_by_detector(q_rows, 'q')

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.4))
    ax_tc, ax_q = axes
    ax_tc_right = ax_tc.twinx()
    ax_q_right = ax_q.twinx()

    for detector in DETECTOR_ORDER:
        style = DETECTOR_STYLE[detector]
        color = style['color']

        rows = tc_grouped.get(detector, [])
        if rows:
            x = np.array([row['T_tc_year'] for row in rows])
            y_constraint = np.array([row['d_sqrt_alpha_dCS_km'] for row in rows])
            ax_tc.plot(x, y_constraint, linewidth=2.1, markersize=5.6, label=detector, **style)

        rows = q_grouped.get(detector, [])
        if rows:
            x = np.array([row['q'] for row in rows])
            y_constraint = np.array([row['d_sqrt_alpha_dCS_km'] for row in rows])
            y_snr = np.array([row['SNR'] for row in rows])
            ax_q.plot(x, y_constraint, linewidth=2.1, markersize=3.0, markevery=10, label=detector, **style)
            ax_q_right.plot(x, y_snr, color=color, linestyle=':', linewidth=2.0, alpha=0.32)

    plot_observation_time_flow_groups(ax_tc_right, tc_grouped)

    ax_tc.set_xlabel(r'$T_{\mathrm{obs}}~[\mathrm{yr}]$', fontsize=16)
    ax_tc.set_ylabel(r'$\Delta\sqrt{\overline{\alpha}_{\mathrm{dCS}}}~[\mathrm{km}]$', fontsize=16)
    ax_tc_right.set_ylabel(r'$f_{\mathrm{low}}~[\mathrm{Hz}]$', fontsize=14, color='dimgray', labelpad=8)
    ax_tc.set_xscale('log')
    ax_tc.set_yscale('log')
    ax_tc_right.set_yscale('log')
    ax_tc.set_xticks(T_OBS_TICKS)
    ax_tc.set_xticklabels(T_OBS_LABELS, rotation=30, ha='right')
    ax_tc.set_title(r'Observation time', fontsize=18, pad=12)
    ax_tc.text(0.03, 0.94, r'$(a)$', transform=ax_tc.transAxes, fontsize=16)
    style_axes(ax_tc)
    ax_tc_right.tick_params(axis='y', labelsize=14, colors='dimgray')

    ax_q.set_xlabel(r'$q=m_1/m_2$', fontsize=16)
    ax_q.set_ylabel(r'$\Delta\sqrt{\overline{\alpha}_{\mathrm{dCS}}}~[\mathrm{km}]$', fontsize=16)
    ax_q_right.set_ylabel(r'$\mathrm{SNR}$', fontsize=14, color='dimgray', labelpad=8)
    ax_q.set_yscale('log')
    ax_q_right.set_yscale('log')
    ax_q.set_xlim(1.2, 8.9)
    ax_q.set_title(r'Mass ratio', fontsize=18, pad=12)
    ax_q.text(0.03, 0.94, r'$(b)$', transform=ax_q.transAxes, fontsize=16)
    style_axes(ax_q)
    ax_q_right.tick_params(axis='y', labelsize=14, colors='dimgray')

    make_combined_legends(ax_q, ax_tc)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.18, wspace=0.52)

    output_dir = f'{results_root}/{date}/fig/'
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, 'sensitivity_combined.pdf'), dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(output_dir, 'sensitivity_combined.png'), dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)








def plot_spin_component_sensitivity(date, results_root="results", data_dir=None):
    """绘制每个探测器的分量自旋 x1-x2 约束热图。"""
    component_path = f'{data_dir or str(results_root) + "/" + str(date) + "/data"}/spin_component_sensitivity.csv'
    rows = read_csv_rows(component_path)
    x1_values = np.unique(np.round(np.array([row['x1'] for row in rows]), 12))
    x2_values = np.unique(np.round(np.array([row['x2'] for row in rows]), 12))
    finite_constraints = np.array([
        row['d_sqrt_alpha_dCS_km'] for row in rows
        if np.isfinite(row['d_sqrt_alpha_dCS_km']) and row['d_sqrt_alpha_dCS_km'] > 0
    ])
    if not len(finite_constraints):
        raise ValueError('Spin sensitivity CSV contains no positive finite constraints.')
    norm = LogNorm(vmin=finite_constraints.min(), vmax=finite_constraints.max())

    fig, axes = plt.subplots(2, 4, figsize=(14, 7.2), sharex=True, sharey=True, constrained_layout=True)
    image = None
    for ax, detector in zip(axes.flat, DETECTOR_ORDER):
        grid = np.full((len(x1_values), len(x2_values)), np.nan)
        for row in rows:
            if row['detector'] != detector:
                continue
            i = np.searchsorted(x1_values, round(row['x1'], 12))
            j = np.searchsorted(x2_values, round(row['x2'], 12))
            value = row['d_sqrt_alpha_dCS_km']
            if np.isfinite(value) and value > 0:
                grid[i, j] = value
        ax.set_facecolor('#eeeeee')
        image = ax.pcolormesh(x1_values, x2_values, grid.T, shading='nearest', cmap='viridis_r', norm=norm)
        spin_contours(ax, x1_values, x2_values, grid, detector)
        ax.add_patch(Rectangle((-0.045, -0.045), 0.09, 0.09, fill=False, hatch='///', edgecolor='black', linewidth=0.7, zorder=3))
        ax.plot([-0.9, 0.9], [-0.9, 0.9], color='white', linestyle=':', linewidth=0.8, alpha=0.55)
        ax.plot(0.4, 0.2, marker='*', markersize=9, markerfacecolor='white', markeredgecolor='black', markeredgewidth=0.9, linestyle='None', zorder=6)
        ax.set_title(detector, fontsize=14)
        ax.set_aspect('equal')
        style_axes(ax)

    for ax in axes[-1, :]:
        ax.set_xlabel(r'$x_1$', fontsize=13)
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$x_2$', fontsize=13)
    colorbar = fig.colorbar(image, ax=axes, shrink=0.92, pad=0.02)
    colorbar.set_label(r'$\Delta\sqrt{\overline{\alpha}_{\mathrm{dCS}}}~[\mathrm{km}]$', fontsize=14)

    output_dir = f'{results_root}/{date}/fig/'
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, 'spin_component_sensitivity.pdf'), dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(output_dir, 'spin_component_sensitivity.png'), dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_inclination_sensitivity(date, results_root="results", data_dir=None):
    """绘制倾角下的 dCS 约束与 SNR 曲线。"""
    rows = read_csv_rows(f'{data_dir or str(results_root) + "/" + str(date) + "/data"}/inclination_sensitivity.csv')
    grouped = group_by_detector(rows, 'iota')
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax_snr = ax.twinx()

    for detector in DETECTOR_ORDER:
        detector_rows = grouped.get(detector, [])
        if not detector_rows:
            continue
        style = DETECTOR_STYLE[detector]
        x = np.array([row['iota'] for row in detector_rows])
        constraint = np.array([row['d_sqrt_alpha_dCS_km'] for row in detector_rows])
        snr_values = np.array([row['SNR'] for row in detector_rows])
        ax.plot(x, constraint, linewidth=2.1, markersize=4.2, markevery=4, label=detector, **style)
        ax_snr.plot(x, snr_values, color=style['color'], linestyle=':', linewidth=1.8, alpha=0.32)

    ax.set_xlabel(r'$\iota~[\mathrm{rad}]$', fontsize=16)
    ax.set_ylabel(r'$\Delta\sqrt{\overline{\alpha}_{\mathrm{dCS}}}~[\mathrm{km}]$', fontsize=16)
    ax_snr.set_ylabel(r'$\mathrm{SNR}$', fontsize=15, color='dimgray')
    ax.set_xlim(0, np.pi / 2)
    ax.axvline(0.9, color='black', linestyle='--', linewidth=1.25, zorder=4)
    ax.text(0.9, 0.04, r'$\iota=0.9$', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=11, color='black')
    ax.set_xticks([0, np.pi / 4, np.pi / 2])
    ax.set_xticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$'])
    ax.set_yscale('log')
    ax_snr.set_yscale('log')
    ax.set_title('Inclination', fontsize=18, pad=12)
    style_axes(ax)
    ax_snr.tick_params(axis='y', labelsize=14, colors='dimgray')
    make_legend(ax)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.27, right=0.86)

    output_dir = f'{results_root}/{date}/fig/'
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, 'inclination_sensitivity.pdf'), dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(output_dir, 'inclination_sensitivity.png'), dpi=600, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_all(date, results_root="results", data_dir=None):
    """生成当前日期下已有输入数据对应的敏感性图。"""
    output_dir = f'{data_dir or str(results_root) + "/" + str(date) + "/data"}/'
    tc_path = os.path.join(output_dir, 'tc_sensitivity.csv')
    q_path = os.path.join(output_dir, 'q_sensitivity.csv')
    if os.path.exists(tc_path) and os.path.exists(q_path):
        plot_combined_sensitivity(date, results_root, data_dir)
    if os.path.exists(os.path.join(output_dir, 'spin_component_sensitivity.csv')):
        plot_spin_component_sensitivity(date, results_root, data_dir)
    if os.path.exists(os.path.join(output_dir, 'inclination_sensitivity.csv')):
        plot_inclination_sensitivity(date, results_root, data_dir)
