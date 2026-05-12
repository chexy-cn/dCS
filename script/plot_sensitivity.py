"""绘制 T_tc 与 q 敏感性分析结果图。"""

from pathlib import Path
import argparse
from datetime import date
import os
import sys

_PROJECT_ROOT = next(
    path for path in Path(__file__).resolve().parents
    if (path / 'getresults').is_dir() and (path / 'math' / 'fim.py').is_file()
)
for _path in (str(_PROJECT_ROOT), str(_PROJECT_ROOT / 'math')):
    if _path not in sys.path:
        sys.path.insert(0, _path)
os.chdir(_PROJECT_ROOT)

from getresults import plot_sensitivity


def current_date_tag():
    """返回当前日期标签，格式为 YYYYMMDD。"""
    return date.today().strftime('%Y%m%d')


def parse_args():
    """解析绘图日期参数。"""
    parser = argparse.ArgumentParser(description='Plot T_tc and q sensitivity results.')
    parser.add_argument(
        '--date',
        default=current_date_tag(),
        help='Result date directory under data/, e.g. 20260506. Defaults to today.',
    )
    return parser.parse_args()


def main():
    """生成论文主图和诊断图。"""
    args = parse_args()
    plot_sensitivity.plot_all(args.date)


if __name__ == '__main__':
    main()
