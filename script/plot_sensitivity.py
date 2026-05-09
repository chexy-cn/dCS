"""绘制 T_tc 与 q 敏感性分析结果图。"""

from pathlib import Path
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


def main():
    """生成论文主图和诊断图。"""
    date = 20260506
    plot_sensitivity.plot_all(date)


if __name__ == '__main__':
    main()
