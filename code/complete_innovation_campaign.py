"""Finish the already-running final campaign and produce its evidence report.

The three L512 source jobs run separately. This job waits for complete ZIP
archives, generates their binary controls in parallel, then verifies the full
campaign and exports results. Existing validated L32--256 binary controls are
preserved with their original source metadata.
"""
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import zipfile

import numpy as np

from innovation_controls import ROOT
from binary_covariance_control import run_binary_cell
from innovation_analysis import analyze

SOURCE = ROOT/'data/raw/innovation_20261031_v3'
BINARY = ROOT/'data/raw/innovation_binary_20261031_v3'
OUT = ROOT/'docs/innovation_20261031/results'


def wait_complete(path):
    while True:
        if path.exists():
            try:
                with zipfile.ZipFile(path) as archive:
                    if archive.testzip() is None:
                        return
            except (OSError, zipfile.BadZipFile):
                pass
        time.sleep(10)


def complete_binary(theta):
    path = SOURCE/f'theta{theta:.2f}_L512.npz'
    wait_complete(path)
    print(f'Source ready: {path.name}', flush=True)
    run_binary_cell(path, BINARY, workers=4)


def main():
    BINARY.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(complete_binary, [.1, .25, .4]))
    for size in [32, 64, 128, 256, 512]:
        for theta in [.1, .25, .4]:
            for folder in [SOURCE, BINARY]:
                wait_complete(folder/f'theta{theta:.2f}_L{size}.npz')
    analyze(SOURCE, BINARY, OUT)
    result = json.loads((OUT/'summary.json').read_text(encoding='utf-8'))
    assert result['status'] == 'COMPLETE' and len(result['cells']) == 15
    assert all(len(c['widths']) == 6 and np.all(np.asarray(c['widths']) > 0)
               for c in result['cells'])
    gates = [c['binary_covariance_gate']['passed'] for c in result['cells']]
    failed = ', '.join('θ={}, L={}'.format(c['theta'], c['L'])
                       for c in result['cells'] if not c['binary_covariance_gate']['passed'])
    lines = ['# 完成状态与物理结论边界', '',
             '全量数据已完成：3个θ×5个尺寸×7种输入×64次界面模拟，共6720次；每次保留T和2T两个观测时长。', '',
             f'独立二值协方差检验通过{sum(gates)}/15个单元。未通过单元：'+failed+'。', '',
             '具体有限尺寸斜率差、宽度比、四点统计及95%置信区间见report_zh.md。界面差异不能在匹配检验失败的单元上归因于高阶噪声统计。', '',
             '这些实验把原始噪声匹配、实际二值驱动匹配和创新结构的影响分开；是否支持新的物理结论应以完整表格及其误差范围为准。有限窗口结果不构成渐近指数、普适类或临界点结论。']
    (OUT/'completion_note_zh.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    command = [sys.executable, '-m', 'pytest', '-q',
               '--basetemp=docs/innovation_20261031/pytest-completed',
               '--junitxml=docs/innovation_20261031/tests_completed.xml']
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    (OUT/'test_output.txt').write_text(completed.stdout+completed.stderr, encoding='utf-8')
    if completed.returncode:
        raise RuntimeError('Final regression checks failed; see test_output.txt')
    outputs = [p for p in OUT.iterdir() if p.is_file()]
    record = dict(status='COMPLETE',
                  completed_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
                  covariance_gates_passed=sum(gates), covariance_gates_total=15,
                  simulations=6720, outputs=[dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in outputs])
    (OUT/'completion.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    print('COMPLETE: full report, figures, provenance and regression results saved.', flush=True)


if __name__ == '__main__':
    main()
