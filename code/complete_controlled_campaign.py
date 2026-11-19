"""Resume the frozen L1024 producer with four scheduling workers.

Only executor concurrency changes; the numerical producer, run identifiers,
seeds, parameters and per-run cache checks remain unchanged. This launcher is
separately archived as execution provenance.
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path

import controlled_roughening as controls

ROOT = Path(__file__).resolve().parents[1]


if __name__ == '__main__':
    directory = ROOT/'data/raw/controlled_20261117/baseline'
    metadata = dict(workers=4, numerical_producer_unchanged=True,
                    scheduling_only=True, seed=20261117, runs=64,
                    launcher_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    producer_sha256=hashlib.sha256(Path(controls.__file__).read_bytes()).hexdigest())
    (directory/'execution_L1024.json').write_text(json.dumps(metadata, indent=2)+'\n')
    controls.ThreadPoolExecutor = lambda max_workers: ThreadPoolExecutor(max_workers=4)
    controls.run_cell(.25, 1024, 64, 491520, 20261117, 2, directory)
