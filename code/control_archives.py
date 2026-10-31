"""Publish a complete control archive before concurrent readers can see it."""
from pathlib import Path
import os
import tempfile

import numpy as np


def atomic_savez(path, **arrays):
    path = Path(path)
    stream = tempfile.NamedTemporaryFile(dir=path.parent, prefix=f'.{path.stem}-',
                                         suffix='.tmp', delete=False)
    temporary = Path(stream.name)
    try:
        with stream:
            np.savez_compressed(stream, **arrays)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
