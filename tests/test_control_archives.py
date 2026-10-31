import numpy as np
import pytest

from control_archives import atomic_savez


def test_atomic_archive_is_loadable_and_failure_keeps_previous_file(tmp_path, monkeypatch):
    path = tmp_path/'cell.npz'
    atomic_savez(path, value=np.arange(10))
    previous = path.read_bytes()

    def fail(stream, **arrays):
        stream.write(b'incomplete zip')
        assert path.read_bytes() == previous
        raise OSError('simulated write interruption')

    monkeypatch.setattr(np, 'savez_compressed', fail)
    with pytest.raises(OSError, match='write interruption'):
        atomic_savez(path, value=np.arange(20))
    assert path.read_bytes() == previous
    assert list(tmp_path.iterdir()) == [path]
    with np.load(path, allow_pickle=False) as z:
        np.testing.assert_array_equal(z['value'], np.arange(10))
