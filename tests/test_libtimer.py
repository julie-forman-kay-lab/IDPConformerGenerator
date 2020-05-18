from idpconfgen.libs import libtimer

import pytest


def test_ProgressBar_1():
    """Test progress bar."""
    kw = {'prefix': 'hello', 'suffix': 'frames', 'decimals': 2}
    with libtimer.ProgressBar(40, **kw) as PB:
        assert PB.counter == 1
        for i in range(40):
            PB.increment()
            assert PB.counter == i + 2

    with pytest.raises(IndexError):
        PB.increment()


def test_ProgressBar_2():
    """Test progress bar."""
    kw = {
        'prefix': 'hello',
        'suffix': 'frames',
        'decimals': 2,
        'bar_length': 70,
        }
    with libtimer.ProgressBar(30, **kw) as PB:
        assert PB.counter == 1
        for i in range(30):
            PB.increment()
            assert PB.counter == i + 2

    with pytest.raises(IndexError):
        PB.increment()
