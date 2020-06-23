from idpconfgen.libs import libtimer

import pytest


#@pytest.mark.parametrize(
#    'items, cls',
#    [
#        (list(range(10)), libtimer.ProgressBar),
#        (range(10), libtimer.ProgressCounter),
#        ]
#    )
def test_ProgressWatcher():
    """
    Test ProgressWatcher constructor.

    Running tests with `tox` the os.get_terminal_size
    raises IOError, so only ProgressFake is expected.
    """
    result = libtimer.ProgressWatcher(range(10))
    assert isinstance(result, libtimer.ProgressFake)


def test_ProgressBar_1():
    """Test progress bar."""
    kw = {'prefix': 'hello', 'suffix': 'frames', 'decimals': 2}
    with libtimer.ProgressBar(40, **kw) as PB:
        assert PB.counter == 1
        for i in range(40):
            PB.increment()
            assert PB.counter == i + 2


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

def test_ProgressFake():
    """Test ProgressFake interface."""
    with libtimer.ProgressFake() as PF:
        PF.increment()

def test_ProgressCounter():
    """Test ProgressCounter."""
    with libtimer.ProgressCounter() as PC:
        assert PC.counter == 0
        for i in range(10):
            PC.increment()
            assert PC.counter == i + 1
