from blimpy import utils
import numpy as np
import pytest


def test_utils():
    assert utils.db(100) == 20.0
    assert utils.lin(20)  == 100.0
    assert utils.closest(np.array([0,1,2,3,4,5]), 2.2) == 2

def test_rebin():
    # 1D
    a = np.array([1, 1, 1, 1])
    aR = utils.rebin(a, 2)
    assert np.allclose(aR, np.array([1, 1]))

    # 2D
    b = np.array([[1,1,1,1], [2,2,2,2]])
    bR = utils.rebin(b, 1, 2)
    assert np.allclose(bR, [[1,1], [2,2]])
    bR = utils.rebin(b, None, 2)
    assert np.allclose(bR, [[1,1], [2,2]])
    bR = utils.rebin(b, 2, 1)
    assert np.allclose(bR, [1.5, 1.5, 1.5, 1.5])
    bR = utils.rebin(b, 2, None)
    assert np.allclose(bR, [1.5, 1.5, 1.5, 1.5])

    c = np.zeros([10, 10, 10])
    with pytest.raises(RuntimeError):
        utils.rebin(c, 2, 2)


if __name__ == "__main__":
    test_utils()
    test_rebin()