from blimpy.utils import unpack_2to8
import numpy as np


def test_2to8():
    # Create an array that should come out as [0, 1, 2, 3, 3, 2, 1, 0]
    # In binary, this is [00, 01, 10, 11, 11, 10, 01, 00]
    # Convert to 8-bit, this is [0b00011011, 0b11100100]

    a = np.array([0b00011011, 0b11100100], dtype=np.uint8)
    b = np.array([0, 1, 2, 3, 3, 2, 1, 0], dtype=np.uint8)

    c = unpack_2to8(a)

    assert np.allclose(b, c)

if __name__ == "__main__":
    test_2to8()