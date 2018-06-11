from blimpy.utils import unpack_1to8, unpack_2to8, unpack_4to8, unpack
import numpy as np
import pytest

def test_1to8():
    a = np.array([0b01010101, 0b10101010], dtype=np.uint8)
    b = np.array([0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0])
    c = unpack_1to8(a)
    assert np.allclose(b, c)
    print(b)
    print(c)

def test_2to8():
    # Create an array that should come out as [0, 1, 2, 3, 3, 2, 1, 0]
    # In binary, this is [00, 01, 10, 11, 11, 10, 01, 00]
    # Convert to 8-bit, this is [0b00011011, 0b11100100]

    a = np.array([0b00011011, 0b11100100], dtype=np.uint8)
    b = np.array([0, 1, 2, 3, 3, 2, 1, 0], dtype=np.uint8)

    c = unpack_2to8(a)

    assert np.allclose(b, c)

def test_4to8():
    # Create an array that should come out as [0, 1, 2, 3, 3, 2, 1, 0]
    # In binary, this is [00, 01, 10, 11, 11, 10, 01, 00]
    # Convert to 8-bit, this is [0b00011011, 0b11100100]

    # Test 4-bit unpack
    a = np.array([0b00000001, 0b00100011], dtype=np.uint8)
    b = np.array([0, 1, 2, 3], dtype=np.uint8)
    c = unpack_4to8(a)
    assert np.allclose(b, c)

def test_unpack():

    # Test 2-bit unpack
    a = np.array([0b00011011, 0b11100100], dtype=np.uint8)
    b = np.array([0, 1, 2, 3, 3, 2, 1, 0], dtype=np.uint8)
    c = unpack(a, 2)
    assert np.allclose(b, c)

    # Catch exceptions
    with pytest.raises(ValueError):
        unpack(a, 16)  # nbit <= 8 is reqd
    with pytest.raises(ValueError):
        unpack(a, 3)  # nbit must divide 8 (1,2,4 or 8)
    z = np.array([1,2,3], dtype='float32')
    with pytest.raises(TypeError):
        unpack(z, 2) # input data must be 8-bit

    # Test 4-bit unpack
    a = np.array([0b00000001, 0b00100011], dtype=np.uint8)
    b = np.array([0, 1, 2, 3], dtype=np.uint8)
    c = unpack(a, 4)
    print(b)
    print(c)
    assert np.allclose(b, c)

    # Test 1-bit unpack
    a = np.array([0b01010101, 0b10101010], dtype=np.uint8)
    b = np.array([0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0])
    c = unpack(a, 1)
    print(b)
    print(c)
    assert np.allclose(b, c)

    # Test 8-bit!
    c = unpack(a, 8)
    assert np.allclose(a, c)

if __name__ == "__main__":
    test_1to8()
    test_2to8()
    test_4to8()
    test_unpack()