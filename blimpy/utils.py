"""
# utils.py
useful helper functions for common data manipulation tasks
"""
import numpy as np


def db(x):
    """ Convert linear to dB """
    return 10 * np.log10(x)


def lin(x):
    """ Convert dB to linear """
    return 10.0 ** (x / 10.0)


def closest(xarr, val):
    """ Return the index of the closest in xarr to value val """
    idx_closest = np.argmin(np.abs(np.array(xarr) - val))
    return idx_closest


def rebin(d, n_x, n_y=None):
    """ Rebin data by averaging bins together

    Args:
    d (np.array): data
    n_x (int): number of bins in x dir to rebin into one
    n_y (int): number of bins in y dir to rebin into one

    Returns:
    d: rebinned data with shape (n_x, n_y)
    """

    if d.ndim == 2:
        if n_y is None:
            n_y = 1
        if n_x is None:
            n_x = 1
        d = d[:int(d.shape[0] // n_x) * n_x, :int(d.shape[1] // n_y) * n_y]
        d = d.reshape((d.shape[0] // n_x, n_x, d.shape[1] // n_y, n_y))
        d = d.mean(axis=3)
        d = d.mean(axis=1)
    elif d.ndim == 1:
        d = d[:int(d.shape[0] // n_x) * n_x]
        d = d.reshape((d.shape[0] // n_x, n_x))
        d = d.mean(axis=1)
    else:
        raise RuntimeError("Only NDIM <= 2 supported")
    return d


def unpack(data, nbit):
    """upgrade data from nbits to np.int8

    Notes: nbit=4 and nbit=1 are untested!
    """
    if data.dtype not in (np.uint8, np.int8):
        raise TypeError("unpack: dtype must be 8-bit")
    if nbit == 8:
        return data.view(np.int8)
    elif nbit == 4:
        data = unpack_4to8(data)
        return data
    elif nbit == 2:
        data = unpack_2to8(data)
        return data
    elif nbit == 1:
        data = unpack_1to8(data)
        return data
    else:
        raise ValueError("unpack: nbit must be 1, 2, 4 or 8")

def unpack_complex64(data, nbit):
    """upgrade data from nbits to np.complex64

    Notes: nbit=4 and nbit=1 are untested!
    """
    data = unpack(data,nbit)
    return data.astype(np.float32).view(np.complex64)

def unpack_1to8(data):
    """ Promote 1-bit unsigned data into np.int8

    Args:
        data: Numpy array with dtype == uint8
    """
    return np.unpackbits(data.view(np.uint8)).astype(np.int8)


def unpack_2to8(data):
    """ Promote 2-bit unsigned data into np.int8

    Args:
        data: Numpy array with dtype == uint8

    Notes:
        DATA MUST BE LOADED as np.array() with dtype='uint8'.

        This works with some clever shifting and AND / OR operations.
        Data is LOADED as 8-bit, then promoted to 32-bits:
        /ABCD EFGH/ (8 bits of data)
        /0000 0000/0000 0000/0000 0000/ABCD EFGH/ (8 bits of data as a 32-bit word)

        Once promoted, we can do some shifting, AND and OR operations:
        /0000 0000/0000 ABCD/EFGH 0000/0000 0000/ (shifted << 12)
        /0000 0000/0000 ABCD/EFGH 0000/ABCD EFGH/ (bitwise OR of previous two lines)
        /0000 0000/0000 ABCD/0000 0000/0000 EFGH/ (bitwise AND with mask 0xF000F)
        /0000 00AB/CD00 0000/0000 00EF/GH00 0000/ (prev. line shifted << 6)
        /0000 00AB/CD00 ABCD/0000 00EF/GH00 EFGH/ (bitwise OR of previous two lines)
        /0000 00AB/0000 00CD/0000 00EF/0000 00GH/ (bitwise AND with 0x3030303)

        Then we change the view of the data to interpret it as 4x8 bit:
        [000000AB, 000000CD, 000000EF, 000000GH]  (change view from 32-bit to 4x8-bit)

        The converted bits are then mapped to values in the range [-40, 40] according to a lookup chart.
        The mapping is based on specifications in the breakthough docs:
        https://github.com/UCBerkeleySETI/breakthrough/blob/master/doc/RAW-File-Format.md

    """
    two_eight_lookup = {0: 40,
                        1: 12,
                        2: -12,
                        3: -40}

    tmp = data.view(np.uint8).astype(np.uint32)
    tmp = (tmp | (tmp << 12)) & 0xF000F
    tmp = (tmp | (tmp << 6)) & 0x3030303
    tmp = tmp.byteswap()
    tmp = tmp.view('uint8')
    mapped = np.array(tmp, dtype=np.int8)
    for k, v in two_eight_lookup.items():
        mapped[tmp == k] = v
    return mapped


def unpack_4to8(data):
    """ Promote 4-bit unsigned data into np.int8

    Args:
        data: Numpy array with dtype == uint8

    Notes:
        # The process is this:
        # ABCDEFGH [Bits of one 4+4-bit value]
        # 00000000ABCDEFGH [astype(uint16)]
        # 0000ABCDEFGH0000 [<< 4]
        # 0000ABCDXXXXEFGH [bitwise 'or' of previous two lines]
        # 0000111100001111 [0x0F0F]
        # 0000ABCD0000EFGH [bitwise 'and' of previous two lines]
        # ABCD0000EFGH0000 [<< 4]
        # which effectively pads the two 4-bit values with zeros on the right
        # Note: This technique assumes LSB-first ordering
    """

    tmpdata = data.view(uint8).astype(np.int16)  # np.empty(upshape, dtype=np.int16)
    tmpdata = (tmpdata | (tmpdata << 4)) & 0x0F0F
    # tmpdata = tmpdata << 4 # Shift into high bits to avoid needing to sign extend
    updata = tmpdata.byteswap()
    return updata.view(np.int8)
