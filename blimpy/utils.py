"""
# utils.py
useful helper functions for common data manipulation tasks
"""
import numpy as np

def db(x):
    """ Convert linear to dB """
    return 10*np.log10(x)

def lin(x):
    """ Convert dB to linear """
    return 10.0**(x / 10.0)

def closest(xarr, val):
    """ Return the index of the closest in xarr to value val """
    idx_closest = np.argmin(np.abs(xarr - val))
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
    """upgrade data from nbits to 8bits

    Notes: Pretty sure this function is a little broken!
    """
    if nbit > 8:
        raise ValueError("unpack: nbit must be <= 8")
    if 8 % nbit != 0:
        raise ValueError("unpack: nbit must divide into 8")
    if data.dtype not in (np.uint8, np.int8):
        raise TypeError("unpack: dtype must be 8-bit")
    if nbit == 8:
        return data
    elif nbit == 4:
        #The process is this:
        #ABCDEFGH [Bits of one 4+4-bit value]
        #00000000ABCDEFGH [astype(uint16)]
        #0000ABCDEFGH0000 [<< 4]
        #0000ABCDXXXXEFGH [bitwise 'or' of previous two lines]
        #0000111100001111 [0x0F0F]
        #0000ABCD0000EFGH [bitwise 'and' of previous two lines]
        #ABCD0000EFGH0000 [<< 4]
        #which effectively pads the two 4-bit values with zeros on the right
        # Note: This technique assumes LSB-first ordering
        tmpdata = data.astype(np.int16)#np.empty(upshape, dtype=np.int16)
        tmpdata = (tmpdata | (tmpdata <<  8)) & 0x0F0F
        tmpdata = tmpdata << 4 # Shift into high bits to avoid needing to sign extend
        updata = tmpdata
        return updata.view(data.dtype)
    elif nbit == 2:
        data = unpack_2to8(data)
        return data
    elif nbit == 1:
        tmpdata = data.astype(np.int64)#np.empty(upshape, dtype=np.int16)
        tmpdata = (tmpdata | (tmpdata << 32)) & 0x0000000F0000000F
        tmpdata = (tmpdata | (tmpdata << 16)) & 0x0003000300030003
        tmpdata = (tmpdata | (tmpdata <<  8)) & 0x0101010101010101
        tmpdata = tmpdata << 7 # Shift into high bits to avoid needing to sign extend
        updata = tmpdata
        return updata.view(data.dtype)



def unpack_2to8(data):
    """ Promote 2-bit unisgned data into 8-bit unsigned data.

    Args:
        data: Numpy array with dtype == uint8

    Notes:
        DATA MUST BE LOADED as np.array() with dtype='uint8'.

        This works with some clever shifting and AND / OR operations.
        Data is LOADED as 8-bit, then promoted to 32-bits:
                                      |ABCD EFGH| (8 bits of data)
        |0000 0000|0000 0000|0000 0000|ABCD EFGH| (8 bits of data as a 32-bit word)

        Once promoted, we can do some shifting, AND and OR operations:
        |0000 0000|0000 ABCD|EFGH 0000|0000 0000| (shifted << 12)
        |0000 0000|0000 ABCD|EFGH 0000|ABCD EFGH| (bitwise OR of previous two lines)
        |0000 0000|0000 ABCD|0000 0000|0000 EFGH| (bitwise AND with mask 0xF000F)
        |0000 00AB|CD00 0000|0000 00EF|GH00 0000| (prev. line shifted << 6)
        |0000 00AB|CD00 ABCD|0000 00EF|GH00 EFGH| (bitwise OR of previous two lines)
        |0000 00AB|0000 00CD|0000 00EF|0000 00GH| (bitwise AND with 0x3030303)

        Then we change the view of the data to interpret it as 4x8 bit:
        [000000AB, 000000CD, 000000EF, 000000GH]  (change view from 32-bit to 4x8-bit)

    """
    tmp = data.astype(np.uint32)
    tmp = (tmp | (tmp << 12)) & 0xF000F
    tmp = (tmp | (tmp << 6))  & 0x3030303
    tmp = tmp.byteswap()
    return tmp.view('uint8')
