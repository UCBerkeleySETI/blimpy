import struct
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
import os
import six

try:
    import h5py
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False


###
# Header parsing
###

# Dictionary of allowed keywords and their types
# Here are the keywordss that a filter bank file may
# contain.  Items marked with "[*]" are not yet # supported.  See docs for
# indivisuabl attribtues for more detailed info.
#
#   * telescope_id (int): 0=fake data; 1=Arecibo; 2=Ooty... others to be added
#   * machine_id (int): 0=FAKE; 1=PSPM; 2=WAPP; 3=OOTY... others to be added
#   * data_type (int): 1=blimpy; 2=time series... others to be added
#   * rawdatafile (string): the name of the original data file
#   * source_name (string): the name of the source being observed by the telescope
#   * barycentric (int): equals 1 if data are barycentric or 0 otherwise
#   * pulsarcentric (int): equals 1 if data are pulsarcentric or 0 otherwise
#   * az_start (double): telescope azimuth at start of scan (degrees)
#   * za_start (double): telescope zenith angle at start of scan (degrees)
#   * src_raj (double): right ascension (J2000) of source (hours, converted from hhmmss.s)
#   * src_dej (double): declination (J2000) of source (degrees, converted from ddmmss.s)
#   * tstart (double): time stamp (MJD) of first sample
#   * tsamp (double): time interval between samples (s)
#   * nbits (int): number of bits per time sample
#   * nsamples (int): number of time samples in the data file (rarely used any more)
#   * fch1 (double): centre frequency (MHz) of first blimpy channel
#   * foff (double): blimpy channel bandwidth (MHz)
#   * FREQUENCY_START [*] (character): start of frequency table (see below for explanation)
#   * fchannel [*] (double): frequency channel value (MHz)
#   * FREQUENCY_END [*] (character): end of frequency table (see below for explanation)
#   * nchans (int): number of blimpy channels
#   * nifs (int): number of seperate IF channels
#   * refdm (double): reference dispersion measure (pc/cm**3)
#   * period (double): folding period (s)
#   * nbeams (int):total number of beams (?)
#   * ibeam (int): number of the beam in this file (?)

header_keyword_types = {
    b'telescope_id' : b'<l',
    b'machine_id'   : b'<l',
    b'data_type'    : b'<l',
    b'barycentric'  : b'<l',
    b'pulsarcentric': b'<l',
    b'nbits'        : b'<l',
    b'nsamples'     : b'<l',
    b'nchans'       : b'<l',
    b'nifs'         : b'<l',
    b'nbeams'       : b'<l',
    b'ibeam'        : b'<l',
    b'rawdatafile'  : b'str',
    b'source_name'  : b'str',
    b'az_start'     : b'<d',
    b'za_start'     : b'<d',
    b'tstart'       : b'<d',
    b'tsamp'        : b'<d',
    b'fch1'         : b'<d',
    b'foff'         : b'<d',
    b'refdm'        : b'<d',
    b'period'       : b'<d',
    b'src_raj'      : b'angle',
    b'src_dej'      : b'angle',
    }


def len_header(filename):
    """ Return the length of the blimpy header, in bytes

    Args:
        filename (str): name of file to open

    Returns:
        idx_end (int): length of header, in bytes
    """
    with  open(filename, 'rb') as f:
        header_sub_count = 0
        eoh_found = False
        while not eoh_found:
            header_sub = f.read(512)
            header_sub_count += 1
            if b'HEADER_END' in header_sub:
                idx_end = header_sub.index(b'HEADER_END') + len(b'HEADER_END')
                eoh_found = True
                break

        idx_end = (header_sub_count -1) * 512 + idx_end
    return idx_end


def read_next_header_keyword(fh):
    """

    Args:
        fh (file): file handler

    Returns:
    """
    n_bytes = np.fromstring(fh.read(4), dtype='uint32')[0]
    #print n_bytes

    if n_bytes > 255:
        n_bytes = 16

    keyword = fh.read(n_bytes)

    #print keyword

    if keyword == b'HEADER_START' or keyword == b'HEADER_END':
        return keyword, 0, fh.tell()
    else:
        dtype = header_keyword_types[keyword]
        #print dtype
        idx = fh.tell()
        if dtype == b'<l':
            val = struct.unpack(dtype, fh.read(4))[0]
        if dtype == b'<d':
            val = struct.unpack(dtype, fh.read(8))[0]
        if dtype == b'str':
            str_len = np.fromstring(fh.read(4), dtype='uint32')[0]
            val = fh.read(str_len)
        if dtype == b'angle':
            val = struct.unpack('<d', fh.read(8))[0]
            val = fil_double_to_angle(val)
            if keyword == b'src_raj':
                val = Angle(val, unit=u.hour)
            else:
                val = Angle(val, unit=u.deg)
        return keyword, val, idx


def is_filterbank(filename):
    """ Open file and confirm if it is a filterbank file or not. """
    with open(filename, 'rb') as fh:
        is_fil = True

        # Check this is a blimpy file
        try:
            keyword, value, idx = read_next_header_keyword(fh)
            try:
                assert keyword == b'HEADER_START'
            except AssertionError:
                is_fil = False
        except KeyError:
            is_fil = False
        return is_fil


def read_header(filename, return_idxs=False):
    """ Read blimpy header and return a Python dictionary of key:value pairs

    Args:
        filename (str): name of file to open

    Optional args:
        return_idxs (bool): Default False. If true, returns the file offset indexes
                            for values

    returns

    """
    with open(filename, 'rb') as fh:
        header_dict = {}
        header_idxs = {}

        # Check this is a blimpy file
        keyword, value, idx = read_next_header_keyword(fh)

        try:
            assert keyword == b'HEADER_START'
        except AssertionError:
            raise RuntimeError("Not a valid blimpy file.")

        while True:
            keyword, value, idx = read_next_header_keyword(fh)
            if keyword == b'HEADER_END':
                break
            else:
                header_dict[keyword] = value
                header_idxs[keyword] = idx

    if return_idxs:
        return header_idxs
    else:
        return header_dict

def fix_header(filename, keyword, new_value):
    """ Apply a quick patch-up to a Filterbank header by overwriting a header value


    Args:
        filename (str): name of file to open and fix. WILL BE MODIFIED.
        keyword (stt):  header keyword to update
        new_value (long, double, angle or string): New value to write.

    Notes:
        This will overwrite the current value of the blimpy with a desired
        'fixed' version. Note that this has limited support for patching
        string-type values - if the length of the string changes, all hell will
        break loose.

    """

    # Read header data and return indexes of data offsets in file
    hd = read_header(filename)
    hi = read_header(filename, return_idxs=True)
    idx = hi[keyword]

    # Find out the datatype for the given keyword
    dtype = header_keyword_types[keyword]
    dtype_to_type = {b'<l'  : np.int32,
                     b'str' : bytes,
                     b'<d'  : np.float64,
                     b'angle' : to_sigproc_angle}
    value_dtype = dtype_to_type[dtype]

    # Generate the new string
    if isinstance(value_dtype, bytes):
        if len(hd[keyword]) == len(new_value):
            val_str = np.int32(len(new_value)).tostring() + new_value
        else:
            raise RuntimeError("String size mismatch. Cannot update without rewriting entire file.")
    else:
        val_str = value_dtype(new_value).tostring()

    # Write the new string to file
    with open(filename, 'rb+') as fh:
        fh.seek(idx)
        fh.write(val_str)

def fil_double_to_angle(angle):
      """ Reads a little-endian double in ddmmss.s (or hhmmss.s) format and then
      converts to Float degrees (or hours).  This is primarily used to read
      src_raj and src_dej header values. """

      negative = (angle < 0.0)
      angle = np.abs(angle)

      dd = np.floor((angle / 10000))
      angle -= 10000 * dd
      mm = np.floor((angle / 100))
      ss = angle - 100 * mm
      dd += mm/60.0 + ss/3600.0

      if negative:
          dd *= -1

      return dd

###
# sigproc writing functions
###

def to_sigproc_keyword(keyword, value=None):
    """ Generate a serialized string for a sigproc keyword:value pair

    If value=None, just the keyword will be written with no payload.
    Data type is inferred by keyword name (via a lookup table)

    Args:
        keyword (str): Keyword to write
        value (None, float, str, double or angle): value to write to file

    Returns:
        value_str (str): serialized string to write to file.
    """

    keyword = bytes(keyword)

    if value is None:
        return np.int32(len(keyword)).tostring() + keyword
    else:
        dtype = header_keyword_types[keyword]

        dtype_to_type = {b'<l'  : np.int32,
                         b'str' : str,
                         b'<d'  : np.float64,
                         b'angle' : to_sigproc_angle}

        value_dtype = dtype_to_type[dtype]

        if value_dtype is str:
            return np.int32(len(keyword)).tostring() + keyword + np.int32(len(value)).tostring() + value
        else:
            return np.int32(len(keyword)).tostring() + keyword + value_dtype(value).tostring()

def generate_sigproc_header(f):
    """ Generate a serialzed sigproc header which can be written to disk.

    Args:
        f (Filterbank object): Filterbank object for which to generate header

    Returns:
        header_str (str): Serialized string corresponding to header
    """

    header_string = b''
    header_string += to_sigproc_keyword(b'HEADER_START')

    for keyword in f.header.keys():
        if keyword == b'src_raj':
            header_string += to_sigproc_keyword(b'src_raj')  + to_sigproc_angle(f.header[b'src_raj'])
        elif keyword == b'src_dej':
            header_string += to_sigproc_keyword(b'src_dej')  + to_sigproc_angle(f.header[b'src_dej'])
        elif keyword == b'az_start' or keyword == b'za_start':
            header_string += to_sigproc_keyword(keyword)  + np.float64(f.header[keyword]).tostring()
        elif keyword not in header_keyword_types.keys():
            pass
        else:
            header_string += to_sigproc_keyword(keyword, f.header[keyword])

    header_string += to_sigproc_keyword(b'HEADER_END')
    return header_string


def to_sigproc_angle(angle_val):
    """ Convert an astropy.Angle to the ridiculous sigproc angle format string. """
    x = str(angle_val)

    if '.' in x:
        if 'h' in x:
                d, m, s, ss = int(x[0:x.index('h')]), int(x[x.index('h')+1:x.index('m')]), \
                int(x[x.index('m')+1:x.index('.')]), float(x[x.index('.'):x.index('s')])
        if 'd' in x:
            d, m, s, ss = int(x[0:x.index('d')]), int(x[x.index('d')+1:x.index('m')]), \
            int(x[x.index('m')+1:x.index('.')]), float(x[x.index('.'):x.index('s')])
    else:
        if 'h' in x:
            d, m, s = int(x[0:x.index('h')]), int(x[x.index('h')+1:x.index('m')]), \
            int(x[x.index('m')+1:x.index('s')])
        if 'd' in x:
            d, m, s = int(x[0:x.index('d')]), int(x[x.index('d')+1:x.index('m')]), \
            int(x[x.index('m')+1:x.index('s')])
        ss = 0
    num = str(d).zfill(2) + str(m).zfill(2) + str(s).zfill(2)+ '.' + str(ss).split(".")[-1]
    return np.float64(num).tostring()


def calc_n_ints_in_file(filename):
    """ Calculate number of integrations in a given file """
    # Load binary data
    h = read_header(filename)
    n_bytes  = int(h[b'nbits'] / 8)

    n_chans = h[b'nchans']
    n_ifs   = h[b'nifs']
    idx_data = len_header(filename)
    f = open(filename, 'rb')
    f.seek(idx_data)
    filesize = os.path.getsize(filename)
    n_bytes_data = filesize - idx_data

    if h[b'nbits'] == 2:
        n_ints = int(4 * n_bytes_data / (n_chans * n_ifs))
    else:
        n_ints = int(n_bytes_data / (n_bytes * n_chans * n_ifs))

    return n_ints
