import h5py
import hdf5plugin
import numpy
import tempfile


def test_is_h5py_correctly_installed():
    """
    If this test fails you probably need to install h5py from source manually:

    $ pip install --no-binary=h5py h5py
    """
    f = h5py.File(tempfile.gettempdir() + '/h5testfile', "w")
    block_size = 0
    dataset = f.create_dataset(
        "data",
        (100, 100, 100),
        dtype='float32',
        **hdf5plugin.Bitshuffle(nelems=0, lz4=True)
    )

    array = numpy.random.rand(100, 100, 100)
    array = array.astype('float32')
    dataset[:] = array
    f.close()

