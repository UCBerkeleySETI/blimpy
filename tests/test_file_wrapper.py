# the lack of some implementations causes a maximum of
# 18 guaranteed misses in codecov

import blimpy as bl
import pytest
from tests.data import voyager_fil, voyager_h5, here

def test_read_fns():
    """ These read functions are currently not implemented. """
    """ L: if that is true, why do the first two lines
        NOT raise NotImplementedError's?"""
    a = bl.Waterfall(voyager_fil)
    b = bl.Waterfall(voyager_h5)
    with pytest.raises(NotImplementedError):
        a.container.read_all()
        a.container.read_row(0)
        a.container.read_rows(0, 2)

        b.container.read_all()
        b.container.read_row(0)
        b.container.read_rows(0, 2)

def test_file_wrapper_open_h5_file():
    from blimpy.io.file_wrapper import open_file
    h5r = open_file(voyager_h5)

    # check all if branches
    h5r.populate_timestamps(update_header=False)
    h5r.populate_timestamps(update_header=True)
  
    with pytest.raises(ValueError):
        h5r.read_blob((300, 300, 300, 300), -1)        

    # testing to make sure all internal functions are already used
    h5r._setup_selection_range(
        f_start=3, f_stop=4, t_stop=3, t_start=4
    )
    h5r._setup_selection_range()    

    h5r._init_empty_selection()
    
    # perform reading again and we make sure that
    # the user-facing functions survive

    h5r.populate_freqs()
    
    # just arbitrary values
    h5r.calc_n_blobs(300)
    h5r.read_blob((300, 300, 300, 300))

    h5r.calc_n_coarse_chan()
    h5r.calc_n_coarse_chan(1)
    # we have at least 2 if branches that we are not exploring

def test_file_wrapper_open_fil_file():
    from blimpy.io.file_wrapper import open_file
    filr = open_file(voyager_fil)
    
    filr.calc_n_blobs(300)
    filr.read_blob((300, 300, 300, 300))

def test_invalid_files():
    from blimpy.io.file_wrapper import open_file
    
    # script should error-out if file does not exist
    with pytest.raises(IOError):
        open_file(here + 'file_does_not_exist.h5')

    # script should error-out if file is not a
    # valid blimpy data file
    with pytest.raises(NotImplementedError):
        open_file(here + '/run_tests.sh')

if __name__ == "__main__":
    test_read_fns()
    test_file_wrapper_open_file()
