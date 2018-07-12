import blimpy as bl
import pytest

def test_read_fns():
    """ These read functions are currently not implemented. """
    a = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.fil')
    b = bl.Waterfall('Voyager_data/Voyager1.single_coarse.fine_res.h5')
    with pytest.raises(NotImplementedError):
        a.container.read_all()
        a.container.read_row(0)
        a.container.read_rows(0, 2)

        b.container.read_all()
        b.container.read_row(0)
        b.container.read_rows(0, 2)

def test_file_wrapper_open_file():
    from blimpy.file_wrapper import open_file
    open_file('Voyager_data/Voyager1.single_coarse.fine_res.h5')
    open_file('Voyager_data/Voyager1.single_coarse.fine_res.fil')

    with pytest.raises(NotImplementedError):
        open_file('run_tests.sh')

if __name__ == "__main__":
    test_read_fns()
    test_file_wrapper_open_file()