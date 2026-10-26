r""" Tests for blimpy package version detection """
import subprocess
import sys

try:
    from importlib.metadata import version
except ImportError:
    from importlib_metadata import version

import blimpy


def test_version_matches_metadata():
    assert blimpy.__version__ == version('blimpy')
    assert not blimpy.__version__.startswith('0.0.0')

def test_import_has_no_pkg_resources_warning():
    # Use a fresh interpreter: blimpy is already imported in this process.
    result = subprocess.run(
        [sys.executable, '-W', 'default', '-c', 'import blimpy'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert result.returncode == 0
    assert 'pkg_resources' not in result.stderr
