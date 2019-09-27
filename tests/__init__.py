import subprocess
import sys
from os import path, listdir

print("Reached blimpy.tests init!")
here = path.dirname(path.abspath(__file__))
print("Running tests from {}".format(here))
if "test_data" not in listdir(here):
    print("Test data has not yet been downloaded. Downloading Data...")
    if sys.version_info >= (3, 0):
        subprocess.run(["sh", "download_data.sh"])
    else:
        subprocess.call(["sh", "download_data.sh"])
