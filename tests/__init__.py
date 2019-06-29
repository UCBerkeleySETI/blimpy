import subprocess
import sys
from os import path, listdir

print("Reached Init!")
here = path.dirname(path.abspath(__file__))
print(listdir(here))
if "Voyager_data" not in listdir(here):
    print("Downloading Data!")
    if sys.version_info >= (3, 0):
        subprocess.run(["sh", "download_data.sh"])
    else:
        subprocess.call(["sh", "download_data.sh"])