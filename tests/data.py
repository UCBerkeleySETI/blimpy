from os import path, listdir
import subprocess
import sys

subprocess.run(["sh", "download_data.sh"])
here = path.dirname(path.abspath(__file__))
print(listdir(here))
if "Voyager_data" not in listdir(here):
    if sys.version_info >= (3, 0):
        subprocess.run(["sh", "download_data.sh"])
    else:
        subprocess.call(["sh", "download_data.sh"])
voyager_fil = path.join(here, 'Voyager_data/Voyager1.single_coarse.fine_res.fil')
voyager_h5 = path.join(here, 'Voyager_data/Voyager1.single_coarse.fine_res.h5')