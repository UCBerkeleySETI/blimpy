from os.path import dirname
from blimpy.stix import cmd_tool
from tests.data import voyager_fil


PLOT_DIR = dirname(voyager_fil)


def execute_command(args):
    print("\ntest_stix: args:", args)
    cmd_tool(args)


def test_stix():

    args = [voyager_fil, "16", "--plot_dir",  PLOT_DIR]
    execute_command(args)

    args = [voyager_fil, "4", "--plot_dir",  PLOT_DIR, "-s", "v"]
    execute_command(args)

    args = [voyager_fil, "4", "-p",  PLOT_DIR, "--stitch", "h"]
    execute_command(args)

    args = [voyager_fil, "4", "--plot_dir",  PLOT_DIR, "--stitch", "n",
            "--dpi", "100", "--width", "8", "--height", "6"]
    execute_command(args)

    args = [voyager_fil, "4", "-p",  PLOT_DIR, "-s", "n", "-d", "100", "-w", "8", "-t", "6"]
    execute_command(args)


if __name__ == "__main__":
    test_stix()
