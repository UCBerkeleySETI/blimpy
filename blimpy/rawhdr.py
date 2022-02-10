r"""
Read the specified raw file.
Examine & print the required fields.
If verbose, print every header field value.
"""

import sys
from argparse import ArgumentParser
import blimpy


DEBUGGING = False


def check_int_field(header, key, valid_values, required=True):
    """
    Check an integer header field for validity.

    Parameters
    ----------
    header : dict
        Header of the .raw file.
    key : str
        Field's key value.
    valid_values : tuple
       The list of valid values or None.
    required : boolean, optional
        Required? The default is True.

    Returns
    -------
    int
        0 : valid value; 1 : invalid or missing (and required).

    """
    if key in header.keys():
        try:
            field = int(header[key])
            if valid_values is None:
                print(F"\t{key} = {field}")
                return 0
            if field in valid_values:
                print(F"\t{key} = {field}")
                return 0
            print(F"\t*** ERROR VALUE *** {key} = {field}")
            return 1
        except:
            print(F"\t*** NOT int *** {key} = {header[key]}")
            return 1

    if required:
        print(F"\t*** MISSING *** {key}")
        return 1
    print(F"\t{key} is not present but not required")
    return 0


def check_float_field(header, key):
    """
    Check a float header field for validity.

    Parameters
    ----------
    header : dict
        Header of the .raw file.
    key : str
        Field's key value.

    Returns
    -------
    int
        0 : valid value; 1 : invalid.

    """
    if key in header.keys():
        try:
            field = float(header[key])
            print(F"\t{key} = {field}")
            return 0
        except:
            print(F"\t*** NOT float *** {key} = {header[key]}")
            return 1

    print(F"\t*** MISSING *** {key}")
    return 1


def examine_header(filepath):
    """
    Examine the critical .raw file header fields.

    Parameters
    ----------
    filepath : str
        Input .raw file path.

    Returns
    -------
    rc : int
        0 : no errors; n>0 : at least one error.

    """
    if DEBUGGING:
        print("DEBUG calling GuppiRaw")
    gr = blimpy.GuppiRaw(filepath)
    if DEBUGGING:
        print("DEBUG called GuppiRaw")
    header, _ = gr.read_header()
    if DEBUGGING:
        print("DEBUG header =", header)
    rc = 0
    rc += check_int_field(header, "OBSNCHAN", None)
    rc += check_int_field(header, "NPOL", [1, 2, 4])
    rc += check_int_field(header, "NBITS", [2, 4, 8, 16])
    rc += check_int_field(header, "BLOCSIZE", None)
    rc += check_int_field(header, "PKTIDX", None)
    rc += check_int_field(header, "DIRECTIO", None, required=False)
    rc += check_int_field(header, "BEAM_ID", None, required=False)
    rc += check_int_field(header, "NBEAM", None, required=False)
    rc += check_int_field(header, "NANTS", None, required=False)
    rc += check_float_field(header, "TBIN")
    rc += check_float_field(header, "OBSFREQ")
    rc += check_float_field(header, "OBSBW")

    return rc


def cmd_tool(args=None):
    """
    rawhdr command line entry point

    Parameters
    ----------
    args : ArgParse, optional
       Command line arguments. The default is None.

    Returns
    -------
    rc : int
        0 : no errors; n>0 : at least one error.

    """
    p = ArgumentParser(description="Show the individual fields of the first header for a given raw file.")
    p.add_argument("filepath", type=str, help="Name of raw guppi file path to access")
    p.add_argument("--verbose", "-v", dest="verbose", action="store_true",
                   help="Show all of the first header fields.")
    if args is None:
        args = p.parse_args()
    else:
        args = p.parse_args(args)

    if args.verbose:
        print("rawhdr: All fields of the first header .....")
        with open(args.filepath, "rb") as fh:
            while True:
                buffer = fh.read(80).decode("utf-8").strip()
                print("\t", buffer)
                if buffer[0:3] == "END":
                    break

    print("rawhdr: Critical rawspec fields .....")
    rc = examine_header(args.filepath)
    if rc != 0:
        print("*** At least one required raw header field is missing or invalid!")
        return rc

    print("rawhdr: No errors found.")
    return rc

if __name__ == "__main__":
    sys.exit(cmd_tool())
