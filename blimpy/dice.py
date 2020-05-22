#!/usr/bin/env python
"""
    Script to dice data to course channel level. From BL FIL of HDF5 files, and outputs HDF5 with '_diced' appended to the file name.

    ..author: Greg Hellbourg (gregory.hellbourg@berkeley.edu)

    March 2018
"""


from .waterfall import Waterfall
import argparse
import math
import sys

#------
# Logging set up
import logging
logger = logging.getLogger(__name__)

level_log = logging.INFO

if level_log == logging.INFO:
    stream = sys.stdout
    format = '%(name)-15s %(levelname)-8s %(message)s'
else:
    stream =  sys.stderr
    format = '%%(relativeCreated)5d (name)-15s %(levelname)-8s %(message)s'

logging.basicConfig(format=format,stream=stream,level = level_log)
#------

def cmd_tool(args=None):
    """ Dices (extracts frequency range) hdf5 or fil files to new file.

    optional arguments:
      -h, --help            show this help message and exit
      -f IN_FNAME, --input_filename IN_FNAME
                            Name of file to write from (HDF5 or FIL)
      -b F_START            Start frequency in MHz
      -e F_STOP             Stop frequency in MHz
      -x OUT_FORMAT, --output_file OUT_FORMAT
                            Output file format [.h5 or .fil].
      -o OUT_FNAME, --output_filename OUT_FNAME
                            Ouput file name to write (to HDF5 or FIL).
      -l MAX_LOAD           Maximum data limit to load. Default:1GB
    """

    parser = argparse.ArgumentParser(description='Dices (extracts frequency range)  hdf5 or fil files and writes to hdf5 or fil.')
    parser.add_argument('-f', '--input_filename', action='store', default=None, dest='in_fname', type=str, help='Name of file to write from (HDF5 or FIL)')
    parser.add_argument('-b', action='store', default=None, dest='f_start', type=float, help='Start frequency in MHz')
    parser.add_argument('-e', action='store', default=None, dest='f_stop', type=float, help='Stop frequency in MHz')
    parser.add_argument('-x', '--output_file', action='store', default=None, dest='out_format', type=str, help='Output file format [.h5 or .fil].')
    parser.add_argument('-o', '--output_filename', action='store', default=None, dest='out_fname', type=str, help='Ouput file name to write (to HDF5 or FIL).')
    parser.add_argument('-l', action='store', default=None, dest='max_load', type=float,help='Maximum data limit to load. Default:1GB')

    if args is None:
        args = sys.argv[1:]

        if len(sys.argv) == 1:
            logger.error('Indicate file name and start and stop frequencies')
            sys.exit()

    args = parser.parse_args(args)

    if args.in_fname is None:
            logger.error('Need to indicate input file name')
            sys.exit()

    if args.out_fname is None:
        if (args.out_format is None) or (args.out_format == 'h5'):
            if args.in_fname[len(args.in_fname)-4:] == '.fil':
                args.out_fname = args.in_fname
                args.out_fname = args.out_fname.replace('.fil','_diced.h5')
            elif args.in_fname[len(args.in_fname)-3:] == '.h5':
                args.out_fname = args.in_fname
                args.out_fname = args.out_fname.replace('.h5','_diced.h5')
            else:
                logger.error('Input file not recognized')
                sys.exit()
        elif args.out_format == 'fil':
            if args.in_fname[len(args.in_fname)-4:] == '.fil':
                args.out_fname = args.in_fname
                args.out_fname = args.out_fname.replace('.fil','_diced.fil')
            elif args.in_fname[len(args.in_fname)-3:] == '.h5':
                args.out_fname = args.in_fname
                args.out_fname = args.out_fname.replace('.h5','_diced.fil')
            else:
                logger.error('input file not recognized.')
                sys.exit()

        else:
            logger.error('Must indicate either output file name or valid output file extension.')
            sys.exit()

    elif (args.out_fname[len(args.out_fname)-4:] == '.fil') and (args.out_format == 'h5'):
        logger.error('Output file extension does not match output file name')
        sys.exit()

    elif (args.out_fname[len(args.out_fname)-3:] == '.h5') and (args.out_format == 'fil'):
        logger.error('Output file extension does not match output file name.')
        sys.exit()

    if (args.out_fname[len(args.out_fname)-3:] != '.h5') and (args.out_fname[len(args.out_fname)-4:] != '.fil'):
        logger.error('Indicate output file name with extension, or simply output file extension.')
        sys.exit()

    if args.f_start == None and args.f_stop == None:
        logger.error('Please give either start and/or end frequencies. Otherwise use fil2h5 or h52fil functions.')
        sys.exit()

    #Read start frequency and bandwidth from data set
    file_big = Waterfall(args.in_fname, max_load = args.max_load)
    f_min_file = file_big.header['fch1']
    f_max_file = file_big.header['fch1'] + file_big.header['nchans'] * file_big.header['foff']

    if args.f_start == None:
        logger.warning('Lower frequency not given, setting to ' + str(f_min_file) + ' MHz to match file.')

    if args.f_stop == None:
        logger.warning('Higher frequency not given, setting to ' + str(f_max_file) + ' MHz to match file.')


    if f_max_file < f_min_file:
            f_max_file,f_min_file = f_min_file,f_max_file

    FreqBWFile = f_max_file-f_min_file
    stdDF = FreqBWFile / float(file_big.calc_n_coarse_chan())   #stdDF = 2.9296875

    if args.f_stop < args.f_start:
            args.f_stop,args.f_start = args.f_start,args.f_stop

    if args.f_start < f_max_file and args.f_start > f_min_file and args.f_stop > f_max_file:
            args.f_stop = f_max_file
            logger.warning('Higher frequency set to ' + str(f_max_file) + ' MHz to match file.')

    if args.f_stop < f_max_file and args.f_stop > f_min_file and args.f_start < f_min_file:
            args.f_start = f_min_file
            logger.warning('Lower frequency set to ' + str(f_min_file) + ' MHz to match file.')

    if args.f_start < f_min_file and args.f_stop > f_max_file:
            args.f_start = f_min_file
            args.f_stop = f_max_file
            logger.warning('Lower frequency set to ' + str(f_min_file) + ' MHz and higher frequency set to ' + str(f_max_file) + ' MHz to match file.')
            # print '\nindicated frequencies include file frequency span - no need to dice\n'
            # sys.exit()

    if min(args.f_start,args.f_stop) < f_min_file or max(args.f_start,args.f_stop) > f_max_file:
            logger.error('Bandwidth to extract must be within ' + str(f_min_file) + ' MHz and ' + str(f_max_file) + ' MHz.')
            sys.exit()

    # calculate real coarse channel begin and end freqs
    f_start_real = math.floor((min(args.f_start,args.f_stop) - f_min_file)/stdDF)*stdDF + f_min_file
    f_stop_real = f_max_file - math.floor((f_max_file - max(args.f_start,args.f_stop))/stdDF)*stdDF

    # print "true start frequency is " + str(f_start_real)
    # print "true stop frequency is " + str(f_stop_real)

    logger.info('Writing to ' + args.out_fname)
    logger.info('Extacting from ' + str(f_start_real) + ' MHz to ' + str(f_stop_real) + ' MHz.')

    # create waterfall object
    file_small = Waterfall(args.in_fname, f_start = f_start_real, f_stop = f_stop_real, max_load = args.max_load)

    # write waterfall object
    if args.out_fname[len(args.out_fname)-4:] == '.fil':
        file_small.write_to_fil(args.out_fname)
    elif args.out_fname[len(args.out_fname)-3:] == '.h5':
        file_small.write_to_hdf5(args.out_fname)
    else:
        logger.error('Error in output file creation : verify output file name and extension.')
        sys.exit()


if __name__ == "__main__":
    cmd_tool()
