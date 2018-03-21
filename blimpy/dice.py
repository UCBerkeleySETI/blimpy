#!/usr/bin/env python
'''
    Script to dice data to course channel level. From BL FIL of HDF5 files, and outputs HDF5 with '_diced' appended to the file name.

    ..author: Greg Hellbourg (gregory.hellbourg@berkeley.edu)

    March 2018
'''


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

def cmd_tool():
    '''read input and output frequency, and output file name
    '''

    parser = argparse.ArgumentParser(description='Dices hdf5 or fil files and writes to hdf5')
    parser.add_argument('-f', '--input_filename', action='store', default=None, dest='in_fname', type=str, help='Name of file to write from (HDF5 or FIL)')
    parser.add_argument('-b', action='store', default=None, dest='f_start', type=float, help='Start frequency in MHz')
    parser.add_argument('-e', action='store', default=None, dest='f_stop', type=float, help='Stop frequency in MHz')
    parser.add_argument('-o', '--output_filename', action='store', default=None, dest='out_fname', type=str, help='Name of file to write to (HDF5)')

    args = parser.parse_args()

    if len(sys.argv) == 1:
            logger.info('indicate file name and start and stop frequencies\n')
            sys.exit()

    if args.in_fname == None:
            logger.info('need to indicate input file name\n'0
            sys.exit()

    if args.out_fname == None:
            if args.in_fname[len(args.in_fname)-4:] == '.fil':
                    args.out_fname = args.in_fname
                    args.out_fname = args.out_fname.replace('.fil','_diced.h5')
            if args.in_fname[len(args.in_fname)-3:] == '.h5':
                    args.out_fname = args.in_fname
                    args.out_fname = args.out_fname.replace('.h5','_diced.h5')

    if args.f_start == None or args.f_stop == None:
            logger.info('start and end frequencies must be informed\n')
            sys.exit()

    # read start frequency and bandwidth from data set

    file_big = Waterfall(args.in_fname)
    f_min_file = file_big.header['fch1']
    f_max_file = file_big.header['fch1'] + file_big.header['nchans']*file_big.header['foff']

    if f_max_file < f_min_file:
            f_max_file,f_min_file = f_min_file,f_max_file

    FreqBWFile = f_max_file-f_min_file
    stdDF = FreqBWFile / float(file_big.calc_n_coarse_chan())   #stdDF = 2.9296875

    if args.f_stop < args.f_start:
            args.f_stop,args.f_start = args.f_start,args.f_stop

    if args.f_start < f_max_file and args.f_start > f_min_file and args.f_stop > f_max_file:
            args.f_stop = f_max_file
            logger.info('\nWarning : higher frequency set to ' + str(f_max_file) + ' MHz to match file\n')

    if args.f_stop < f_max_file and args.f_stop > f_min_file and args.f_start < f_min_file:
            args.f_start = f_min_file
            logger.info('\nWarning : lower frequency set to ' + str(f_min_file) + ' MHz to match file\n')

    if args.f_start < f_min_file and args.f_stop > f_max_file:
            args.f_start = f_min_file
            args.f_stop = f_max_file
            logger.info('\nWarning : lower frequency set to ' + str(f_min_file) + ' MHz\nand higher frequency set to ' + str(f_max_file) + ' MHz to match file\n')
            # print '\nindicated frequencies include file frequency span - no need to dice\n'
            # sys.exit()

    if min(args.f_start,args.f_stop) < f_min_file or max(args.f_start,args.f_stop) > f_max_file:
            logger.info('\nbandwidth to extract must be within ' + str(f_min_file) + ' MHz and ' + str(f_max_file) + ' MHz\n')
            sys.exit()

    # calculate real coarse channel begin and end freqs
    f_start_real = math.floor((min(args.f_start,args.f_stop) - f_min_file)/stdDF)*stdDF + f_min_file
    f_stop_real = f_max_file - math.floor((f_max_file - max(args.f_start,args.f_stop))/stdDF)*stdDF

    # print
    # print "true start frequency is " + str(f_start_real)
    # print "true stop frequency is " + str(f_stop_real)

    logger.info('\nwriting to ' + args.out_fname + ' extacting from ' + str(f_start_real) + ' MHz to ' + str(f_stop_real) + ' MHz\n')

    # create waterfall object
    file_small = Waterfall(args.in_fname, f_start = f_start_real, f_stop = f_stop_real)

    # write waterfall object
    file_small.write_to_hdf5(args.out_fname)



if __name__ == "__main__":

    cmd_tool()
