import blimpy
import argparse
import math
import sys

stdDF = 2.9296875


# read input and output frequency, and output file name

parser = argparse.ArgumentParser(description='Dices hdf5 or fil files and writes to hdf5')
parser.add_argument('-f', '--input_filename', action='store', default=None, dest='in_fname', type=str, help='Name of file to write from (HDF5 or FIL)')
parser.add_argument('-b', action='store', default=None, dest='f_start', type=float, help='Start frequency in MHz')
parser.add_argument('-e', action='store', default=None, dest='f_stop', type=float, help='Stop frequency in MHz')
parser.add_argument('-o', '--output_filename', action='store', default=None, dest='out_fname', type=str, help='Name of file to write to (HDF5)')

args = parser.parse_args()

if len(sys.argv) == 1:
        print '\nindicate file name and start and stop frequencies\n'
        sys.exit()

if args.in_fname == None:
        print '\nneed to indicate input file name\n'
        sys.exit()

if args.out_fname == None:
        if args.in_fname[len(args.in_fname)-4:] == '.fil':
                args.out_fname = args.in_fname
                args.out_fname = args.out_fname.replace('.fil','_diced.h5')
        if args.in_fname[len(args.in_fname)-3:] == '.h5':
                args.out_fname = args.in_fname
                args.out_fname = args.out_fname.replace('.h5','_diced.h5')

if args.f_start == None or args.f_stop == None:
        print '\nstart and end frequencies must be informed\n'
        sys.exit()

# read start frequency and bandwidth from data set

file_big = blimpy.Waterfall(args.in_fname)
f_min_file = file_big.header['fch1']
f_max_file = file_big.header['fch1'] + file_big.header['nchans']*file_big.header['foff']
if f_max_file < f_min_file:
        f_max_file,f_min_file = f_min_file,f_max_file

if args.f_stop < args.f_start:
        args.f_stop,args.f_start = args.f_start,args.f_stop

if args.f_start < f_max_file and args.f_start > f_min_file and args.f_stop > f_max_file:
        args.f_stop = f_max_file
        print '\nWarning : higher frequency set to ' + str(f_max_file) + ' MHz to match file\n'

if args.f_stop < f_max_file and args.f_stop > f_min_file and args.f_start < f_min_file:
        args.f_start = f_min_file
        print '\nWarning : lower frequency set to ' + str(f_min_file) + ' MHz to match file\n'
		
if args.f_start < f_min_file and args.f_stop > f_max_file:
        args.f_start = f_min_file
        args.f_stop = f_max_file
        print '\nWarning : lower frequency set to ' + str(f_min_file) + ' MHz\nand higher frequency set to ' + str(f_max_file) + ' MHz to match file\n'
        # print '\nindicated frequencies include file frequency span - no need to dice\n'
        # sys.exit()

if min(args.f_start,args.f_stop) < f_min_file or max(args.f_start,args.f_stop) > f_max_file:
        print '\nbandwidth to extract must be within ' + str(f_min_file) + ' MHz and ' + str(f_max_file) + ' MHz\n'
        sys.exit()

# calculate real coarse channel begin and end freqs

f_start_real = math.floor((min(args.f_start,args.f_stop) - f_min_file)/stdDF)*stdDF + f_min_file
f_stop_real = f_max_file - math.floor((f_max_file - max(args.f_start,args.f_stop))/stdDF)*stdDF

# print
# print "true start frequency is " + str(f_start_real)
# print "true stop frequency is " + str(f_stop_real)

print '\nwriting to ' + args.out_fname + ' extacting from ' + str(f_start_real) + ' MHz to ' + str(f_stop_real) + ' MHz\n'

# create waterfall object

file_small = blimpy.Waterfall(args.in_fname, f_start = f_start_real, f_stop = f_stop_real)

# write waterfall object

file_small.write_to_hdf5(args.out_fname)