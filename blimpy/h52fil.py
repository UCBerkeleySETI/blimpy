#!/usr/bin/env python
'''
    Simple script for quicly making an .fil file from a .h5.

    ..author: Emilio Enriquez (jeenriquez@gmail.com)

    July 28th 2017
'''

try:
    from .waterfall import Waterfall
except:
    from waterfall import Waterfall

from optparse import OptionParser
import sys
import os

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


def make_fil_file(filename,out_dir='./', new_filename=None, max_load = None):
    ''' Converts file to Sigproc filterbank (.fil) format.  Default saves output in current dir.
    '''

    fil_file = Waterfall(filename, max_load = max_load)
    if not new_filename:
        new_filename = out_dir+filename.replace('.h5','.fil').split('/')[-1]

    if '.fil' not in new_filename:
        new_filename = new_filename+'.fil'

    fil_file.write_to_fil(new_filename)

def cmd_tool():
    '''
    '''

    p = OptionParser()
    p.set_usage('Command line utility for converting HDF5 (.h5) to Sigproc filterbank (.fil) format \n >>h52fil <FULL_PATH_TO_FIL_FILE> [options]')
    p.add_option('-o', '--out_dir', dest='out_dir', type='str', default='./', help='Location for output files. Default: local dir. ')
    p.add_option('-n', '--new_filename', dest='new_filename', type='str', default='', help='New filename. Default: replaces extention to .fil')
    p.add_option('-d', '--delete_input', dest='delete_input', action='store_true', default=False, help='This option deletes the input file after conversion.')
    p.add_option('-l', action='store', default=None, dest='max_load', type=float,help='Maximum data limit to load. Default:1GB')

    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        logger.info('Please specify a file name \nExiting.')
        sys.exit()
    else:
        filename = args[0]

    make_fil_file(filename, out_dir = opts.out_dir, new_filename=opts.new_filename, max_load = opts.max_load)

    if opts.delete_input:
        logger.info("'Deleting input file: %s"%(filename))
        os.remove(filename)

if __name__ == "__main__":

    cmd_tool()
