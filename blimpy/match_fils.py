#!/usr/bin/env python

from optparse import OptionParser
import socket
import subprocess
import sys
import os
from .waterfall import Waterfall

#import pdb #pdb.set_trace()

local_host = socket.gethostname()

if 'bl' in local_host:
    header_loc = '/usr/local/sigproc/bin/header' #Current location of header command in GBT.
else:
    raise IOError('Spcript only able to run in BL systems.')

def reset_outs():
    '''
    '''

    return None,None

def find_header_size(filename):
    ''' Script to find the header size of a filterbank file'''

    # open datafile
    filfile=open(filename,'rb')
    # go to the start of the file
    filfile.seek(0)
    #read some region larger than the header.
    round1 = filfile.read(1000)
    headersize = round1.find('HEADER_END')+len('HEADER_END')

    return headersize

def cmd_tool(args=None):
    """ Command line tool to make a md5sum comparison of two .fil files. """

    p = OptionParser()
    p.set_usage('compare_fils <FIL_FILE1> <FIL_FILE2>')
    opts, args = p.parse_args(sys.argv[1:])

    file1 = args[0]
    file2 = args[1]

    #------------------------------------
    #First checksum

    headersize1 = find_header_size(file1)
    file_size1 = os.path.getsize(file1)

    command=['tail','-c',str(file_size1-headersize1),file1,'|','md5sum']
    print '[matchfils] '+' '.join(command)

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    check_sum1 = out.split()[0]

    out,err = reset_outs()

    command=[header_loc,file1]
    print '[matchfils] Header information:'

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    header1 = out
    print header1

    #------------------------------------
    #Second checksum

    out,err = reset_outs()

    headersize2 = find_header_size(file2)
    file_size2 = os.path.getsize(file2)

    command=['tail','-c',str(file_size2-headersize2),file2,'|','md5sum']
    print '[matchfils] '+' '.join(command)

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    check_sum2 = out.split()[0]

    out,err = reset_outs()

    command=[header_loc,file2]
    print '[matchfils] Header information:'

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    header2 = out
    print header2

    #------------------------------------
    #check the checksums

    if check_sum1 != check_sum2:
        print '[matchfils] Booo! Checksum does not match between files.'
    else:
        print '[matchfils] Hooray! Checksum matches between files.'

if __name__ == "__main__":
    cmd_tool()

