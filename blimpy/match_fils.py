#!/usr/bin/env python

from optparse import OptionParser
import socket
import subprocess
import sys
import os
try:
    from .waterfall import Waterfall
except:
    from blimpy import Waterfall

#import pdb #pdb.set_trace()

local_host = socket.gethostname()

def reset_outs():
    '''
    '''

    return None,None

def make_batch_script():
    """ This creates a batch script for the heavy lifting"""

    script_text  = '#! /bin/bash\n# This script calculates the md5sum for the last number of bits (nbits).# Usage :\n#     ./tail_sum file bnits\n#\n\nfile_name=$1\nbit_num=$2\n\ntail -c $bit_num $file_name | md5sum\n'

    with open('tail_sum.sh', 'w') as batch_script:
        batch_script.write(script_text)

    os.chmod('tail_sum.sh', 0o775)

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

    if 'bl' in local_host:
        header_loc = '/usr/local/sigproc/bin/header' #Current location of header command in GBT.
    else:
        raise IOError('Script only able to run in BL systems.')

    p = OptionParser()
    p.set_usage('matchfils <FIL_FILE1> <FIL_FILE2>')
    opts, args = p.parse_args(sys.argv[1:])

    file1 = args[0]
    file2 = args[1]

    #------------------------------------
    #Create batch script

    make_batch_script()

    #------------------------------------
    #First checksum

    headersize1 = find_header_size(file1)
    file_size1 = os.path.getsize(file1)

    #Strip header from file, and calculate the md5sum of the rest.
    #command=['tail','-c',str(file_size1-headersize1),file1,'|','md5sum']
    command=['./tail_sum.sh',file1,str(file_size1-headersize1)]
    print('[matchfils] '+' '.join(command))

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    check_sum1 = out.split()[0]
    print('[matchfils] Checksum is:', check_sum1)

    if err:
        raise Error('There is an error.')

    #---
    out,err = reset_outs()

    command=[header_loc,file1]
    print('[matchfils] Header information:')

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    header1 = out
    print(header1)

    #------------------------------------
    #Second checksum

    out,err = reset_outs()

    headersize2 = find_header_size(file2)
    file_size2 = os.path.getsize(file2)

    #Strip header from file, and calculate the md5sum of the rest.
    command=['./tail_sum.sh',file2,str(file_size2-headersize2)]
    print('[matchfils] '+' '.join(command))

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    check_sum2 = out.split()[0]
    print('[matchfils] Checksum is:', check_sum2)

    if err:
        raise Error('There is an error.')

    #---
    out,err = reset_outs()

    command=[header_loc,file2]
    print('[matchfils] Header information:')

    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()

    header2 = out
    print(header2)

    #------------------------------------
    #check the checksums

    if check_sum1 != check_sum2:
        print('[matchfils] Booo! Checksum does not match between files.')
    else:
        print('[matchfils] Hooray! Checksum matches between files.')

    #------------------------------------
    #Remove batch script

    os.remove('tail_sum.sh')

if __name__ == "__main__":
    cmd_tool()

