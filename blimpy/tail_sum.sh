#! /bin/bash

# This script calculates the md5sum for the last number of bits (nbits).
# Usage :
#     ./tail_sum file bnits
#

file_name=$1
bit_num=$2

tail -c $bit_num $file_name | md5sum
