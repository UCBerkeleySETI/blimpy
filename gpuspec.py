#!/usr/bin/env python
"""
# gpuspec.py

Python script to produce spectra from GuppiRaw data files.

This script requires PyCuda and scikit-cuda to be installed and a GPU card
to go along with it.

"""

import numpy as np
import os
import sys

try:
	import pycuda.driver as cuda
	import pycuda.autoinit
	import pycuda.gpuarray as gpuarray
except ImportError:
	raise RuntimeError("could not load pycuda. Please check your install.")

try:
	from skcuda.fft import Plan
	from skcuda.fft import fft as cufft
except ImportError:
	raise RuntimeError("Could not load scikit-cuda. Please check your install")

import pprint
import pylab as plt
import time

from guppi import GuppiRaw

def gpuspec(raw, n_int, f_avg, blank_dc_bin):
	header = raw.read_first_header()
	d_avg = np.zeros(header["BLOCSIZE"] / dec_fac / 4)

	fft_plan = None
	df_gpu   = None

	t00 = time.time()
	for ii in range(n_int):
		print "\nIntegration ID:	 %i" % ii

		t1 = time.time()
		header, data = raw.read_next_data_block()
		t2 = time.time()
		print "(Data load:		 %2.2fs)" % (t2 - t1)

		d_xx = data[..., 0]
		d_yy = data[..., 1]

		if not fft_plan:
			t1 = time.time()
			#		  Plan(N_fft,				 input_dtype,  output_dtype, batch_size
			fft_plan = Plan(d_xx.shape[1]/dec_fac, np.complex64, np.complex64, d_xx.shape[0])
			t2 = time.time()
			print "FFT Plan:		   %2.2fs" % (t2 - t1)

		t1 = time.time()

		#print d_xx.dtype, d_xx.shape

		# Malloc on GPU
		if not df_gpu:
			print d_xx.shape[0], d_xx.shape[1]/dec_fac
			df_gpu = gpuarray.empty((d_xx.shape[0], d_xx.shape[1]/dec_fac), np.complex64)

		## XX POL
		d_gpu  = gpuarray.to_gpu(d_xx)
		cufft(d_gpu, df_gpu, fft_plan)
		d_xx_fft  = df_gpu.get()

		## YY POL
		d_gpu  = gpuarray.to_gpu(d_yy)
		cufft(d_gpu, df_gpu, fft_plan)
		d_yy_fft  = df_gpu.get()

		t2 = time.time()
		print "cuFFT:			  %2.2fs" % (t2 - t1)

		t1 = time.time()
		d_pow = np.abs(d_xx_fft)  + np.abs(d_yy_fft)

		if blank_dc_bin:
			d_pow[:, 0] = (d_pow[:, 1] + d_pow[:, -1]) / 2   # Blank DC bin

		d_pow = np.fft.fftshift(d_pow)
		d_avg += d_pow.flatten()

		t2 = time.time()
		print "Square + integrate: %2.2fs" % (t2 - t1)

	t01 = time.time()
	print "\nTotal time:		 %2.2fs" % (t01 - t00)

	d_avg = d_avg.reshape(d_avg.shape[0] / f_avg, f_avg).mean(axis=1)

	# A roll is required for some reason, probably to do with the FFTshift
	d_avg = np.roll(d_avg, d_avg.shape[0]/2)
	return d_avg

if __name__ == "__main__":

	from argparse import ArgumentParser

	parser = ArgumentParser(description="Command line utility for creating \
										 spectra from GuppiRaw files.")

	parser.add_argument('filename', type=str,
						help='Name of file to read')
	parser.add_argument('-n', action='store', default=None, dest='n_int', type=int,
						help='Number of blocks to read from raw file')
	parser.add_argument('-f', action='store', default=None, dest='f_avg', type=int,
						help='Number of channels to average together after FFT')
	parser.add_argument('-b', action='store_false', default=True, dest='blank_dc_bin',
						help='Turn off blanking DC bins of coarse channels')

	args = parser.parse_args()

	# Open filterbank data
	filename = args.filename
	raw = GuppiRaw(filename)

	if args.n_int:
		n_int = args.n_int
	else:
		n_int = raw.n_blocks

	if args.f_avg:
		f_avg = args.f_avg
	else:
		f_avg = 1

	blank_dc_bin = args.blank_dc_bin

	pprint.pprint(raw.read_header()[0])

	print "Num. blocks: %s " % raw.n_blocks

	dec_fac = 1  		# This may be bad to use

	print "Decimation factor: %i" % dec_fac
	print "Num. integrations: %i" % n_int

	# Compute spectrum from raw file
	d_spec = gpuspec(raw, n_int, f_avg, blank_dc_bin)

	# plot data
	print "Plotting..."
	plt.plot(10*np.log10(d_spec))
	plt.show()
