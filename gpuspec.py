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

import h5py

try:
	import pycuda.driver as cuda
	import pycuda.autoinit
	import pycuda.gpuarray as gpuarray
except ImportError:
	raise RuntimeError("could not load pycuda. Please check your install.")

try:
	from skcuda.fft import Plan
	from skcuda.fft import fft as cufft
	from skcuda.misc import multiply as cumultiply
except ImportError:
	raise RuntimeError("Could not load scikit-cuda. Please check your install")

import pprint
import pylab as plt
import time

from guppi import GuppiRaw
from filterbank import Filterbank

def gpuspec(raw, n_int, f_avg, blank_dc_bin):
	header = raw.read_first_header()
	f_xx_avg = np.zeros(header["BLOCSIZE"] / dec_fac / 4)
	f_yy_avg = np.zeros(header["BLOCSIZE"] / dec_fac / 4)
	f_xy_avg = np.zeros(header["BLOCSIZE"] / dec_fac / 4, dtype='complex64')

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
		f_xx  = cumultiply(df_gpu, df_gpu).view('float32')[::2].get() 		

		## YY POL
		d_gpu  = gpuarray.to_gpu(d_yy)
		cufft(d_gpu, df_gpu, fft_plan)
		d_yy_fft  = df_gpu.get()
                f_yy = cumultiply(df_gpu, df_gpu).view('float32')[::2].get()

		## XY CROSS POL	
		# Reuse d_gpu
		d_gpu = gpuarray.to_gpu(d_xx_fft)
		f_xy = cumultiply(d_gpu, df_gpu.conj()).get()
		
		t2 = time.time()
		print "cuFFT:			  %2.2fs" % (t2 - t1)

		t1 = time.time()

	        	
		#f_xx = np.abs(d_xx_fft)
		#f_yy = np.abs(d_yy_fft)
		#f_xy = d_xx_fft * d_yy_fft.conj()
		#print np.allclose(f_xy, f_xy2)

		if blank_dc_bin:
			print "BLANKING"
			f_xx[:, 0] = (f_xx[:, 1] + f_xx[:, -1]) / 2
			f_yy[:, 0] = (f_yy[:, 1] + f_yy[:, -1]) / 2
			f_xy[:, 0] = (f_xy[:, 1] + f_xy[:, -1]) / 2

		f_xx = np.fft.fftshift(f_xx)
		f_xx_avg += f_xx.flatten()

		f_yy = np.fft.fftshift(f_yy)
		f_yy_avg += f_yy.flatten()

		f_xy = np.fft.fftshift(f_xy)
		f_xy_avg += f_xy.flatten()

		t2 = time.time()
		print "Square + integrate: %2.2fs" % (t2 - t1)

	t01 = time.time()
	print "\nTotal time:		 %2.2fs" % (t01 - t00)

	f_xx_avg = f_xx_avg.reshape(f_xx_avg.shape[0] / f_avg, f_avg).mean(axis=1)
	f_yy_avg = f_yy_avg.reshape(f_yy_avg.shape[0] / f_avg, f_avg).mean(axis=1)
	f_xy_avg = f_xy_avg.reshape(f_xy_avg.shape[0] / f_avg, f_avg).mean(axis=1)

	# A roll is required for some reason, probably to do with the FFTshift
	f_xx_avg = np.roll(f_xx_avg, f_xx_avg.shape[0]/2)
	f_yy_avg = np.roll(f_yy_avg, f_yy_avg.shape[0]/2)
	f_xy_avg = np.roll(f_xy_avg, f_xy_avg.shape[0]/2)
	return (f_xx_avg, f_yy_avg, f_xy_avg)


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
	(xx, yy, xy) = gpuspec(raw, n_int, f_avg, blank_dc_bin)
	print xx.shape, yy.shape, xy.shape

	fil_header = raw.generate_filterbank_header(nchans=xx.shape[0])
	fil_data = np.row_stack((xx, yy)).reshape(2, 1, xx.shape[0])

	#fb = Filterbank(filename='test.h5', header_dict=fil_header, data_array=fil_data)
	#fb.write_to_hdf5('test_gpuspec.h5')

	# plot data
	print "Plotting..."
	plt.plot(10*np.log10(xx))
	plt.plot(10*np.log10(yy))

	plt.show()
