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

def gpuspec(raw, n_win, n_int, f_avg, blank_dc_bin):
	header = raw.read_first_header()
	
	n_spec   = n_win / n_int
	print "Number output spectra: %s" % n_spec
	f_xx_avg = np.zeros((n_spec, header["BLOCSIZE"] / 4 / f_avg), dtype='float32')
	f_yy_avg = np.zeros((n_spec, header["BLOCSIZE"] / 4 / f_avg), dtype='float32')
	f_xy_avg = np.zeros((n_spec, header["BLOCSIZE"] / 4 / f_avg), dtype='complex64')

	fft_plan = None
	df_gpu   = None
	

	t00 = time.time()
	for ii in range(n_spec):
		for jj in range(n_int):
			print "\nIntegration ID:	 %i (%i of %i)" % (ii, jj+1, n_int)

			t1 = time.time()
			header, data = raw.read_next_data_block()
			t2 = time.time()
			print "(Data load:		 %2.2fs)" % (t2 - t1)

			d_xx = data[..., 0]
			d_yy = data[..., 1]

			if not fft_plan:
				t1 = time.time()
				#		  Plan(N_fft,				 input_dtype,  output_dtype, batch_size
				fft_plan = Plan(d_xx.shape[1], np.complex64, np.complex64, d_xx.shape[0])
				t2 = time.time()
				print "FFT Plan:		   %2.2fs" % (t2 - t1)

			t1 = time.time()

			#print d_xx.dtype, d_xx.shape

			# Malloc on GPU
			if not df_gpu:
				df_gpu = gpuarray.empty((d_xx.shape[0], d_xx.shape[1]), np.complex64)

			## XX POL
			t00 = time.time()
			d_gpu  = gpuarray.to_gpu(d_xx)
			t01 = time.time()
			cufft(d_gpu, df_gpu, fft_plan)
			t02 = time.time()
			d_xx_fft  = df_gpu.get()
			t03 = time.time()
			f_xx  = cumultiply(df_gpu, df_gpu.conj()).get().real
			t04 = time.time()
			
			print "    ..memcopy to gpu:     %2.2fs" % (t01 - t00)
			print "    ..cuFFT:              %2.2fs" % (t02 - t01)
			print "    ..memcopy from gpu:   %2.2fs" % (t03 - t02)
			print "    ..cuMultiply:         %2.2fs" % (t04 - t03)
			
			## YY POL
			d_gpu  = gpuarray.to_gpu(d_yy)
			cufft(d_gpu, df_gpu, fft_plan)
			d_yy_fft  = df_gpu.get()
                	f_yy = cumultiply(df_gpu, df_gpu.conj()).get().real

			## XY CROSS POL	
			# Reuse d_gpu
			d_gpu = gpuarray.to_gpu(d_xx_fft)
			f_xy = cumultiply(d_gpu, df_gpu.conj()).get()
		
			t2 = time.time()
			print "cuFFT:			  %2.2fs" % (t2 - t1)


			if blank_dc_bin:
				t1 = time.time()
				f_xx[:, 0] = (f_xx[:, 1] + f_xx[:, -1]) / 2
				f_yy[:, 0] = (f_yy[:, 1] + f_yy[:, -1]) / 2
				f_xy[:, 0] = (f_xy[:, 1] + f_xy[:, -1]) / 2
				t2 = time.time()
				print "DC blanking:              %2.2fs" % (t2 - t1)
		
			t1 = time.time()
			f_xx = np.fft.fftshift(f_xx).ravel()
			f_yy = np.fft.fftshift(f_yy).ravel()
			f_xy = np.fft.fftshift(f_xy).ravel()
			t2 = time.time()

			print "FFT shift:                %2.2fs" % (t2 - t1)
			t1 = time.time()
			fs0 = f_xx.shape[0] / f_avg
			f_xx_avg[ii] += np.roll(f_xx.reshape(fs0, f_avg).mean(axis=1), fs0/2)
			f_yy_avg[ii] += np.roll(f_yy.reshape(fs0, f_avg).mean(axis=1), fs0/2)
			f_xy_avg[ii] += np.roll(f_xy.reshape(fs0, f_avg).mean(axis=1), fs0/2)
			t2 = time.time()
			print "Accumulate:               %2.2fs" % (t2 - t1)



		t01 = time.time()
		print "\nTotal time:		 %2.2fs" % (t01 - t00)


		print f_xx_avg.shape
	return (f_xx_avg, f_yy_avg, f_xy_avg)


if __name__ == "__main__":

	from argparse import ArgumentParser

	parser = ArgumentParser(description="Command line utility for creating \
										 spectra from GuppiRaw files.")

	parser.add_argument('filename', type=str,
						help='Name of file to read')
	parser.add_argument('-n', action='store', default=None, dest='n_win', type=int,
						help='Number of blocks to read from raw file')
	parser.add_argument('-f', action='store', default=None, dest='f_avg', type=int,
						help='Number of channels to average together after FFT')
	parser.add_argument('-b', action='store_false', default=True, dest='blank_dc_bin',
						help='Turn off blanking DC bins of coarse channels')
	parser.add_argument('-N', action='store', default=1, dest='n_int', type=int,
						help='number of integrations per dump')
	parser.add_argument('-w', action='store_true', default=False, dest='plot_waterfall',
						help='Plot waterfall instead of spectrum')
	parser.add_argument('-s', action='store', default='', dest='outfile', type=str,
						help='Save data to file with given filename')
	args = parser.parse_args()

	# Open filterbank data
	filename = args.filename
	raw = GuppiRaw(filename)

	if args.n_win:
		n_win = args.n_win
	else:
		n_win = raw.n_blocks


	if args.f_avg:
		f_avg = args.f_avg
	else:
		f_avg = 1

	blank_dc_bin = args.blank_dc_bin
	n_int = args.n_int

	pprint.pprint(raw.read_header()[0])

	print "Num. blocks: %s " % raw.n_blocks
	print "Num. windows:      %i" % n_win
	print "Num. int per dump: %i" % n_int

	# Compute spectrum from raw file
	(xx, yy, xy) = gpuspec(raw, n_win, n_int, f_avg, blank_dc_bin)
	print xx.shape, yy.shape, xy.shape

	if args.outfile != '':
		print "Saving to %s" % args.outfile
		fil_header = raw.generate_filterbank_header(nchans=xx.shape[0])
		fil_data = np.row_stack((xx, yy, xy.real, xy.imag))#.reshape(2, 1, xx.shape[0])
			
		import h5py
		h5 = h5py.File(args.outfile, 'w')
		h5.create_dataset('data', data=fil_data)
		for key, val in fil_header.items():
			h5.attrs[key] = val
		#fb = Filterbank(filename='test.h5', header_dict=fil_header, data_array=fil_data)
		#fb.write_to_hdf5('test_gpuspec.h5')

	# plot data
	print np.max(xx[0]), np.min(xx[0])
	print "Plotting..."
	if args.plot_waterfall:
		plt.imshow(10*np.log10(xx), aspect='auto', cmap='viridis', interpolation='nearest')
		plt.colorbar()
	else:
		plt.plot(10*np.log10(xx[0]))#,# aspect='auto', cmap='viridis')
		plt.plot(10*np.log10(yy[0]))#, aspect='auto', cmap='viridis')

	plt.show()
