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
	import pycuda.cumath as cumath
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
	n_coarse = int(header["OBSNCHAN"])
	blk_size = int(header["BLOCSIZE"])
	f_xx_avg = np.zeros((n_spec, n_coarse, blk_size / 4 / n_coarse / f_avg), dtype='float32')
	f_yy_avg = np.zeros((n_spec, n_coarse, blk_size / 4 / n_coarse / f_avg), dtype='float32')
	f_xy_avg = np.zeros((n_spec, n_coarse, blk_size / 4 / n_coarse / f_avg), dtype='complex64')

	fft_plan_init = False
	df_gpu_init   = False
	d_gpu_init    = False

	d_xx_init = False
	d_yy_init = False


	t000 = time.time()
	for ii in range(n_spec):
		for jj in range(n_int):
			print "\nIntegration ID:	 %i (%i of %i)" % (ii, jj+1, n_int)

			t1 = time.time()
			header, data = raw.read_next_data_block()
			t2 = time.time()
			print "(Data load:		 %2.2fs)" % (t2 - t1)

			# Memcopy to GPU
			if not d_xx_init:
				d_xx_init = True
				d_xx = np.ascontiguousarray(data[..., 0])
				cuda.register_host_memory(d_xx)		# pagelock memory
			else:
				d_xx[:] = data[..., 0]

			if not d_yy_init:
				d_yy_init = True
				d_yy = np.ascontiguousarray(data[..., 1])
				cuda.register_host_memory(d_yy)		# Pagelock memory
			else:
				d_yy[:] = data[..., 1]

			if not fft_plan_init:
				t1 = time.time()
				fft_plan_init = True
				#		  Plan(N_fft,				 input_dtype,  output_dtype, batch_size
				fft_plan = Plan(d_xx.shape[1], np.complex64, np.complex64, d_xx.shape[0])
				t2 = time.time()
				print "FFT Plan:		   %2.2fs" % (t2 - t1)

			t_gpu = time.time()

			#print d_xx.dtype, d_xx.shape

			# Malloc on GPU
			if not df_gpu_init:
				t1 = time.time()
				df_gpu_init = True
				df_gpu = gpuarray.empty((d_xx.shape[0], d_xx.shape[1]), np.complex64)
				df_gpu2 = gpuarray.empty((d_xx.shape[0], d_xx.shape[1]), np.complex64)
				f_xx_gpu_avg = gpuarray.zeros((d_xx.shape[0], d_xx.shape[1]), np.float32)
				f_yy_gpu_avg = gpuarray.zeros((d_xx.shape[0], d_xx.shape[1]), np.float32)
				f_xy_gpu_avg = gpuarray.zeros((d_xx.shape[0], d_xx.shape[1]), np.complex64)
				t2 = time.time()
				print "MALLOC TIME:                %2.2fs" % (t2 - t1)

			## XX POL
			print "XX POL"
			t00 = time.time()
			d_gpu  = gpuarray.to_gpu(d_xx)
			t01 = time.time()
			cufft(d_gpu, df_gpu, fft_plan)
			t02 = time.time()
			t03 = time.time()
			f_xx  = cumultiply(df_gpu, df_gpu.conj()).real
			t04 = time.time()

			print "    ..memcopy to gpu:     %2.2fs" % (t01 - t00)
			print "    ..cuFFT:              %2.2fs" % (t02 - t01)
			#print "    ..memcopy from gpu:   %2.2fs" % (t03 - t02)
			print "    ..cuMultiply:         %2.2fs" % (t04 - t03)

			## YY POL
			d_gpu  = gpuarray.to_gpu(d_yy)
			cufft(d_gpu, df_gpu2, fft_plan)
			f_yy = cumultiply(df_gpu2, df_gpu2.conj()).real

			## XY CROSS POL
			f_xy = cumultiply(df_gpu, df_gpu2.conj())

			t01 = time.time()
			f_xx_gpu_avg += f_xx
			f_yy_gpu_avg += f_yy
			f_xy_gpu_avg += f_xy
			t02 = time.time()
			print "Accumulate:               %2.2fs" % (t02 - t01)
			
			t_gpu2 = time.time()
			print "GPU total:		  %2.2fs" % (t_gpu2 - t_gpu)


		t1 = time.time()	
		f_xx_avg[ii] = f_xx_gpu_avg.get()
		f_yy_avg[ii] = f_yy_gpu_avg.get()
		f_xy_avg[ii] = f_xy_gpu_avg.get()
		t2 = time.time()
		print "GPU -> CPU memcopy:       %2.2fs" % (t2 - t1)
		t01 = time.time()
		print "BLOCK TIME:		 %2.2fs" % (t01 - t00)


	print "\n### CPU POST PROCESSING ###"
	if blank_dc_bin:
		t1 = time.time()
		f_xx_avg[:, :, 0] = (f_xx_avg[:, :, 1] + f_xx_avg[:, :, -1]) / 2
		f_yy_avg[:, :, 0] = (f_yy_avg[:, :, 1] + f_yy_avg[:, :, -1]) / 2
		f_xy_avg[:, :, 0] = (f_xy_avg[:, :, 1] + f_xy_avg[:, :, -1]) / 2
		t2 = time.time()
		print "DC blanking:              %2.2fs" % (t2 - t1)
			
	t1 = time.time()
	s0, s1, s2 = f_xx_avg.shape
	f_xx_avg = np.fft.fftshift(f_xx_avg, axes=2).reshape(s0, s1 * s2)
	f_yy_avg = np.fft.fftshift(f_yy_avg, axes=2).reshape(s0, s1 * s2)
	f_xy_avg = np.fft.fftshift(f_xy_avg, axes=2).reshape(s0, s1 * s2)
	t2 = time.time()
	print "FFT shift:                %2.2fs" % (t2 - t1)

	t001 = time.time()
	print "\nTotal gpuspec time:       %2.2fs" % (t001 - t000)
	return (f_xx_avg, f_yy_avg, f_xy_avg)


if __name__ == "__main__":

	from argparse import ArgumentParser

	parser = ArgumentParser(description="Command line utility for creating spectra from GuppiRaw files.")

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
	#print xx.shape, yy.shape, xy.shape

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
	print "Plotting..."
	if args.plot_waterfall:
		plt.figure('xx')
		plt.imshow(10*np.log10(xx), aspect='auto', cmap='viridis', interpolation='nearest')
		plt.colorbar()
		plt.figure('yy')
		plt.imshow(10*np.log10(yy), aspect='auto', cmap='viridis', interpolation='nearest')
		plt.colorbar()
	else:
		if xx[0].shape[0] > 1e6:
			xx = xx[0].reshape(xx[0].shape[0] / 1024, 1024).mean(axis=1)
			yy = yy[0].reshape(yy[0].shape[0] / 1024, 1024).mean(axis=1)
		else:
			xx, yy = xx[0], yy[0]
		plt.plot(10*np.log10(xx))#,# aspect='auto', cmap='viridis')
		plt.plot(10*np.log10(yy))#, aspect='auto', cmap='viridis')

	plt.show()
