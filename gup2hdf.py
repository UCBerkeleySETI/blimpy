from guppi import GuppiRaw
import h5py
import bitshuffle.h5
import time
import os
import glob

parser = ArgumentParser(description="Command line utility for creating HDF5 Raw files.")
parser.add_argument('dirname', type=str, help='Name of directory to read')
args = parser.parse_args()

filelist = glob.glob(os.path.join(args.dirname, '*.raw'))

for filename in filelist:
    if not os.path.exists(filename + '.h5'):
        
        t0 = time.time()
        print "\nReading %s header..." % filename
        r = GuppiRaw(filename)
        h5 = h5py.File(filename + '.h5', 'w')
        h5.attrs['CLASS'] = 'GUPPIRAW'

        for ii in range(r.n_blocks):
            header, data = r.read_next_data_block()

            print "Creating block %i of %i" % (ii+1, r.n_blocks)
            
            block_size = 0      # This is chunk block size
            dset = h5.create_dataset('block_%02i' % ii,
                          data=data,
                          compression=bitshuffle.h5.H5FILTER,
                          compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
                          dtype=data.dtype)

            # Copy over header information as attributes
            for key, value in header.items():
                dset.attrs[key] = value

        h5.close()

        t1 = time.time()
        print "Conversion time: %2.2fs" % (t1- t0)
