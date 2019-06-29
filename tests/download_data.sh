#!/usr/bin/env bash
apt-get update; apt-get install curl
echo "Works!"
mkdir data
curl --url "http://blpd0.ssl.berkeley.edu/Voyager_data/Voyager1_block1.npy" -o ./data/Voyager1_block1.npy
curl --url "http://blpd0.ssl.berkeley.edu/Voyager_data/blc3_2bit_guppi_57396_VOYAGER1_0006.0013.raw" -o ./data/blc3_2bit_guppi_57396_VOYAGER1_0006.0013.raw
curl --url "http://blpd0.ssl.berkeley.edu/Voyager_data/Voyager1.single_coarse.fine_res.fil" -o ./data/Voyager1.single_coarse.fine_res.fil
curl --url "http://blpd0.ssl.berkeley.edu/Voyager_data/Voyager1.single_coarse.fine_res.h5" -o ./data/Voyager1.single_coarse.fine_res.h5
ls ./data