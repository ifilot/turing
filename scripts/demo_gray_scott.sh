#!/bin/bash

#
# Simple demo simulation for reaction diffusion model
#

# Perform time integration and write data to data.bin

../build/turing --Da 2e-5 --Db 1e-5 --dx 0.005 --dt 0.1 --width 256 --height 256 \
--steps 20 --tsteps 1000 --outfile "data.bin" --reaction gray-scott \
--parameters "f=0.06;k=0.0609"

# Execute Python script to generate graphs
python vis.py data.bin 0.0 1.0 0.0 1.0

# create a gif file from the series of png files
convert {0000..0015}.png -resize 25% -loop 0 gray-scott.gif
