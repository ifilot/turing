#!/bin/bash

#
# Simple demo simulation for reaction diffusion model
#

# Perform time integration and write data to data.bin

../build/turing --Da 5.0 --Db 0.0 --dx 1.0 --dt 0.001 --width 128 --height 128 \
--steps 20 --tsteps 10000 --outfile "data.bin" --reaction barkley \
--parameters "alpha=0.75;beta=0.06;epsilon=13.0"

# Execute Python script to generate graphs
python vis.py data.bin 0.0 1.0 0.0 0.5

# create a gif file from the series of png files
convert {0000..0020}.png -resize 25% -loop 0 barkley.gif
