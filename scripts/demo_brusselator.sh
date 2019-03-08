#!/bin/bash

#
# Simple demo simulation for reaction diffusion model
#

# Perform time integration and write data to data.bin

../build/turing --Da 2 --Db 16 --dx 1.0 --dt 0.005 --width 256 --height 256 \
--steps 20 --tsteps 1000 --outfile "data.bin" --reaction brusselator \
--parameters "alpha=4.5;beta=7.50" --pbc

# Execute Python script to generate graphs
python vis.py data.bin 2 6.0 1 2.0

# create a gif file from the series of png files
convert {0000..0020}.png -resize 25% -loop 0 brusselator.gif
