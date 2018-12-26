#!/bin/bash

#
# Simple demo simulation for reaction diffusion model
#

# Perform time integration and write data to data.bin

../build/turing --Da 1 --Db 100 --dx 1.0 --dt 0.001 --width 100 --height 100 \
--steps 20 --tsteps 1000 --outfile "data.bin" --reaction fitzhugh-nagumo \
--parameters "alpha=-0.005;beta=10.0"

# Execute Python script to generate graphs
python vis.py data.bin -0.8 0.8 -0.35 0.35

# create a gif file from the series of png files
convert {0000..0020}.png -resize 25% -loop 0 fitzhugh-nagumochr.gif
