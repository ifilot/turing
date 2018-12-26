#!/bin/bash

#
# Simple demo simulation for reaction diffusion model
#

# Perform time integration and write data to data.bin
../build/turing --Da 2e-5 --Db 1e-5 --dx 0.005 --dt 0.01 --width 256 --height 256 \
--steps 100 --tsteps 100 --outfile "data.bin" --reaction lotka-volterra \
--parameters "alpha=2.3333;beta=2.6666;gamma=1.0;delta=1.0"

# Execute Python script to generate graphs
python vis.py

# create a gif file from the series of png files
convert {0080..0100}.png -resize 25% -loop 0 lotka-volterra.gif
