#!/bin/bash

time ../build/turing --Da 2e-5 --Db 1e-5 --dx 0.005 --dt 0.1 --width 256 --height 256 --steps 1000 --tsteps 1 --outfile "data.bin"
python vis.py
