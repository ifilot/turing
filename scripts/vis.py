#!/usr/bin/env python3

import struct
import numpy as np
import matplotlib.pyplot as plt
import sys

# set boundaries for cmap
vmin1 = float(sys.argv[2])
vmax1 = float(sys.argv[3])
vmin2 = float(sys.argv[4])
vmax2 = float(sys.argv[5])

with open(sys.argv[1], "rb") as f:
    width = struct.unpack('i', f.read(4))[0]
    height = struct.unpack('i', f.read(4))[0]
    steps = struct.unpack('i', f.read(4))[0]

    for i in range(0, steps+1):
        a = np.fromfile(f, dtype=np.dtype('d'), count=width * height)
        b = np.fromfile(f, dtype=np.dtype('d'), count=width * height)

        ap = a.reshape((height, width))
        bp = b.reshape((height, width))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,4))

        im1 = ax1.imshow(ap, origin='lower', interpolation='bicubic', vmin=vmin1, vmax=vmax1)
        plt.colorbar(im1, ax=ax1)

        im2 = ax2.imshow(bp, origin='lower', interpolation='bicubic', vmin=vmin2, vmax=vmax2, cmap='PiYG')
        plt.colorbar(im2, ax=ax2)

        ax1.set_title('Concentration A')
        ax2.set_title('Concentration B')
        filename = '%04i.png' % i
        print("Writing image: %s" % filename)
        plt.savefig(filename, dpi=72)
        plt.close()
