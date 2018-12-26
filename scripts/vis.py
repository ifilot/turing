#!/usr/bin/env python3

import struct
import numpy as np
import matplotlib.pyplot as plt

with open("data.bin", "rb") as f:
    width = struct.unpack('i', f.read(4))[0]
    height = struct.unpack('i', f.read(4))[0]
    steps = struct.unpack('i', f.read(4))[0]

    for i in range(0, steps+1):
        a = np.fromfile(f, dtype=np.dtype('d'), count=width * height)
        b = np.fromfile(f, dtype=np.dtype('d'), count=width * height)

        # if i < steps - 1:
        #     continue

        ap = a.reshape((height, width))
        bp = b.reshape((height, width))

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,4))

        im1 = ax1.imshow(ap, origin='lower', interpolation='bicubic', vmin=0, vmax=1)
	plt.colorbar(im1, ax=ax1)

        im2 = ax2.imshow(bp, origin='lower', interpolation='bicubic', vmin=0, vmax=1, cmap='PiYG')
	plt.colorbar(im2, ax=ax2)

        ax1.set_title('Concentration A')
        ax2.set_title('Concentration B')
        filename = '%04i.png' % i
        print("Writing image: %s" % filename)
        plt.savefig(filename, dpi=72)
        plt.close()
