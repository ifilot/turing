#!/usr/bin/env python3

import struct
import numpy as np
import matplotlib.pyplot as plt

with open("data.bin", "rb") as f:
    width = struct.unpack('i', f.read(4))[0]
    height = struct.unpack('i', f.read(4))[0]
    steps = struct.unpack('i', f.read(4))[0]

    print(width, height, steps)

    for i in range(0, steps+1):
        a = np.fromfile(f, dtype=np.dtype('d'), count=width * height)
        b = np.fromfile(f, dtype=np.dtype('d'), count=width * height)

        # if i < steps - 1:
        #     continue

        ap = a.reshape((height, width))
        bp = b.reshape((height, width))

        plt.figure(dpi=240)
        im = plt.imshow(ap, origin='lower', interpolation='nearest', vmin=0, vmax=1)
        plt.colorbar(im)
        plt.title('%04i' % i)
        filename = '%04i.png' % i
        print("Writing image: %s" % filename)
        plt.savefig(filename)
        plt.close()
