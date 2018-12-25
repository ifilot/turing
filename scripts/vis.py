#!/usr/bin/env python3

import struct
import numpy as np
import matplotlib.pyplot as plt

with open("data.bin", "rb") as f:
    width = struct.unpack('i', f.read(4))[0]
    height = struct.unpack('i', f.read(4))[0]
    steps = struct.unpack('i', f.read(4))[0]

    print(width, height, steps)

    for i in range(0, steps):
        a = np.fromfile(f, dtype=np.dtype('d'), count=width * height)
        b = np.fromfile(f, dtype=np.dtype('d'), count=width * height)

        ap = a.reshape((height, width))
        bp = b.reshape((height, width))

        plt.figure()
        plt.imshow(ap, origin='lower', interpolation='bicubic')
        plt.title('%03i' % i)
        filename = '%03i.png' % i
        print("Writing image: %s" % filename)
        plt.savefig(filename)
        plt.close()
