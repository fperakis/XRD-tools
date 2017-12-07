import numpy as np
import sys
import os
import matplotlib.pyplot as plt
def readcbf(fname,nrows,ncols):
        f = open(fname, 'r')
        while True:
                c = np.fromfile(f, dtype='uint8', count=1)
                if c == 12:
                        c = np.fromfile(f, dtype='uint8', count=1)
                        if c == 26:
                                c = np.fromfile(f, dtype='uint8', count=1)
                                if c == 4:
                                        c = np.fromfile(f, dtype='uint8', count=1)
                                        if c == 213:
                                                break
        data = np.zeros(nrows*ncols, dtype='int32')
        j = 1
        total = 0
        maxval = -100
        minval = 10000
        while True:
                c = np.fromfile(f, dtype='int8', count=1)
                if c == -128:
                        c = np.fromfile(f, dtype='int16', count=1)
                        if c == -32768:
                                c = np.fromfile(f, dtype='int32', count=1)
                data[j] = data[j-1] + c[0]
                total += data[j]
                if data[j] > maxval:
                        maxval = data[j]
                if data[j] < minval:
                        minval = data[j]
                j += 1
                if j == nrows*ncols:
                        break
        f.close()
        #avg=total / (nrows*1.0) / ncols
        #print "## ", avg, maxval, minval
        #        print data.shape
        return data#(data, avg)

