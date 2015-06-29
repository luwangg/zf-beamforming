#!/usr/bin/env python

import scipy as sp
import matplotlib.pyplot as plt
from cmath import phase, sqrt
from math import pi
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--sync-index", type="int",
                  dest="start")
parser.add_option("-n", "--num-samples", type="int",
                  dest="samples")
parser.add_option("-N", "--num-tx-samples", type="int",
                  dest="txsamples")

(options, args) = parser.parse_args()

tx1 = sp.fromfile(open("/tmp/tx1.dat"), dtype=sp.complex64)
tx2 = sp.fromfile(open("/tmp/tx2.dat"), dtype=sp.complex64)
rx1 = sp.fromfile(open("/tmp/rx1.dat"), dtype=sp.complex64)
rx2 = sp.fromfile(open("/tmp/rx2.dat"), dtype=sp.complex64)

if(options.start != None):
  start = options.start
else:
  start = 0

if(options.samples != None):
  end = start + options.samples
else:
  end = len(rx1)

if(options.txsamples != None):
  txend = 0 + options.txsamples
else:
  txend = len(tx1)

plt.figure(1)
axes = plt.gca()
title = "TX Signal"
plt.title(title)
plt.plot([x.real for x in tx1[:txend]], 'r', label="Ch1-Real")
plt.plot([x.imag for x in tx1[:txend]], 'b', label="Ch1-Imag")
plt.plot([x.real for x in tx2[:txend]], 'g', label="Ch2-Real")
plt.plot([x.imag for x in tx2[:txend]], 'y', label="Ch2-Imag")
plt.legend(loc=4)

plt.figure(2)
axes = plt.gca()
title = "RX Signal"
plt.title(title)
plt.plot([x.real for x in rx1[start:end]], 'r', label="Ch1-Real")
plt.plot([x.imag for x in rx1[start:end]], 'b', label="Ch1-Imag")
#plt.plot([x.real for x in rx2], 'g', label="Ch2-Real")
#plt.plot([x.imag for x in rx2], 'y', label="Ch2-Imag")
plt.legend(loc=4)
plt.show()
