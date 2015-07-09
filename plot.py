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
parser.add_option("-p", "--plateau-index", type="int",
                  dest="plateau_start")

(options, args) = parser.parse_args()

tx1 = sp.fromfile(open("/tmp/tx1.dat"), dtype=sp.complex64)
tx2 = sp.fromfile(open("/tmp/tx2.dat"), dtype=sp.complex64)
rx1 = sp.fromfile(open("/tmp/rx1.dat"), dtype=sp.complex64)
rx2 = sp.fromfile(open("/tmp/rx2.dat"), dtype=sp.complex64)
sco = sp.fromfile(open("/tmp/f_sc_out"), dtype=sp.float32)
s0c = sp.fromfile(open("/tmp/f_s0_corr"), dtype=sp.float32)
ac11 = sp.fromfile(open("/tmp/f_ac11"), dtype=sp.float32)
ac12 = sp.fromfile(open("/tmp/f_ac12"), dtype=sp.float32)
ac13 = sp.fromfile(open("/tmp/f_ac13"), dtype=sp.float32)
ac21 = sp.fromfile(open("/tmp/f_ac21"), dtype=sp.float32)
ac22 = sp.fromfile(open("/tmp/f_ac22"), dtype=sp.float32)
ac23 = sp.fromfile(open("/tmp/f_ac23"), dtype=sp.float32)
f_s0c = open("f_s0c.txt", "w")
f_ac11 = open("f_ac11.txt", "w")
f_ac12 = open("f_ac12.txt", "w")
f_ac13 = open("f_ac13.txt", "w")
f_ac21 = open("f_ac21.txt", "w")
f_ac22 = open("f_ac22.txt", "w")
f_ac23 = open("f_ac23.txt", "w")

if(options.start != None):
  start = options.start
elif(options.plateau_start != None):
  start = options.plateau_start - 200
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

f, axarr = plt.subplots(2, sharex=True)
title = "TX Signal"
axarr[0].set_title(title)
axarr[0].plot([x.real for x in tx1[:txend]], 'r', label="Ch1-Real")
axarr[0].plot([x.imag for x in tx1[:txend]], 'b', label="Ch1-Imag")
axarr[1].plot([x.real for x in tx2[:txend]], 'g', label="Ch2-Real")
axarr[1].plot([x.imag for x in tx2[:txend]], 'y', label="Ch2-Imag")
axarr[0].legend(loc=4)
axarr[1].legend(loc=4)

plt.figure(2)
axes = plt.gca()
title = "RX Signal"
plt.title(title)
plt.plot([x.real for x in rx1[start:end]], 'r', label="Real")
plt.plot([x.imag for x in rx1[start:end]], 'b', label="Imag")
#plt.plot([x.real for x in rx2], 'g', label="Ch2-Real")
#plt.plot([x.imag for x in rx2], 'y', label="Ch2-Imag")
plt.legend(loc=4)

if(options.plateau_start != None):
  plt.figure(3)
  axes = plt.gca()
  title = "SC Algorithm Output"
  plt.title(title)
  plt.plot(sco[options.plateau_start - 2*64:options.plateau_start + 2*64], 'r', label="SC Out")
  plt.legend(loc=4)

if(len(s0c)):
  plt.figure(4)
  axes = plt.gca()
  title = "Correlation"
  plt.title(title)
  plt.plot(s0c, label="S0")
  plt.plot(ac11, label="AC11")
  plt.plot(ac12, label="AC12")
  plt.plot(ac13, label="AC13")
  plt.plot(ac21, label="AC21")
  plt.plot(ac22, label="AC22")
  plt.plot(ac23, label="AC23")
  plt.legend(loc=4)
  for i in range(len(s0c)):
    f_s0c.write(str(s0c[i]) + "\n")
    f_ac11.write(str(ac11[i]) + "\n")
    f_ac12.write(str(ac12[i]) + "\n")
    f_ac13.write(str(ac13[i]) + "\n")
    f_ac21.write(str(ac21[i]) + "\n")
    f_ac22.write(str(ac22[i]) + "\n")
    f_ac23.write(str(ac23[i]) + "\n")

f_s0c.close() 
f_ac11.close()
f_ac11.close()
f_ac13.close()
f_ac21.close()
f_ac22.close()
f_ac23.close()

plt.show()
