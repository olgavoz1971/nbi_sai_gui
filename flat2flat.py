#!/usr/bin/env python

import sys
if len(sys.argv) != 3:
  sys.exit("Usage: ./flat2flat.py file1.fits file2.fits")

fn1 = sys.argv[1]
fn2 = sys.argv[2]

import pysao
import pyfits

d=pysao.ds9(wait_time=180)
d.set('scale histequ')
d.set('scale mode zscale')
d.view(pyfits.getdata(fn1).astype(float)/pyfits.getdata(fn2).astype(float))
raw_input("Press Enter to continue")
