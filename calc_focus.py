#!/usr/bin/python3

# Script precesses focus-sequences image. The first and last star image must me marked on ds9 window. Uses sex ldac files.
# Invoked from foc_server.tcl

from jkclass import DS9_NAME

import os
import sys
import argparse
import numpy as np
import math
from skimage.measure import moments, moments_central

from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.stats import mad_std

import pyds9
from pyds9 import DS9 as ds9
from pygrace import grace
import pyregion
import prf_module as prf

def round2int(s):
  try:
    num = int(np.rint(float(s)))
  except ValueError:
    num = -1
  return num  

def round2float(s, ndigits):
#  if (type(s) is float) or (type(s) is int):
    try:
      num = round(float(s), ndigits)
    except ValueError:
      num = -1
    return num
#  else: return -1
    
def circ(x, y, xc, yc, rad):
    return (x-xc)**2 + (y-yc)**2 >= rad**2

#import matplotlib.pyplot as plt
#from astropy.visualization import SqrtStretch, simple_norm

def calc_fwhm(xx, yy, rad, img):
  impix = img[yy-rad:yy+rad,xx-rad:xx+rad].copy()
  image_centre_coord, fwhm_empir, com, r, fl = \
       prf.profile(impix, xx-rad, yy-rad, 10, prf.center_method_list[0])
  ind = r.argsort()
  xsort = r[ind]
  ysort = fl[ind]
  xrevers = xsort[: :-1]
  yrevers = ysort[: :-1]
  xdouble = np.hstack((-1.*xrevers,xsort))
  ydouble = np.hstack((yrevers,ysort))
 
  fitmodel = prf.fitmodel_default
  fwhm_model, gfit = prf.fitmodel_dict[fitmodel](xdouble, ydouble, 10000., 0., 7.)
  return fwhm_empir, fwhm_model

def calcmom(xx, yy, rad, img):
  impix = img[yy-rad:yy+rad,xx-rad:xx+rad].copy()

  x,y = np.indices(impix.shape)
  circ_mask = circ(x, y, rad-0.5, rad-0.5, rad)  # around new centre
  
  impix_ravel = impix.ravel()
  idx = int(round(0.1 * len(impix_ravel)))  # 10%
  sorted_idx = np.argsort(impix_ravel)
  i = sorted_idx[idx]
  backgr = impix_ravel[i] # 10% quantile of pixel distribution
#  print('Background =', backgr)
  impix -= backgr
  impix[circ_mask] = 0

  thresh = impix.max() * 0.05
  impix[impix <= thresh] = 0

#  norm = simple_norm(impix, 'sqrt', percent=99)
#  plt.imshow(impix, norm=norm, origin='lower', cmap='Greys_r',
#           interpolation='nearest')
#  plt.show()

  Mc = moments_central(impix)
  M =  moments(impix)
  sum_ = np.sum(impix)
  mxx = Mc[0,2]/Mc[0,0]
  myy = Mc[2,0]/Mc[0,0]

  xc = M[0,1]/M[0,0] + xx - rad
  yc = M[1,0]/M[0,0] + yy - rad

  return xc,yc,mxx,myy

def find_nearest_data(point, x, y, param):
# find nearest from sextractor data
   r2min = -1
   point_nearest = (-1,-1,-1)
   for p in zip(x, y, param):
     r2 = (point[0]-p[0])**2 + (point[1]-p[1])**2
     if (r2min < 0) or (r2 < r2min):
       r2min = r2
       point_nearest = p
#     print(p, point, r2min, r2)
#     print(point_nearest)
   return point_nearest

# Calculate minimum
#   foci xc yc mxx myy

def calc_minimum(x, y):
#  print('\n',x,y)
  poly = np.poly1d(np.polyfit(x,y, args.degree))
  roots = poly.deriv(1).r
#  max or min?
  test_poly = poly.deriv(2)(roots)
  x_min_poly = roots[test_poly>0]
  x_max_poly = roots[test_poly<0]
  
  if len(x_min_poly) < 1.0: 
    focus = -1
  else:
    focus = round2float(x_min_poly[0], 2)
  return focus, poly

def draw(g, x, y, poly):
  g.restart()
  g.x = x
  g.y = y
  g('s0 symbol color "green"')
  g('s0 line type 0')
  g('s0 symbol 1')
  g('s0 symbol fill color "green"')
  g('s0 symbol fill pattern 1')
  g('s0 symbol size 0.5')
  g('plot(x,y)')

#  x.sort()
  g.hold(1)
  nsteps = 100.0
  step = (x[-1] - x[0]) / nsteps
  x_detailed = np.arange(x[0],x[-1]+step,step)
  g.plot(x_detailed, poly(x_detailed))
#  input("Press Enter to continue...")
  return

# Mark foci on ds9
#def mark_on_ds9(data, d, x_col, y_col, text_col):
def mark_on_ds9(d, x, y, text):
#  for row in data[:]:
  for i in range(len(x)):   
    foci_reg = 'physical; point('+str(x[i])+','+str(y[i])+') # text={'+str(round(text[i],2))+'} point=cross color=red'
    d.set('regions', foci_reg)

def get_region(d):
  x0,y0,rad = 0,0,0
  r = d.get('region')
  try:
   regions_all = pyregion.parse(r)
  except ValueError as e:
   input('Please, select first and last star image and press Enter')
   return -1,-1,-1

# Clean regins (I use only circles)
  regions = []
  for reg in regions_all:
    if reg.name == 'circle':
     regions.append(reg)

  num_reg = len(regions)
  if num_reg < 2:
    key = input('Please, select two circle regions and press Enter')
    return -1,-1,-1
  
# Use only two last circles
  xc_beg = regions[num_reg-2].coord_list[0]
  yc_beg = regions[num_reg-2].coord_list[1]
  xc_end = regions[num_reg-1].coord_list[0]
  yc_end = regions[num_reg-1].coord_list[1]

# Tolerance is the biggest circle radius
  tol_rad = max(regions[num_reg-2].coord_list[2], regions[num_reg-1].coord_list[2])
   
  return (xc_beg, yc_beg), (xc_end, yc_end), tol_rad

# rotate ortobasis:
def fxy(x,y,alp):
  x1 = x * math.cos(alp) - y * math.sin(alp)
  y1 = 1.0 * x * math.sin(alp) + y * math.cos(alp)
  return x1,y1

# --- MAIN: -------------------------------------
def main_job(args):

  ds9_list = pyds9.ds9_targets()
  if ds9_list is None:
    sys.exit("Please, open ds9 and select first and last star image of the focus secuence")
  if(len(ds9_list) > 1):
    d = ds9('DS9:'+DS9_NAME)
  else:
    d = ds9()

#  hdulist = fits.open(args.infile)

  fits_image_full_name = args.infile
  hdulist = fits.open(fits_image_full_name, uint=False)  # working with float data
  nhdu = len(hdulist)
  img = hdulist[nhdu-1].data

# TMP
  hdusexlist = fits.open(args.catfile)
  nhdusex    = len(hdusexlist)
  print('nhdusex =', nhdusex)
  data2sex   = hdusexlist[nhdusex-1].data
  headsex    = hdusexlist[0].header
  xsex       = data2sex['X_IMAGE']
  ysex       = data2sex['Y_IMAGE']
# Second moments:
#  x2sex      = data2sex['X2_IMAGE']
#  y2sex      = data2sex['Y2_IMAGE']
#  xysex      = data2sex['XY_IMAGE']
  fwhmsex    = data2sex['FWHM_IMAGE']*args.scale # pix to srgsec
# ----

  g = grace()

  doit = True
  while doit:

    moments_all = []
    foci = np.linspace(args.focbeg+(args.focnum-1)*args.focstep, args.focbeg, args.focnum)
    print('foci =', foci)

    point_beg, point_end, tol_rad = get_region(d)
    print('point_beg, point_end, tol_rad =', point_beg, point_end, tol_rad)
    if tol_rad < 0: break
#   Calculate rotaion angle (Altitude (along focus-sequence) to matrix-column)
    alp = math.atan2((point_end[0]-point_beg[0]),(point_end[1]-point_beg[1]))

# !!! calc baricenter here to correct x,y
    xb,yb = point_beg
    xe,ye = point_end
    xb, yb, mxx, myy = calcmom(round(xb), round(yb), round(tol_rad), img)
    xe, ye, mxx, myy = calcmom(round(xe), round(ye), round(tol_rad), img)

    ymax = yb > ye and yb or ye
    ymin = yb < ye and yb or ye
    xmax = xb > xe and xb or xe
    xmin = xb < xe and xb or xe

# Check foci direction:
# Check dero-angle, search upper/bottom side
    deroa = Angle(args.dero*u.deg)
    pos0=0 # !!! TMP Instang?
    zenith = Angle(pos0*u.deg) + deroa
    alpa = Angle(alp*u.rad)
    angle_difference = abs((alpa-zenith).wrap_at(180*u.deg).degree)
    if angle_difference > 90:
      print('reverse foci, dangle =',angle_difference)
      foci = foci[::-1] 
      print(foci)

    step_x = (xe-xb) / (len(foci) - 1)
    step_y = (ye-yb) / (len(foci) - 1)
#    print('Step x,y =', step_x, step_y)

    imsize_array = np.zeros([len(foci)])
    x_array  = np.zeros([len(foci)])
    y_array  = np.zeros([len(foci)])

#    rad = 32 
    rad = round(tol_rad)
    for step in range(0, len(foci)):
      xc = xb + step_x * step
      yc = yb + step_y * step
      xc_corr, yc_corr, mxx, myy = calcmom(int(round(xc)), int(round(yc)), rad, img)

      scale = 2.385875519959664 # for a4
# if moments
      a4 = scale * np.sqrt( (mxx + myy) * 0.5 ) / 1.102
# if psf:
      fwhm_empir_px, fwhm_model_px  = calc_fwhm(int(round(xc)), int(round(yc)), rad, img)
      fwhm_empir = fwhm_empir_px*args.scale
      fwhm_model = fwhm_model_px*args.scale
      xc_sex,yc_sex,fwhm_sex = find_nearest_data((xc_corr,yc_corr), xsex, ysex, fwhmsex)
      print('focus: %7.3f %7.3f %7.3f %7.3f %7.3f' % (foci[step], a4, fwhm_empir, fwhm_model, fwhm_sex))
      imsize_array[step] = fwhm_empir
      x_array[step] = xc_corr
      y_array[step] = yc_corr

    focus, poly = calc_minimum(foci, imsize_array )
    mark_on_ds9(d, x_array, y_array, foci) # Mark foci on ds9
    draw(g, foci, imsize_array, poly)

    print('\n\nfocus, poly =', focus, poly)


    key = input('Enter Q(uit) or any other key to continue ...')
    print(key)
    if key.upper() == 'Q': break

#raw_input("Press 'Return' to continue...")    #python 2


parser = argparse.ArgumentParser()
parser.add_argument('-f','--infile', required=True, nargs='?', default='multifocus.fits', help='input fits-file')
parser.add_argument('-r','--dero',   nargs='?', default=180.0, type=float, help='DERO-angle in degrees')
parser.add_argument('-d','--degree', nargs='?', default=3, type=int, help='degree of polynomials')
parser.add_argument('-p','--param',  nargs='?', default='FWHM_IMAGE', help='parameter (fwhm_image, elongation, x2_image etc)')
parser.add_argument('--focbeg',  nargs='?', required=True, type=float, help='first focus position')
parser.add_argument('--focstep', nargs='?', required=True, type=float, help='focus step')

parser.add_argument('--focnum',  nargs='?', required=True, type=int,   help='number of focus positions')
parser.add_argument('--enlargement', nargs='?', default=True, type=bool,  help='Is the focus increasing?')
parser.add_argument('--scale', nargs='?', default=0.155, type=float,  help='Image scale in \"/px')
# !!! TMP  --------
parser.add_argument('-c','--catfile', required=False, nargs='?', default='work_ldac.cat', help='input sex-ldac-file')
args = parser.parse_args()

main_job(args)
