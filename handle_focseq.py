#!/usr/bin/python3

# Script precesses focus-sequences image. The first and last star image must me marked on ds9 window. Uses sex ldac files.
# Invoked from foc_server.tcl

DS9_NAME = 'KGO'

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
def calcmom(xx,yy, img):
  rad = 32
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

def find_nearest_data(point, data, x_col, y_col, param_col, x2_col, y2_col):
# x_col, y_col, param_col, x2_col, y2_col - culumn numbers in data-file
#def find_nearest(point, point_list):
   r2min = -1
   point_nearest = (-1,-1,-1,-1,-1)
   for p in zip(data[:,x_col],data[:,y_col],data[:,param_col],data[:,x2_col],data[:,y2_col]):
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

def draw(g,x, y,poly):
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
  g.plot(x, poly(x))
#  input("Press Enter to continue...")
  return

# Mark foci on ds9
def mark_on_ds9(data, d, x_col, y_col, text_col):
  for row in data[:]:
    foci_reg = 'physical; point('+str(row[x_col])+','+str(row[y_col])+') # text={'+str(row[text_col])+'} point=cross color=red'
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

  hdulist = fits.open(args.infile)
  nhdu = len(hdulist)
  data2 = hdulist[nhdu-1].data
  head = hdulist[0].header

  x = data2['X_IMAGE']
  y = data2['Y_IMAGE']

# Second moments:
  x2 = data2['X2_IMAGE']
  y2 = data2['Y2_IMAGE']
  xy = data2['XY_IMAGE']

  param = data2[args.param]
  if args.param == 'FWHM_IMAGE':
    param = param*args.scale

  fits_image_full_name = args.infile.replace('_ldac.cat','.fits')
  hdulist = fits.open(fits_image_full_name, uint=False)  # working with float data
  nhdu = len(hdulist)
  img = hdulist[nhdu-1].data

  g = grace()

  doit = True
  while doit:

    moments_all = []
    foci = np.linspace(args.focbeg+(args.focnum-1)*args.focstep, args.focbeg, args.focnum)
    print('foci =', foci)

    point_beg, point_end, tol_rad = get_region(d)
    if tol_rad < 0: continue
#   Calculate rotaion angle (Altitude (along focus-sequence) to matrix-column)
    alp = math.atan2((point_end[0]-point_beg[0]),(point_end[1]-point_beg[1]))
    xb,yb = point_beg
    xe,ye = point_end
    ymax = yb > ye and yb or ye
    ymin = yb < ye and yb or ye
    xmax = xb > xe and xb or xe
    xmin = xb < xe and xb or xe

# Check foci direction:
# Check dero-angle, search upper/bottom side
    deroa = Angle(args.dero*u.deg)
    pos0=0
    zenith = Angle(pos0*u.deg) + deroa
    alpa = Angle(alp*u.rad)
    angle_difference = abs((alpa-zenith).wrap_at(180*u.deg).degree)
    if angle_difference > 90:
      print('reverse foci, dangle =',angle_difference)
      foci = foci[::-1] 
      print(foci)

# Cut sequence-table from total sex-table:
    poly_seq = np.poly1d(np.polyfit([point_beg[0], point_end[0]], [point_beg[1],  point_end[1]],1))
    mask = (abs(y-poly_seq(x)) <= tol_rad) & (y > ymin-tol_rad) & (y < ymax+tol_rad) & (x > xmin-tol_rad) & (x < xmax+tol_rad)
#    print('xmask ymask:')
#    for xy in zip(x[mask],y[mask]): print(xy)

# focus, x, y, x2, y2,  along_alt_coord, fwhm
#   0    1  2   3   4        5            6
    data_cut = np.zeros([len(x[mask]),7])
    data_cut[:,1] = x[mask]
    data_cut[:,2] = y[mask]
    data_cut[:,6] = param[mask]   #fwhm
    data_cut[:,3] = x2[mask] # sex xx2
    data_cut[:,4] = y2[mask] # sex yy2
    x2_col = 3; y2_col = 4; alt_col = 5
    x_col = 1; y_col = 2
#    print(data_cut)

#### to new ortobasis
#    x1,y1 = fxy(x, y, alp)
    y1, x1 = fxy(data_cut[:,x_col], data_cut[:,y_col], alp)
    data_cut[:,alt_col] = x1
# Sort along altitude-column (y1)
    data_cut_sorted = data_cut[data_cut[:,alt_col].argsort()]

    step_x = (xe-xb) / (len(foci) - 1)
    step_y = (ye-yb) / (len(foci) - 1)
    print('Step =', step_x, step_y)

#   foci xc yc mxx myy  a4 fwhm
#     0   1  2  3   4   5   6
    da = np.zeros([len(foci), 7])
    x_col = 1;    y_col = 2; fwhm_col = 6
    foci_col = 0; a4_col = 5
    for step in range(0, len(foci)):
      xc = xb + step_x * step
      yc = yb + step_y * step
      xc_corr, yc_corr, mxx, myy = calcmom(int(round(xc)), int(round(yc)), img)
      xc_sex , yc_sex, fwhm, xx2, yy2 = find_nearest_data((xc_corr,yc_corr), data_cut_sorted, x_col, y_col, fwhm_col, x2_col, y2_col)
      print('%7.3f %7.3f\n%7.3f %7.3f\n%7.3f %7.3f\n' % (xc, yc, xc_corr, yc_corr, xc_sex , yc_sex))
      scale = 2.385875519959664
      a4_1 = scale * np.sqrt( (mxx + myy) * 0.5 ) / 1.102
      a4_2 = scale * np.sqrt( (xx2 + yy2) * 0.5 ) / 1.102
      a4 = a4_1
      da[step] = (foci[step], xc_corr, yc_corr, mxx, myy, a4, fwhm)
      print('foci[step] =', foci[step], fwhm, xx2, yy2)

#    xb1,yb1 = fxy(point_beg[0], point_beg[1], alp)
#    xe1,ye1 = fxy(point_end[0], point_end[1], alp)
#    mask1 = (x1>(xb1-tol_rad)) & (x1<(xe1+tol_rad)) & (y1>(yb1-tol_rad)) & (y1<(ye1+tol_rad))

#    print(da)

    focus, poly = calc_minimum(da[:, foci_col], da[:, fwhm_col])
#    focus, poly = calc_minimum(da[:,0], da[:,6])
#    draw(g, da[:,foci_col], da[:,fwhm_col], poly)
    mark_on_ds9(da, d, x_col, y_col, foci_col) # Mark foci on ds9
    draw(g, da[:,0], da[:,6], poly)

    print('\n\nfocus, poly =', focus, poly)

#    for xx,yy,xx2,yy2 in zip(x[mask],y[mask],x2[mask],y2[mask]): 
#      xc,yc,mxx,myy = calcmom(int(round(xx)), int(round(yy)), img)
#      scale = 2.385875519959664
#      a4_1 = scale * np.sqrt( (mxx + myy) * 0.5 ) / 1.102
#      a4_2 = scale * np.sqrt( (xx2 + yy2) * 0.5 ) / 1.102
#      moments_all.append((xx, yy, mxx, myy, xx2, yy2, a4_1, a4_2))

    i_focus1 = 0
    i_focus2 = len(foci) - 1
#    f_minus = foci[i_minus]
#    f_plus  = foci[i_plus]
#    A4_minus = moments_all[i_minus][7]
#    A4_plus  = moments_all[i_plus][7]
#lll
    
    if da[i_focus1, foci_col] < da[i_focus2, foci_col]:
      i_minus = i_focus1
      i_plus  = i_focus2
    else:
      i_minus = i_focus2
      i_plus  = i_focus1

    f_minus  = da[i_minus, foci_col]
    f_plus   = da[i_plus,  foci_col]
    A4_minus = da[i_minus, a4_col]
    A4_plus  = da[i_plus , a4_col]
    print('i_minus, i_plus, f_minus, f_plus, A4_minus, A4_plus')
    for ent in (i_minus, i_plus, f_minus, f_plus, A4_minus, A4_plus):
      print( '%7.3f' % ent, end='')
    print()

    if focus > 0:
      F0 = focus
      alpha = (A4_plus - A4_minus) / (f_plus + f_minus - 2*F0)
      print('alpha =', alpha)
#    F0 = (f_plus + f_minus)/2.0 + (A4_minus - A4_plus)/(2.0*alpha)

    for i in range(0, len(foci)):
      print('-- %7.3f %7.3f %7.3f  %7.3f' % \
          (da[i,foci_col], da[i,foci_col]-focus, focus, da[i, a4_col]), end='')
#      for moment in moments_all[i]:
#          print(' %7.3f' % moment, end='')
      print()

    key = input('Enter Q(uit) or any other key to continue ...')
    print(key)
    if key.upper() == 'Q': break

#raw_input("Press 'Return' to continue...")    #python 2


parser = argparse.ArgumentParser()
parser.add_argument('-f','--infile', required=True, nargs='?', default='multifocus.fits', help='input sextractor-cat-file in ldac format (*ldac.cat)')
parser.add_argument('-r','--dero',   nargs='?', default=180.0, type=float, help='DERO-angle in degrees')
parser.add_argument('-d','--degree', nargs='?', default=3, type=int, help='degree of polynomials')
parser.add_argument('-p','--param',  nargs='?', default='FWHM_IMAGE', help='parameter (fwhm_image, elongation, x2_image etc)')
parser.add_argument('--focbeg',  nargs='?', required=True, type=float, help='first focus position')
parser.add_argument('--focstep', nargs='?', required=True, type=float, help='focus step')
parser.add_argument('--focnum',  nargs='?', required=True, type=int,   help='number of focus positions')
parser.add_argument('--enlargement', nargs='?', default=True, type=bool,  help='Is the focus increasing?')
parser.add_argument('--scale', nargs='?', default=0.155, type=float,  help='Image scale in \"/px')
args = parser.parse_args()

main_job(args)
