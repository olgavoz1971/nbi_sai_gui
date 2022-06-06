#!/usr/bin/python3

#from pygrace import grace

import os
import sys
import numpy as np
import math
from skimage.measure import moments, moments_central
from time import sleep

from jkclass import printd, printdd, printe, round2float

from astropy.io import fits
from astropy.coordinates import Angle
from astropy import units as u
from astropy.stats import mad_std

import pyregion
import prf_module as prf

def circ(x, y, xc, yc, rad):
    return (x-xc)**2 + (y-yc)**2 >= rad**2

#import matplotlib.pyplot as plt
#from astropy.visualization import SqrtStretch, simple_norm

def calc_fwhm(xx, yy, rad, img):
  impix = img[yy-rad:yy+rad,xx-rad:xx+rad].copy()
  image_centre_coord, fwhm_empir, com, r, fl = \
       prf.profile(impix, xx-rad, yy-rad, 10, prf.center_method_list[0])
  return fwhm_empir
#  ind = r.argsort()
#  xsort = r[ind]
#  ysort = fl[ind]
#  xrevers = xsort[: :-1]
#  yrevers = ysort[: :-1]
#  xdouble = np.hstack((-1.*xrevers,xsort))
#  ydouble = np.hstack((yrevers,ysort))
#  fitmodel = prf.fitmodel_default
#  fwhm_model, gfit = prf.fitmodel_dict[fitmodel](xdouble, ydouble, 10000., 0., 7.)
#  return fwhm_empir, fwhm_model

def calc1mom(xx, yy, rad, img):
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

#  Mc = moments_central(impix)
  M =  moments(impix)
#  sum_ = np.sum(impix)
#  mxx = Mc[0,2]/Mc[0,0]
#  myy = Mc[2,0]/Mc[0,0]

  xc = M[0,1]/M[0,0] + xx - rad
  yc = M[1,0]/M[0,0] + yy - rad

  return xc,yc

# Calculate minimum
#   foci xc yc mxx myy

def calc_minimum(x, y, degree=3):
#  print('\n',x,y)
  poly = np.poly1d(np.polyfit(x,y, degree))
  roots = poly.deriv(1).r
#  max or min?
  test_poly = poly.deriv(2)(roots)
  x_min_poly = roots[test_poly>0]
  x_max_poly = roots[test_poly<0]
  
  if len(x_min_poly) < 1.0: 
    focus = -1
  else:
    focus = round2float(x_min_poly[0], 3)
  return focus, poly

#import matplotlib.pyplot as plt
#from matplotlib.pyplot import plot, draw, show
#def draw_matplotlib(x, y, poly):
##  fig, ax = plt.subplots()
#  plt.clf()
#  plt.scatter(x, y)
##  ax.set(xlabel='foci, m2_pos', ylabel='psf, asec', title='Looking for the best focus')
#  nsteps = 100.0
#  step = (x[-1] - x[0]) / nsteps
#  x_detailed = np.arange(x[0],x[-1]+step,step)
#  plt.plot(x_detailed, poly(x_detailed))
##  ax.grid()
#  plt.draw()
#  plt.pause(0.001)

#def draw_grace(g, x, y, poly):
#  g.restart()
#  g.x = x
#  g.y = y
#  g('s0 symbol color "green"')
#  g('s0 line type 0')
#  g('s0 symbol 1')
#  g('s0 symbol fill color "green"')
#  g('s0 symbol fill pattern 1')
#  g('s0 symbol size 0.5')
#  g('plot(x,y)')
#
##  x.sort()
#  g.hold(1)
#  nsteps = 100.0
#  step = (x[-1] - x[0]) / nsteps
#  x_detailed = np.arange(x[0],x[-1]+step,step)
#  g.plot(x_detailed, poly(x_detailed))
##  input("Press Enter to continue...")
#  return

# Mark foci on ds9
#def mark_on_ds9(data, d, x_col, y_col, text_col):
def mark_on_ds9(d, x, y, text):
#  for row in data[:]:
  for i in range(len(x)):   
    foci_reg = 'physical; point('+str(x[i])+','+str(y[i])+') # text={'+str(round(text[i],2))+'} point=cross color=red'
    d.set('regions', foci_reg)

def get_region(d, error_tkvar):
  x0,y0,rad = 0,0,0
  r = d.get('region')
  try:
   regions_all = pyregion.parse(r)
  except ValueError as e:
   error_tkvar.set('Please, select first and last star image and try again')
   return None

# Clean regins (I use only circles)
  regions = []
  for reg in regions_all:
    if reg.name == 'circle':
     regions.append(reg)

  num_reg = len(regions)
  if num_reg < 2:
    error_tkvar.set('Please, select two circle regions and press Enter')
    return None
  
# Use only two last circles
  xc_beg = regions[num_reg-2].coord_list[0]
  yc_beg = regions[num_reg-2].coord_list[1]
  xc_end = regions[num_reg-1].coord_list[0]
  yc_end = regions[num_reg-1].coord_list[1]

# Tolerance is the biggest circle radius
  tol_rad = max(regions[num_reg-2].coord_list[2], regions[num_reg-1].coord_list[2])
   
  return (xc_beg, yc_beg), (xc_end, yc_end), tol_rad

# rotate ortobasis:
#def fxy(x,y,alp):
#  x1 = x * math.cos(alp) - y * math.sin(alp)
#  y1 = 1.0 * x * math.sin(alp) + y * math.cos(alp)
#  return x1,y1

def precalc(ds9, fullpath2focusfile, focbeg, focstep, focnum, error_tkvar):
# Check for file with bias:
    try:
      hdulist = fits.open(fullpath2focusfile, uint=False)  # working with float data
      nhdu = len(hdulist)
      img = hdulist[nhdu-1].data
      head = fits.header.Header()
      for hdu in hdulist:
          head = head + hdu.header
    except Exception as e:
      error_tkvar.set('There are problems with the file '+fullpath2focusfile+' '+str(e))
      printe('There are problems with the file '+fullpath2focusfile, e)
      return None

    foci = np.linspace(focbeg+(focnum-1)*focstep, focbeg, focnum)
    printd('foci =', foci)
    res = get_region(ds9, error_tkvar)
    if res is None: return None
    point_beg, point_end, tol_rad = res
    printd('point_beg, point_end, tol_rad =', point_beg, point_end, tol_rad)
#    if tol_rad < 0: break
#   Calculate rotaion angle (Altitude (along focus-sequence) to matrix-column)
    alp = math.atan2((point_end[0]-point_beg[0]),(point_end[1]-point_beg[1]))
    print('alp =', alp)
# !!! calc baricenter here to correct x,y
    xb,yb = point_beg
    xe,ye = point_end
    xb, yb = calc1mom(round(xb), round(yb), round(tol_rad), img)
    xe, ye = calc1mom(round(xe), round(ye), round(tol_rad), img)
    print(xb,yb,xe,ye)

    ymax = yb > ye and yb or ye
    ymin = yb < ye and yb or ye
    xmax = xb > xe and xb or xe
    xmin = xb < xe and xb or xe

# Check foci direction:
# Check dero-angle, search upper/bottom side
    custom_seq = False
    try:
      parang = float(head.get('PARANG', None))
      posang = float(head.get('POSANG', None))
      if parang is None or posang is None:
          printe('PARANG or POSANG is None. Use custom donut sequence')
          custom_seq = True
#      parang = float(ds9.get('fits header keyword ' + 'PARANG'))
#      posang = float(ds9.get('fits header keyword ' + 'POSANG'))
      printd(parang, posang, type(parang), type(posang))
    except Exception as e:
      printe('Invalid PARANG or POSANG value:', e)
      printe('Use custom donut sequence')
      custom_seq = True

    if not custom_seq:
      zenith = Angle(parang*u.deg) - Angle(posang*u.deg)
      alpa = Angle(alp*u.rad)
      angle_difference = abs((alpa-zenith).wrap_at(180*u.deg).degree)
      printd('alpa =', alpa, 'zenith =', zenith, 'angle_difference =', angle_difference)
      if angle_difference > 90:
        print('reverse foci, dangle =',angle_difference)
        foci = foci[::-1] 
        printd(foci)

    step_x = (xe-xb) / (len(foci) - 1)
    step_y = (ye-yb) / (len(foci) - 1)
    rad = round(tol_rad)
    print(step_x, step_y, rad)
    return img, foci, xb, yb, step_x, step_y, rad

def calc_focus_psf(ds9, fullpath2focusfile, focbeg, focstep, focnum, scale, error_tkvar):
    res = precalc(ds9, fullpath2focusfile, focbeg, focstep, focnum, error_tkvar)
    if res is None: return None
    img, foci, xb, yb, step_x, step_y, rad = res
    fwhm_array = np.zeros([len(foci)])
    x_array  = np.zeros([len(foci)]) # for marking on ds9
    y_array  = np.zeros([len(foci)]) # for marking on ds9
    for step in range(0, len(foci)):
      xc = xb + step_x * step
      yc = yb + step_y * step
#      xc_corr, yc_corr, mxx, myy = calcmom(int(round(xc)), int(round(yc)), rad, img)
      xc_corr, yc_corr = calc1mom(int(round(xc)), int(round(yc)), rad, img)
      fwhm_empir_px = calc_fwhm(int(round(xc)), int(round(yc)), rad, img)
      fwhm_empir = fwhm_empir_px * scale
      fwhm_array[step] = fwhm_empir
      x_array[step] = xc_corr
      y_array[step] = yc_corr
      print('focus fwhm: %7.3f %7.3f' % (foci[step], fwhm_empir))

#    g = grace()
    focus_best, poly = calc_minimum(foci, fwhm_array)
    mark_on_ds9(ds9, x_array, y_array, foci) # Mark foci on ds9
#    draw(g, foci, fwhm_array, poly)
#    draw_matplotlib(foci, fwhm_array, poly)
    fwhm_best = poly(focus_best)
    printd(focus_best, fwhm_best)
#    key = input('Press Return to continue ...')
    return focus_best, fwhm_best, foci, fwhm_array, poly

# -------------- Donuts part
import don11
import common

def extract(img_tot, xc, yc):
    if img_tot is None: return
    print(xc, yc, common.fovpix)
    impix = don11.extract(img_tot.copy(), int(round(xc)), int(round(yc)), common.fovpix)

    thresh = impix.max() * common.donpar.thresh
    impix[impix <= thresh] = 0
    Mc = moments_central(impix)
    mxx = Mc[0,2]
    myy = Mc[2,0]
    scale = common.npixperpix / (common.ngrid / common.Rpix)
    printd(' ------------------------------ scale =', scale)
    a4 = scale * np.sqrt( (mxx + myy) * 0.5 ) / 1.102
    printd('a4 =', a4)
    return impix

def fit(impix, efoc):
    common.donpar.efoc   = efoc # 1 or -1
    printd('Fit it!')
    zres, impix, chi2 = don11.fit(impix, display=None)
    defocus_sai = zres[3]*(-2)*8**2
    print('chi2 =', chi2, 'defocus =', zres[3], 'defocus, mkm =', defocus_sai)
    seeing = zres[0]
    return defocus_sai, seeing

def calc_focus_donuts(ds9, fullpath2focusfile, focbeg, focstep, focnum, error_tkvar):
    filename_params = 'donut.par'
    don11.init(filename_params)

    img,foci,xb,yb,step_x,step_y,rad = \
    res = precalc(ds9, fullpath2focusfile, focbeg, focstep, focnum, error_tkvar)
    if res is None: return None
    img,foci,xb,yb,step_x,step_y,rad = res

    if foci[0] < foci[-1]: efoc = 1
    else: efoc = -1
    defocus = np.zeros([2])
    seeing  = np.zeros([2])
    for step in [0, 1]: # first and last
      if step == 0: # first donut
        xc = xb
        yc = yb
      else:         # last donut
        xc = xb + step_x * (len(foci) - 1)
        yc = yb + step_y * (len(foci) - 1)
        printd('donut: xc_last, yc_last =', xc, yc)
      printd('--------- DONUTS -----', step, xc, yc)
      xc_corr, yc_corr = calc1mom(int(round(xc)), int(round(yc)), 32, img)
      impix = extract(img, xc, yc)
      defocus_sai, see = fit(impix, efoc)
      defocus[step] = defocus_sai
      seeing[step]  = see
      efoc *= -1

    alpha_donuts = (defocus[0] - defocus[1]) / (foci[0] - foci[-1])
    print('alpha_donuts =', alpha_donuts)

    print('Best donuts focus 1 is', foci[0]  - defocus[0]/abs(alpha_donuts))
    print('Best donuts focus 2 is', foci[-1] - defocus[1]/abs(alpha_donuts))
    print('Seeing =', seeing[0], seeing[1])

    m2pos1 = foci[0]  - defocus[0]/abs(alpha_donuts)
    m2pos2 = foci[-1] - defocus[1]/abs(alpha_donuts)
    tol = 0.001
    if abs(m2pos1-m2pos2)<tol:
       printd('Yes')
       return m2pos1, np.mean(seeing)
    return -1, np.mean(seeing)
      
if __name__ == '__main__':
    import pyds9
    from pyds9 import DS9 as ds9

    DS9_NAME = 'KGO'

    def redraw():
      global canvas, ax, fig
      ax.clear()
#      fig = Figure(figsize=(8,4), dpi=100, facecolor='black')
#      fig, ax = plt.subplots()

      x = [0,1 ,2 ,3 ,4 ,5]
      y = [2,12,22,27,25,29]
      ax.scatter(x,y)
      canvas.draw() 

    import tkinter as Tk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure 

    fig = Figure(figsize=(8,4), dpi=100, facecolor='black')
    x = [0,1 ,2 ,3 ,4 ,5]
    y = [1,11,21,23,25,27]
    fig, ax = plt.subplots()
#    plt.clf()
    ax.scatter(x,y)
    root= Tk.Tk()

    button = Tk.Button(root, text='Draw', command = redraw)
    button.pack()
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack()

    root.mainloop() 

    exit(0)    
    
    bar1 = FigureCanvasTkAgg(figure1, root)
    bar1.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)

    figure = plt.Figure(figsize=(6,5), dpi=100)
    ax = figure.add_subplot(111)
    ax.plot(x,y)
    ax.set_title('Title')
    input('aaa')
    exit(0)


#    class Tkvar:
#      var = ''
#      def set(self, value):
#        self.var = value    

    ds9_list = pyds9.ds9_targets()
    if ds9_list is None:
      sys.exit("Please, open ds9 and select first and last star image of the focus secuence")
    if(len(ds9_list) > 1):
      d = ds9('DS9:' + DS9_NAME)
    else:
      d = ds9()
    
    fullpath2focusfile = 'work_foc_2021-03-04T23:32:49.fits'
#    fullpath2focusfile = d.get('file').replace('[im1]','')
    printd('fullpath2focusfile =', fullpath2focusfile)
    focbeg  = 10.75
    focstep = 0.05
    focnum  = 7
 
    error_tkvar = Tk.StringVar(value='')
    focus_best, fwhm_best = calc_focus_psf(d, fullpath2focusfile, focbeg, focstep, focnum, 0.155, error_tkvar)
#    focus_best, seeing = calc_focus_donuts(d, fullpath2focusfile, focbeg, focstep, focnum, error_tkvar)
#    fwhm_best = seeing
    printd(focus_best, fwhm_best)
    sleep(1)
    focus_best, fwhm_best = calc_focus_psf(d, fullpath2focusfile, focbeg, focstep, focnum, 0.3, error_tkvar)
   
    input("Press Enter to finsh...")
