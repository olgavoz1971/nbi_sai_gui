#!/usr/bin/python3

# Module for nbi camera
# Maintains geometry methods and calculations
# Contains methods from prf1.py
SCALE = 0.155 # '' / px
gain  = 0.28
read_out_noise = 6.4 # e-

# !!! TMP !!!! use INSTANGS for it
NBI_TDS_DIFFERENCE = 5.0 + 85.1    # instang_tds - instang_nbi
# NBI_POSANG = TDS_POSANG + 90.1

import debug
from jkclass import printd, printdd, printe, PV, sec2deg, pol2cart, cart2pol, calc_noise
#DS9_NAME = debug.ds9_name

import sys

import tkinter as Tk
import numpy as np
#from astropy.modeling import models, fitting
from astropy.coordinates import Angle
import astropy.units as u
#import parse as prs

from math import radians, sin, cos, isnan, degrees, atan

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import prf_module as prf

#import pyds9
#from pyds9 import DS9 as ds9

SAI25_ICS       = PV('SAI25:ICS')

if __name__ == '__main__':
  root = Tk.Tk()
  root.wm_title("Geometry and profile")

# -------------------------------------------------------

def do_nothing(*args):
  return

from photutils import DAOStarFinder
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from photutils import  (datasets, aperture_photometry,
                       CircularAperture, CircularAnnulus, IRAFStarFinder)

class Geometry():

#  def __init__(self, nbi_cut=do_nothing, nbi_fres=do_nothing, ics_shiftEquM=do_nothing):
  def __init__(self, ds9nbi, nbi, ics_shiftEquM=do_nothing):
    self.nbi    = nbi
    self.ds9nbi = ds9nbi
#    self.nbi_cut = nbi_cut
#    self.nbi_fres = nbi_fres
    self.ics_shiftEquM = ics_shiftEquM


# slit position TDS assistant
    self.slit_mode= Tk.BooleanVar(value=False)    # TDS assistant or not mode
    self.x_slit_c = Tk.DoubleVar(value=2147)
#    self.y_slit_c = Tk.DoubleVar(value=1980)
    self.y_slit_c = Tk.DoubleVar(value=1958)
    self.dy_image4TDS = Tk.DoubleVar()
    self.dy_image4TDS.set(1160)
    self.dero_xc = Tk.DoubleVar(value=2141)   #derotation centre
    self.dero_yc = Tk.DoubleVar(value=1972)

# move telescope to the selected point

    self.xfrom = Tk.DoubleVar()
    self.x_to  = Tk.DoubleVar()
    self.yfrom = Tk.DoubleVar()
    self.y_to  = Tk.DoubleVar()


    self.posang    = Tk.DoubleVar()
    self.dram_dsec = Tk.DoubleVar()
    self.ddec_dsec = Tk.DoubleVar()
    self.ddero_deg = Tk.DoubleVar()

    self.xbeg  = Tk.IntVar()
    self.ybeg  = Tk.IntVar()
    self.xsize = Tk.IntVar()
    self.ysize = Tk.IntVar()
    self.xbeg.set(1)
    self.ybeg.set(1)
#    self.xsize.set(NBI_XTOT.get()*2)
#    self.ysize.set(NBI_YTOT.get())
    self.xsize.set(self.nbi.get_xtot()*2)
    self.ysize.set(self.nbi.get_ytot())

    self.xbeg_float = 0.0
    self.ybeg_float = 0.0
    self.image_centre_coord = 0.0,0.0   # normal image's, not python's x,y (x - first)

    self.fwhm_px   = Tk.DoubleVar()
    self.fwhm_dsec = Tk.DoubleVar()
    self.fwhm_model_px   = Tk.DoubleVar()
    self.fwhm_model_dsec = Tk.DoubleVar()
    self.fwhm_xc   = Tk.DoubleVar()
    self.fwhm_yc   = Tk.DoubleVar()

    self.img = np.ones((1,1))
#    self.ds9 = -1
    
    self.fitmodel = Tk.StringVar()
    self.fitmodel.set(prf.fitmodel_default)
    self.center_method = Tk.StringVar()
    self.center_method.set(prf.center_method_list[0])
#    self.center_method.set("maximum_count")
    self.Factor = Tk.IntVar()
    self.Factor.set(10)

    self.y0 = 0; self.x0  = 0; self.rad = 0


  def initExchange(self):
    pass

  def profile(self):
     self.image_centre_coord, fwhm_empir, com, r, fl = \
     prf.profile(self.img, self.xbeg.get(), self.ybeg.get(), self.Factor.get(), self.center_method.get())
     print('self.xbeg.get(), self.ybeg.get()', self.xbeg.get(), self.ybeg.get())
     print('com =', com)
     s2n = self.phot()
     return fwhm_empir, com, r, fl, s2n

  def phot(self):
#      r_app = 10
      r_app = self.rad
      r_in  = r_app*2
      r_out = r_app*3

      print('self.image_centre_coord =', self.image_centre_coord)
      print('self.xbeg self.ybeg =', self.xbeg.get(), self.ybeg.get())
      y,x =  self.image_centre_coord
#      print('self.x0, self.y0, self.rad', self.x0, self.y0, self.rad)
#      print('self.xbeg.get(), self.ybeg.get()', self.xbeg.get(), self.ybeg.get())
#      img = self.ds9nbi.get_img()[0].data[xbeg:xbeg+xsize, ybeg:ybeg+ysize]
#      img = self.ds9nbi.get_img()[0].data[int(self.x0-self.rad):int(self.x0+self.rad), int(self.y0-self.rad):int(self.y0+self.rad)]
      x0 = int(x - r_out - 1)
      x1 = int(x + r_out + 1)
      y0 = int(y - r_out - 1)
      y1 = int(y + r_out + 1)
#      img = self.ds9nbi.get_img()[0].data[int(x):int(x+2*self.rad), int(self.y0-self.rad):int(self.y0+self.rad)]
      img = self.ds9nbi.get_img()[0].data[x0:x1, y0:y1]
      print('img.shape =', img.shape)
      
      x = x - x0
      y = y - y0
      position = (x,y)
      print('position =', position)
      aperture = CircularAperture(position, r=r_app)
      npix     = aperture.area

      annulus_aperture = CircularAnnulus(position, r_in=r_in, r_out=r_out)
      annulus_mask     = annulus_aperture.to_mask(method='center')
      n_b              = annulus_aperture.area
      print('n_b =', n_b)
      annulus_data     = annulus_mask.multiply(img)
      annulus_data_1d  = annulus_data[annulus_mask.data > 0]
      mean, median_sigclip, std = sigma_clipped_stats(annulus_data_1d)
      adu_bkg_per_px   = median_sigclip

      phot_table = aperture_photometry(img-adu_bkg_per_px, aperture)
      print(phot_table)

      adu_star = phot_table['aperture_sum'][0]

      noise          = calc_noise(adu_star, adu_bkg_per_px, npix, gain, read_out_noise)/gain
      print('signal[ADU] =', adu_star, 'noise[ADU] =', noise, 's2n =', adu_star/noise)
      return adu_star/noise


  def set_posang_from_fits(self):
#     self.find_ds9()
#     p = self.get_fits_keyword('POSANG')
     p = self.ds9nbi.get_fits_keyword('POSANG')
     try:
       self.posang.set(round(self.correct_posang(float(p)),2))
     except Exception as e:
       p_string = '' 
       self.posang.set(float('NaN'))


  def get_two_regions(self):
#     xc, yc, rad = self.get_circle_region(0) # the last selected region
     xc, yc, rad = self.ds9nbi.get_circle_region(0) # the last selected region
     self.xfrom.set(round(xc,1)) 
     self.yfrom.set(round(yc,1))
#     xc, yc, rad = self.get_circle_region(1) # the last selected region
     xc, yc, rad = self.ds9nbi.get_circle_region(1) # the last selected region
     self.x_to.set(round(xc,1)) 
     self.y_to.set(round(yc,1))
     self.set_posang_from_fits()

  def get_region(self):
     x0,y0,rad = 0,0,0
     y0,x0,rad = self.ds9nbi.get_circle_region(-1) # HERE ARE PYTHON's order y and x !!!
     self.y0  = y0
     self.x0  = x0
     self.rad = rad

     self.xbeg_float = x0 - rad if x0 >= rad else 0.0
     self.ybeg_float = y0 - rad if y0 >= rad else 0.0

     xbeg = int(round(self.xbeg_float,0))
     xend = int(x0 + rad)

     ybeg = int(round(self.ybeg_float,0))
     yend = int(y0 + rad)

     printd('---- xbeg, xend, ybeg, yend:', xbeg, xend, ybeg, yend)           # python's x y 
     printd('---- xbeg_float, ybeg_float:', self.xbeg_float, self.ybeg_float) # python's x y
     self.xbeg.set(xbeg)
     self.xsize.set(xend - xbeg)
     self.ybeg.set(ybeg)
     self.ysize.set(yend - ybeg)


  def cut_region(self):
#    return self.nbi_cut(self.xbeg.get(), self.ybeg.get(), self.xsize.get(), self.ysize.get())
# change python's order of x,y to the normal order of x,y:
    printd('cut_region:', self.ybeg.get(), self.xbeg.get(), self.ysize.get(), self.xsize.get())
    return self.nbi.cut(self.ybeg.get(), self.xbeg.get(), self.ysize.get(), self.xsize.get())

  def format_reset(self):
    self.nbi.fres()

  def fit_profile(self):
#    self.find_ds9()
#    self.img = self.ds9.get_pyfits()[0].data[self.xbeg.get():self.xbeg.get()+self.xsize.get(),\
#                                             self.ybeg.get():self.ybeg.get()+self.ysize.get()]

    self.img = self.ds9nbi.get_img()[0].data[self.xbeg.get():self.xbeg.get()+self.xsize.get(),\
                                             self.ybeg.get():self.ybeg.get()+self.ysize.get()]
    if (self.img.size < 4): return

# lll
    fwhm_empir, com, x, y, s2n = self.profile()
    ind = x.argsort()
    xsort = x[ind]
    ysort = y[ind]
    xrevers = xsort[: :-1]
    yrevers = ysort[: :-1]
    xdouble = np.hstack((-1.*xrevers,xsort))
    ydouble = np.hstack((yrevers,ysort))

    fwhm_model, gfit = prf.fitmodel_dict[self.fitmodel.get()](xdouble, ydouble, 10000., 0., 7.)
#    fwhm_model, gfit = select_fitmodel(self.fitmodel.get())(xdouble, ydouble, 10000., 0., 7.)
    printd('MODEL fwhm', fwhm_model, 'model details:', gfit)

# ----------- Plot it :
    self.a.clear()
    self.b.clear()
    self.a.set_aspect('auto')
    self.a.plot(xsort, ysort, 'bo')
    self.a.plot(xsort, gfit(xsort), 'r-', lw=2, label='Gaussian')
    self.b.contour(self.img)
#    self.a.imshow(self.img, vmin=self.img.min(), vmax=self.img.max(), norm=matplotlib.colors.LogNorm(), cmap ='gray', origin='bottom')
# -------------

    printd('FWHM='+str(round(fwhm_empir*SCALE,2))+'\'\'')
    self.fwhm_px.set(round(fwhm_empir,1))
    self.fwhm_dsec.set(round(fwhm_empir*SCALE,2))
    
    self.fwhm_model_px.set(round(fwhm_model,1))
    self.fwhm_model_dsec.set(round(fwhm_model*SCALE,2))

    self.fwhm_xc.set(round(self.image_centre_coord[0],1))  # here NOT python's x,y, here are normal image x,y
    self.fwhm_yc.set(round(self.image_centre_coord[1],1))

#         verticalalignment='top', horizontalalignment='right',

#    self.a.text(0.5, 0.5, "Profile",
#         verticalalignment='center', horizontalalignment='center',
#         transform = self.a.transAxes,
#         color='green', fontsize=20)

    y0_text = 0.97
    dy_text = 0.065
    self.a.text(0.9, y0_text, str(self.fwhm_model_dsec.get()) + "\'\'",
         verticalalignment='top', horizontalalignment='right',
         transform = self.a.transAxes,
         color='black', fontsize=10)
    y0_text -= dy_text

    self.a.text(0.9, y0_text, str(self.fwhm_dsec.get()) + "\'\'",
         verticalalignment='top', horizontalalignment='right',
         transform = self.a.transAxes,
         color='black', fontsize=10)
    y0_text -= dy_text

    s = 's/n:' + "{:3.0f}".format(s2n)
    self.a.text(0.9, y0_text, s,
         verticalalignment='top', horizontalalignment='right',
         transform = self.a.transAxes,
         color='black', fontsize=10)

    self.canvas.draw()

  def calc_equ_shift(self):
    self.set_posang_from_fits()
    for var in [self.x_to, self.y_to, self.xfrom, self.yfrom]:
      self.check_and_correct_tkfloat_var(var)
    if isnan(self.posang.get()): return -1
    posang_radians = radians(self.posang.get())
    cos_posang = cos(posang_radians)
    sin_posang = sin(posang_radians)
  
  # Here! Check signs !
    dx = (self.x_to.get() - self.xfrom.get())
    dy = (self.y_to.get() - self.yfrom.get())

    printdd('calc_equ_shift')
    printdd('dx =',dx,'dy =', dy, 'r*scale =', np.sqrt( dx**2 + dy**2)*SCALE )

    ddec_dx = dx * SCALE * sin_posang
    dram_dx = -1 * dx * SCALE * cos_posang

    ddec_dy = -1 * dy * SCALE * cos_posang
    dram_dy = -1 * dy * SCALE * sin_posang
    
    printd('calc_equ_shift')
    printd('dram_dx = ', dram_dx, 'dram_dy =', dram_dy, 'dram =', dram_dx + dram_dy)
    printd('ddec_dx = ', ddec_dx, 'ddec_dy =', ddec_dy, 'ddec =', ddec_dx + ddec_dy)

    self.dram_dsec.set(round(dram_dx + dram_dy,2))
    self.ddec_dsec.set(round(ddec_dx + ddec_dy,2))
    self.ddero_deg.set(0.0)

  def do_equ_shift(self):
#    self.calc_equ_shift()
    for var in [self.dram_dsec, self.ddec_dsec, self.ddero_deg]:
      self.check_and_correct_tkfloat_var(var)
    dram_deg  = sec2deg(self.dram_dsec.get())
    ddec_deg  = sec2deg(self.ddec_dsec.get())
    ddero_deg = self.ddero_deg.get()
    return self.ics_shiftEquM(dram_deg, ddec_deg, ddero_deg)

  def correct_posang(self, posang):
    print (' !!!!!!!!!!!!!!!!!!!!!!! ----- correct_posang', posang)
    try:
       posang_float = float(posang)
    except Exception as e:
       posang_float = float('NaN')
#    if not self.slit_mode.get(): return posang
    if SAI25_ICS.get().upper() == 'TDS':
       print ('TDS', posang_float + NBI_TDS_DIFFERENCE)
       return posang_float + NBI_TDS_DIFFERENCE       
    else:
       return posang_float

  def show_slit(self):
   slit_len_dminuts = 3.0
#   self.find_ds9()  # in set_posang
   self.set_posang_from_fits()

   posang = self.posang.get()
   try:
      if not isnan(posang):
          vangle = Angle(90,u.degree) + Angle(posang,u.degree)
          vangle.wrap_at('360d', inplace=True)
          vlen = 500
#          vxstart = float(self.ds9.get("fits width"))/2
#          vystart = float(self.ds9.get("fits height"))/2
          vxstart = float(self.ds9nbi.get_fits_keyword("fits width"))/2
          vystart = float(self.ds9nbi.get_fits_keyword("fits height"))/2
          vtext = 'N'
#          north_arrow = 'physical; # vector('+str(vxstart)+\
#           ','+str(vystart)+','+str(vlen)+','+str(vangle.degree)\
#           +') vector=1 width=3 color=green font=\"helvetica 14 bold roman\" text={N}'
          north_arrow = self.ds9nbi.norh_arrow_region(vxstart, vystart, vlen, vangle.degree, color='green')
##          self.ds9.set('regions', north_arrow)
          self.ds9nbi.set_region(north_arrow)

##  DETWIN1 = '[1:4296, 1401:2561]'                 [x0,yend, y0, yend]                                
      s = self.ds9nbi.get_fits_keyword('DETWIN1')
#      s = self.get_fits_keyword('DETWIN1')
      detwin_list = s.replace('[','').replace(']','').replace(':',',').split(',')

      slit_len = round(60 * slit_len_dminuts / SCALE,1)
      xstart = self.x_slit_c.get() - slit_len / 2.0
      xend   = self.x_slit_c.get() + slit_len / 2.0 

      y0 = float(detwin_list[2])
      ystart = self.y_slit_c.get() - y0
      yend   = ystart
      
      print ('y0, NBI_YTOT, ystart\n',\
           y0, NBI_YTOT.get(), ystart)

##      slit = 'physical; line('+str(xstart)+','+str(ystart)+','+str(xend)+','+str(yend)+') # line=1 width=3 color=green select=1 highlite=1 dash=0 fixed=0 edit=1'
##      slit = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\nline(1567,1980,2727,1980) # line=1 1 width=3'
#      slit = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\nline('+str(xstart)+','+str(ystart)+','+str(xend)+','+str(yend)+') # line=1 0 width=3'
      slit = self.ds9nbi.slit_region(self, xstart, ystart, xend, yend, color='green')
      print ('|' + slit + '|')
#      self.ds9.set('regions', slit)
      self.ds9nbi.set_region(slit)
   except Exception as e:
          printe('show_slit() ds9: Exception occured: ', e)

  def cut_and_rotate(self):
#    res = self.nbi_cut(self.ybeg.get(), self.xbeg.get(), self.ysize.get(), self.xsize.get())
#    ybeg = int(round(IMAGE_YFULL / 2.0 - self.dy_image4TDS.get() / 2.0))
    ybeg = int(self.y_slit_c.get() - self.dy_image4TDS.get() / 2.0) + 1   # ds9 counts pixels from 1
    printd('cut_and_rotate, ybeg =', ybeg)
#    res = self.nbi_cut(1, ybeg, NBI_XTOT.get(), self.dy_image4TDS.get())
    res = self.nbi.cut(1, ybeg, self.nbi.get_xtot(), self.dy_image4TDS.get())
    d = self.ds9nbi.find_ds9()
    d.set('rotate to 270')
##    self.slit_mode = True
# rotate 270
    return res

  def restore_image(self):
    self.format_reset()
    d = self.ds9nbi.find_ds9()
    d.set('rotate to 0')
#    self.find_ds9()
#    self.ds9.set('rotate to 0')
##    self.slit_mode = False
    return

  def check_and_correct_tkfloat_var(self, tkfloat_var):
    try:
      val_float = float(tkfloat_var.get())
      return 0
    except ValueError:
      tkfloat_var.set(0)
      return -1

  def calc_equ_shift_project_to_the_slit(self):
    self.set_posang_from_fits()
    for var in [self.x_to, self.y_to, self.xfrom, self.yfrom]:
      self.check_and_correct_tkfloat_var(var)
    dx = (self.x_to.get() - self.xfrom.get())
    dy = (self.y_to.get() - self.yfrom.get())
    printdd(dx,dy)
    if dx == 0:
      ddero_radians = np.pi/2.0
      self.ddero_deg.set(90) # check sign !!!
    else: 
      printdd(dy/dx, atan(dy/dx), degrees(atan(dy/dx)))
      ddero_radians = atan(dy/dx)
      self.ddero_deg.set(round(degrees(ddero_radians),2))

    printd('calc_equ_shift_project')
    printd('dx =',dx,'dy =', dy, 'ddero_deg=', self.ddero_deg.get())

# calc predicted new x,y:
# And take into account the NEW DERO posang_new = posang - ddero (?)  !!!!
#    xc,yc = 2147,1980
    x = self.xfrom.get() - self.dero_xc.get()
    y = self.yfrom.get() - self.dero_yc.get()
    rho, phi = cart2pol(x,y)
#    print 'FROM x,y for pol =',x,y, 'rho, phi, phi_deg =',rho, phi, degrees(phi)
    phi1 = phi - ddero_radians
    x1,y1 = pol2cart(rho, phi1)
#    print 'FROM phi1_deg = ',degrees(phi)-degrees(ddero_radians), degrees(phi),'-',degrees(ddero_radians)

    x = self.x_to.get() - self.dero_xc.get()
    y = self.y_to.get() - self.dero_yc.get()
    rho, phi = cart2pol(x,y)
#    print 'TO x,y for pol =',x,y, 'rho, phi, phi_deg =',rho, phi, degrees(phi)
    phi1 = phi - ddero_radians
    x2,y2 = pol2cart(rho, phi1)
 #   print 'TO phi1_deg = ',degrees(phi)-degrees(ddero_radians), degrees(phi),'-',degrees(ddero_radians)

 #   printdd('rho, phi, phi-ddero_radians',rho, phi,phi-ddero_radians)
    printdd('New predicted x,y  are', x1+self.dero_xc.get(),y1+self.dero_yc.get(),'or', x2+self.dero_xc.get(),y2+self.dero_yc.get())
    self.dram_dsec.set(0.0)
    self.ddec_dsec.set(0.0)

# Move to the slit along y:
    y = y1 + self.dero_yc.get()
    dx = 0.0
    if isnan(self.posang.get()): return -1
# Take into account the new predicted posang
    posang_new = self.posang.get() - self.ddero_deg.get()
    self.posang.set(posang_new)

    posang_radians = radians(posang_new)
    cos_posang = cos(posang_radians)
    sin_posang = sin(posang_radians)
  
    dy = (self.y_slit_c.get() - y)

    printdd('calc_equ_shift')
    printdd('dx =',dx,'dy =', dy, 'r*scale =', np.sqrt( dx**2 + dy**2)*SCALE, 'new posang =',posang_new )

    ddec_dx = dx * SCALE * sin_posang
    dram_dx = -1 * dx * SCALE * cos_posang

    ddec_dy = -1 * dy * SCALE * cos_posang
    dram_dy = -1 * dy * SCALE * sin_posang
    
    printd('calc_equ_shift')
    printd('dram_dx = ', dram_dx, 'dram_dy =', dram_dy, 'dram =', dram_dx + dram_dy)
    printd('ddec_dx = ', ddec_dx, 'ddec_dy =', ddec_dy, 'ddec =', ddec_dx + ddec_dy)

    self.dram_dsec.set(round(dram_dx + dram_dy,2))
    self.ddec_dsec.set(round(ddec_dx + ddec_dy,2))
#    self.ddero_deg.set(0.0)


# -----------------------  GUI -----------------------

  def gwindow(self, parent):
    width = 10
    pdx = 2
    pdy = 1
    pdy_frame = 1

# ----------------- slit position -------
    frame01 = Tk.LabelFrame(parent, text='TDS assistant', padx=10, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame01.pack(fill=Tk.X)

    Tk.Checkbutton(frame01, text = 'TDS',  fg='darkgreen', variable=self.slit_mode).grid(row=1, column=0, padx=pdx, pady = pdy, sticky=Tk.W)

    row = 0; column = 1; width = 10
    for text in ['x_slit_centre', 'y_slit_centre', 'image height', 'x_dero_centre','y_dero_centre']:
      Tk.Label(frame01, text=text, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
      column = column + 1

    row = row+1; column = 1
    for var in [self.x_slit_c, self.y_slit_c, self.dy_image4TDS, self.dero_xc, self.dero_yc]:
      Tk.Entry(frame01, width=width, textvariable=var).grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
      column = column + 1

    row = row+1; column = 1; width = 7
    Tk.Button(frame01, text='Show slit',width=width, command=self.show_slit).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + 1
    Tk.Button(frame01, text='Cut, turn',width=width, command=self.cut_and_rotate).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + 1
    Tk.Button(frame01, text='Restore'  ,width=width, bg='lightgreen', command=self.restore_image).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)

# ----------------------- regions for shitft --------------------------

    frame00 = Tk.LabelFrame(parent, text='Place into the slit', padx=10, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame00.pack(fill=Tk.X)

    row = 0; column = 0; columnspan = 4
    Tk.Label(frame00,text='Select two regions (line)').\
       grid(row=row, column=column, columnspan=columnspan, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + columnspan
    Tk.Button(frame00, text='Calc turn and shift', width=19, command=self.calc_equ_shift_project_to_the_slit).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)

    row = row+1; column = 0; columnspan = 4
    Tk.Label(frame00,text='Select two regions (\'from\' and \'to\')').\
          grid(row=row, column=column, columnspan=columnspan, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + columnspan
    Tk.Button(frame00, text='Calc shift from 1 to 2', width=19, command=self.calc_equ_shift).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)

    row = row+1; column = 0
    for text in ['x1', 'y1', 'x2', 'y2']:
      Tk.Label(frame00, text=text, width=width, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1
    row = row+1; column = 0
    for var in [self.xfrom, self.yfrom, self.x_to, self.y_to]:
      Tk.Entry(frame00, width=width, textvariable=var).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1

    Tk.Button(frame00, text='Get two regs',command=self.get_two_regions).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)

#    Tk.Button(frame00, text='Get from',  command=self.get_region_from).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
#    column = column + 1
#    Tk.Button(frame00, text='Get to    ',command=self.get_region_to).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)

    row = row+1; column = 0
    for text in ['DRAM\"', 'DDEC\"', 'DDERO', 'POSANG']:
      Tk.Label(frame00, text=text, width=width, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1
    row = row+1; column = 0
    for var in [self.dram_dsec, self.ddec_dsec, self.ddero_deg, self.posang]:
      Tk.Entry(frame00, width=width, textvariable=var).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1
    Tk.Button(frame00, text='Go!', command=self.do_equ_shift).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)

# ---------------------------- Cut region and calc profile --------

    frame0 = Tk.LabelFrame(parent, text='Select and Cut region', padx=10, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame0.pack(fill=Tk.X)

    row = 0
    column = 0
    for text in ['xbeg', 'ybeg', 'xsize', 'ysize']:
      Tk.Label(frame0, text=text, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1

    row = row+1
    column = 0
    for var in [self.ybeg, self.xbeg, self.ysize, self.xsize]:    # Change x<->y to show normal image coordinates
      Tk.Entry(frame0, width=width, textvariable=var).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1

    row = row+1
    column = 0
    Tk.Button(frame0, text='Get reg ', command=self.get_region).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + 1
    Tk.Button(frame0, text='Cut reg ', command=self.cut_region).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + 1
    Tk.Button(frame0, text='Profile ', command=self.fit_profile).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + 1
    Tk.Button(frame0, text='F.Reset ',  bg='lightgreen', command=self.format_reset).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
    column = column + 1
    
# -------------- Profile   RESULTS ----------------------------------

    frame1 = Tk.LabelFrame(parent, text='Profile', padx=10, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame1.pack(fill=Tk.X)

    row = 0
    column = 0
    for text in ['fwhm px', 'fwhm \'\'', 'fwhm fit px', 'fwhm fit \'\'', 'star xc', 'star yc']:
      Tk.Label(frame1, text=text, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1

    row = row + 1
    column = 0
#    for var in [self.fwhm_px, self.fwhm_dsec, self.fwhm_xc, self.fwhm_yc]:
    for var in [self.fwhm_px, self.fwhm_dsec, self.fwhm_model_px, self.fwhm_model_dsec, self.fwhm_xc, self.fwhm_yc]:
      Tk.Label(frame1, textvariable=var, fg='blue').grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1

#    row = row + 1
#    column = 0
#    Tk.Label(frame1, width=width, text='star xc:').grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
#    column = column + 1
#    Tk.Label(frame1, width=width, textvariable=self.fwhm_xc, fg='blue').grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
#    column = column + 1
#    Tk.Label(frame1, width=width, text='star yc:').grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
#    column = column + 1
#    Tk.Label(frame1, width=width, textvariable=self.fwhm_yc, fg='blue').grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)

# ---------------------------- Profile Canvas  --------

    frame3 = Tk.Frame(parent, padx=20, borderwidth=3, relief=Tk.RIDGE)
    frame3.pack(fill=Tk.X)

    figure = Figure(figsize=(5,2.5), dpi=100)
#    self.a = figure.add_subplot(121, aspect='auto')
    self.a = figure.add_subplot(121, aspect='equal')
    self.b = figure.add_subplot(122, aspect='equal')
    self.a.text(0.5, 0.5, "Profile",
         verticalalignment='center', horizontalalignment='center',
         transform = self.a.transAxes,
         color='green', fontsize=20)
    self.b.text(0.5, 0.5, "Contour",
         verticalalignment='center', horizontalalignment='center',
         transform = self.b.transAxes,
         color='green', fontsize=20)

    self.canvas = FigureCanvasTkAgg(figure, master=frame3)
    self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, padx=0, pady=pdy, expand=1)
    self.canvas.draw()
##    self.toolbar = NavigationToolbar2TkAgg(self.canvas, frame3)
#    self.toolbar = NavigationToolbar2Tk(self.canvas, frame3)
#    self.toolbar.update()
# -------------------------- Profile Options ---------------------------------

    frame2 = Tk.Frame(parent, padx=10, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame2.pack(fill=Tk.X)

    frame_fitmodels = Tk.Frame(frame2, padx=20, borderwidth=3)
    frame_fitmodels.pack(side=Tk.LEFT, fill=Tk.Y)
    for model_name in prf.fitmodel_dict:
        Tk.Radiobutton(frame_fitmodels, text=model_name, variable=self.fitmodel, value=model_name).pack(anchor=Tk.W)

    frame_center = Tk.Frame(frame2, padx=20, borderwidth=3)
    frame_center.pack(side=Tk.LEFT, fill=Tk.Y)

    for method in prf.center_method_list:
      Tk.Radiobutton(frame_center, text=method,        variable=self.center_method, value=method).pack(anchor=Tk.W)
#    Tk.Radiobutton(frame_center, text='Maximum count',     variable=self.center_method, value="maximum_count").pack(anchor=Tk.W)
#    Tk.Radiobutton(frame_center, text='Maximum count filt',variable=self.center_method, value="maximum_count_filt").pack(anchor=Tk.W)
    Tk.Entry(frame_center,       textvariable=self.Factor, width=5).pack(anchor=Tk.W)

# ------------------- End of GUI part -----------------------



#
if __name__ == '__main__':
  geometry = Geometry()
  geometry.gwindow(root)
  root.mainloop()
