from __future__ import print_function

import tkinter as Tk
import sys
import debug
epics = debug.epics
debug_log = debug.debug_log
import os.path
import numpy as np
import parse as prs
import pyregion
from time import sleep

MAX_PV_LENGTH = 39
SCALE = 0.155 # asec / px
DS9_NAME = 'KGO'

if epics:
  from epics import PV

# --- auxiliary functions ---
from astropy.coordinates import Angle
from astropy import units as u
from math import isnan
from datetime import datetime

def printd(*args):
  if debug_log:
     print(*args)
#    for arg in argv: print(arg),   # print without \n

def printdd(*args):
#  if debug_log:
     print(*args)

def printe(*args, **kwargs):
    print('Error!', *args, file=sys.stderr, **kwargs)


# ------------ DS9 -------------------
import pyds9  
from pyds9 import DS9 as ds9  
class DS9NBI:

  def find_ds9(self):
     try:
       ds9_list = pyds9.ds9_targets()
       if(len(ds9_list) > 1): d = ds9('DS9:'+DS9_NAME)
       else: d = ds9()
     except Exception as e:
       printe('DS9NBI:find_ds9 exception ',e)
       return None
     return d

  def get_path_to_file(self):     # only local path :(
       print('DS9NBI: get_path_to_file')
       fullfilename = None
       d = self.find_ds9()
       if d is None: return fullfilename
       fullfilename = d.get('file').replace('[im1]','').replace('[1]','').replace('[1]','')
       print('fullfilename =', fullfilename)
       return fullfilename

  def get_fits_keyword(self):
       d = self.find_ds9()
       if d is None: return None
       keyvalue = d.get('fits header keyword ' + keyword)
       return keyvalue

  def find_all_circle_regions(self, regions_string):
    reg = regions_string
    circ_list = []
    pos = reg.find('circle')
    while pos >= 0:
      circ = prs.search('circle({:g},{:g},{:g}', reg, pos=pos)
      circ_list.append(circ)
      pos = reg.find('circle',pos+1)
    return circ_list

  def set_region(self, region):
     d = self.find_ds9()
     if d is None: return None
     d.set('regions', region)

  def get_circle_region(self, which=0):   # number in the list, m.b 0,1, ... -1,-2
     d = self.find_ds9()
     if d is None: return 0,0,0
     regs = d.get('region').replace('\"','')
#     printd('regions:', type(regs), '|'+regs+'|')
     circ_list = self.find_all_circle_regions(regs)
     if which in range(-len(circ_list), len(circ_list)):
       circ = circ_list[which]
       return circ[0], circ[1], circ[2] # ordinary image x_centre, y_centre, radius
     return 0,0,0

  def get_region(self):
    d = self.find_ds9()
    if d is None: return None
    x0,y0,rad = 0,0,0
    r = d.get('region')
    try:
     regions_all = pyregion.parse(r)
    except ValueError as e:
     return None

  def get_img(self):
    d = self.find_ds9()
    if d is None: return None
    img = d.get_pyfits()
    return img

  def norh_arrow_region(self, vxstart, vystart, vlen, vangle_degree, color='orange'):
    north_arrow = 'physical; # vector('+str(vxstart)+\
      ','+str(vystart)+','+str(vlen)+','+str(vangle_degree)\
      +') vector=1 width=3 color='+color+' font=\"helvetica 14 bold roman\" text={N}'
    return north_arrow

  def slit_region(self, xstart, ystart, xend, yend, color='green'):
      slit = '# Region file format: DS9 version 4.1\nglobal color='\
      +color+' dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\nline('\
      +str(xstart)+','+str(ystart)+','+str(xend)+','+str(yend)+') # line=1 0 width=3'


  def view(self, fullpath, draw_arrow=True):
       d = self.find_ds9()
       if d is None: return None       
       try:
            d.set("file " + fullpath)
            if not draw_arrow: return
            posang = d.get('fits header keyword ' + 'posang')
#            posang = geometry.correct_posang(posang)
            print('posang=|'+str(posang)+'|', type(posang))
            if posang == '': posang = float('NaN')
            else:
              posang = float(posang)
#            posang = POSANG.get()
            if not isnan(posang):
              vangle = Angle(90,u.degree) + Angle(posang,u.degree)
              vangle.wrap_at('360d', inplace=True)
              vlen = 500
#              printd('+++++++++++++++++++ d.get("fits width"), d.get("fits height")', d.get("fits width"), d.get("fit
              vxstart = float(d.get("fits width"))/2
              vystart = float(d.get("fits height"))/2
#              vtext = 'N'
              north_arrow = self.norh_arrow_region(vxstart, vystart, vlen, vangle.degree)
#              
##              north_arrow_ics_corrected = 'physical; # vector('+str(vxstart)+\
##               ','+str(vystart)+','+str(vlen)+','+str(vangle.degree)\
##               +') vector=1 width=3 color=green font=\"helvetica 14 bold roman\" text={N}'

              d.set('regions', north_arrow)            
       except ValueError as e:
            printe('ds9 view: Exception occured: ', e)


# ------------------------------------
def calc_noise(adu_star, adu_bkg_per_px, npix, gain, read_out_noise):  # in ADU  
# Classic formula
   n_dark = 0.0
   n_star = adu_star * gain
   n_bkg  = adu_bkg_per_px * gain
   noise = np.sqrt(n_star + npix * (n_bkg + n_dark + read_out_noise**2))
   print(n_star, npix, n_bkg, n_dark, read_out_noise)
   print('noise =', noise)
   return noise

def composeSAIoutfile():
    return datetime.utcnow().strftime("%Y%m%d%H%M%S")

def sec2deg(dsec):
  return dsec/60.0/60.0

def deg2dms(degree):
    try:
      if isnan(degree) : return ''
      a = Angle(degree, u.degree)
      ret = a.to_string(unit=u.degree, sep=':')
    except Exception as e:
      printe('deg2dms ', degree, ':', e)
      return ''
    return ret

def deg2hms(degree):
    try:
      if isnan(degree) : return ''
      a = Angle(degree, u.degree)
      ret = a.to_string(unit=u.hour, sep=':')
    except Exception as e:
      printe('deg2hms ', degree, ':', e)
      return ''
    return ret

def dms2deg(d,m,s):
   a = Angle((d, m, s), unit=u.deg)  # (d, m, s)
   return a.degree

def hms2deg(h,m,s):
   a = Angle((h, m, s), unit='hourangle')  # (h, m, s)
   return a.degree

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def round2int(s):
  try:
    num = int(np.rint(float(s)))
  except ValueError:
    num = -1
  return num  

def round2float(s,ndigits):   # if ndigits < 0 -> no convertion to float/int
  try:
    if ndigits > 0:
      num = round(float(s), ndigits)
    else:
      num = int(round(float(s)))    
  except Exception as e:
    printe('round2float', e)
    num = -1
  return num  

def correct_outfile(filepath_, fname_):
    printd( filepath_)
    printd( fname_)
    if (fname_ == 'autofocusing') or (fname_ == 'focusing'):
      fullpath = os.path.join(filepath_, fname_+'.fits')
      exist = os.path.exists(fullpath)
      if exist:
        print('file ', fname_, 'exist. Removing...')
        os.remove(fullpath)
        sleep(1)        
      return fname_
    i = 0
    fname = fname_
    while i < 1000:
      fullpath = os.path.join(filepath_, fname+'.fits')
      exist = os.path.exists(fullpath)
      if not exist:
        return fname
      fname = fname_ + str(i)
      i = i+1
    return fname


# ------------- PV epics based classes

class JkDoubleVar(Tk.Variable):
        _default = 0
        def __init__(self, master=None, value=-1, ndigits=-1):
                self.ndigits = ndigits
                Tk.Variable.__init__(self, master, value)
        def get(self):
                varch = self._tk.globalgetvar(self._name)
                if varch == 'None' :
                  return None
                else:
                  return self._tk.globalgetvar(self._name)
        def set(self, value):
                value_corr = value
                if self.ndigits >= 0:
                   value_corr = round2float(value,self.ndigits)
                self._tk.globalsetvar(self._name, value_corr)


# --- Epics ----------
if not epics:
  class PV:
    def __init__(self, epics_name):
      self.value = 0
      self.pvname = epics_name
      printd( '+++++++++++++++++++++++++++++++++++++++++++++++++++++')

    def put(self, val):
      self.value = val

    def get(self):
      return self.value
      
    def add_callback(self, *args, **kwargs):
      pass

    def connect(self,timeout):
      return True


empty_pv = PV('Empy_PV')
empty_pv.put(None)

class PVTkbase():
  def __init__(self, parent, pv=None, ndigits=-1):
    printd( self.__class__.__name__, 'PVTkbase.__init__')
    self.parent  = parent

  def readPV(self):
     pass

  def initExchange(self):
     pass

class PVTk(PVTkbase):
  def __init__(self, parent, pv=None, ndigits=-1):
    PVTkbase.__init__(self,parent)
    printd( self.__class__.__name__, '__init__')
    self.parent  = parent
    self.pv = pv
    self.tkvar = JkDoubleVar(ndigits=ndigits)

  def readPV(self):
    printd( self.__class__.__name__, 'readPV')
    if epics:
       self.tkvar.set(self.pv.value)

  def pvcallback(self,pvname=None, value=-1, char_value=-1, **kw):
    if epics:
       self.tkvar.set(value)

  def initExchange(self):
    printd( self.__class__.__name__,'initExchange')
    if epics:
      printd( self.pv.pvname)
      try:
         self.pv.add_callback(self.pvcallback, with_ctrlvars=False)
      except Exception as e:
         printe( self.pv.pvname, ' add_callback exception:', e)
    else:
      self.tkvar.set(-1)
    printd( self.__class__.__name__, 'end of initExchange')



# ---- Widgets -------;

color_status = {'READY':'lightgreen', 'LOOP':'yellow', 'NOSTAR':'orange', 'ERROR':'red', 'BUSY':'cyan', 'OFF':'gray', 'NOHOME':'yellow', 'MANUAL':'magenta', 'PARK':'blue', 'OFF':'gray', -1:'white', None:'white',-1:'white'}
enum_onoff   = {1:'ON', 0:'OFF'}
color_status_onoff = {'ON':'lightgreen', 'OFF':'gray', 1:'lightgreen', 0:'gray'}


class JKLabel(Tk.Label):
  def __init__(self, parent, pv=None, *args, **kwargs):
    self.pv = pv
    self.parent = parent
    relief = Tk.RIDGE
    self.tkvar = JkDoubleVar()
    self.tkvar.set('None')
#    self.help_text = 'LLL'
    Tk.Label.__init__(self, parent, textvariable=self.tkvar, relief=relief, *args, **kwargs)
    self.readPV()
    self.initExchange()

  def readPV(self):
    self.tkvar.set(self.pv.get())

  def pvcallback(self, pvname=None, value=-1, char_value=-1, **kw):
    self.tkvar.set(value)

  def initExchange(self):
    try:
         self.pv.add_callback(self.pvcallback,    with_ctrlvars=False)
    except Exception as e:
         printd( self.__class__.__name__, ' add_callback exception:', e)


class IndicatorStatus(Tk.Label): # PV
  def __init__(self, parent, pv=empty_pv, *args, **kwargs):
    self.pv = pv
    self.parent = parent
    self.color_code = color_status
    relief = Tk.RIDGE
    Tk.Label.__init__(self, parent, relief=relief, *args, **kwargs)
    self.text = kwargs.get('text')
    
    self.bind('<Enter>', self.on_enter)
    self.bind('<Leave>', self.on_leave)

    self.readPV()
    self.initExchange()
    
  def pvcallback(self,pvname=None, value=-1, char_value=-1, **kw):
    self.configure(bg=self.color_code.get(value,'gray'))

  def readPV(self):
     self.configure(bg=self.color_code.get(self.pv.value,'gray'))

  def initExchange(self):
#    printd( self.__class__.__name__,"initExchange")
    try:
           self.pv.add_callback(self.pvcallback, with_ctrlvars=False)
    except Exception as e:
           printe( self.pv.pvname, " add_callback exception:",e)
#    print self.__class__.__name__, 'end of initExchange'

  def on_enter(self, event):
    self.configure(text=self.pv.get())

  def on_leave(self, event):
    self.configure(text=self.text)

class IndicatorStatusTkvar(Tk.Label):
  def __init__(self, parent, tkvar, color_code=color_status_onoff, enum=enum_onoff, *args, **kwargs):
    self.tkvar = tkvar
    self.parent = parent
    self.color_code = color_code
    self.enum = enum
    relief = Tk.RIDGE
    Tk.Label.__init__(self, parent, relief=relief, *args, **kwargs, bg=self.color_code.get(self.tkvar.get(),'blue'))
    self.text = kwargs.get('text')
    
    self.bind('<Enter>', self.on_enter)
    self.bind('<Leave>', self.on_leave)
    self.tkvar.trace('w', self.tkcallback)

  def tkcallback(self, *args, **kw):
    value = self.tkvar.get()
    self.configure(bg=self.color_code.get(value,'blue'))

  def on_enter(self, event):
    if self.enum is not None:
      self.configure(text=self.enum.get(self.tkvar.get(),'Undefined'))
    else:
      self.configure(text=self.tkvar.get())

  def on_leave(self, event):
    self.configure(text=self.text)

