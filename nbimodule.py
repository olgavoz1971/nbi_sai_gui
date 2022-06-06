#!/usr/bin/env python3
# Version 3.0
# Exchange via epics
# Now NO EPICS - trace tkvar instead, and threading instead of many callbacks

LOG_PM = True
LOG_FNAME = '20211204_pm.log'

import tkinter as Tk
import debug
ivydebug = debug.ivydebug
#epics = debug.epics
ABSOLUTE_PATH = '/var/spool/ccd3/images/4kx4k/'
#ABSOLUTE_PATH = '.'
MIN_EXPTIME = 0.02
PM_LOG = False

#from jkclass import PV
from jkclass import round2float
from jkclass import printd, printdd, printe
from jkclass import deg2hms
from jkclass import deg2dms, MAX_PV_LENGTH
#maxn = MAX_PV_LENGTH

from time import sleep

from epics import PV
# TMP   for V.A.Senik
NBI_TEMPERL     = PV('NBI:TEMPERL')
NBI_TEMPERA     = PV('NBI:TEMPERA')
NBI_PRESSURE    = PV('NBI:PRESSURE')

ERROR_MESSAGE_tkvar	= Tk.StringVar(value='')
EXPOSURE_tkvar		= Tk.IntVar(value=0)
FILEISREADY_tkvar	= Tk.IntVar(value=0)
EXPPROGRESS_tkvar	= Tk.DoubleVar(value=0)
DRVPROGRESS_tkvar	= Tk.DoubleVar(value=0)
TEMPERL_tkvar		= Tk.StringVar(value='')
TEMPERA_tkvar		= Tk.StringVar(value='')
PRESSURE_tkvar		= Tk.StringVar(value='')
FILENAME_tkvar		= Tk.StringVar(value='')
# !!!!!!!!!!!!!!!!!!!!!! What is it ?????????
RESPONSE_tkvar		= Tk.StringVar(value='')
OBSERVER_tkvar		= Tk.StringVar(value='Somebody')

import sys
import numpy as np
import os
from datetime import datetime

'''
def send2Ivy(message):
   if ivydebug:
     return
   sleep(0.3)
   IvyStart('')
   sleep(0.3)
   print "nbi: send2Ivy",message
   IvySendMsg(message)
   sleep(0.3)
   IvyStop()
'''

import os
# PM-model C1 CORR logging

SAI25_RA            = PV('SAI25:RA')
SAI25_DEC           = PV('SAI25:DEC')
TCS2_AZ_POS         = PV('SAI25:AZ')
TCS2_ALT_POS        = PV('SAI25:ALT')
TCS2_C1_POS         = PV('SAI25:C1')

#TCS2_TRACK_PM_AZ    = PV('TCS2:TRACK_PM_AZ')
#TCS2_TRACK_PM_ALT   = PV('TCS2:TRACK_PM_ALT')
#TCS2_TRACK_PM_DERO  = PV('TCS2:TRACK_PM_DERO')

#TCS2_AIR_TEMPERATURE = PV('TCS2:AIR_TEMPERATURE')
#TCS2_AIR_PRESSURE    = PV('TCS2:AIR_PRESSURE')
ECS_AIR_TEMPERATURE = PV('ECS:TAMB')
ECS_AIR_PRESSURE    = PV('ECS:PRES')
ECS_RH               = PV('ECS:RH')

def log_pm_corr(is_beg):
  if not PM_LOG: return
  fname = LOG_FNAME
  with open(fname,'a') as f:
    if is_beg: prefix = 'B: '
    else: prefix = 'E: '
    lt = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    filename = FILENAME_tkvar.get()

    format = '%15.7f '
    data = ''
#    pvlist = [SAI25_RA, SAI25_DEC, TCS2_AZ_POS, TCS2_ALT_POS, TCS2_C1_POS, TCS2_TRACK_PM_AZ, TCS2_TRACK_PM_ALT, TCS2_TRACK_PM_DERO, TCS2_AIR_TEMPERATURE, TCS2_AIR_PRESSURE, ECS_RH]
    pvlist = [SAI25_RA, SAI25_DEC, TCS2_AZ_POS, TCS2_ALT_POS, TCS2_C1_POS, ECS_AIR_TEMPERATURE, ECS_AIR_PRESSURE, ECS_RH]

    format_all = ''
    data = ''
    for pv in pvlist:
       data += format % pv.get()
#    data = format_all + ' % '
#    for pv in pvlist:
#       data += str(pv.get()) + ','
#    print('format =', format_all)
    print('data =', data)
    message = filename + ' ' + prefix + lt + ' ' + data + '\n'
    print(message)

    f.write(message)

# Not used now
def check_process(myveryown_pid, pstring):
  for line in os.popen("ps ax | grep " + pstring + " | grep -v grep"):
    pid = int(line.split()[0])
    if pid != myveryown_pid:
#      os.kill(int(pid), signal.SIGKILL)
      printe( "own pid", myveryown_pid)
      printe( line)
      printe("nbimodule.py:", pstring, "process pid=",pid,"is running already. Exitted...")
      return 1
  return 0


if not ivydebug:
    from ivy.ivy import IvyApplicationDisconnected
    from ivy.std_api import *
    from ivy.ivy import ivylogger
    import logging


class NBI():

  def get_filename(self):
      global FILENAME_tkvar
      filename = FILENAME_tkvar.get()
      printdd ('nbimodule.get_filename %%%%%%%%%%%%%%%%%%%%%% ', filename)
      return filename

  def get_filepath(self):
      return self.filepath_current

# status variables
  ccd3comm_initialized		= 0  # == 1 if initialized
  filepath_current		= ABSOLUTE_PATH

  ivy_is_running = False
#  nbi_exposure   = 0    # doubler of NBI_EXPOSURE
  suspended      = 0    # suspend exposure (@hold 1 or @hold 0)

  ivybind_id_list = []  # remember to unbinding

  min_exptime = MIN_EXPTIME
  in_onfileclose = False
  readout = 1
  time_ms = 1000

  def __init__(self):
    self.in_onfileclose = False
    self.exposure_lasts = False
    printd( 'NBI constructor')
    if ivydebug:
      return
    IVYAPPNAME = 'nbi.py'
    IvyInit(IVYAPPNAME,"hi",0)

    self.xsiz_current	= 0
    self.ysiz_current	= 0
    self.xtot		= 0
    self.ytot		= 0
    self.xbin		= 0
    self.ybin		= 0
    self.vbha_current	= 0
    self.vbha_right	= 24.250
    
#    self.in_onexposureend = False

  def get_xtot(self):
    return self.xtot

  def get_ytot(self):
    return self.ytot

  def isCCD3init(self):
    print('isCCD3init:', self.vbha_current, self.vbha_right)
    return self.vbha_current == self.vbha_right

  def init(self):
# move to the constructor later
    print('NBI-INIT')
    if self.ivy_is_running:
      printe('Ivy is running already, skipped')
      return
    
#    IvyStart('')   # ??? Here ???
    delay_time = 0.5
    self.ivybind_id_list.append(IvyBindMsg(self.onxtot,   '(ccd3comm.cam.reply !xtot (.*))'))
    sleep(delay_time)
    print(1)
    self.ivybind_id_list.append(IvyBindMsg(self.onfileclose, '(ccd3comm.file.close)'))
    sleep(delay_time)
    print(2)
    self.ivybind_id_list.append(IvyBindMsg(self.onfilename,  '(ccd3comm.file.name.* (.*))'))
    sleep(delay_time)
    print(3)
    self.ivybind_id_list.append(IvyBindMsg(self.onfilepath,  '(ccd3comm.file.path.* (.*))'))
    sleep(delay_time)
    print(4)
    self.ivybind_id_list.append(IvyBindMsg(self.onexposureend, '(ccd3comm.exposure.end)')) # Doesn't work stable ? Check in the new version ccd3comm!
    sleep(delay_time)
    print(5)
    self.ivybind_id_list.append(IvyBindMsg(self.onexposurestart, '(ccd3comm.exposure.start)'))
    sleep(delay_time)
##    self.ivybind_id_list.append(IvyBindMsg(self.onexposureend, '(ccd3comm.transfer.end)')
    print(6)
    self.ivybind_id_list.append(IvyBindMsg(self.ondrvprogress, '(ccd3comm.cam.reply !drv progress (.*))'))
    sleep(delay_time)
    print(7)
    self.ivybind_id_list.append(IvyBindMsg(self.onexpprogress, '(ccd3comm.cam.reply !tima (.*))'))
    sleep(delay_time)
    print(8)
    self.ivybind_id_list.append(IvyBindMsg(self.onresponse,    '(ccd3comm.con.response (.*))'))
    sleep(delay_time)
    print(9)
    self.ivybind_id_list.append(IvyBindMsg(self.onxbin,   '(ccd3comm.cam.reply !xbin (.*))'))
    sleep(delay_time)
    print(13)
    self.ivybind_id_list.append(IvyBindMsg(self.onybin,   '(ccd3comm.cam.reply !ybin (.*))'))
    sleep(delay_time)
    print(14)
    self.ivybind_id_list.append(IvyBindMsg(self.onytot,   '(ccd3comm.cam.reply !ytot (.*))'))
    sleep(delay_time)
    print(15)
    self.ivybind_id_list.append(IvyBindMsg(self.onxsiz,   '(ccd3comm.cam.reply !xsiz (.*))'))
    sleep(delay_time)
    print(16)
    self.ivybind_id_list.append(IvyBindMsg(self.onysiz,   '(ccd3comm.cam.reply !ysiz (.*))'))
    sleep(delay_time)
    print(17)
    self.ivybind_id_list.append(IvyBindMsg(self.onfres,   '(ccd3comm.cam.reply !fres (.*))'))
    sleep(delay_time)
    print(18)
    self.ivybind_id_list.append(IvyBindMsg(self.onvbha,   '(ccd3comm.cam.reply !vbha (.*))'))
    sleep(delay_time)
    print(19)
    self.ivybind_id_list.append(IvyBindMsg(self.onhold,   '(ccd3comm.cam.reply !hold (.*))'))
    sleep(delay_time)

    print(20)
    self.ivybind_id_list.append(IvyBindMsg(self.on_application_stop,   '(ccd3comm.application.stop)'))
    sleep(delay_time)
    print(21)
    self.ivybind_id_list.append(IvyBindMsg(self.on_application_start,  '(ccd3comm.application.start)'))
    sleep(delay_time)
    print(22)
    self.ivybind_id_list.append(IvyBindMsg(self.on_application_error,  '(ccd3comm.application.error)'))
    sleep(delay_time)
    print(23)

#    self.ivybind_id_list.append(IvyBindMsg(self.oncamreply, '(ccd3comm.cam.reply !tmpa (.*))'))
    self.ivybind_id_list.append(IvyBindMsg(self.oncamreply, '(ccd3comm.cam.reply (.*))'))
    sleep(delay_time)
    print(24)

#    self.ivybind_id_list.append(IvyBindMsg(self.ontemperL, '(ccd3comm.cam.reply !tmpl (.*))'))
#    sleep(delay_time)
#    print('tmpa')
#    self.ivybind_id_list.append(IvyBindMsg(self.onpress,   '(ccd3comm.cam.reply !pres (.*))'))
#    sleep(delay_time)
#    print('pres')
#    sleep(delay_time)

    ivylogger.setLevel(logging.WARN)
    IvyStart('')
    self.ivy_is_running = True

#    sleep(5)
#    self.ivybind_id_list.append(IvyBindMsg(self.ontemperA, '(ccd3comm.cam.reply !tmpa (.*))'))
    print('request params:')
    self.request_params()
    print('End of nbi.init()')

  def write_empty_fits_keyword(self, key, comment):
    IvySendMsg('ccd3comm.con.command keyword 0 TSTRING '+key+' \" \" \"'+comment+'\"')

  def write_fits_keyword(self, key, value, comment, vtype='string'):
#   https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html
    dict = {'STRING' : 'TSTRING', 'INT' : 'TINT', 'DOUBLE' : 'TDOUBLE', 'HMS' : 'TSTRING', 'DMS' : 'TSTRING'}
    T = dict.get(vtype.upper(), 'TSTRING')
    if vtype.upper() == 'HMS': value = deg2hms(value)
    if vtype.upper() == 'DMS': value = deg2dms(value)
    IvySendMsg('ccd3comm.con.command keyword 0 '+T+' '+key+' \"'+str(value)+'\" \"'+comment+'\"')
#    printd( 'ccd3comm.con.command keyword 0 '+T+' '+key+' \"'+str(value)+'\" \"'+comment+'\"')

  def write_fits_keyword254(self, key, value, comment, vtype='string'):
#   https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html
    dict = {'STRING' : 'TSTRING', 'INT' : 'TINT', 'DOUBLE' : 'TDOUBLE', 'HMS' : 'TSTRING', 'DMS' : 'TSTRING'}
    T = dict.get(vtype.upper(), 'TSTRING')
    if vtype.upper() == 'HMS': value = deg2hms(value)
    if vtype.upper() == 'DMS': value = deg2dms(value)
    IvySendMsg('ccd3comm.con.command keyword 254 '+T+' '+key+' \"'+str(value)+'\" \"'+comment+'\"')
#    printd( 'ccd3comm.con.command keyword 254 '+T+' '+key+' \"'+str(value)+'\" \"'+comment+'\"')

  def write_nbi_header(self):
     global OBSERVER_tkvar

     self.write_fits_keyword('ORIGIN', 'SAI MSU', 'Institution name', 'string')
     self.write_fits_keyword('OBSERVAT', 'CMO SAI', 'Observatory name','string')
     self.write_fits_keyword('TELESCOP', 'SAI-2.5', 'Telescope name', 'string')
     self.write_fits_keyword('LONGITUD', 42.6675,  '[deg] telescope geodetic longitude, eastward', 'double')
     self.write_fits_keyword('LATITUD',  43.73611,'[deg] telescope geodetic latitude', 'double')	
     self.write_fits_keyword('ELEVAT', 2112.0, '[m] telescope elevation above sea level', 'double')
     self.write_fits_keyword('INSTRUME',  'NBI4K+FSU', 'Instrument used', 'string')
     self.write_fits_keyword('EQUINOX', 2000.0, 'Standard ICRS', 'string')

     self.write_fits_keyword('OBSERVER', OBSERVER_tkvar.get(), 'Camera operator name', 'string')
     self.write_fits_keyword('IMAGETYP', self.args.imagetype, '', 'string')
     self.write_fits_keyword('FILTER'  , self.args.filter, '', 'string')
     self.write_fits_keyword('FILTERID', self.args.filterid, '', 'string')
     self.write_fits_keyword('PROGID',   self.args.progid, 'Observation program identifier','string')

     self.write_fits_keyword254('CRVAL1', 1, 'Physical value of the reference pixel', 'int')
     self.write_fits_keyword254('CRVAL2', 1, 'Physical value of the reference pixel', 'int')
     self.write_fits_keyword254('CDELT1', 1, 'Size projected into a detector pixel in axis1', 'int')
     self.write_fits_keyword254('CDELT2', 1, 'Size projected into a detector pixel in axis2', 'int')
 
  def calc_exposure_and_reading_time_sec(self):  # relevant after self.load only!

    reading_rate = 400000.0 # px/sec
    timelag = 10.0 # sec
    res = self.xsiz_current * self.ysiz_current / reading_rate + float(self.args.time) + timelag
    return int(round(res))

  def load(self, argstr):
    global ERROR_MESSAGE_tkvar, EXPOSURE_tkvar
# loading exposure params
## f.e.  nbimodule.load(['--fname','autofocusing','--time',str(exptime.get()),'--view','1','--preclear','0','--readout','1']) 
#
    from argparse import ArgumentParser
    parser = ArgumentParser()
#    parser.add_argument('-o','--observer', nargs='?', default='', help='OBSERVER key')
    parser.add_argument('-p','--progid',   nargs='?', default='UNK', help='PROGID key')
    parser.add_argument('-f','--filter',   nargs='?', default='', help='FILTER key')
    parser.add_argument('-l','--filterid', nargs='?', default='', help='FILTERID key')
    parser.add_argument('-n','--fname',    nargs='?', default='tmp', help='output fits file name without \"fits\" extension')
    parser.add_argument('-t','--time',     nargs='?', default=0.002, type=float, help='exposition time in sec')
    parser.add_argument('-i','--imagetype',nargs='?', default='light', help='bias or light')
    parser.add_argument('-e','--fres',     nargs='?', default=0, type=int, help='reset frame format (0 or 1)')

    parser.add_argument('-x','--xbeg',     nargs='?', const=-1, default=-1, type=int, help='xbeg (if -1 - don\'t change)')
    parser.add_argument('-y','--ybeg',     nargs='?', const=-1, default=-1, type=int, help='ybeg (if -1 - don\'t change)')
    parser.add_argument('-w','--xsiz',     nargs='?', const=-1, default=-1, type=int, help='xsiz (if -1 - don\'t change)')
    parser.add_argument('-z','--ysiz',     nargs='?', const=-1, default=-1, type=int, help='ysiz (if -1 - don\'t change)')
    parser.add_argument('-c','--preclear', nargs='?',  const=1,  default=1, type=int, help='preclear image 1 or 0')
    parser.add_argument('-g','--readout',  nargs='?',  const=1,  default=1, type=int, help='readout image after exposition, 1 or 0')

    self.args = parser.parse_args(argstr)
#    printd( "nbimodule: observer:", self.args.observer, "imagetype:", self.args.imagetype)   # observer from epics
    
    self.time_ms = int(round(self.args.time*1000,0))
    self.readout = self.args.readout
    printd( "nbimodule.py: time_ms=", self.time_ms)
    printd( "nbimodule.py: filter:", self.args.filter)

    if EXPOSURE_tkvar.get() != 0:
        str_error = 'Exposure is running already, skipped'
        printe( datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'nbimodule.load: exposure is running already, skipped')
        ERROR_MESSAGE_tkvar.set(str_error)
        return -1

    IvySendMsg('ccd3comm.con.command file %s.fits' % self.args.fname)
    if str.lower(self.args.imagetype) == "dark":
        shut=0
    elif str.lower(self.args.imagetype) == "bias":
        shut=0
        self.time_ms=2
    else:
        shut=1

    IvySendMsg('ccd3comm.con.command @shut %i' % shut)
    
    if self.args.xbeg > 0:
       IvySendMsg('ccd3comm.con.command @xbeg %i' % self.args.xbeg)
    if self.args.ybeg > 0:
       IvySendMsg('ccd3comm.con.command @ybeg %i' % self.args.ybeg)
    if self.args.xsiz > 0:
       IvySendMsg('ccd3comm.con.command @xsiz %i' % self.args.xsiz)
    if self.args.ysiz > 0:
       IvySendMsg('ccd3comm.con.command @ysiz %i' % self.args.ysiz)
    if self.args.fres == 1:
       print('++++++++++++++++++++++ FRES ++++++++++++++++++++++++')
       self.xsiz_current = self.xtot
       self.ysiz_current = self.ytot
       IvySendMsg("ccd3comm.con.command @fres")

#    self.request_params()

    IvySendMsg('ccd3comm.con.command @pclr %i' % self.args.preclear) # 0 or 1, For Normal exposition preclear=1
    IvySendMsg('ccd3comm.con.command ?pclr')
    IvySendMsg('ccd3comm.con.command @rexp %i' % self.readout) # 0 or 1, For Normal exposition redout=1
    IvySendMsg('ccd3comm.con.command ?rexp')
    IvySendMsg('ccd3comm.con.command @time %i' % self.time_ms)
    sleep(1) # !!! TMP
    return 0

  def suspend(self, val):    # 1 - suspend exposure;  0 - resume exposure
    if (val != 0) and (val != 1): return -1
    IvySendMsg('ccd3comm.con.command @hold %i' % val)
    return 0

  def request_params(self):
    print('request_params', 1)
    delay_time = 0.05
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?xtot")
    print('request_params', 2)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?ytot")
    print('request_params', 3)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?xsiz")
    print('request_params', 4)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?ysiz")
    print('request_params', 5)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?tmpl")
    print('request_params', 6)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?tmpa")
    print('request_params', 7)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?pres")
    print('request_params', 8)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?xbin")
    print('request_params', 9)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?ybin")
    print('request_params', 10)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?hold")
    print('request_params', 11)
    sleep(delay_time)
    IvySendMsg("ccd3comm.con.command ?vbha")  # to CCD3INIT
#    sleep(delay_time)
    print('request_params', 12, 'end')

  def start_expose(self):
    global ERROR_MESSAGE_tkvar, EXPOSURE_tkvar, FILEISREADY_tkvar, EXPPROGRESS_tkvar, DRVPROGRESS_tkvar
    if EXPOSURE_tkvar.get() != 0:
        str_error = 'Exposure is running already, skipped'
        printe(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'nbimodule.start_expose: exposure is running already, skipped')
        ERROR_MESSAGE_tkvar.set(str_error)
        return -1
    EXPOSURE_tkvar.set(1)
    FILEISREADY_tkvar.set(0)
    EXPPROGRESS_tkvar.set(0)
    DRVPROGRESS_tkvar.set(0)

    if ivydebug: return 0

##    IvySendMsg("ccd3comm.con.command exp")  Does not work with readout=0 !
##    sleep(0.5)
  
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), '--------------------------- @sint')
#    self.in_onexposureend = False
    self.exposure_lasts = True
    IvySendMsg("ccd3comm.con.command @sint")
    return 0

  def fres(self):
    IvySendMsg("ccd3comm.con.command @fres")

  def cut(self, xbeg, ybeg, xsize, ysize):
# Check if this point is on the right chip and col correct coords
#    global XBIN, YBIN
    printd('nbimodule:cut', xbeg, ybeg, xsize, ysize)
    xbin = self.xbin
    ybin = self.ybin
    if xbin <= 0: xbin = 1
    if ybin <= 0: ybin = 1
# !!! Change here to self.xtot self.ytot ????
    xtot = 4296/xbin
    ytot = 4102/ybin

    wmin = 10
    xhalf = xtot/2
    printd('xhalf =', xhalf)
    if (xbeg > xhalf):  # Right detector
      xbeg = xbeg - xhalf
# And a few additional checks:
    if xbeg < 1: xbeg = 1
    if ybeg < 1: ybeg = 1
    if xbeg > xhalf - wmin: xbeg = xhalf - wmin
    if xbeg + xsize > xhalf+1: xsize = xhalf - xbeg + 1
    if ybeg > ytot - wmin: ybeg = ytot - wmin
    if ybeg + ysize > ytot+1: ysize = ytot - ybeg + 1

    printd('nbi.cut: xbeg, ybeg, xsize, ysize:', xbeg, ybeg, xsize, ysize)

#  set xhalf [expr int(round($::xtot / 2))] 
#  puts "$xhalf $::ytot"
#  set wmin 10
#  if {$::xbeg_desired > $xhalf} {
## Right detector
#    set ::xbeg_desired [expr $::xbeg_desired-$xhalf]
#  }
## Else left detector
#  if {$::xbeg_desired < 1} {
#    set ::xbeg_desired 1
#  }
#  if {$::ybeg_desired < 1} {
#    set ::ybeg_desired 1
#  }
#
#  if {$::xbeg_desired >  [expr $xhalf-$wmin]} {
#    set ::xbeg_desired [expr $xhalf-$wmin]
#  }
#  if {[expr $::xbeg_desired+$::xsize_desired] > [expr $xhalf+1]} {
#    set ::xsize_desired [expr $xhalf-$::xbeg_desired]
#  }
#
#  if {$::ybeg_desired >  [expr $::ytot-$wmin]} {
#    set ::ybeg_desired [expr $::ytot-$wmin]
#  }
#  if {[expr $::ybeg_desired+$::ysize_desired] > [expr $::ytot+1]} {
#    set ::ysize_desired [expr $::ytot-$::ybeg_desired]
#  }

    IvySendMsg('ccd3comm.con.command @xsiz %i' % xsize)
    IvySendMsg('ccd3comm.con.command @ysiz %i' % ysize)
    IvySendMsg('ccd3comm.con.command @xbeg %i' % xbeg)
    IvySendMsg('ccd3comm.con.command @ybeg %i' % ybeg)

  def on_application_error(self, agent, *larg):
    global ERROR_MESSAGE_tkvar
    printe('on_application_error', agent,larg[0])
    msg = larg[0]
    ERROR_MESSAGE_tkvar.set(str_error)

  def on_application_stop(self, agent, *larg):
#    printe('on_application_stop', agent,larg[0])
    self.ccd3comm_initialized = 0

  def on_application_start(self, agent, *larg):
#    printe('on_application_start', agent,larg[0])
    self.ccd3comm_initialized = 1

  def initCCD3(self):
     IvySendMsg('ccd3comm.con.command on_appstart')

  def quitCCD3(self):
    IvySendMsg('ccd3comm.con.command quit')

  def stop(self):
    IvySendMsg('ccd3comm.con.command @timr 0')
    printe('User Stop')
    self.in_onfileclose = False
    global EXPOSURE_tkvar; EXPOSURE_tkvar.set(0)

  def onresponse(self, agent, *larg):
#    printd( 'onresponse', larg[0])
    RESPONSE_tkvar.set(larg[1])

# lll
  def oncamreply(self, agent, *larg):
    global TEMPERL_tkvar, TEMPERA_tkvar, PRESSURE_tkvar
    resp = larg[0].split(' ')
#    print(resp)
    if len(resp) < 2:
      printe( 'strange oncamreply callback')
      return
    if resp[1] == '!tmpl':
        NBI_TEMPERL.put(resp[2])
        TEMPERL_tkvar.set(resp[2])
    elif resp[1] == '!tmpa':
        NBI_TEMPERA.put(resp[2])
        TEMPERA_tkvar.set(resp[2])
    elif resp[1] == '!pres':
        NBI_PRESSURE.put(resp[2])
        PRESSURE_tkvar.set(resp[2])

  def ondrvprogress(self, agent, *larg):
    global DRVPROGRESS_tkvar
#    printd( 'ondrvrogress',larg[0])
    if len(larg) < 2:
      printe('strange ondrvprogress callback')
      return
    DRVPROGRESS_tkvar.set(larg[1])

  def onexposurestart(self, agent, *larg):
    print('on exposurestart')
    log_pm_corr(True)
  
  def onexposureend(self, agent, *larg):
    global EXPOSURE_tkvar
    print('onexposureend, self.exposure_lasts =', self.exposure_lasts)
    if not self.exposure_lasts:
      return
    self.exposure_lasts = False
#    if self.in_onexposureend:
#      return
#    self.in_onexposureend = True
    log_pm_corr(False)
    printd(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), ' ------------------ onexposureend', larg[0])
    if self.readout == 0:
      print('nbimodule: -------- end of exposure without readout')
      sleep(0.5)
      EXPOSURE_tkvar.set(0)
#      self.nbi_exposure = 0
#    self.in_onexposureend = False


  def onfileclose(self, agent, *larg):
#    printd( 'onfileclose',agent,larg[0])
    if self.in_onfileclose:
      printe('nbimodule.onfileclose: repetead enter, skipped')
      return
    self.in_onfileclose = True
#    printd(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), ' ------------------ onefileclose', larg[0])
    global EXPOSURE_tkvar; EXPOSURE_tkvar.set(0)
#    self.nbi_exposure = 0
    global FILEISREADY_tkvar; FILEISREADY_tkvar.set(1)
#    print('Stop on fileclose')
    self.in_onfileclose = False

  def onexpprogress(self, agent, *larg):
    if len(larg) < 2:
      printe('strange onexpprogress callback')
      return
    ms = float(larg[1])
    progress_sec = round(ms / 1000)
    global EXPPROGRESS_tkvar; EXPPROGRESS_tkvar.set(progress_sec)

  def onfilename(self, agent, *larg):
    if len(larg) < 2:
      printe( 'strange onfilename callback')
      return
    global FILENAME_tkvar; FILENAME_tkvar.set(larg[1])

  def onfilepath(self, agent, *larg):
    if len(larg) < 2:
      printe( 'strange onfilepath callback')
      return
    self.filepath_current = larg[1]   # Need we that class variable?

  def ontemperL(self, agent, *larg):
#    print('ontemperL')
    global TEMPERL_tkvar
    if len(larg) < 2:
      printe( 'strange ontemperL callback')
      return
    TEMPERL_tkvar.set(larg[1])
    NBI_TEMPERL.put(larg[1])
    self.ccd3comm_initialized = 1

  def ontemperA(self, agent, *larg):
    global TEMPERA_tkvar
    if len(larg) < 2:
      printe( 'strange ontemperA callback')
      return
    NBI_TEMPERA.put(larg[1])
    TEMPERA_tkvar.set(larg[1])

  def onhold(self, agent, *larg):
#    print('onhold')
    if len(larg) < 2:
      printe('strange onhold callback')
      return
    self.suspended = int(larg[1])

  def onpress(self, agent, *larg):
    global PRESSURE_tkvar
    if len(larg) < 2:
      printe('strange onpress callback')
      return
    PRESSURE_tkvar.set(larg[1])    
    NBI_PRESSURE.put(larg[1])

  def onxbin(self, agent, *larg):
    if len(larg) < 2:
      printe( 'strange onxbin callback')
      return
    self.xbin = int(larg[1])

  def onybin(self, agent, *larg):
    if len(larg) < 2:
      printe( 'strange onybin callback')
      return
    self.ybin = int(larg[1])

  def onxtot(self, agent, *larg):
#    print('onxtot')
    if len(larg) < 2:
      printe('strange onxtot callback')
      return
    self.xtot = int(larg[1])

  def onytot(self, agent, *larg):
#    print('onxtot')
    if len(larg) < 2:
      printe('strange onytot callback')
      return
    self.ytot = int(larg[1])
#    print('onytot')

  def onxsiz(self, agent, *larg):
    if len(larg) < 2:
      printe('strange onxsiz callback')
      return
    self.xsiz_current = int(larg[1])
#    print('onxsiz')

  def onysiz(self, agent, *larg):
    if len(larg) < 2:
      printe( 'strange onysiz callback')
      return
    self.ysiz_current = int(larg[1])
#    print('onysiz')

  def onfres(self, agent, *larg):
    self.xsiz_current = self.xtot
    self.ysiz_current = self.ytot

  def onvbha(self, agent, *larg):
    if len(larg) < 2:
      printe( 'strange onvbha callback', larg)
      return
    li = larg[1].split('  ')  # !!! two spaces here !!!
    if len(li) < 2:
      printe( 'strange onvbha callback, args =', larg)
      return
    self.vbha_current = float(li[1])
#    print('self.vbha_current =', self.vbha_current)

  def stop_ivy(self):
    if not self.ivy_is_running:
      printe('Ivy is not running, skip STOP IVY')
      return
    for id in self.ivybind_id_list:
      IvyUnBindMsg(id)
    self.ivybind_id_list = []

    IvyStop()
    printd('Ivy has stopped')
    self.ivy_is_running = False


#if __name__ == '__main__':
#    init()
#    expose(sys.argv[1:])
#    sys.exit(0)

