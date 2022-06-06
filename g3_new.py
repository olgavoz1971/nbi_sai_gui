#!/usr/bin/python3

myICS = 'NBI'
debug_focus = False

import sys
import subprocess
import psycopg2

# For scenarios
import threading
import queue
# ----

import debug
#import imageserver_adapter as imageserver
import upload2bd
import BDmeasurements
from iq2bd import iq
import upload_iq

import focus_module
import numpy as np
from jkclass import SCALE
from jkclass import PV, correct_outfile, composeSAIoutfile
from jkclass import MAX_PV_LENGTH
from jkclass import PVTkbase, IndicatorStatus, IndicatorStatusTkvar, JKLabel
from jkclass import printd, printdd, printe
from jkclass import deg2hms, deg2dms, sec2deg
from jkclass import round2float
from jkclass import DS9NBI


from astropy.coordinates import Angle
import astropy.units as u
from math import isnan
from os import system

from datetime import datetime

#import tkFileDialog
from tkinter import filedialog as tkFileDialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

global mainApp

epics = debug.epics
debug_nbi = debug.debug_nbi

SAI25_STATUS    = PV('SAI25:STATUS')
#SAI25_STATUS    = PV('NBI:XBIN')

FSU_FILTER      = PV('FSU:FILTER')
FSU_SHUTTER     = PV('FSU:SHUTTER')
FSU_MIRROR      = PV('FSU:MIRROR')
FSU_LAMP        = PV('FSU:LAMP')
FSU_STATUS      = PV('FSU:STATUS')
FSU_CMD_RUN     = PV('FSU:CMD_RUN')
SAI25_OPERATOR  = PV('SAI25:OPERATOR')
SAI25_PROGID    = PV('SAI25:PROGID')
SAI25_ICS       = PV('SAI25:ICS')
AGU_C1_TRAN     = PV('SAI25:AGUTRAN')
AGSTATUS        = PV('SAI25:AGSTATUS')

M2_POS          = PV('SAI25:M2')
#DERO            = PV('TCS2:C1_POS')
#PARANG          = PV('TCS2:PARANG')
DERO            = PV('SAI25:DERO')
PARANG          = PV('SAI25:PARANG')
POSANG          = PV('SAI25:POSANG')

class FitsHeader():

  keys = (\
  ('OBJECT'    , PV('SAI25:OBJECT'  ),	'Target as given by the user', 'string'),\
#  ('PROGID'    , PV('SAI25:PROGID'  ),	'Observation program identifier','string')\
  ('TARRA'     , PV('SAI25:TARRA'   ),	'[hour] target RA2000 at julian day start','hms'),\
  ('TARDEC'    , PV('SAI25:TARDEC'  ),  '[deg] target DEC2000 at julian day start', 'dms'),\
  ('TARPMRA'   , PV('SAI25:TARPMRA' ),  '[\"/s] target proper motion set in RA', 'double'),\
  ('TARPMDEC'  , PV('SAI25:TARPMDEC'),  '[\"/s] target proper motion set in Dec', 'double'),\
  ('AZ'        , PV('SAI25:AZ'      ),  '[deg] Mount azimuth at exposure start', 'double'),\
  ('ALT'       , PV('SAI25:ALT'     ),  '[deg] Mount elevation at exposure start', 'double'),\
  ('DERO'      , PV('SAI25:DERO'    ),  '[deg] Active rotator angle', 'double'),\
  ('POSANG'    , PV('SAI25:POSANG'  ),  '[deg] Instrument position angle', 'double'),\
  ('INSTANG'   , PV('SAI25:INSTANG' ),  '[deg] Rotator angle for vertical instrum./slit', 'double'),\
  ('PARANG'    , PV('SAI25:PARANG'  ),  '[deg] Parallactic angle at exposure start', 'double'),\
  ('RA'        , PV('SAI25:RA'      ),  '[hour] Telescope pointing at exp start', 'hms'),\
  ('DEC'       , PV('SAI25:DEC'     ),  '[deg] Telescope pointing at exp start', 'dms'),\
  ('TELSTATE'  , PV('SAI25:STATUS'  ),	'Telescope control system status','string'),\
  ('M2'        , PV('SAI25:M2'      ),	'[mm] M2 mirror position (focus)','double',),\
  ('M3'        , PV('SAI25:M3'      ),	'Nasmyth mirror position identifier','string'),\
  ('M1COVERF'  , PV('SAI25:COVER'   ),	'M1 mirror opening fraction','double'),\
  ('AGSTATUS'  , PV('SAI25:AGSTATUS'),  'Autoguiding status (guiding==LOOP)', 'string'),\
  ('AGUGSID'   , PV('SAI25:AGUGSID' ),  'Guiding star identifier', 'string'),\
  ('AGUGSFX'   , PV('SAI25:AGUFLUX' ),  '[ph/s] Guiding star median', 'flux'),\
  ('AGUROT'    , PV('SAI25:AGUROT'  ),	'[deg] Autoguider rotator angle','double'),\
  ('AGUTRAN'   , PV('SAI25:AGUTRAN' ),  '[mm] Autoguider translator shift wrt FoV centre','double'),\
  ('DOMEAZ'    , PV('SAI25:DOMEAZ'  ),	'[deg] Dome slit azimuth, from south westward','double'),\
  ('DOMESLIT'  , PV('SAI25:DOMESLIT'),	'Dome slit state','double'),\
  ('TM1'       , PV('SAI25:TM1'	    ),	'[deg C] Primary mirror temperature','double'),\
  ('TTUBE'     , PV('SAI25:TTUBE'   ),	'[deg C] Telescope tube temperature','double'),\
  ('TOUT'      , PV('SAI25:TOUT'    ),	'[deg C] Outside temperature','double'),\
  ('ATMPRESS'  , PV('SAI25:ATMPRESS'),	'[mbar] Atmospheric pressure','double'),\
  ('RELHUM'    , PV('SAI25:RELHUM'  ),	'[%] Outside relative humidity','double'),\
  ('DEWPOINT'  , PV('SAI25:DEWPOINT'),	'[deg C] Outside dew point','double'),\
  ('WIND'      , PV('SAI25:WIND'    ),	'[m/s] Outside average wind speed (Vaisala)','double'),\
  ('WINDDIR'   , PV('SAI25:WINDDIR' ),	'[deg] Outside average wind direction, geodetic','double'),\
  ('OPERATOR'  , PV('SAI25:OPERATOR'),  'Telescope operator name', 'string'),\
  )

  def write_all_keys(self):
    print('write_all_keys')
    for key,pv,comment,vtype in self.keys:
#      printd( pv.pvname, comment)
      self.write_key(key, pv, comment, vtype)

  def write_key(self, key, pv, comment, vtype):
    value = pv.get()

    if (pv.type.find('str') >= 0):
       if (pv.get().upper() == 'NAN') or (pv.get().upper() == ''):
         return
    elif isnan(pv.get()):
       return
#      nbi.write_empty_fits_keyword(key, comment)
    nbi.write_fits_keyword(key, value, comment, vtype)

#import pyds9  
#from pyds9 import DS9 as ds9  
#class DS9NBI:
#
#  def find_ds9(self):
#     try:
#       ds9_list = pyds9.ds9_targets()
#       if(len(ds9_list) > 1): d = ds9('DS9:'+DS9_NAME)
#       else: d = ds9()
#     except Exception as e:
#       printe('DS9NBI:find_ds9 exception ',e)
#       return None
#     return d
#
#  def get_path_to_file(self):     # only local path :(
#       fullfilename = None
#       d = self.find_ds9()
#       if d is None: return fullfilename
#       fullfilename = d.get('file').replace('[im1]','')
#       return fullfilename
#
#  def view(self, fullpath):
#       d = self.find_ds9()
#       try:
#            d.set("file " + fullpath)
#            posang = d.get('fits header keyword ' + 'posang')
#            posang = geometry.correct_posang(posang)
#            print('posang=|'+str(posang)+'|', type(posang))
#            if posang == '': posang = float('NaN')
#            else:
#              posang = float(posang)
##            posang = POSANG.get()
#            if not isnan(posang):
#              vangle = Angle(90,u.degree) + Angle(posang,u.degree)
#              vangle.wrap_at('360d', inplace=True)
#              vlen = 500
##              printd('+++++++++++++++++++ d.get("fits width"), d.get("fits height")', d.get("fits width"), d.get("fits height"))
#              vxstart = float(d.get("fits width"))/2
#              vystart = float(d.get("fits height"))/2
#              vtext = 'N'
#              north_arrow = 'physical; # vector('+str(vxstart)+','+str(vystart)+','+str(vlen)+','+str(vangle.degree)+') vector=1 width=3 color=orange font=\"helvetica 14 bold roman\" text={N}'
##              
##              north_arrow_ics_corrected = 'physical; # vector('+str(vxstart)+\
##               ','+str(vystart)+','+str(vlen)+','+str(vangle.degree)\
##               +') vector=1 width=3 color=green font=\"helvetica 14 bold roman\" text={N}'
#
#              d.set('regions', north_arrow)            
#       except ValueError as e:
#            printe('g1:view() ds9: Exception occured: ', e)
##     except ValueError as e:
##        printe('g1:view() ds9: Exception occured: ', e)


#dedug = True

# NBIserver (nbiserver.py) is running on the binary (nbi-machine)
# Use for active command to OCS
IP_NBIserver = '127.0.0.1'
port_NBIserver = 33001
#  image_server (image_server.tcl) is running on the tower
# Use it for long calculations
#IP_image_server = '192.168.15.51'
#port_image_server = 33000
from time import sleep
import os.path

import tkinter as Tk
pdx=1
pdy=1

# ---------------- Queue buttons stuff ----------
ready_color        = 'SeaGreen3' #'green3'
ready_color_active = 'SeaGreen1' #'green2'
busy_color         = 'dark turquoise'
busy_color_active  = 'turquoise'

from collections import namedtuple
Task = namedtuple('Task', ['on_start', 'on_end', 'arglist'], verbose=False)
Task.__new__.__defaults__ = (None,None,[])

lwidth = 10
class ButtonQueue(Tk.Button):
  def __init__(self, parent, queue, text='', queue_command=None, on_end_command=None, queue_tkvar_list=[], *args, **kwargs):
# on_end_command executes in the Main GUI thread
    Tk.Button.__init__(self, parent, text=text, activebackground=ready_color_active,\
          bg=ready_color, command=self.run_task, *args, **kwargs)
    self.queue = queue
    self.queue_command    = queue_command
    self.queue_tkvar_list = queue_tkvar_list
    self.on_end_command   = on_end_command

  def run_task(self): # in the worker thread
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'run task in the thread start')
# Read tkvars:
    arglist = []
    for tkvar in self.queue_tkvar_list:
         print('name =', str(tkvar), tkvar.get())
         arglist.append(tkvar.get())
    print('arglist=', arglist)
    task = Task(on_start=self.queue_command, on_end=self.on_end, arglist=arglist)
    self.configure(bg=busy_color, activebackground=busy_color_active)
    self.queue.put(task)
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'run task in the thread end')

  def on_end(self): # in the main GUI thread
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'on end start')
    if self.on_end_command is not None:
      try:
        self.on_end_command()
      except Exception as e:
        printe('ButtonQueue:on_end', e)
    self.configure(bg=ready_color, activebackground=ready_color_active)
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'on end end')

class ProcessQueue(): # like ButtonQueue but not a button
  def __init__(self, queue, queue_command=None, on_end_command=None, queue_tkvar_list=[], *args, **kwargs):
# on_end_command executes in the Main GUI thread
    self.queue = queue
    self.queue_command    = queue_command
    self.queue_tkvar_list = queue_tkvar_list
    self.on_end_command   = on_end_command

  def run_task(self): # in the worker thread
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'run task in the thread start')
# Read tkvars:
    arglist = []
    for tkvar in self.queue_tkvar_list:
         print('name =', str(tkvar), tkvar.get())
         arglist.append(tkvar.get())
    print('arglist=', arglist)
    task = Task(on_start=self.queue_command, on_end=self.on_end, arglist=arglist)
#    self.configure(bg=busy_color, activebackground=busy_color_active)
    self.queue.put(task)
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'run task in the thread end')

  def on_end(self): # in the main GUI thread
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'on end start')
    if self.on_end_command is not None:
      try:
        self.on_end_command()
      except Exception as e:
        printe('ProcessQueue:on_end', e)
#    self.configure(bg=ready_color, activebackground=ready_color_active)
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'on end end')

# ------------------------------------

def unfold1(s):   # work with []   f.e. '2 * [L,300 + 2*[B+D]+F,B+L, 1+2*[F,B+L,1]]'
    o = s.rfind('[')
    c = s.find(']',o)
    if c < 0:
      str_error = 'Wrong scenario '+s
      mainApp.messages.error_message.set(str_error)
      printe(s)
      return ''
    if s[o-1] != '*':
      str_error = 'Wrong scenario '+s
      mainApp.messages.error_message.set(str_error)
      printe( s)
      return ''
    sub = s[o+1:c]
#    print 'sub =',sub,o,c
    p = s[:o].rfind('+')
    try:
      num = int(s[p+1:o-1])
    except Exception as e:
      str_error = 'Scenario exception '+s
      mainApp.messages.error_message.set(str_error)
      printe('Scenario exception ', s, e)
      return ''
    newsub = ''
    for i in range(0,num): newsub = newsub + sub + '+'
    newsub = newsub[:-1]
#    print 'newsub =', newsub
    sub2replace = s[p+1:c+1]
#    print 'sub2replace=', sub2replace
    return s.replace(sub2replace, newsub)

def unfold(stri):   # f.e. '2 * [L,300 + 2*[B+D]+F,B+L, 1+2*[F,B+L,1]]'
    stri = stri.replace(' ','')
    while True:
      printd(stri)
      if stri.find('[') < 0: break
      stri = unfold1(stri)
      if stri == '':  return stri
    stri = stri.replace('++','+')
    return stri


class ScenarioStep():
  def __init__(self, name, method):
     self.name = name
     self.method = method

class Scenario:
  def __init__(self):

    self.fitsheader = FitsHeader()

    self.iq_dict = {}
    self.processQueue = ProcessQueue(queue=nbi_command_queue, queue_command=self.do_it_in_the_thread_when_file_is_ready, \
                      on_end_command=self.set_iq_tkvars)

    self.Exposition = ScenarioStep('Expose',   self.expose)
    self.Filter     = ScenarioStep('Filter',   self.set_filter)
    self.Focus      = ScenarioStep('Focus',    self.move_focus)
    self.ShiftEqu   = ScenarioStep('ShiftEqu', self.shift_equ)
    self.ShiftHor   = ScenarioStep('ShiftHor', self.shift_hor)

    self.scenario_is_running = False
    self.long_scenario_is_running = False   # Real ling scenario - BD or user scn-file scenarios
    self.exposure_and_reading_time_ms = 0
    self.auto_target_done = False
    self.auto_start       = False

    self.scenario_step_dict = {\
      'Exp'     : self.Exposition,\
      'Filter'  : self.Filter,\
      'Focus'   : self.Focus,\
      'ShiftEqu': self.ShiftEqu,\
      'ShiftHor': self.ShiftHor\
      
    }
 
    self.step_collection = [self.Exposition, self.Filter, self.Focus, self.ShiftEqu, self.ShiftHor]
    self.siterator  = 0  # help callbacks to go through the scenario
    self.scenario   = []   #    of(Step, argument) paires
    self.step_name  = 'MysteriousStep'
    self.stop_tkvar = Tk.BooleanVar(value=False)

# Smart sleep

    self.exposure_wait_variable		= nbimodule.EXPOSURE_tkvar
    self.fileready_wait_variable	= nbimodule.FILEISREADY_tkvar
    self.fsu_wait_variable              = Tk.IntVar(value=0)

    self.current_wait_variable     = -1 # must be equal to one of the upper wait_variables
    self.timer_id                  = -1   # one for all steps

    nbimodule.ERROR_MESSAGE_tkvar.set('')
  
  def initExchange(self):
    FSU_CMD_RUN.add_callback(self.pvcallback_filter, with_ctrlvars=False)

  def load(self, scenario_string_wrapped):
    self.siterator = 0
    self.scenario = []
    scenario_string = unfold(scenario_string_wrapped)
    printdd(scenario_string)
    if scenario_string == '': return -1
    for step_string in scenario_string.split('+'):
      coma = step_string.find(',')
      step = self.scenario_step_dict.get(step_string[:coma], -1)
      printd( 'step_string =', step_string, step_string[:coma])
      if step == -1:
        mainApp.messages.error_message.set('Wrong scenario '+scenario_string_wrapped)
        printe( 'Wrong scenario '+scenario_string_wrapped)
        return -1
      arg = step_string[coma+1:]
      self.scenario.append((step, arg))

#    self.scenario = [(self.Filter, '0'), (self.Exposition, (Light,100.0),(self.Filter, '0'), (self.Exposition, 300.0)]
#    self.scenario = [(self.Exposition, 1.0)]

  def run(self, auto_target_done=False, auto_start=False):
    printdd('Scenario:Run')
    self.auto_target_done = auto_target_done
    self.auto_start = auto_start
    if self.scenario_is_running:
      mainApp.messages.error_message.set('Scenario is running already')
      printe('Scenario:run Scenario is running already')
      return

    self.scenario_is_running = True

    mainApp.messages.error_message.set('')
    mainApp.messages.current_message.set('Ready')
    self.stop_tkvar.set(False)

    thread_scenario = threading.Thread(target=self.execute_it_in_the_thread)
    thread_scenario.start()
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Scenario thread start')

  def execute_it_in_the_thread(self):    # Scenario
    # put scenario into The Thread

    my_error = False
    if self.auto_start: sleep(10)
    for ent in self.scenario:
       if self.stop_tkvar.get(): break        #! check and handle ! Stop
       if nbi.ccd3comm_initialized != 1:
          str_error = 'ccd3comm does not work, Start it first'
          printe(str_error)
          mainApp.messages.error_message.set(str_error)
          my_error = True
          break

       step = ent[0]
       argument = ent[1]
       self.step_name = step.name
       method = step.method
       print('Scenario:Run: name, method, argument', self.step_name, method, argument, 'error =', my_error)
       try:
         res = method(argument)   # Wait in each method
         print('execute_it_in_the_thread try: res =',res, type(res), 'error =', my_error)
#         print('\n\n\n 2222 Set my erreor True \n\n\n')
         if res != 0: my_error = True # m.b. 0 - ok, -1 - error or errcode > 1
         print('res =', res, res != 0, my_error)
# lll
         if self.step_name.upper() == 'EXPOSE' and scenario.long_scenario_is_running:
# TMP !!! Check case of individuals and multi expositions !!!
            gui_queue.queue.put(mainApp.expose_frame.sci_scenario_table.next)   # to the main GUI thread
#             mainApp.expose_frame.sci_scenario_table.next()
       except Exception as e:
         print('Exception:', e, 'res =', res, type(res))
         printe('\n\n\n Scenario in the thread exception: ', e, 'Set Error True\n\n\n')
         my_error = True
    self.end_of_scenario()
    print('error = ', my_error, 'self.auto_target_done: ', self.auto_target_done)
    if self.auto_target_done and not my_error:  icsmodule.targetDone()

    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Scenario thread done')
#   completed
 

  def suspend(self, value):  # Scenario class method
    if value == 1:
      nbi.suspend(1)
      # fix exposure elapsed time here
      self.exp_progress_on_suspend = nbimodule.EXPPROGRESS_tkvar.get()
      self.datetime_suspended = datetime.now()
      # stop timer here!
      mainApp.after_cancel(self.timer_id)
    else:
      nbi.suspend(0)
      # recalculate self.expose_timeout
      self.exposure_and_reading_time_ms -= int(round(self.exp_progress_on_suspend * 1000)) # must be integer!
      print('-------------- new exposure_and_reading_time_ms = ', self.exposure_and_reading_time_ms)
      # release timer here!
      self.timer_id = mainApp.after(self.exposure_and_reading_time_ms, self.expose_timeout)


  def stop(self, err):
    printd( 'Stop Scenario err=',err, type(err), '\n')
    self.siterator = len(self.scenario)
    self.stop_tkvar.set(True)
    if err < 0:
      printe( 'Stop scenario. An error occured, code=', err, mainApp.messages.error_message.get(), \
           nbimodule.ERROR_MESSAGE_tkvar.get())
#      self.end_of_scenario()
    return err
    
  def end_of_scenario(self):
    printd( 'End_of_scenario')
    mainApp.messages.current_message.set('Done')
    self.scenario_is_running = False
    self.long_scenario_is_running = False

  def step_timeout(self):
    printe('--------- !!!!!!!!!! -----', self.step_name, 'timeout')
    mainApp.messages.error_message.set(self.step_name+' timeout')
    mainApp.after_cancel(self.timer_id)
    self.current_wait_variable.set(-1)
    self.stop(-1)

  def do_it_when_file_is_ready(self):   # Scenario    in the thread !
    value = nbimodule.FILEISREADY_tkvar.get()
    printd( 'Scenario:file_is_ready', value)
    if value != 1: return
    printd( 'Exposure File is ready\n')
    
# Fix filename here BEFORE new exposure start !!!!!!!!!!!!!!!!
    self.fixed_file_name = nbimodule.FILENAME_tkvar.get()
    self.fixed_file_path = nbi.get_filepath()
    
    sleep(5)
    self.processQueue.run_task()  # New thread for save exposure/telescope time
# And wait a while. Oterwise CCD may overflow itself somehow

  def do_it_in_the_thread_when_file_is_ready(self):   # What to do in thread to save telescope time (in the Process task)
    fullpath = os.path.join(self.fixed_file_path, self.fixed_file_name)
    exist = os.path.exists(fullpath)
    if not exist:
        err_mess = 'File '+fullpath+' does not exist'
        mainApp.messages.error_message.set(err_mess)
        printe( 'File', fullpath, 'does not exist')
        return

    if mainApp.ocs.ds9view.get() == 1:
         printd ('\nds9view.get()=', mainApp.ocs.ds9view.get())
         ds9nbi.view(fullpath)
    if mainApp.ocs.bd.get() == 1:
        upload2bd.upload(self.fixed_file_name, self.fixed_file_path)

    if mainApp.ocs.iq.get() == 1:
       print('self.calc_iq(fullpath), fullpath =', fullpath)
       self.calc_iq(fullpath)
       print('After  self.calc_iq')

    print('+++++++++++End of do_it_in_the_thread_when_file_is_ready++++++++++++')

  def calc_iq(self, fullname):  # in the worker thread
     print('\n\n ------------- CALC IQ ----------\n\n')
# lll
     res = upload_iq.check_fits(fullname)
     if res < 0: return res
   
     print('CALC IQ: check_fits: OK')
     res = iq.do(fullname)
     print('CALC IQ: iq.do res =', res)
     if type(res) is int:
        print(fullname+' : iq calculation error')
        return res
     
     self.iq_dict = res       
     
     path, filename = os.path.split(fullname)
     print('CALC IQ: upload_iq...', filename, self.iq_dict)
     upload_iq.upload_iq(filename, self.iq_dict)
     print('CALC IQ: uploaded')

     return 0
      
  def set_iq_tkvars(self):  # In the main gui thread
     print('set_iq_tkvars')
     print('self.iq_dict:', self.iq_dict)
     mainApp.expose_frame.iq_window.set_tkvars(self.iq_dict)
#     self.fwhma_tkvar.set(self.iq_dict['fwhma'])
# lll

  def expose(self, argstring):
#  Strong arguments order:
#      1 - exptime(obligatory)
#      2 - imagetype      default = 'Light'
#      3 - userfilename   default = '' (SAI-style filename f.e.'20200517121418'). *tmp* files will not be put into BD, if 'focusing' - will be  rewrited
#      4 - preclear (0/1) default = 1
#      5 - readout (0/1)  default = 1

    self.exposure_and_reading_time_ms = 0
    printd( 'Scenario:expose')
    imagetype_valid_list = ['LIGHT', 'BIAS', 'FLAT', 'DARK']
    arglist = argstring.split(',')
    try:
      exptime = float(arglist[0])
    except Exception as e:
      mainApp.messages.error_message.set('Exposition time '+arglist[0]+ ' error')
      printe( 'Exposition time '+arglist[0]+ ' error ',e)
      self.stop(-1)
      return -1

    imagetype = 'light'
    userfilename = ''
    preclear = 1
    readout  = 1

    if len(arglist) > 1:
      imagetype = arglist[1]
      if imagetype.upper() not in imagetype_valid_list:
        mainApp.messages.error_message.set('Wrong imagetype '+imagetype)
        printe( 'Wrong imagetype '+imagetype)
        return -1

    if len(arglist) > 2:
      userfilename = arglist[2]

    if len(arglist) > 3:
      preclear = int(arglist[3])
      if preclear not in [1, 0]:
        printe('Wrong preclear value=', preclear)
        return -1

    if len(arglist) > 4:
      readout = int(arglist[4])
      if readout not in [1, 0]:
        printe('Wrong readout value=', readout)
        return -1

    if nbimodule.EXPOSURE_tkvar.get() != 0:
      printe( 'Exposure is running already, skiped')
      mainApp.messages.error_message.set('Exposure is running already')
      self.stop(-1)
      return -1

    if userfilename == '':
      printd( "SAI style")
      userfilename = composeSAIoutfile()
    else:
#      if userfilename != 'focusing':
      userfilename = correct_outfile(nbi.get_filepath(), userfilename)
      printd ('Not SAI style', userfilename)

#    progid = SAI25_PROGID.get()
    progid = mainApp.ocs.progid.get()
    if not mainApp.ocs.change_progid.get():
      if (imagetype.lower() == 'dark') or (imagetype.lower() == 'bias') or (imagetype.lower() == 'flat'): 
        progid = 'calibration'
    if progid == '': progid = 'Unknown'

    if AGU_C1_TRAN.get() < agutran_min.get():
#      printe('AGU_C1_TRAN.get() agutran_min.get()', AGU_C1_TRAN.get(), agutran_min.get())
      if (imagetype.lower() != 'dark') and (imagetype.lower() != 'bias') and (SAI25_ICS.get() == myICS):
        str_error = 'Take off AGU paw!'
        printe(str_error)
        mainApp.messages.error_message.set(str_error)
        self.stop(-1)
        return -1
    arg_list = ['--fname', userfilename,  '--time', str(exptime),  '--preclear',str(preclear),  '--readout', str(readout),\
    '--imagetype', imagetype, '--progid', progid, '--filter', fsu.get_filter_name(),\
    '--filterid', FSU_FILTER.get(), '--fres','0']

    printd( 'arg_list =', arg_list)
    nbi.load(arg_list)
    if not nbi.isCCD3init():
      str_error = 'ccd3comm has not initiated, Press InitCCD3 first'
      printe( str_error)
      mainApp.messages.error_message.set(str_error)
      self.stop(-1)
      return -1
    nbi.write_nbi_header()
    timelag = 1.0
    readtime = 40

    self.fitsheader.write_all_keys()

    try:
#     must be integer
      self.exposure_and_reading_time_ms = nbi.calc_exposure_and_reading_time_sec() * 1000  # in ms
#      expose_and_reading_time = int((float(exptime) + timelag + readtime)*1000) # in ms
    except Exception as e:
      printe(e)
      printe(exptime, timelag, readtime)
      mainApp.messages.error_message.set(str(exptime)+'+'+str(timelag)+'+'+str(readtime))
      self.stop(-1)
      return -1

    print('-------------- exposure_and_reading_time_ms = ', self.exposure_and_reading_time_ms)
    nbi.start_expose()
    if readout == 0:
      current_wait_variable = self.exposure_wait_variable
    else: current_wait_variable = self.fileready_wait_variable
    ret = self.wait_step(self.exposure_and_reading_time_ms, current_wait_variable, timeout_callback=self.expose_timeout)
    if ret < 0: return ret
# Put it to the OTHER thread !
# lll
    self.do_it_when_file_is_ready()  # in the thread !

    return 0

#  Charge The end of exposure is callback     

  def pvcallback_filter(self, pvname=None, value=-1, char_value=-1, **kw):
     printd( 'Scenario: pvcallback_filter', value)
     if value == 0:
       printdd('Callback fsu: Unlock')
       self.fsu_wait_variable.set(0)
     printd( 'Callback: filter is ready value=', value)


  def expose_timeout(self):
    print('expose_timeout ------------- ')
    self.step_timeout()
    nbimodule.EXPOSURE_tkvar.set(0)

  def set_filter(self, filter):
    printd( 'Scenario : set_filter and mirror', filter)

    res = fsu.set_mirror_target(mainApp.filters.mirror_target.get())
    if res < 0: return self.stop(res)

    res = fsu.set_lamp_target(mainApp.filters.lamp_target.get())
    if res < 0: return self.stop(res)

    res = fsu.set_shutter_target(mainApp.filters.shutter_target.get())
    if res < 0: return self.stop(res)

    print('fsu.set_filter_target')
    res = fsu.set_filter_target(filter)
    print('Yes fsu.set_filter_target, res =', res)
    if res < 0: return self.stop(res)
    print('Move fsu')
    res = fsu.move_fsu()
    print('Yes! res =', res)
    if res < 0: return self.stop(res)
    set_time = fsu.max_move_time_sec*1000 # in ms
    print('set_time =', set_time, ' ms')
    return self.wait_step(set_time, self.fsu_wait_variable)

  def wait_step(self, timeout_time, step_wait_variable, timeout_callback=None):
#    step_timeout - method, timeout callback
    if timeout_callback is None: timeout_cb = self.step_timeout
    else: timeout_cb = timeout_callback

    self.timer_id = mainApp.after(timeout_time, timeout_cb) # timeout_time in usec
    self.current_wait_variable = step_wait_variable # for timeout cancel
    print('Wait step ', self.step_name)
#   Wait here until step will be ready
    mainApp.messages.current_message.set('Wait '+ self.step_name + '...')

    mainApp.wait_variable(step_wait_variable)
    mainApp.after_cancel(self.timer_id)
    print('Wait step: Unlocked, step_wait_variable =', step_wait_variable.get())
    return int(step_wait_variable.get())

  def move_focus(self, m2pos):
    printd( 'Scenario:move_focus', m2pos)
    res = icsmodule.setFocus(m2pos)
    printd ('Scenario: move_focus before wait res =', res, type(res))
    if res < 0:
      return self.stop(res)      

    set_time = 60000 # in ms    for timeouot
    res = self.wait_step(set_time, icsmodule.OCS_FOCUS_RUN_tkvar)
    print('Scenario: move_focus after wait: res=', res, type(res))
    if res != 0: 
      str_error = 'Move focus: error code' + str(res)
      mainApp.messages.error_message.set(str_error)
      self.stop(res)
    return res

  def shift_equ(self, argstring):   # args in degrees
    printd( 'Scenario:shift_equ', argstring)
    arglist = argstring.split(',')
    try:
      dra  = float(arglist[0])
      ddec = float(arglist[1])
    except Exception as e:
      mainApp.messages.error_message.set('ShiftEqu '+arglist+ ' error')
      printe('Exposition time '+arglist + ' error ',e)
      self.stop(-1)
      return -1
    res = icsmodule.shiftEqu(dra, ddec, 0)
    printd('shift_equ res=', res)
    if res < 0:
      return self.stop(res)     
    set_time = 60000 # in ms
    res = self.wait_step(set_time, icsmodule.OCS_EQU_RUN_tkvar)
    if res != 0: 
      str_error = 'Shift Equ: error code ' + str(res)
      mainApp.messages.error_message.set(str_error)
      self.stop(res)
    return res

  def shift_hor(self, argstring):  # args in degrees
    printd( 'Scenario:shift_hor', argstring)
    arglist = argstring.split(',')
    try:
      daz  = float(arglist[0])
      dalt = float(arglist[1])
    except Exception as e:
      mainApp.messages.error_message.set('ShiftHor '+arglist+ ' error')
      printe('Exposition time '+arglist + ' error ',e)
      self.stop(-1)
      return -1
    res = icsmodule.shiftHor(daz, dalt)
    printd ('shift_hor res=', res)
    if res < 0:
      return self.stop(res)
    set_time = 60000 # in ms
    res = self.wait_step(set_time, icsmodule.OCS_HOR_RUN_tkvar)
    print('shift_hor  after wait res =', res, type(res))
    if res != 0:
      str_error = 'Shift Hor: error code ' + str(res)
      printe(str_error) 
      mainApp.messages.error_message.set(str_error)
      self.stop(res)
    return res


# ---------------------- GUI ----------------------

class IQwindowFrame(Tk.Frame):
# lll
  def __init__(self, parent, *args, **kwargs):
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    self.parent     = parent
    self.fwhma      = Tk.DoubleVar(value=0.0)
    self.fwhmb      = Tk.DoubleVar(value=0.0)
    self.background = Tk.DoubleVar(value=0.0)
    self.nstars     = Tk.DoubleVar(value=0.0)

    self.tkvar_dict = {'fwhma':self.fwhma, 'fwhmb': self.fwhmb, \
                    'background':self.background, 'nstars':self.nstars}
    column = 0
    for key in self.tkvar_dict:
      txt = key
      tkvar = self.tkvar_dict[key]
      label_txt   = Tk.Label(self, text=txt,  width=13) #, anchor=Tk.E)
      label_tkvar = Tk.Label(self, textvariable=tkvar, fg='blue') #, anchor=Tk.E)
      label_txt.grid  (row=0, column=column, sticky=Tk.EW)
      label_tkvar.grid(row=1, column=column, sticky=Tk.EW)
      column += 1

  def set_tkvars(self, iq_dict):
    print('IQwindowFrame:set_tkvars', iq_dict)
    for key in self.tkvar_dict:
        print(key, iq_dict.get(key, -1.0))
        self.tkvar_dict[key].set(iq_dict.get(key, 0.0))

class ExposeFrame(Tk.Frame, PVTkbase):
  def __init__(self, parent, *args, **kwargs):
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    PVTkbase.__init__(self, parent)
    self.parent = parent

    self.sci_scenario_table = sciscenario_window.SciScenarioTable()    

    self.wlist4disable = []

    self.measurements = BDmeasurements.Measurements()

    self.imagetype    = Tk.StringVar()
    self.exptime      = Tk.StringVar() # exposure time in sec
    self.numexp       = Tk.StringVar()
    self.userfilename = Tk.StringVar()  # user's specific file name
    self.filename     = nbimodule.FILENAME_tkvar  # real fits filename
    self.user_scenario_string = ''
    self.user_scenario_name   = Tk.StringVar()
    self.auto_start           = Tk.BooleanVar(value=False)
    self.sai25_status_old     = ''


    self.exposure     = Tk.StringVar()  # 1 - if exposure continues, 0 - in the other case
    self.expprogress  = nbimodule.EXPPROGRESS_tkvar
    self.drvprogress  = nbimodule.DRVPROGRESS_tkvar

    self.tmpa         = nbimodule.TEMPERA_tkvar
    self.tmpl         = nbimodule.TEMPERL_tkvar
    self.press        = nbimodule.PRESSURE_tkvar
    self.response     = nbimodule.RESPONSE_tkvar
    
    self.cb_auto_index = None

    self.timer_id = 0

    self.imagetype.set('light')
    self.exptime.set('1.0')
    self.numexp.set(1)
    self.userfilename.set('')   
    self.user_scenario_name.set('')
    
    radio = Tk.Frame(self)
    
    pdy = 5; pdx = 5 
    for imt in ['bias','dark','light','flat']:
      rb = Tk.Radiobutton(radio, text=imt, variable=self.imagetype, value=imt).pack(anchor=Tk.W, side=Tk.LEFT, fill=Tk.BOTH, pady=pdy, padx=12)

    AGU_lim_label = Tk.Label(radio, text='AGU tr.limit:')
    AGU_lim_label.pack(side=Tk.LEFT, anchor=Tk.E, fill=Tk.BOTH, pady=pdy, padx=pdx)
    AGU_lim_entry = Tk.Entry(radio, width=7, textvariable=agutran_min)
    AGU_lim_entry.pack(side=Tk.LEFT, anchor=Tk.W, fill=Tk.BOTH, pady=pdy, padx=pdx)
    exp_frame = Tk.Frame(self)

    anch = Tk.W
    exp_label          = Tk.Label(exp_frame, text='Exp.time (s):',  anchor=anch)
    numexp_label       = Tk.Label(exp_frame, text='num mexp:',      anchor=anch)
    userfilename_label = Tk.Label(exp_frame, text='Outfile:',       anchor=anch)

    e_width = 10
    pdy = 2
    pdx = 2
    exp_entry     = Tk.Entry(exp_frame, width=e_width, textvariable=self.exptime)
    numexp_entry  = Tk.Entry(exp_frame, width=e_width, textvariable=self.numexp)
    userfilename_entry = Tk.Entry(exp_frame, width=e_width, textvariable=self.userfilename)
 
    bv = 9
    expose_button = Tk.Button(exp_frame, text='Expose',    width=bv,  pady=pdy, highlightbackground='lightgreen',  command = self.run_expose)
    self.wlist4disable.append(expose_button)
    multi_button  = Tk.Button(exp_frame, text='MultiExp',  width=bv,  pady=pdy, highlightbackground='lightgreen',  command = self.run_multi)
    self.wlist4disable.append(multi_button)
    probe_button       = Tk.Button(exp_frame, text='Probe',      width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.run_probe)
    self.wlist4disable.append(probe_button)
    stop_multi_button   = Tk.Button(exp_frame, text='Stop multi',width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.stop_multi)
    scenario_button  = Tk.Button(exp_frame, text='Run scen',     width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.run_expose_scenario)
    self.wlist4disable.append(scenario_button)
    stop_scenario_button   = Tk.Button(exp_frame,      text='Stop scen',    width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.stop_scenario)
    auto_start_check       = Tk.Checkbutton(exp_frame, text='Auto start',   fg='darkgreen', variable=self.auto_start, command=self.refresh_auto_start)
    
    load_scenario_button   = Tk.Button(exp_frame, text='Load scen',   width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.load_scenario)
    load_bdscenario_button = Tk.Button(exp_frame, text='Load DBscen', width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.load_BDscenario)
    open_sci_scenario_window_button  = Tk.Button(exp_frame, text='Open' , width=bv,  pady=pdy, highlightbackground='lightgreen', command = self.open_sci_scenario_window)

    stop_button            = Tk.Button(exp_frame, text='Stop',        width=bv,  pady=pdy, highlightbackground='red',        command = self.stop_expose)
    self.suspend_button    = Tk.Button(exp_frame, text='Suspend',     width=bv,  pady=pdy, highlightbackground='red',        command = self.suspend_expose)
#    self.wlist4disable.append(self.suspend_button)

    scenario_name_label =  Tk.Label (exp_frame, text='Scenario:',  anchor=anch)
    scenario_label      =  Tk.Label (exp_frame, textvariable=self.user_scenario_name,  anchor=anch)

    expprogress_name_label =  Tk.Label (exp_frame, text='lasts',  anchor=anch)
    drvprogress_name_label =  Tk.Label (exp_frame, text='read',   anchor=anch)
    expprogress_label =  Tk.Label (exp_frame, textvariable=self.expprogress,  anchor=anch)
    drvprogress_label =  Tk.Label (exp_frame, textvariable=self.drvprogress,  anchor=anch)

    tmpa_label  =  Tk.Label (exp_frame, textvariable=self.tmpa,   width=7, anchor=anch)
    tmpl_label  =  Tk.Label (exp_frame, textvariable=self.tmpl,   width=7, anchor=anch)
    press_label =  Tk.Label (exp_frame, textvariable=self.press,  width=7, anchor=anch)
    response_label = Tk.Label (exp_frame, textvariable=self.response,  anchor=anch)
    
    filename_label = Tk.Label(exp_frame, textvariable=self.filename, anchor=anch)

    row = 0
    for e in [exp_entry, numexp_entry, userfilename_entry]:
      e.grid(row=row, column=1, sticky=Tk.W)
      row = row+1

    row = 0
    for l in [exp_label, numexp_label, userfilename_label]:
      l.grid(row=row, column=0, sticky=Tk.E)
      row = row+1

    row = 0
    for b in [expose_button, multi_button, scenario_button]:
      b.grid(row=row, column=2, sticky=Tk.W, pady=3, padx=5)
      row = row+1

    row = 0
    for b in [probe_button, stop_multi_button, stop_scenario_button, load_bdscenario_button]:
#      b.grid(row=row, column=3, columnspan=2, sticky=Tk.W)
      b.grid(row=row, column=3, sticky=Tk.W)
      row = row+1
      
    stop_button.grid                    (row=0, column=4, padx=5, sticky=Tk.W)
    self.suspend_button.grid            (row=1, column=4, padx=5, sticky=Tk.W)
    auto_start_check.grid               (row=2, column=4, padx=5, sticky=Tk.W)
    open_sci_scenario_window_button.grid(row=3, column=4, padx=5, sticky=Tk.W)

    column = 0; row = 3
    for l in [scenario_name_label, scenario_label, load_scenario_button]:
       l.grid(row=row, column=column, sticky=Tk.W, pady=3, padx=5)
       column = column + 1

    row = row + 1; column = 0
    for l in [expprogress_name_label, expprogress_label, drvprogress_name_label]:
       l.grid(row=row, column=column, sticky=Tk.W, pady=3, padx=5)
       column = column + 1
    drvprogress_label.grid(row=row, column=3, columnspan=2, sticky=Tk.W, pady=3, padx=5)

    row = row + 1
    filename_label.grid (row=row, column=0, columnspan=2, sticky=Tk.W)
    tmpa_label.grid     (row=row, column=2, sticky=Tk.W)
    tmpl_label.grid     (row=row, column=3, sticky=Tk.W)
    press_label.grid    (row=row, column=4, sticky=Tk.W)

    row = row + 1
    response_label.grid (row=row, column=0, columnspan=5, sticky=Tk.W)

    radio.pack    (anchor=Tk.W, side=Tk.TOP, pady=pdy)
    exp_frame.pack(anchor=Tk.W, side=Tk.TOP, pady=pdy)

    self.iq_window = IQwindowFrame(self)
    self.iq_window.pack(anchor=Tk.W, side=Tk.TOP, pady=pdy)

    nbimodule.EXPOSURE_tkvar.trace('w', self.tkcallback_exposure)

#    self.initPVcallbacks()
# SAI25_STATUS

  def initExchange(self):
    pass
#    SAI25_STATUS.add_callback(self.pvcallback_auto_start, with_ctrlvars=False)
#    pass
#    self.initPVcallbacks()

  def refresh_auto_start(self):
    if self.auto_start.get():
       self.cb_auto_index = SAI25_STATUS.add_callback(self.pvcallback_auto_start, with_ctrlvars=False)
    else:
       SAI25_STATUS.remove_callback(index=self.cb_auto_index)
    

  def prepare_exposure(self):
#    self.parent.ocs.writePV_observer()
    exp_str = self.exptime.get()
    filter_name = self.parent.filters.filter_target.get()
    imagetype = self.imagetype.get()
    userfilename = self.userfilename.get()

    printd(filter_name, exp_str, imagetype, userfilename)
    return exp_str, filter_name, imagetype, userfilename
  

  def run_multi(self):
    printdd( 'Run multi', self.numexp.get(), 'exp of', self.exptime.get(), 'sec')
    try:
      numexp = int(self.numexp.get())
    except Exception as e:
      mainApp.messages.error_message.set('Wrong exposition number '+self.numexp.get())
      printe( 'Wrong exposition number '+self.numexp.get())
      return -1
    exp_str, filter_name, imagetype, userfilename = self.prepare_exposure()
    scenario_string = 'Filter,'+filter_name+'+'+str(numexp)+'*[Exp,'+exp_str+','+imagetype+','+userfilename+']'
    printd(scenario_string)
    if scenario.load(scenario_string) == -1: return -1
    return scenario.run()

  def run_expose(self):
    printdd('Scenario:run_expose')
    return self.run_expose_base(1, None)

  def run_probe(self):
    printd( 'Run probe', self.exptime.get(), 'sec')
    return self.run_expose_base(1,'focusing')

  def pvcallback_auto_start(self, pvname=None, value=-1, char_value=-1, **kw):
     if value == self.sai25_status_old: return
     self.sai25_status_old = value

     if value == 'TRACK' and self.auto_start.get() and not scenario.scenario_is_running:
       if SAI25_ICS.get() == myICS:
         printd( '\n\n\n -------------- Scenario: pvcallback_auto_start', value, '-----\n\n\n')
         self.load_BDscenario()
         self.run_expose_scenario()

  def run_expose_scenario(self):
    printdd('Scenario:run_expose_scenario')
    if self.user_scenario_string == '': return -1
    if scenario.load(self.user_scenario_string) == -1: return -1
    scenario.long_scenario_is_running = True
    return scenario.run(auto_target_done=self.auto_start.get(), auto_start=self.auto_start.get())

  def run_expose_base(self, numexp, userfilename_): # if userfilename_ is None compose SAIstyle filename
    exp_str, filter_name, imagetype, userfilename = self.prepare_exposure()
    if userfilename_ is not None: userfilename = userfilename_
#    scenario_string = str(numexp)+'*[Exp,'+exp_str+','+imagetype+','+userfilename+']'
    scenario_string = str(numexp)+'*[Filter,'+filter_name+'+Exp,'+exp_str+','+imagetype+','+userfilename+']'
    printd(scenario_string)
    if scenario.load(scenario_string) == -1: return -1
    return scenario.run()

  def stop_multi(self):
    printd( 'Stop multi')
    scenario.stop(0)

  def stop_scenario(self):
    printd( 'Stop scenario')
    scenario.stop(0)

  def stop_expose(self):
    printd( 'Stop expose')
    nbi.stop()

  def suspend_expose(self):
    printd( 'Suspend expose')
    if nbimodule.EXPOSURE_tkvar.get() != 1: return
#    if nbi.nbi_exposure != 1: return
    if nbi.suspended == 0:
      scenario.suspend(1)
      self.suspend_button.configure(text='Continue')
#      self.datetime_suspended = datetime.now()
#      mainApp.after_cancel(scenario.timer_id)
      # stop timer here!
    else:
      scenario.suspend(0)
      self.suspend_button.configure(text='Suspend')
# recalucale self.expose_timeout
#      self.timer_id = mainApp.after(self.exposure_and_reading_time_ms, self.expose_timeout) # if something go wron
      # release timer here!

  def load_BDscenario(self):
    print('load_BDscenario')
    self.user_scenario_string = self.measurements.get_scenario_string()
    print('user_scenario_string =', self.user_scenario_string)
    self.sci_scenario_table.load(unfold(self.user_scenario_string))

  def open_sci_scenario_window(self):
    name = 'scenario'
# Any MY Tk.Toplevel of CommonButtons MUST have name !
    for widget in self.winfo_children():
      if isinstance(widget, Tk.Toplevel) and widget.name == name:
        widget.destroy()
    self.scenario_window = Tk.Toplevel(self)
    self.scenario_window.name = name
    self.scenario_window.wm_title('Scenario')
    self.sci_scenario_table.gwindow(self.scenario_window)

  def load_scenario(self):
    scenario_fname =  tkFileDialog.askopenfilename(initialdir = ".",\
      title = "Select file",filetypes = (("scenario files","*.scn"),("all files","*.*")))

#    scenario_fname = 'RW_Aur.scn'


    self.scenario_name = scenario_fname.replace('.scn','')
    userfilename = ''

    with open(scenario_fname,'r') as f:
      lines = f.readlines()

    ncycles = '1'
    for line in lines:
      tu = line.split()
      if len(tu) < 2: continue
      if tu[0].upper() == 'NCYCLES':
        ncycles = tu[1]
        break;

    self.user_scenario_string = str(ncycles) + '*['

    for line in lines:
      if line[0] == '#': continue
      tu = line.split()
      if len(tu) < 2: continue
      if tu[0] == 'NCYCLES' or tu[0] == 'OBSERVER' or tu[0] == 'PROGID' or tu[0] == 'OBJECT' or tu[0] == 'FWHM': continue
      filter,exp = tu[0],tu[1]
      if not fsu.filter_name_in_fsu(filter):
#      if filter not in fsu.FILTER_NAMES:
        printe('Load scenario: Wrong filter name '+filter+' skipped')
        continue
      self.user_scenario_string = self.user_scenario_string + 'Filter,'+filter + '+Exp,' + str(exp) + ',Light,' + userfilename + '+'  

    self.user_scenario_string = self.user_scenario_string[:-1] + ']'
    print(self.user_scenario_string)
    self.sci_scenario_table.load(unfold(self.user_scenario_string))
    return 0
  # Check filters Ask FSU 
  # save filename as scenario_name


  def disable_due_to_exposure(self):
#    printd( self.__class__.__name__, 'disable_due_to_exposure')
    for w in self.wlist4disable:
      w.configure(state='disabled')

  def enable_after_exposure(self):
#    printd( self.__class__.__name__, 'enable_after_exposure')
# restore Continue to Suspend if it hung
    self.suspend_button.configure(text='Suspend')
    for w in self.wlist4disable:
      w.configure(state='normal')

  def end_of_exposure_callback(self): #  empty for charging
    pass

#  def initPVcallbacks(self):
#    pass

  def tkcallback_exposure(self, *args):
    value = nbimodule.EXPOSURE_tkvar.get()
    printd( 'ExposeFrame tkcallback_exposure', value)
    if value == 1:
      mainApp.disable_due_to_exposure()
    elif value == 0:
      mainApp.enable_after_exposure()

class OCS(Tk.Frame):

  def __init__(self, parent, *args, **kwargs):
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    self.parent = parent

    self.progid      = Tk.StringVar()
#    self.observer = Tk.StringVar()
    self.observer    = nbimodule.OBSERVER_tkvar
    self.object      = Tk.StringVar()
    self.ra          = Tk.StringVar()
    self.dec         = Tk.StringVar()
#    self.ics_current = Tk.StringVar()
    self.newics_target  = Tk.StringVar(value=SAI25_ICS.get())
    
    self.iq      = Tk.IntVar()
    self.ds9view = Tk.IntVar()
    self.change_progid = Tk.IntVar()
    self.bd = Tk.IntVar()
    for v in [self.iq, self.ds9view, self.bd]: v.set(1)
    self.change_progid.set(0)


    e_width = 12
    self.observer.set(SAI25_OPERATOR.get())
#    self.ra.set('06:06:15.00')
    observer_entry  = Tk.Entry (self, width=e_width, textvariable=self.observer)  # NBI:OBSERVER
#    observer_entry.bind('<Return>', self.writePV_observer)

    progid_entry    = Tk.Entry(self, width=e_width, textvariable=self.progid)    # SAI25:PROGID
    object_entry    = Tk.Entry(self, width=e_width, textvariable=self.object)    # SAI25:OBJECT
    ra_entry        = Tk.Entry(self, width=e_width, textvariable=self.ra)        # SAI25:RA
    dec_entry       = Tk.Entry(self, width=e_width, textvariable=self.dec)       # SAI25:DEC

    anch = Tk.W
    observer_label = Tk.Label(self, text='Observer',  anchor=anch)
    progid_label   = Tk.Label(self, text='ProgID',    anchor=anch)
    object_label   = Tk.Label(self, text='Object',    anchor=anch)
    ra_label       = Tk.Label(self, text='RA:',       anchor=anch)
    dec_label      = Tk.Label(self, text='Dec:',      anchor=anch)

    ocs_connection_frame = Tk.Frame(self)
#    ocs_conn_label       = Tk.Label (ocs_connection_frame, text='OCS:', width=5,     anchor=Tk.E)
    ocs_connection       = IndicatorStatusTkvar(ocs_connection_frame, text='OCS', width=7,  tkvar=icsmodule.ocs_client_lbl)
    target_done_button   = Tk.Button(ocs_connection_frame,  width=8, text='Next', command=icsmodule.targetDone)
    for w in [ocs_connection, target_done_button]:
      w.pack(side=Tk.LEFT, padx=3)
#    ocs_conn_label.pack    (side=Tk.LEFT, padx=3)
#    ocs_connection.pack    (side=Tk.LEFT, padx=3)
#    target_done_button.pack(side=Tk.LEFT, padx=3)

    agu_frame  = Tk.Frame(self)
    self.agu_status = IndicatorStatus(agu_frame, AGSTATUS, text='AGU', width=7, anchor=Tk.CENTER)
    self.agu_on_button  = Tk.Button(agu_frame, width=2, text='ON',  command=icsmodule.aguOn)
    self.agu_off_button = Tk.Button(agu_frame, width=2, text='OFF', command=icsmodule.aguOff)
    for w in [self.agu_status, self.agu_on_button, self.agu_off_button]:
      w.pack(side=Tk.LEFT, padx=3)
#
    iq_check      = Tk.Checkbutton(self, text = 'IQ',      fg='darkgreen', variable=self.iq)
    ds9_check     = Tk.Checkbutton(self, text = 'DS9',     fg='darkgreen', variable=self.ds9view)
    bd_check      = Tk.Checkbutton(self, text = 'BD',      fg='darkgreen', variable=self.bd)
    progid_check  = Tk.Checkbutton(self, text = 'Fix progid', fg='darkgreen', variable=self.change_progid, command=self.refresh_progid)

#    ocs_host_label = Tk.Label (self, text='OCS:',      anchor=anch) 
#    self.ocs_host_target = Tk.StringVar(value=icsmodule.ocsd_default_host)
#    ocs_host_menu = Tk.OptionMenu(self,  self.ocs_host_target, *icsmodule.ocs_host_choices)
#    self.ocs_host_target.trace('w', self.ocs_host_tkcallback)

    ics_frame = Tk.Frame(self)
    ics_menu           = Tk.OptionMenu(ics_frame, self.newics_target, *icsmodule.ics_choices)
    ics_menu.config(width=2)
    newics_button      = Tk.Button(ics_frame, text='New ICS', width=8, command=self.newics)
    for w in [ics_menu, newics_button]:
       w.pack(side=Tk.LEFT, padx=3)
       

#    icon_scale = 11
#    self.drink = Tk.IntVar(value=1)
#    drunk_path = 'drunk.png'
#    sober_path = 'sober.png'
#    exist = os.path.exists(drunk_path)
#    if not os.path.exists(drunk_path):
#        printe( 'File', drunk_path, 'does not exist')
#        return
#    if not os.path.exists(sober_path):
#        printe( 'File', sober_path, 'does not exist')
#        return
#
#    self.drunk_image = Tk.PhotoImage(file = drunk_path).subsample(icon_scale)
#    self.sober_image = Tk.PhotoImage(file = sober_path).subsample(icon_scale)
#    drink_check = Tk.Checkbutton(self, image=self.sober_image, selectimage=self.drunk_image, indicatoron=False,
#                     onvalue=1, offvalue=0, variable=self.drink)

    row = 0
    for label in [observer_label, progid_label, object_label, ra_label, dec_label]:
      label.grid(row=row, column=0, sticky=Tk.E)    
      row = row + 1

    row = 0
    for w in [observer_entry, progid_entry, object_entry, ra_entry, dec_entry]:
      w.grid(row=row, column=1, columnspan=2, sticky=Tk.W, pady=3, padx=5)    
      row = row + 1

    ocs_connection_frame.grid(row=row, column=0, columnspan=3, sticky=Tk.EW)
    row = row + 1
    ics_frame.grid           (row=row, column=0, columnspan=3, sticky=Tk.EW)
    row = row + 1
    agu_frame.grid           (row=row, column=0, columnspan=3, sticky=Tk.EW)

    row = row+1
    iq_check.grid(row=row, column=0, sticky=Tk.W)
    ds9_check.grid(row=row,column=1, sticky=Tk.W)
#    drink_check.grid(row=row, rowspan=2, column=2, sticky=Tk.W)
    row = row+1
    bd_check.grid(row=row, column=0, sticky=Tk.W)
    progid_check.grid(row=row, column=1, sticky=Tk.W)
#    row += 1
#    ocs_host_label.grid(row=row, column=0, sticky=Tk.E)
#    ocs_host_menu.grid( row=row, column=1, columnspan=2, sticky=Tk.W)

    self.pvra = PV('SAI25:TARRA')
    self.pvdec = PV('SAI25:TARDEC')
    self.pvobject = PV('SAI25:OBJECT')
    self.pvprogid = PV('SAI25:PROGID')
#    self.initExchange()

  def newics(self):
    icsmodule.newICS(self.newics_target.get())

  def ocs_host_tkcallback(self, *args):
    print('new ocs host IP =', self.ocs_host_target.get())
    icsmodule.OCS_IP = self.ocs_host_target.get()

#  def writePV_observer(self, *kwargs):
#    global NBI_OBSERVER
#    printd( 'writePV_observer, self.observer =', self.observer.get())
#    NBI_OBSERVER.put(self.observer.get())
#    printd( 'NBI_OBSERVER =', NBI_OBSERVER)
#    if not debug_nbi:
#      NBI_OBSERVER.put(self.observer.get())

  def disable_due_to_exposure(self):
    pass
#    printd( self.__class__.__name__, 'disable_due_to_exposure')

  def enable_after_exposure(self):
    pass
#    printd( self.__class__.__name__, 'enable_after_exposure')

  def pvcallback_ra(self, pvname=None, value=-1, char_value=-1, **kw):
    ra_string = deg2hms(value)
    self.ra.set(ra_string)

  def pvcallback_dec(self, pvname=None, value=-1, char_value=-1, **kw):
    dec_string = deg2dms(value)
    self.dec.set(dec_string)

  def pvcallback_object(self, pvname=None, value=-1, char_value=-1, **kw):
    self.object.set(value)

  def pvcallback_progid(self, pvname=None, value=-1, char_value=-1, **kw):
    if not self.change_progid.get():
      self.progid.set(value)

  def pvcallback_ics(self, pvname=None, value=-1, char_value=-1, **kw):
#     self.ics_current.set(value)
     self.newics_target.set(value=value)

  def refresh_progid(self):
    if not self.change_progid.get():
      self.progid.set(self.pvprogid.get())

  def readPV(self):
    self.progid.set(self.pvprogid.get())
    self.object.set(self.pvobject.get())
    ra_string = deg2hms(self.pvra.get())
    self.ra.set(ra_string)
    dec_string = deg2dms(self.pvdec.get())
    self.dec.set(dec_string)

  def initExchange(self):
    self.readPV()
    try:
#      self.agu_status.initExchange()
      self.pvra.add_callback    (self.pvcallback_ra,     with_ctrlvars=False)
      self.pvdec.add_callback   (self.pvcallback_dec,    with_ctrlvars=False)
      self.pvobject.add_callback(self.pvcallback_object, with_ctrlvars=False)
      self.pvprogid.add_callback(self.pvcallback_progid, with_ctrlvars=False)
      SAI25_ICS.add_callback    (self.pvcallback_ics,    with_ctrlvars=False)

    except Exception as e:
      printe('OCS add_callback exception:', e)

    
class Filters(Tk.Frame, PVTkbase):

  def __init__(self, parent, *args, **kwargs):

#    self.method_to_run_after_filter_setting = empty_method

    self.wheelA = fsu.get_filter_names_in_wheel('A')
    self.wheelB = fsu.get_filter_names_in_wheel('B')
    
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    PVTkbase.__init__(self, parent)

    radio = Tk.Frame(self)
    self.filter_current = Tk.StringVar()
    self.filter_current.set(fsu.get_filter_name())

    self.mirror_current = Tk.StringVar()
    self.mirror_current.set(fsu.get_mirror())

    self.lamp_current = Tk.StringVar()
    self.lamp_current.set(fsu.get_lamp())

    self.filter_target = Tk.StringVar()
    self.filter_target.set(fsu.get_filter_name())

    self.mirror_target = Tk.StringVar()
    self.mirror_target.set(fsu.get_mirror())

    self.lamp_target = Tk.StringVar()
    self.lamp_target.set(fsu.get_lamp())

    self.shutter_target = Tk.StringVar()
    self.shutter_target.set(fsu.get_shutter())

    self.wlist4disable = []

    width = 9

    self.status          = IndicatorStatus(radio, FSU_STATUS, text='FSU', width=width, anchor=Tk.W)
    filter_current_label = Tk.Label(radio, width=width, textvariable=self.filter_current, bg='yellow', anchor=Tk.W, relief = Tk.RIDGE)
#    self.cmd_run_label   = JKLabel (radio, width=width, pv=FSU_CMD_RUN, anchor=Tk.W)
    shutter_label        = JKLabel (radio, width=width, pv=FSU_SHUTTER, anchor=Tk.W)
    mirror_current_label = Tk.Label(radio, width=width, textvariable=self.mirror_current, bg='yellow', anchor=Tk.W, relief = Tk.RIDGE)
    lamp_current_label   = Tk.Label(radio, width=width, textvariable=self.lamp_current, bg='yellow', anchor=Tk.W, relief = Tk.RIDGE)

    row = 0
    self.status.grid         (row=row, column=0, sticky=Tk.W)
    filter_current_label.grid(row=row, column=1, sticky=Tk.W)
#    self.cmd_run_label.grid  (row=row, column=2, sticky=Tk.W)
    shutter_label     .grid  (row=row, column=2, sticky=Tk.W)
    mirror_current_label.grid(row=row, column=3, sticky=Tk.W)
    lamp_current_label.grid  (row=row, column=4, sticky=Tk.W)

    row = row + 1
    for wheel in [self.wheelA, self.wheelB]:
      printd( wheel[0])
      col = 0
      for filter in wheel:         
         rb = Tk.Radiobutton(radio, text=filter, variable=self.filter_target, value=filter) #.pack(anchor=Tk.W, side=Tk.LEFT, pady=pdy)
         rb.grid(row=row, column=col, sticky=Tk.W)
         self.wlist4disable.append(rb)
         col = col + 1
      row = row + 1

    manual = Tk.Frame(self)

    entry_target  = Tk.Entry (manual, width=17, textvariable=self.filter_target)
#    entry_target  = Tk.Entry (manual, textvariable=self.filter_target)
    self.wlist4disable.append(entry_target)
    move   = Tk.Button(manual, text='Move', width=6,   command = self.move_fsu ) #.pack(anchor=Tk.W, side=Tk.LEFT, pady=pdy)
    self.wlist4disable.append(move)
    
    mirror_target_menu = Tk.OptionMenu(manual, self.mirror_target, *fsu.mirror_target_choices)
    self.wlist4disable.append(mirror_target_menu)

    lamp_target_menu = Tk.OptionMenu(manual, self.lamp_target, *fsu.lamp_target_choices)
    self.wlist4disable.append(lamp_target_menu)

    shutter_target_menu = Tk.OptionMenu(manual, self.shutter_target, *fsu.shutter_target_choices)
    self.wlist4disable.append(shutter_target_menu)

    row=0
    entry_target.grid (row=row, column=0, columnspan=2, sticky=Tk.W)
    row = row+1
    mirror_target_menu.grid (row=row,   column=0, sticky=Tk.W)
    lamp_target_menu.grid   (row=row+1, column=0, sticky=Tk.W)
    move.grid               (row=row,   column=1, sticky=Tk.W)
    shutter_target_menu.grid(row=row+1, column=1, sticky=Tk.W)

    radio.pack(side = Tk.LEFT, padx = 3, anchor=Tk.W)
    manual.pack(side = Tk.LEFT, padx = 3, anchor=Tk.W)
    
#    self.init_epics()

#  def set_after_method(self, method):
#    self.method_to_run_after_filter_setting = method

  def pvcallback_cmd_run(self, pvname=None, value=-1, char_value=-1, **kw):
    if value == 0 :
#      self.method_to_run_after_filter_setting()
      self.enable_fsu()
    else:
      self.disable_fsu()

  def disable_fsu(self):
    printd( self.__class__.__name__, 'disable')
    for w in self.wlist4disable:
       w.configure(state='disabled')

  def enable_fsu(self):
    printd( self.__class__.__name__, 'enable')
    for w in self.wlist4disable:
       w.configure(state='normal')

  def move_fsu(self):
    printd('Filters:move_fsu')
    fsu.set_filter_target(self.filter_target.get())
    fsu.set_mirror_target(self.mirror_target.get())
    fsu.set_lamp_target(self.lamp_target.get())
    fsu.set_shutter_target(self.shutter_target.get())
    fsu.move_fsu()
    printd('End of Filters:move_fsu')
#  self.filter may be 'U' or F001UBES or 'A0' or 'A0,B6' or 'CA1'(for change filter in the A1 hole)

  def disable_due_to_exposure(self):
    pass
#    printd( self.__class__.__name__, 'disable_due_to_exposure')

  def enable_after_exposure(self):
    pass
#    printd( self.__class__.__name__, 'enable_after_exposure')

  def readPV(self):
    self.filter_current.set(fsu.get_filter_name())
    self.mirror_current.set(fsu.get_mirror())
    self.lamp_current.set(fsu.get_lamp())
  
  def pvcallback_filter(self, pvname=None, value=-1, char_value=-1, **kw):
#    printd('Filters:pvcallback_filter',value)
    self.filter_current.set(fsu.filterid2name(value))

  def pvcallback_mirror(self, pvname=None, value=-1, char_value=-1, **kw):
    self.mirror_current.set(value)

  def pvcallback_lamp(self, pvname=None, value=-1, char_value=-1, **kw):
    self.lamp_current.set(value)

  def initExchange(self):
    printd( self.__class__.__name__,'initExchange')
    self.readPV()
    try:
         FSU_FILTER.add_callback (self.pvcallback_filter,  with_ctrlvars=False)
         FSU_MIRROR.add_callback (self.pvcallback_mirror,  with_ctrlvars=False)
         FSU_LAMP.add_callback   (self.pvcallback_lamp,    with_ctrlvars=False)
         FSU_CMD_RUN.add_callback(self.pvcallback_cmd_run, with_ctrlvars=False)
    except Exception as e:
         printe( 'FSU_FILTER add_callback exception:', e)
    printd( self.__class__.__name__, 'end of initExchange')

class Messages(Tk.Frame):
  def __init__(self, parent, *args, **kwargs):
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    self.parent = parent
    self.error_message   = nbimodule.ERROR_MESSAGE_tkvar
    self.current_message = Tk.StringVar()

    error_lbl   = Tk.Label(self, textvariable=self.error_message, bg='orange',  anchor=Tk.W)
    current_lbl = Tk.Label(self, textvariable=self.current_message, anchor=Tk.W)
    current_lbl.grid(row=0, sticky=Tk.EW)
    error_lbl.grid  (row=1, sticky=Tk.EW)

  def disable_due_to_exposure(self):
    pass

  def enable_after_exposure(self):
    pass

  def initExchange(self):
    pass

class Focuser():
  def __init__(self):
    self.altstep_sec = Tk.DoubleVar();
    self.altstep_sec.set(10.0) # arcsec
    self.nsteps = Tk.IntVar()
    self.nsteps.set(7)
    self.foc_step = Tk.DoubleVar()
    self.foc_step.set(0.05)
    self.foc_start = Tk.DoubleVar()
    self.set_foc_start()
    self.focus_best_tkvar = Tk.DoubleVar(value=M2_POS.get())
    self.focus_best = None
    self.fwhm_tkvar       = Tk.DoubleVar(value=-1.0)
    self.fwhm             = -1.0
    self.exptime = Tk.DoubleVar()
    self.exptime.set(5.0) # sec
    self.m2_pos = Tk.DoubleVar()
    self.m2_pos.set(M2_POS.get())
    self.ds9 = None
#    self.foc_step.trace('w', self.set_foc_start)
#    self.nsteps.trace('w', self.set_foc_start)
    self.initFocusPSF()
    self.dero_start = None
    self.parang_start = None

    self.fullpath2focusfile = None

  def initFocusPSF(self):
    self.foci = []
    self.fwhm_array = []
    self.poly = None


  def set_foc_start(self, *args):
    self.foc_start.set(round2float(self.calc_foc_start(),3))

  def calc_foc_start(self):
    return M2_POS.get() - self.nsteps.get()/2.0*self.foc_step.get()

  def initExchange(self):
    M2_POS.add_callback(self.pvcallback_m2, with_ctrlvars=False)

  def focus_window(self, parent):
    frame0 = Tk.Frame(parent, padx=10, pady=10, borderwidth=3, relief=Tk.RIDGE)
    frame0.pack(fill=Tk.X)
    frame1 = Tk.Frame(parent, padx=10, pady=10, borderwidth=3, relief=Tk.RIDGE)
    frame1.pack(fill=Tk.X)

    width = 7

    row = 0
    column = 0
    Tk.Label(frame0, text='M2 start', anchor=Tk.W ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Entry(frame0, width=width, textvariable=self.foc_start).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Label(frame0, text='M2 step ', anchor=Tk.W ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    step_entry = Tk.Entry(frame0, width=width, textvariable = self.foc_step)
    step_entry.bind('<Return>', self.set_foc_start)
    step_entry.grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Label(frame0, text='N steps ', anchor=Tk.W ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    nsteps_entry = Tk.Entry(frame0, width=width, textvariable = self.nsteps   )
    nsteps_entry.bind('<Return>', self.set_foc_start)
    nsteps_entry.grid(row=row, column=column, sticky=Tk.W)
    row += 1; column = 0
    Tk.Label(frame0, text='M2 now', anchor=Tk.W                        ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Label(frame0, width=width, textvariable=self.m2_pos, fg='blue', anchor=Tk.W).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Label(frame0, text='Expose', anchor=Tk.W                        ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Entry(frame0, width=width, textvariable=self.exptime            ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Label(frame0, text='Alt step', anchor=Tk.W                      ).grid(row=row, column=column, sticky=Tk.W)
    column = column + 1
    Tk.Entry(frame0, width=width, textvariable = self.altstep_sec      ).grid(row=row, column=column, sticky=Tk.W)


    bw = 7
    row += 1; column = 0
    Tk.Button(frame0, text='Go', fg="blue", command=self.run_focusing).grid(row=row, column=column, sticky=Tk.EW)
    column += 1
#    Tk.Button(frame0, text='PSF', fg="blue", width=bw, command=self.calc_focus_psf).grid(row=row, column=column, sticky=Tk.W) 
    psf_button = ButtonQueue(frame0, queue=nbi_command_queue, text='PSF', queue_command=self.calc_focus_psf, on_end_command=self.psf_on_end) #width=7
    psf_button.grid(row=row, column=column, sticky=Tk.EW)
    column += 1
#    Tk.Button(frame0, text='Donuts', fg="blue", width=bw, command=self.calc_focus_donuts).grid(row=row, column=column, sticky=Tk.W) 
    donuts_button = ButtonQueue(frame0, queue=nbi_command_queue, text='Donuts', queue_command=self.calc_focus_donuts, on_end_command=self.donuts_on_end)#, width=10)
    donuts_button.grid(row=row, column=column, sticky=Tk.EW)
    column += 1
    setM2_button  = Tk.Button(frame0, text='set Best', fg="blue", command=self.set_best)
    setM2_button.grid(row=row, column=column, sticky=Tk.EW)

    column += 1
    upload2bdM2_button  = Tk.Button(frame0, text='save M2', fg="blue", command=self.upload2BD_M2)
    upload2bdM2_button.grid(row=row, column=column, sticky=Tk.EW)

    row += 1; column = 0
    Tk.Label(frame0, text='Best focus', anchor=Tk.W                      ).grid(row=row, column=column, sticky=Tk.W)
    column += 1
    Tk.Entry(frame0, width=width, textvariable=self.focus_best_tkvar, fg='blue').grid(row=row, column=column, sticky=Tk.W)
    column += 1
    Tk.Label(frame0, text='fwhm', anchor=Tk.W                      ).grid(row=row, column=column, sticky=Tk.W)
    column += 1
    Tk.Label(frame0, width=width, textvariable=self.fwhm_tkvar, fg='blue', anchor=Tk.W).grid(row=row, column=column, sticky=Tk.W)

#    fig = Figure(figsize=(8,4), dpi=100, facecolor='black')
    fig, self.ax = plt.subplots(figsize=(5, 3))  # inch
    self.canvas = FigureCanvasTkAgg(fig, master=frame1)
    self.canvas.get_tk_widget().pack()

  def upload2BD_M2(self):
    try:
       m2_pos = self.focus_best_tkvar.get()
       req = 'update instruments set M2_POS0=' + str(m2_pos) + ' where ics=\'NBI\''
       conn = psycopg2.connect("dbname='sai2p5' user='ocsuser' host='192.168.10.87' password='?ocsuser=' ")
       cursor = conn.cursor()
       cursor.execute(req)
       conn.commit()
    except Exception as e:
      printe(e)

  def set_best(self):
     icsmodule.setFocus(self.focus_best_tkvar.get())

  def pvcallback_m2(self, pvname=None, value=-1, char_value=-1, **kw):
     self.m2_pos.set(value)

  def calc_focus_psf(self):
    printd('calc_focus_psf')
    try:
      self.initFocusPSF()
      self.calc_focus('psf')
    except Exception as e:
      printe(e)
    printd('end of calc_focus_psf')

  def calc_focus_donuts(self):
    try:
      self.calc_focus('donuts')
    except Exception as e:
      printe(e)

  def calc_focus(self, method):
    print('calc_focus method =', method)

    self.ds9 = ds9nbi.find_ds9()
    if self.ds9 is None:
      str_error = 'Open ds9 window, load fits file and try again'
      mainApp.messages.error_message.set(str_error)
      printe(str_error)
      return
    self.fullpath2focusfile = ds9nbi.get_path_to_file()
    if self.fullpath2focusfile is None:
       str_error = 'Open DS9, load fits file and try again'
       mainApp.messages.error_message.set(str_error)
       return
    if method.upper() == 'DONUTS':
      target = self.execute_donuts
    else:
      target = self.execute_psf
    target()

  def psf_on_end(self):
    if self.poly is not None:
      self.focus_best_tkvar.set(round2float(self.focus_best, 3))
      self.fwhm_tkvar.set(round2float(self.fwhm,2))
      print(' ------------ ', round2float(self.fwhm,2))
      self.draw_on_canvas(self.foci, self.fwhm_array, self.poly)

  def donuts_on_end(self):
    self.focus_best_tkvar.set(round2float(self.focus_best, 3))
    self.fwhm_tkvar.set(round2float(self.fwhm,2))
    print(' ------------ ', round2float(self.fwhm,2))

  def draw_on_canvas(self, x, y, poly):
    self.ax.clear()
    self.ax.scatter(x, y)
    nsteps = 100.0
    step = (x[-1] - x[0]) / nsteps
    x_detailed = np.arange(x[0],x[-1]+step,step)
    if poly is not None:
       self.ax.plot(x_detailed, poly(x_detailed))
    self.ax.set(xlabel='foci, m2_pos', ylabel='psf, asec', title='Looking for the best focus')
    self.canvas.draw()

  def execute_psf(self):
    print('execute_psf')
    res = focus_module.calc_focus_psf(self.ds9, \
       self.fullpath2focusfile, self.foc_start.get(), self.foc_step.get(), self.nsteps.get(), SCALE, mainApp.messages.error_message)
    print('res =', res)
    if res is None: return
    self.focus_best, self.fwhm, self.foci, self.fwhm_array, self.poly  = res
#    print('focus_best, fwhm_best, foci, fwhm_array, poly =', self.focus_best, self.fwhm, self.foci, self.fwhm_array, self.poly)
#    self.draw_on_canvas(foci, fwhm_array, poly)
    print(self.focus_best, self.fwhm)
#    input('press something...')
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'calcfocus psf end')

  def execute_donuts(self):
    res = focus_module.calc_focus_donuts(self.ds9, \
       self.fullpath2focusfile, self.foc_start.get(), self.foc_step.get(), self.nsteps.get(), mainApp.messages.error_message)
    if res is None: return
    self.focus_best, self.fwhm = res
    print(self.focus_best, self.fwhm)
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'calcfocus donuts end')

  def run_focusing(self):
    printd('Run focusing')
    altstep_deg = sec2deg(self.altstep_sec.get())
    scenario_string = ''
    printd ('------------ FOCI ------------------:')
    self.dero_start = DERO.get()
    self.parang_start = PARANG.get()
    for step in range(0, self.nsteps.get()):
      preclear = 1 if step == 0       else 0
      readout = 0
#      readout  = 1 if step == self.nsteps.get()-1 else 0
      foc = self.foc_start.get() + step * self.foc_step.get()
      printd(foc)
      if debug_focus:
        scenario_string = scenario_string + \
        'Exp,' + str(self.exptime.get()) + ',' + 'Light' + ',' + 'autofocusing' + ',' + str(preclear) + ',' \
        + str(readout) + '+'
      else:
        scenario_string = scenario_string + 'Focus,' + str(foc) + '+' + 'ShiftHor,0.0,' + str(altstep_deg) + \
        '+' + 'Exp,' + str(self.exptime.get()) + ',' + 'Light' + ',' + 'autofocusing' + ',' + str(preclear) + ',' \
        + str(readout) + '+'

# And the control shot:
    readout  = 1
    preclear = 0
    azstep_sec = 10.0
    azstep_deg = sec2deg(azstep_sec)
    if debug_focus:
      scenario_string = scenario_string + \
        '+' + 'Exp,' + str(self.exptime.get()) + ',' + 'Light' + ',' + 'autofocusing' + ',' + str(preclear) + ',' \
        + str(readout)
    else:
      scenario_string = scenario_string + 'ShiftHor,' + str(azstep_deg) + ',' + str(altstep_deg) + \
        '+' + 'Exp,' + str(self.exptime.get()) + ',' + 'Light' + ',' + 'autofocusing' + ',' + str(preclear) + ',' \
        + str(readout) + '+'
# And return:
      scenario_string = scenario_string + 'ShiftHor,' + str(-1.0*azstep_deg) + ',' + str(-1.0*(self.nsteps.get()+1)*altstep_deg)

#    scenario_string = scenario_string[:-1]  # remove extra '+'
    printd(scenario_string)
    if scenario.load(scenario_string) == -1: return -1
    return scenario.run()

class CommonButtons(Tk.Frame):
  def __init__(self, parent, *args, **kwargs):
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    self.parent = parent
    self.pult_window = -1
    self.focuser_window = -1
    self.geometry_window = -1
    butwidth = 10

#anchor=Tk.W, side=Tk.LEFT, fill=Tk.BOTH, pady=pdy, padx=12
    frame1 = Tk.Frame(self)
    frame1.pack(side=Tk.TOP, fill=Tk.BOTH, pady=2)
    frame2 = Tk.Frame(self)
    frame2.pack(side=Tk.TOP, fill=Tk.BOTH, pady=2)
    
    ocs_pult_button = Tk.Button(frame1, text='OCS pult', command=self.openOCSpult    )
    ocs_pult_button.pack(side = Tk.LEFT, padx = 5, pady = 2)
    geometry_buton  = Tk.Button(frame1, text='Geometry',command=self.openGeometry   )
    geometry_buton.pack(side = Tk.LEFT, padx = 5, pady = 2)
    focuser_button  = Tk.Button(frame1, text='Focuser' , command=self.openFocusWindow)
    focuser_button.pack(side = Tk.LEFT, padx = 5, pady = 2)
    flat_button     = Tk.Button(frame1, text='Flat'    , command=self.openFlatWindow )
    flat_button.pack(side = Tk.LEFT, padx = 5, pady = 2)
    wakeupFSU_button     = Tk.Button(frame1, text='FSU wakeup', command=fsu.wakeup   )
    wakeupFSU_button.pack(side = Tk.LEFT, padx = 5, pady = 2)
    shutter_button     = Tk.Button(frame1, text='Reset shutter', command=fsu.shutter_reset   )
    shutter_button.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.quitbutton   = Tk.Button(frame2, text='Quit', highlightbackground='blue', command=self.quit_script)
    self.quitbutton.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.init_exchange_button   = Tk.Button(frame2, text='Init Exch', highlightbackground='blue', command=self.initExchange)
    self.init_exchange_button.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.close_exchange_button   = Tk.Button(frame2, text='Close Exch', highlightbackground='blue', command=self.closeExchange)
    self.close_exchange_button.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.quitCCD3   = Tk.Button(frame2, text="QuitCCD3", highlightbackground='blue', command=self.run_quitCCD3)
    self.quitCCD3.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.initCCD3   = Tk.Button(frame2, text="InitCCD3", highlightbackground='blue', command=self.run_initCCD3)
    self.initCCD3.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.enable   = Tk.Button(frame2, text="Enable", command=self.enable_and_set_exposure_0)
    self.enable.pack(side = Tk.LEFT, padx = 5, pady = 2)

    ds9_start = ButtonQueue(frame2, queue=nbi_command_queue, text='DS9', queue_command=self.start_ds9_window)
#    ds9_start   = Tk.Button(frame2, text="DS9", bg='lightgreen', highlightbackground='green', command=self.start_ds9_window)
    ds9_start.pack(side = Tk.LEFT, padx = 5, pady = 2)

    self.wlist4disable = [self.quitCCD3, self.initCCD3]

    self.wlist4disable_init = [ocs_pult_button, geometry_buton, focuser_button, flat_button, self.close_exchange_button, \
    self.quitCCD3, self.initCCD3, self.enable]

  def start_ds9_window(self):
#    system('./ds9_start.sh &')
    subprocess.call('./ds9_start.sh')

  def openOCSpult(self):
    name = 'pult'
# Any MY Tk.Toplevel of CommonButtons MUST have name !
    for widget in self.winfo_children():
      if isinstance(widget, Tk.Toplevel) and widget.name == name:
        widget.destroy()
    self.pult_window = Tk.Toplevel(self)
    self.pult_window.name = name
    self.pult_window.wm_title('nbi OCS pult')
    icsmodule.pult(self.pult_window, M2_POS.get())

  def openGeometry(self):
    name = 'geometry'
# Any MY Tk.Toplevel of CommonButtons MUST have name !
    for widget in self.winfo_children():
      if isinstance(widget, Tk.Toplevel) and widget.name == name:
        widget.destroy()
    self.geometry_window = Tk.Toplevel(self)
    self.geometry_window.name = name
    self.geometry_window.wm_title('Geometry')
    geometry.gwindow(self.geometry_window)

  def openFlatWindow(self):
    name = 'flat'
# Any MY Tk.Toplevel of CommonButtons MUST have name !
    for widget in self.winfo_children():
      if isinstance(widget, Tk.Toplevel) and widget.name == name:
        widget.destroy()
    self.flat_window = Tk.Toplevel(self)
    self.flat_window.name = name
    self.flat_window.wm_title('Flatmaker')
    flatmaker.gwindow(self.flat_window)

  def openFocusWindow(self):
    name = 'focuser'
    for widget in self.winfo_children():
      if isinstance(widget, Tk.Toplevel) and widget.name == name:
        widget.destroy()
    self.focuser_window = Tk.Toplevel(self)
    self.focuser_window.name = name
    self.focuser_window.wm_title('Focuser')
    focuser.focus_window(self.focuser_window)

  def disable_due_to_exposure_all(self):
    mainApp.disable_due_to_exposure()
  
  def enable_and_set_exposure_0(self):
    self.enable_after_exposure_all()
    nbimodule.EXPOSURE_tkvar.set(0)

  def enable_after_exposure_all(self):
    mainApp.enable_after_exposure()

  def init_it_in_the_thread(self):
    icsmodule.initExchange()
    sleep(1)
    nbi.init()
#    sleep(1)
    scenario.initExchange()
#    sleep(1)
    mainApp.ocs.initExchange()
#    sleep(1)
    mainApp.filters.initExchange()
#    sleep(1)
    mainApp.messages.initExchange()
#    sleep(1)
    focuser.initExchange()
#    sleep(1)
    mainApp.expose_frame.initExchange()
#    sleep(1)
    flatmaker.initExchange()
#    sleep(1)
    geometry.initExchange()
    
    self.enable_after_init()
    print('Usjo')
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Init thread finished')

  def initExchange(self):
#     self.init_it_in_the_thread()  # Testing
    thread_init = threading.Thread(target=self.init_it_in_the_thread)
    thread_init.start()
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Init thread start')
    
  def closeExchange(self):
    thread_close = threading.Thread(target=self.close_it_in_the_thread)
    thread_close.start()
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Close thread start')

  def close_it_in_the_thread(self):
    printd( 'Close Exchange')
    gui_queue.do_periodicCall = False
    nbi.stop_ivy()    
#    scenario.sleep_wait_variable.set(1)
#    scenario.after_cancel(scenario.sleep_timer_id)
    for pv in [\
      FSU_FILTER,\
      FSU_SHUTTER,\
      FSU_MIRROR,\
      FSU_LAMP,\
      FSU_STATUS,\
      FSU_CMD_RUN,\
      M2_POS,\
      SAI25_ICS,\
      AGSTATUS,\
      SAI25_STATUS] :
      pv.clear_callbacks()
    self.disable_before_init()
#    system('xpaset -p KGO quit')
##    system('fuser -k -n tcp 33000')
##    subprocess.call('fuser -k -n tcp 33000')
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Close thread finished')

  def quit_script(self):
    printd( 'Quit')
#    nbi.stop_ivy()
#    self.closeExchange()
    printd( 'mainApp.quit')
    self.quit()

  def run_quitCCD3(self):
    printd ('Quit CCD3')
    nbi.quitCCD3()

  def run_initCCD3(self):
    printd('Init CCD3')
    nbi.initCCD3()

  def disable_due_to_exposure(self):
#    printd( self.__class__.__name__, 'disable_due_to_exposure')
    for w in self.wlist4disable:
      w.configure(state='disabled')

  def enable_after_exposure(self):
#    printd( self.__class__.__name__, 'enable_after_exposure')
    for w in self.wlist4disable:
      w.configure(state='normal')

  def disable_before_init(self):
    printd( self.__class__.__name__, 'disable_before_init')
    for w in self.wlist4disable_init:
      w.configure(state='disabled')
    mainApp.disable_due_to_exposure()
    mainApp.filters.disable_fsu()
    self.init_exchange_button.configure(state='normal')
    self.quitbutton.configure(state='normal')
    

  def enable_after_init(self):
    printd( self.__class__.__name__, 'disable_after_init')
    for w in self.wlist4disable_init:
      w.configure(state='normal')
    self.init_exchange_button.configure(state='disabled')
    self.quitbutton.configure(state='disabled')
    mainApp.enable_after_exposure()
    mainApp.filters.enable_fsu()


class MainApplication(Tk.Frame):
  def __init__(self, parent, *args, **kwargs):
    Tk.Frame.__init__(self, parent, *args, **kwargs)
    self.parent = parent

    self.wlist = []
    
    pdx = 5
    pdy = 5
    self.expose_frame  = ExposeFrame(self, borderwidth=1, relief=Tk.RIDGE, padx=pdx, pady=pdy)
    self.ocs           = OCS(        self, borderwidth=1, relief=Tk.RIDGE, padx=pdx, pady=pdy)


    self.filters  = Filters(self, borderwidth=1, relief=Tk.RIDGE, padx=pdx, pady=pdy)
    self.messages = Messages(self, borderwidth=1, relief=Tk.RIDGE, padx=pdx, pady=pdy) 
    self.buttons  = CommonButtons(self, borderwidth=1, relief=Tk.RIDGE, padx=pdx, pady=pdy)

    self.wlist = [self.expose_frame, self.ocs, self.filters, self.messages, self.buttons]

    self.expose_frame.grid (row=0, column=0, sticky=Tk.NSEW)
    self.ocs.grid     (row=0, column=1, sticky=Tk.NSEW)
    self.filters.grid (row=1, column=0, columnspan=2, sticky=Tk.NSEW)
    self.messages.grid(row=2, column=0, columnspan=2, sticky=Tk.NSEW)
    self.buttons.grid (row=3, column=0, columnspan=2, sticky=Tk.NSEW)
    
    for col in [0,1]:
      self.columnconfigure(col, weight=1)
    for row in [0,1,2]:
      self.rowconfigure(row, weight=1)

  def disable_due_to_exposure(self):
     for w in self.wlist:
       w.disable_due_to_exposure()

  def enable_after_exposure(self):
     for w in self.wlist:
       w.enable_after_exposure()

# ------------------ worker threads -------------------
def worker(myqueue, gui_queue):
    while True:
        task = myqueue.get()
        err_mess = 0
        try:
            res = task.on_start(*task.arglist)
        except Exception as e:
            err_mess = e
        gui_queue.put(task.on_end)
        myqueue.task_done()

class GUIqueue:
    def __init__(self, master):
      self.queue = queue.Queue()
      self.master = master
      self.do_periodicCall = True

    def periodicCall(self):
      while self.queue.qsize():
        try:
          item = self.queue.get(0)
          item()
        except queue.Empty:
          print('Queue.Empty')
          pass
      if self.do_periodicCall:
         self.master.after(200, self.periodicCall)

if __name__ == "__main__":
    global mainApp, nbi, fsu, geometry, ds9nbi, scenario, focuser, flatmaker, agutran_min # , exposew # m b change to expose
    ds9nbi = DS9NBI()
    root = Tk.Tk()

    nbi_command_queue = queue.Queue()
    gui_queue = GUIqueue(root)
    agutran_min = Tk.DoubleVar(value=100.0)
# Workers for long functions:
    n_thread = 3
    for i in range(0, n_thread):
      th = threading.Thread(target=worker, args=(nbi_command_queue, gui_queue.queue))
      th.daemon = True
      th.start()
    # GUI-thread queue polling:
    gui_queue.periodicCall()

    import icsmodule          # Uses Tk.Vars
    import geometry_module    # - " -
    import sciscenario_window # - " -
    import flatmaker_module   # - " -
    import nbimodule          # - " -
    import fsumodule          # - " -
    fsu = fsumodule.FSU()
    focuser = Focuser()       # Uses Tk.Vars -> import after Tk.Tk()

    nbi = nbimodule.NBI()
#    nbi.init()
    scenario = Scenario()
    geometry  = geometry_module.Geometry  (ds9nbi, nbi, icsmodule.shiftEquM)  # get two nbi functions
    flatmaker = flatmaker_module.Flatmaker(ds9nbi, fsu, nbi, scenario.fitsheader, icsmodule.shiftHor, agutran_min)     # can do exposure

    root.resizable(height = None, width = None)
    root.wm_title("nbi GUI")
    #SZH ico
    icon_path = os.path.dirname(os.path.abspath(__file__) )+'/icon_nbi.png'
    root.tk.call('wm', 'iconphoto', root._w, Tk.PhotoImage(file=icon_path) )

    mainApp = MainApplication(root, padx = 1, pady = 1, borderwidth=3, relief=Tk.RIDGE)
    mainApp.pack(side="top", fill="both", expand=True)
    mainApp.buttons.disable_before_init()

    root.mainloop()

