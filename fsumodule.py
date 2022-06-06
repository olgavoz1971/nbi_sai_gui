#!/usr/bin/env python
# FSU manager for NBI camera
# Epics version

import debug
import nbimodule    # ERROR_MESSAGE here

epics = debug.epics
debug_fsu = debug.debug_fsu

import sys
#data = " ".join(sys.argv[1:])
from time import sleep

from jkclass import PV, printd, printdd, printe
#maxn = MAX_PV_LENGTH
#from jkclass import SET_ERROR

FSU_FILTERA_LIST13 = PV('FSU:FILTERA_LIST13')
FSU_FILTERA_LIST46 = PV('FSU:FILTERA_LIST46')
FSU_FILTERB_LIST13 = PV('FSU:FILTERB_LIST13')
FSU_FILTERB_LIST46 = PV('FSU:FILTERB_LIST46')
FSU_FILTER_TARGET  = PV('FSU:FILTER_TARGET')
FSU_FILTER         = PV('FSU:FILTER')
FSU_MIRROR         = PV('FSU:MIRROR')
FSU_MIRROR_TARGET  = PV('FSU:MIRROR_TARGET')
FSU_LAMP           = PV('FSU:LAMP')
FSU_LAMP_TARGET    = PV('FSU:LAMP_TARGET')
FSU_INIT_TARGET    = PV('FSU:INIT_TARGET')

FSU_SHUTTER        = PV('FSU:SHUTTER')
FSU_SHUTTER_TARGET = PV('FSU:SHUTTER_TARGET')

FSU_CMD_RUN        = PV('FSU:CMD_RUN')
FSU_CMD_PARK       = PV('FSU:CMD_PARK')
FSU_CMD_INIT       = PV('FSU:CMD_INIT')
#NBI_ERROR          = PV('NBI:ERROR')
  

class FSU():
#  STATUS = {'LAMP':'?', 'FILTER':'?', 'SHUTTER':'?', 'MIRROR':'?'}
  WHEEL = {'FILTERA':'?', 'FILTERB':'?'}
  FILTER_IDS   = [] # list of filter_id's F001UBES, F002BBES...
#  FILTER_NAMES = [] # list of filter_names U, B,V,R,...

  mirror_target_choices  = [ 'TDS', 'WFI', 'CAL', 'HRS', 'LID']
  lamp_target_choices    = [ 'CON', 'ARC', 'OFF' ]
  shutter_target_choices = [ 'OPEN','CLOSE','AUTO']

  wheel_filter_names_A = []
  wheel_filter_names_B = []

  filter_id_dict = {}
  filter_name_dict = {}
  filter_alias_filename = 'filter_alias.dat'
  
  max_move_time_sec = 100

  def __init__(self):
    try:

      with open(self.filter_alias_filename, 'r') as f:
        txt = f.readlines()
        for line in txt:
          if line[0] == '#':
            continue
          l = line.split()
          if len(l) < 2:
            continue
          self.filter_id_dict[l[0]] = l[1]
          self.filter_name_dict[l[1]] = l[0]
    except IOError:
      sys.exit("NBIFSU: There are some problems with "+self.filter_alias_filename+" file\n")

    self.read_status()


  def read_status(self):
    printd( self.__class__.__name__,'read_status')
    if not epics:
#      self.STATUS['FILTER'] = 'F002BBES'
#      self.WHEEL = {'FILTERA':'F001UBES,F002BBES,F003VBES,F004RBES,F005IBES,F0586430',\
#               'FILTERB':'F0000000,F0596563,F006USDS,F007GSDS,F008RSDS,F009ISDS'}

      FSU_FILTERA_LIST13.put('F001UBES,F002BBES,F003VBES')
      FSU_FILTERA_LIST46.put('F004RBES,F005IBES,F0586430')
      FSU_FILTERB_LIST13.put('F0000000,F0596563,F006USDS')
      FSU_FILTERB_LIST46.put('F007GSDS,F008RSDS,F009ISDS')
      FSU_FILTER.put('F007GSDS')

    try:
      self.WHEEL['FILTERA'] = FSU_FILTERA_LIST13.get() + ',' + FSU_FILTERA_LIST46.get()
      self.WHEEL['FILTERB'] = FSU_FILTERB_LIST13.get() + ',' + FSU_FILTERB_LIST46.get()
    except Exception as e:
      str_error = 'fsumodule exception'
      printe(str_error, e)
      self.WHEEL['FILTERA'] = 'F0000000,F0000000,F0000000,F0000000,F0000000,F0000000'
      self.WHEEL['FILTERB'] = 'F0000000,F0000000,F0000000,F0000000,F0000000,F0000000'


    printd( self.WHEEL['FILTERA'])
    printd( self.WHEEL['FILTERB'])
    self.FILTER_IDS = []
    for key in self.WHEEL:
        self.FILTER_IDS = self.FILTER_IDS + self.WHEEL[key].split(',')

    printd( 'all filters:', self.FILTER_IDS)

#    fill FILTER_NAMES list f.e. ['U','B','V'...]
    self.wheel_filter_names_A = []
    self.wheel_filter_names_B = []

    wheelA_id = self.WHEEL['FILTERA'].split(',')
    wheelB_id = self.WHEEL['FILTERB'].split(',')

    for filter_id in wheelA_id:
      self.wheel_filter_names_A.append(self.filter_name_dict.get(filter_id, filter_id))

    for filter_id in wheelB_id:
      self.wheel_filter_names_B.append(self.filter_name_dict.get(filter_id, filter_id))

  def filterid2name(self, filter_id):
    filter_name = self.filter_name_dict.get(filter_id, filter_id)
    return filter_name

  def filtername2id(self, filter_name):
    filter_id = self.filter_id_dict.get(filter_name, filter_name)
    return filter_id

  def get_filter_name(self):
    return self.filterid2name(FSU_FILTER.get())

#  def get_filter(self):
#    filt_id = self.STATUS['FILTER']
#    return self.filter_name_dict.get(filt_id, filt_id)

#    if filt in self.wheel_filter_names_A + self.wheel_filter_names_B:
#      return self.filter_name_dict[filt]
#    else:
#      return filt

  def is_filter_ready(self):
    return FSU_FILTER_TARGET.get() == FSU_FILTER.get()

  def get_mirror(self):
    return FSU_MIRROR.get()

  def get_lamp(self):
    return FSU_LAMP.get()

  def get_shutter(self):
    return FSU_SHUTTER.get()

  def set_filter_target(self, filter_):
    filter_id = self.filter_id_dict.get(filter_, filter_)
    printd( self.__class__.__name__, 'set_filter_target: filter_id = ', filter_id)

    res = FSU_FILTER_TARGET.put(filter_id) # 'A0'..'B6','CA1'..'CB6' or A5,B1 or F001UBES
    if res < 0:
      str_error = 'FSU FILTER_TARGET error, code='+str(res)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
#      NBI_ERROR.put(str_error[:maxn])
      printe(str_error)
      printd(self.__class__.__name__, 'set_filter_target: result =', res)
    return res

  def set_filter(self, filter_):
    res = self.set_filter_target(filter_)
    if res < 0: 
      printe('FSU:set_filter error res =',res)
      return res
    return self.move_fsu()

  def set_mirror_target(self, mirror):
    res = FSU_MIRROR_TARGET.put(mirror)  # WFI or TDS or CAL or HRS or LID
    printd( self.__class__.__name__, 'set_mirror_target: mirror_target = ', mirror)
    if res < 0:
      str_error = 'FSU MIRROR_TARGET error, code='+str(res)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
#      NBI_ERROR.put(str_error[:maxn])
      printe(str_error)
      printd(self.__class__.__name__, 'set_mirror_target: result =', res)
      return res
    return 0

  def set_lamp_target(self, lamp):
    res = FSU_LAMP_TARGET.put(lamp)  # CON or ARC or OFF
    printd( self.__class__.__name__, 'set_lamp_target: lamp_target = ', lamp)
    if res < 0:
      str_error = 'FSU LAMP_TARGET error, code='+str(res)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
#      NBI_ERROR.put(str_error[:maxn])
      printe(str_error)
      printd(self.__class__.__name__, 'set_lamp_target: result =', res)
      return res
    return 0

  def set_shutter_target(self, shutter):
    res = FSU_SHUTTER_TARGET.put(shutter)  # OPEN CLOSE AUTO
    printd( self.__class__.__name__, 'set_shutter_target: shutter_target = ', shutter)
    if res < 0:
      str_error = 'FSU SHUTTER_TARGET error, code='+str(res)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
#      NBI_ERROR.put(str_error[:maxn])
      printe(str_error)
      printd(self.__class__.__name__, 'set_shutter_target: result =', res)
      return res
    return 0

  def shutter_reset(self):
    FSU_INIT_TARGET.put('SHUTTER')
    FSU_CMD_RUN.put(1)
#    caput FSU:INIT_TARGET SHUTTER
#    caput FSU:CMD_RUN 1

  def wakeup(self):
    FSU_CMD_PARK.put(1)
    sleep(1)
    FSU_CMD_INIT.put(1)

  def move_fsu(self):
    printd('FSU:move_fsu')
    if debug_fsu:
      return 0
    res = FSU_CMD_RUN.put(1)
    if res < 0:
      str_error = 'FSU CMD_RUN error, code='+str(res)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
#      NBI_ERROR.put(str_error[:maxn])
      printe(str_error)     
      return res
    return 0

  def filter_name_in_fsu(self,name):
    return name in self.wheel_filter_names_A + self.wheel_filter_names_B

  def get_filter_names_in_wheel(self, wheel_name):
    if wheel_name.upper() == 'A':
      return self.wheel_filter_names_A
    elif wheel_name.upper() == 'B':
      wheel_full_name = 'FILTERB'
      return self.wheel_filter_names_B
    else:
      sys.exit('NBIFSU: invalid Filter Wheel Name '+wheel_name)

  def get_filter_names(self):
     return self.wheel_filter_names_A + self.wheel_filter_names_B

  def get_filterID(self):
     return FSU_FILTER.get()
