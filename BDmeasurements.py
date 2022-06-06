#!/usr/bin/python3

import psycopg2
import pandas.io.sql as sqlio
import os
import numpy as np
from jkclass import PV

SAI25_TRACKID = PV('SAI25:TRACKID')

class Measurements():
  user_scenario_string = ''
  meases = None
  nrep_min = 0
  ncycles  = 0

  def init_scenario(self):
    self.user_scenario_string = ''
    self.nrep_min = 0
    self.ncycles  = 0

  def add2scenario(self):
    nrep_array = self.meases['repeat_num'].values
    nrep_masked = nrep_array[nrep_array > self.nrep_min]
    if len(nrep_masked) < 1:
      self.user_scenario_string = self.user_scenario_string[:-1]
      return False

    self.nrep_min = nrep_masked.min()
    self.ncycles = self.nrep_min - self.ncycles
    self.user_scenario_string += str(self.ncycles) + '*['

    for numrow, row in self.meases.iterrows():
      if row['repeat_num'] - self.nrep_min >= 0:
        filter,exp,mode,nframes = row['band'], row['exposure'], row['mode'], row['nframes']  # mode == light
        if nframes > 1:
          self.user_scenario_string += 'Filter,' + filter + '+' + str(nframes) + '*[Exp,' + str(exp) + ','+mode+']' + '+'
        else:
          self.user_scenario_string += 'Filter,' + filter + '+Exp,' + str(exp) + ','+mode+',' + '+'
    self.user_scenario_string = self.user_scenario_string[:-1] + ']+'
    return True

  def get_meas(self):
        self.init_scenario()
#        print('get_meas')
        conn = psycopg2.connect("dbname='sai2p5' user='ocsuser' host='192.168.10.87' password='?ocsuser=' ")
#        print('conn =', conn)
# TMP !!!!!!!
        track_id = SAI25_TRACKID.get()
##        track_id = 1381306
#        track_id = 1377944
        print('track_id =', track_id)
##        track_id = 1377529
# ---------------------

        mysql = """SELECT measurements.id, measurements.ord, measurements.ics ,measurements.exposure,measurements.snr,measurements.active,
                                measurements.nframes, measurements.repeat_num, modes.mode, bands.band
                        from measurements 
                        JOIN tracks ON measurements.pointing_id = tracks.pointing_id  
                        JOIN bands ON bands.id = measurements.band_id
                        JOIN modes ON modes.id = measurements.mode_id                           
                        WHERE tracks.id = {} AND measurements.active=true ORDER BY measurements.ord;""" .format(track_id)
        self.meases = sqlio.read_sql_query(mysql, conn)

# !!! TMP
#        self.meases['repeat_num']    = 5
#        self.meases['repeat_num'][0] = 4
#        self.meases['repeat_num'][2] = 4
#        print('self.meases =', self.meases)

        while True:
          if not self.add2scenario(): break

#        print(self.user_scenario_string)
        conn.close()
        return 0

  def get_scenario_string(self):
     self.get_meas()
     return self.user_scenario_string
  

def load_scenario():

    scenario_fname =  'RW_Aur.scn'
    scenario_name = scenario_fname.replace('.scn','')
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

    user_scenario_string = str(ncycles) + '*['

    for line in lines:
      if line[0] == '#': continue
      tu = line.split()
      if len(tu) < 2: continue
      if tu[0] == 'NCYCLES' or tu[0] == 'OBSERVER' or tu[0] == 'PROGID' or tu[0] == 'OBJECT' or tu[0] == 'FWHM': continue
#      print(tu)
      filter,exp = tu[0],tu[1]
#      if not fsu.filter_name_in_fsu(filter):
#      if filter not in fsu.FILTER_NAMES:
#        printe('Load scenario: Wrong filter name '+filter+' skipped')
#        continue
      user_scenario_string = user_scenario_string + 'Filter,'+filter + '+Exp,' + str(exp) + ',Light,' + userfilename + '+'

    user_scenario_string = user_scenario_string[:-1] + ']'
#    print(user_scenario_string)
    return user_scenario_string

def unfold1(s):   # work with []   f.e. '2 * [L,300 + 2*[B+D]+F,B+L, 1+2*[F,B+L,1]]'
    o = s.rfind('[')
    c = s.find(']',o)
    if c < 0:
      str_error = 'Wrong scenario '+s
#      print(s)
      return ''
    if s[o-1] != '*':
      str_error = 'Wrong scenario '+s
#      print(s)
      return ''
    sub = s[o+1:c]
    p = s[:o].rfind('+')
    try:
      num = int(s[p+1:o-1])
    except Exception as e:
      str_error = 'Scenario exception '+s
      print('Scenario exception ', s, e)
      return ''
    newsub = ''
    for i in range(0,num): newsub = newsub + sub + '+'
    newsub = newsub[:-1]
    sub2replace = s[p+1:c+1]
    return s.replace(sub2replace, newsub)

def unfold(stri):   # f.e. '2 * [L,300 + 2*[B+D]+F,B+L, 1+2*[F,B+L,1]]'
#    print('\n\n unfold stri =',stri, '\n')
#    stri = '2*[Filter,U+2*[Exp,300,Light]+Filter,B+Exp,100,Light]'
    stri = stri.replace(' ','')
    while True:
#      print(stri)
      if stri.find('[') < 0: break
      stri = unfold1(stri)
      if stri == '':  return stri
    stri = stri.replace('++','+')
    return stri

#measurements = Measurements()
#scenario_string = measurements.get_scenario_string()
#print('-------------------\n', scenario_string)
#scenario_string = unfold(scenario_string)
#print(scenario_string)

