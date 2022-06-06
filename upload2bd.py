#!/usr/bin/env python

bd = True

import sys
import os
import numpy as np
import datetime
from astropy.io import fits
import psycopg2

from astropy import units as u
from astropy.coordinates import SkyCoord

#object_dict = {' '           : '_',
#               'CQC'         : 'CQ_C',
#               'CXC'         : 'CX_C',
#               '_NBI'        : '',
#               '_SHIFTED'    : '',
#               '_VIS'        : '',
#               'BD+75D325'   : 'BD+75_325',
#               'BD+25D4655'  : 'BD+25_4655',
#               'BD+25-4655'  : 'BD+25_4655',
#               'BD+33D2642'  : 'BD+33_2642',
#               'M81_ULS-1'   : 'M81-ULS1',
#               'IRAS02143'   : 'IRAS02143+5852',
#               'IRAS07253'   : 'IRAS07253-2001',
#               'AGK+81D266'  : 'AG+81_266',
#               'GRW+70D5824' : 'GRW+70_5824',
#               'SEXTANSB'    : 'SEXTANS_B',
#               'J120922+295551' : '[SST2011] J120922.18+295559.7',
#               'VYTAU'          : 'VY_TAU',
#               'RYTAU'          : 'RY_TAU',
#               'RW AUR SHIFTED' : 'RW_AUR',
#               'AICMI'          : 'AI_CMI',
#               'V455AND'        : 'V455_AND',
#               'IRAS19356'      : 'IRAS19356+0754',
#               'RRTAU'          : 'RR_TAU',
#               'V407CYG'        : 'V407_CYG',
#               'V444CYG'        : 'V444_CYG',
#               'V616MON'        : 'V616_MON',
#               'ASASSN_'        : 'ASASSN-',
#               'DF_TAU_JET'     : 'DF_TAU',
#               'RW_AUR_SII_DECENTER'    : 'RW_AUR',
#               'SN_2017MF'      :       'SN2017MF',
#               '_DEROT+180'     : '',
#               'M81_X'          : 'M81-ULS1' }
#
#progid_dict = {'SN_TSVETKOV' : 'SN',
#               'WR_IGOR'     : 'WR',
#               'UXOR_LAMZIN' : 'UXOR',
#               'POSTAGB'     : 'PAGB' }

def handle_file(fullname):
  sys.stderr.write("BD Handle file : "+fullname+"\n")
  global cursor, logfile
  path, filename = os.path.split(fullname)
  if 'FOCUS' in filename.upper() or 'TMP' in filename.upper():
      sys.stderr.write("Temporary file "+filename+" has been skipped\n")
      return 0

  cursor.execute("select \"ID\",\"FILENAME\" from measurements WHERE \"FILENAME\" = \'%s\'" % (filename)) 
  row = cursor.fetchall()
  if row:
    sys.stderr.write('BD contains file '+filename+' already. Skipped...\n')
    return 0



#  global object_dict, progid_dict #progid_object_dict, progid_dict

  try:
    hdulist = fits.open(fullname)
  except IOError:
    sys.stderr.write(fullname+" is not a fits file, skipped\n")
    return 0


  observatory = "CMO"
  telescope = "2.5m"
  coordSrcString = 'T'
  ccdname = "nbi4k"

  try:
    accumtime = hdulist[0].header['EXPTIME']
    modeString = hdulist[0].header['IMAGETYP'].lower()
    filterString = 'UNKNOWN'
    filter = ''
    if 'FILTER' in hdulist['PRIMARY'].header:
        filter = hdulist['PRIMARY'].header['FILTER']
    if 'FILTERID' not in hdulist['PRIMARY'].header:
        if (modeString == 'light') or (modeString == 'flat'):
            sys.stderr.write(fullname+" : No FILTERID key\n")
            return -1
    else:
        filterString = hdulist['PRIMARY'].header['FILTERID']

    if 'A:' in filter and 'B:' in filter:
        filter = 'OUT'
        filterString = 'OUT'

    observerString = hdulist['PRIMARY'].header['OBSERVER'];
    progid = hdulist[0].header['PROGID'] #.upper()
    dateobs = hdulist[0].header['DATE-OBS']
  except KeyError:
    sys.stderr.write(fullname+" : No EXPTIME,FILTERID,OBSERVER,IMAGETYP, PROGID or DATE-OBS key\n")
    return -1


  if 'CURRA' in hdulist[0].header:
     rastring = hdulist[0].header['CURRA'];
  elif 'RA' in hdulist[0].header:
     rastring = hdulist[0].header['RA'];
  elif 'TARRA' in hdulist[0].header:
     rastring = hdulist[0].header['TARRA'];
  else:
    sys.stderr.write(fullname+" : No RA key\n")
    rastring = rastring = '00:00:00'
#    return -1

  if 'CURDEC' in hdulist[0].header:
     decstring = hdulist[0].header['CURDEC'];
  elif 'DEC' in hdulist[0].header:
     decstring = hdulist[0].header['DEC'];
  elif 'TARDEC' in hdulist[0].header:
     decstring = hdulist[0].header['TARDEC'];
  else:
     sys.stderr.write(fullname+" : No DEC key\n")
     decstring = '00:00:00'
#     return -1

  naxis1 = hdulist[1].header['NAXIS1']
  naxis2 = hdulist[1].header['NAXIS2']

  if (modeString == 'flat') or (modeString == 'bias') or (modeString == 'dark'):
      if (rastring == 'hh:mm:ss') or (rastring == 'unknown') or (rastring == ''):
          rastring = '00:00:00'
      if (decstring == 'dd:mm:ss') or (decstring == 'unknown') or (decstring == ''):
          decstring = '00:00:00'

  try:
     radec = SkyCoord(rastring+' '+decstring, unit=(u.hourangle, u.deg))
  except ValueError:
     if (modeString != 'flat') and (modeString != 'bias') and (modeString != 'dark'):
         sys.stderr.write(fullname+" : Invalid ra or dec string, check it"+rastring+" "+decstring+"\n")
         return -1
         
  coordString = ( '(%.8f D,%.8f D)' % (radec.ra.rad, radec.dec.rad) )


  if 'OBJECT' in hdulist[0].header:
      object = hdulist[0].header['OBJECT']#.upper()
  else:
      sys.stderr.write(fullname+" fits key OBJECT is not found\n")
      object='unknown'

# reduce has been removed in python3 M.b. functools.reduce... but I remove it too
#  objectString = reduce(lambda x, y: x.replace(y, object_dict[y]), object_dict,object)
  objectString = object

  if progid == 'CALIBR':
      progid = 'CALIBRATION'
# And here
#  progIDstring = reduce(lambda x, y: x.replace(y, progid_dict[y]), progid_dict, progid)
  progIDstring = progid

# Check path to calibration file:
  filenameString = filename
  path = os.path.normpath(path)
#  lastdir = path.split(os.sep)[-1]
  if (modeString == 'bias') or (modeString == 'flat'):
#    if lastdir.lower() != modeString.lower():
#        sys.stderr.write(fullname+" : Check path to calibration file! IMAGETYP="+modeString+"\n")
    if progIDstring.upper() != 'CALIBRATION':
        sys.stderr.write(fullname+" : imagetyp="+modeString+" : wrong progid="+progIDstring+"\n")
        progIDstring = 'CALIBRATION'

  validV = "True"
  deletedV = "False"
  currentTime = datetime.datetime.strptime(dateobs[0:19], "%Y-%m-%dT%H:%M:%S")
  timeString = dateobs[0:19]



##     consider last observation
  cursor.execute("select \"ID\",\"OBJECT\",\"PROGID\" from observations ORDER BY \"ID\" DESC LIMIT 1")
  obsRecord = cursor.fetchall()
  if obsRecord:
             obsID = obsRecord[0][0]
             cursor.execute("select \"SCRDATE\" from measurements WHERE \"OBSERV_ID\" = %d ORDER BY \"ID\" DESC LIMIT 1" % (obsID))
             timeRecord = cursor.fetchall()
             if timeRecord:
                lastTime = timeRecord[0][0]
                dt = currentTime - lastTime
                delta_seconds = dt.seconds
#               delta_seconds = (dt.microseconds + 0.0 + (dt.seconds + dt.days * 24 * 3600) * 10 ** 6) / 10 ** 6               
                delta_seconds = delta_seconds - accumtime 
             else:
                delta_seconds = 86400

             lastObject = obsRecord[0][1]
             lastProgID = obsRecord[0][2]

  else:
              lastObject = "";
              lastProgID = "";
              delta_seconds = 100000   
##     if last observation has different object or progID or shutter status or ended more than 1 hour ago, create new observation
  if ( objectString != lastObject ) or ( progIDstring.upper() != lastProgID.upper() ) or ( delta_seconds > 3600.0 ):
              cursor.execute("INSERT INTO observations(\"OBJECT\",\"PROGID\",\"OBSERVER\") VALUES ('%s','%s','%s')"
                                                 % (objectString,progIDstring,observerString));
              cursor.execute("select \"ID\",\"OBJECT\" from observations ORDER BY \"ID\" DESC LIMIT 1")
              obsRecord = cursor.fetchall()
              obsID = obsRecord[0][0]        

##     create new measurement

  cursor.execute("INSERT INTO measurements(\"OBSERV_ID\",\"CCDNAME\",\"OBSERVATORY\",\"TELESCOPE\",\"FILENAME\",\"SCRDATE\",\"ACCUMTIME\",\"COORDEQU\",\"COORDEQUSRC\",\"MODE\",\"BAND\",\"VALID\",\"DELETED\")  \
                                         VALUES(%d,         '%s',        '%s',            '%s',         '%s',        '%s',        %f,        '%s',          '%s',        '%s',     '%s',    '%s',      '%s')"
                                          %   (obsID,      ccdname,   observatory,    telescope,  filenameString, timeString, accumtime, coordString, coordSrcString, modeString, filterString, validV, deletedV))
#  cursor.execute("INSERT INTO measurements(\"OBSERV_ID\",\"OBSERVATORY\",\"TELESCOPE\",\"FILENAME\",\"SCRDATE\",\"ACCUMTIME\",\"COORDEQU\",\"COORDEQUSRC\",\"MODE\",\"BAND\",\"VALID\",\"DELETED\")  \
#                                         VALUES(%d,          '%s',            '%s',         '%s',        '%s',        %f,        '%s',          '%s',        '%s',     '%s',    '%s',      '%s')"
#                                          %   (obsID,      observatory,    telescope,  filenameString, timeString, accumtime, coordString, coordSrcString, modeString, filterString, validV, deletedV))
  cursor.execute("select \"ID\" from measurements ORDER BY \"ID\" DESC LIMIT 1")                                                             

  measRecord = cursor.fetchall()

  measID = measRecord[0][0]


##     upload metadata
  for i in range(0,len(hdulist)):
          for key in list(hdulist[i].header.keys()):
              if key.upper() == 'COMMENT':
                  continue; 
              cursor.execute("select \"ID\",\"DATEBEG\",\"DATEEND\" from params WHERE \"NAME\" = '%s'" % key)
              paramRecord = cursor.fetchall()
                                            
              if paramRecord:
                 datebegString = paramRecord[0][1]
                 dateendString = paramRecord[0][2]
                 paramID = paramRecord[0][0]
                 if datebegString > currentTime:
                     datebegString = timeString
                     cursor.execute("UPDATE params SET \"DATEBEG\"=('%s') WHERE \"ID\"=('%s')" % (datebegString, paramID))
                 elif dateendString < currentTime:
                     dateendString = timeString
                     cursor.execute("UPDATE params SET \"DATEEND\"=('%s') WHERE \"ID\"=('%s')" % (dateendString, paramID))
              else:
                
                  datebegString = timeString
                  dateendString = timeString
                  keycomment = hdulist[i].header.comments[key]
                  keytype =  hdulist[i].header[key].__class__.__name__ # O! I know what a trick!
                  if keytype == 'str':
                      keytype = 'string'
                  cursor.execute("INSERT INTO params(\"NAME\",\"COMMENT\",\"DATEBEG\",\"DATEEND\",\"DATATYPE\") VALUES('%s',            '%s',           '%s',         '%s',        '%s')"
                                                  %   (key, keycomment, datebegString, dateendString, keytype))
                
                  cursor.execute("select \"ID\" from params WHERE \"NAME\" = '%s'" % key)
                  paramRecord = cursor.fetchall()
                                            
              paramID = paramRecord[0][0]
              cursor.execute("INSERT INTO metadata(\"P_ID\",\"M_ID\",\"VALUE\") VALUES(%d,%d,'%s')" % (paramID, measID, hdulist[i].header[key]))
                  
#  conn.commit()
  currentTime_now_str = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
  logfile.write(currentTime_now_str+" "+filename+"\n")

  return 1



def handle_dir(dirname):
  sys.stderr.write("handle_dir "+dirname+"\n")
  if not os.path.exists(dirname):
      sys.stderr.write(dirname+" : Invalid path, no such file or directory\n")
      return -1

  for entity in sorted(os.listdir(dirname)):
    fullname = os.path.join(dirname,entity)
    if os.path.isdir(fullname):
        p, shname = os.path.split(fullname)
        if shname.upper() != 'TEST':
            handle_dir(fullname)
    elif os.path.isfile(fullname):
        handle_file(fullname) 
    else:
        sys.stderr.write(fullname+" : I don't know what is it. Please handle it by yourself\n")
        return -1


def upload(dir_subdir_mask, path2data_base):

  global cursor, conn, logfile
 
  subpath, mask = os.path.split(dir_subdir_mask)
  path2data = os.path.join(path2data_base, subpath)
#  path2data = os.path.join('/home/oper/nbi4k/fileserver/', subpath)

  import fnmatch
  dirlist = fnmatch.filter(sorted(os.listdir(path2data)), mask)
  if len(dirlist) < 1:
      sys.stderr.write(path2data+"/"+mask+" : Nothing has been found\n")
      return -1
      

  conn_string = "host='192.168.10.87' dbname='nbidata' user='nbiuser' password='?nbiuser='"
  sys.stderr.write("Connecting to database\n	-> "+conn_string+"\n")
  conn = psycopg2.connect(conn_string)
  cursor = conn.cursor()
  sys.stderr.write("Connected!\n")
  print(conn_string)

  logfile = open('upload_iq.log','a')
  

  for entity in dirlist:
    fullname = os.path.join(path2data,entity)
    if not os.path.exists(fullname):
        sys.exit(fullname+" : Invalid path, no such file or directory")
    if os.path.isdir(fullname):
        handle_dir(fullname)
    else:
        handle_file(fullname)
#        sys.exit( "is a file, not dir. Check!")
  conn.commit()
  logfile.close()
    


if __name__ == '__main__':

        if len(sys.argv) != 3:
           sys.exit("Usage: ./upload2bd.py dirname_mask [f.e. 20181112 or 201801* or 20180112/W* ] path2data [/home/oper/nbi4k/fileserver/]")
        dir_subdir_mask = sys.argv[1]
        path2data_base = sys.argv[2]

        if bd:
          upload(dir_subdir_mask, path2data_base)
