#!/usr/bin/python3
import sys
from collections import OrderedDict
import datetime

import psycopg2
import os
from astropy.io import fits



def upload_metadata(cursor, key, value, measurement_id):
  currentTime = datetime.datetime.utcnow()
  timeString = currentTime.strftime('%Y-%m-%dT%H:%M:%S')

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
     keycomment = "image quality"
     keytype =  'string'
     cursor.execute("INSERT INTO params(\"NAME\",\"COMMENT\",\"DATEBEG\",\"DATEEND\",\"DATATYPE\") VALUES('%s','%s','%s','%s','%s')" \
                 %   (key, keycomment, datebegString, dateendString, keytype))

     cursor.execute("select \"ID\" from params WHERE \"NAME\" = '%s'" % key)
     paramRecord = cursor.fetchall()

  paramID = paramRecord[0][0]

  cursor.execute("select \"ID\",\"M_ID\" from metadata WHERE \"M_ID\" = (%d) and \"P_ID\" = (%d)" % (measurement_id,paramID))
  metadataRecord = cursor.fetchall()

  if not metadataRecord:
      cursor.execute("INSERT INTO metadata(\"P_ID\",\"M_ID\",\"VALUE\") VALUES(%d, %d, '%s')" % (paramID, measurement_id, value))
  else:
      metadataID = metadataRecord[0][0]
      cursor.execute("UPDATE metadata SET \"VALUE\"=('%s') WHERE \"ID\"=('%s')" % (value, metadataID))
      sys.stderr.write("BD contains metadata of M_ID = "+str(measurement_id)+" and P_ID = "+str(paramID)+" already. Updated.\n")
  return 0


def upload_iq(filename, iq_dict):
    print('upload_iq.upload_iq', filename, iq_dict)
    conn_string = "host='192.168.10.87' dbname='nbidata' user='nbiuser' password='?nbiuser='"
    sys.stderr.write("Connecting to database\n       ->"+conn_string+"\n")
    conn = psycopg2.connect(conn_string)
    cursor = conn.cursor()
    sys.stderr.write("Connected!\n")
    logfile = open('upload_iq.log','a')


    cursor.execute("select \"ID\",\"FILENAME\" from measurements WHERE \"FILENAME\" = \'%s\'" % (filename)) 
    measurement_record = cursor.fetchall()
    if not measurement_record:
        sys.stderr.write('BD does not contain file '+filename+' yet. Skipped...\n')
        return -1
    measurement_id = measurement_record[0][0]

    currentTime_str = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
    logfile.write(currentTime_str+" "+filename+"\n")
    for key in iq_dict:
        sys.stdout.write(key+" "+str(iq_dict[key])+" ")
        upload_metadata(cursor, key, str(iq_dict[key]), measurement_id)

    conn.commit()
    logfile.close()
    print('upload_iq.upload_iq: return')

def check_fits(fullname):
        path, filename = os.path.split(fullname)

        if 'FOCUS' in filename.upper() or 'TMP' in filename.upper():
            sys.stderr.write("Temporary file "+filename+" has been skipped\n")
            return -1

# Check if it is not fits file
        try:
            hdulist = fits.open(fullname)
        except IOError:
            sys.stderr.write(fullname+" is not a fits file, skipped\n")
            return -1
    
# Check size
        if len(hdulist) > 1:
           hdunum = 1
        else:
           hdunum = 0

        if ( (hdulist[hdunum].header.get('NAXIS1') != 4296) or (hdulist[hdunum].header.get('NAXIS2') != 4102)):
            sys.stderr.write(fullname+": Unsuitable Naxis1 or(and) Naxis2 = "+str(hdulist[hdunum].header.get('NAXIS1'))+\
            " "+str(hdulist[hdunum].header.get('NAXIS2'))+"\n")
            return -1
        return 0

if __name__ == '__main__':
        from iq2bd import iq

        if len(sys.argv) != 2:
            sys.exit("Usage: ./upload_iq.py fullfilename(*.fits)")
        fullname = sys.argv[1]

        if check_fits(fullname) < 0:
            sys.exit(-1)

        res = iq.do(fullname)
        if type(res) is int:
            sys.stderr.write(fullname+" : uploading has been skipped\n")
            sys.exit(-1)

        iq_dict = res
        path, filename = os.path.split(fullname)
        upload_iq(filename, iq_dict)

