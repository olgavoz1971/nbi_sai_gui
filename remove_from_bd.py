#!/usr/bin/env python3

bd = True

import sys
import os
import datetime
import psycopg2


def remove_from_tables(filename):
  conn_string = "host='192.168.10.87' dbname='nbidata' user='nbiuser' password='?nbiuser='"
  sys.stderr.write("Connecting to database\n	-> "+conn_string+"\n")
  conn = psycopg2.connect(conn_string)
  cursor = conn.cursor()
  sys.stderr.write("Connected!\n")

  cursor.execute("select \"ID\",\"FILENAME\",\"OBSERV_ID\" from measurements WHERE \"FILENAME\" = \'%s\'" % (filename)) 
  meas_row = cursor.fetchall()
  if not meas_row:
    sys.stderr.write('BD does not contain file '+filename+' yet. Skipped...\n')
    return 0
  meas_id = meas_row[0][0]
  observ_id = meas_row[0][2]
  print("meas_id=",meas_id, "observ_id=",observ_id, meas_row[0])

  cursor.execute("select * from metadata WHERE \"M_ID\" = %d" % (meas_id)) 
  meta_row = cursor.fetchall()
#  print meta_row[0]
  
  cursor.execute("select count(*)  from  measurements where \"OBSERV_ID\"=%d" % (observ_id)) 
  row = cursor.fetchall()
#  print "count measurements with observ_ID=", observ_id,"=",row[0][0]
  n_meas_in_observ = row[0][0]


  cursor.execute("select *  from  metadata_pre  WHERE \"M_ID\" = %d" % (meas_id)) 
  metadata_pre_row = cursor.fetchall()
  if metadata_pre_row:
#    print metadata_pre_row[0]
    answ = raw_input("The are some metadata_pre rows linked to our image. Remove anyway? Y/N (default N): ")
    if answ.upper() != 'Y':
      print("Skipped...")
      return 0
    cursor.execute("delete from  metadata_pre where \"M_ID\" = %d" % (meas_id))
    print("records with M_ID=", meas_id, "have been removed from metadata_pre")

  else:
    print("No metadata_pre linked to our image")
  

  delete_meas = True
  delete_observ = True
  if delete_observ:
    if n_meas_in_observ == 1:
      print("delete observation", observ_id)
    else:
      delete_observ = False

  if delete_meas:
     cursor.execute("delete from  metadata where \"M_ID\" = %d" % (meas_id))
     print("records with M_ID=",  meas_id, "have been removed from metadata")
     cursor.execute("delete from  measurements where \"FILENAME\" = \'%s\'" % (filename))
     print("record with ID=", meas_id, "has been removed from measurements")
  if delete_observ:
     cursor.execute("delete from observations where \"ID\" = %d" % (observ_id))
     print("record with ID=",observ_id, "has been removed from observations")
  
  conn.commit()

  return 0
  

if __name__ == '__main__':
        if len(sys.argv) != 2:
           print(sys.argv)
           sys.exit("Usage: ./remove_from_bd dirname_mask[f.e. 20181112 or 201801* or 20180112/W* ]")
        remove_from_tables(sys.argv[1])
