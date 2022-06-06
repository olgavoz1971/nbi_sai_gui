#!/usr/bin/python
import sys
import datetime

tol_hour = 0.03

infile = sys.argv[1]
#str_date = '2021-03-02T21:36:08'
str_date = infile.replace('work_foc_','').replace('.fits','')
then = datetime.datetime.strptime(str_date,'%Y-%m-%dT%H:%M:%S')

with open('autofocus.log') as f: lines = f.readlines()
for line in lines:
  s = line.find(': -f')
  if s < 0: continue
  strl = line[:s]
  logtime = datetime.datetime.strptime(strl,'%Y-%m-%dT%H:%M:%S')

  t_delta = then - logtime
  dt = t_delta.days*24 + t_delta.seconds/3600 - 3
  if dt < 0 and abs(dt) < tol_hour:
    args_line = line[s+1:-1]
    args_line = args_line.replace('-f work_foc_ldac.cat', ' ')
    print(args_line)


