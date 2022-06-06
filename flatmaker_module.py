#!/usr/bin/python

import sys
import numpy as np
import os #, errno
import threading

import nbimodule
import upload2bd

from jkclass import printdd, printe, round2float, correct_outfile, MAX_PV_LENGTH, composeSAIoutfile, sec2deg
#maxn = MAX_PV_LENGTH
from time import sleep
from datetime import datetime

from astropy.io import fits

from epics import PV
import tkinter as Tk
from tkinter import filedialog

FSU_CMD_RUN     = PV('FSU:CMD_RUN')
FSU_FILTER      = PV('FSU:FILTER')

if __name__ == '__main__':
  root = Tk.Tk()
  root.wm_title('Flatmaker')


def do_nothing(*args):
  return

class Flatmaker():
  def __init__(self, ds9nbi, fsu = None, nbi = None, fitsheader = None, ics_shiftHor = do_nothing, agutran_min_tkvar=None):

    self.agutran = PV('SAI25:AGUTRAN')

    self.fsu = fsu
    self.nbi = nbi
    self.ds9nbi = ds9nbi
    self.fitsheader = fitsheader
    self.ics_shiftHor = ics_shiftHor
    self.agutran_min_tkvar = agutran_min_tkvar

    self.count   = Tk.DoubleVar(value=100000.0) # Desired count in ADU
    self.mintime = 0.5      # min exposition time in sec
    self.maxtime = 20       # max exposition time in sec

    self.twrate  = Tk.DoubleVar()
    self.flux    = Tk.DoubleVar(value=0)
    self.sunset  = Tk.BooleanVar(value=True)     # Sunrise(False) or sunset(True)

    self.xbeg = Tk.IntVar(value=3)     # xbeg of right detector for count estimation
    self.ybeg = Tk.IntVar(value=2000)  # ybeg for count estimation
    self.xsiz = Tk.IntVar(value=150)   # xsiz for count estimation'
    self.ysiz = Tk.IntVar(value=150)   # ysiz for count estimation'  # check
    self.dalt_dsec = Tk.IntVar(value=10)       # alt shift in dsec
    self.wait_interval = Tk.IntVar(value=30) # sec

    self.exptime   = Tk.DoubleVar(value=1.0)   # Will be calculated -- current exposition time in sec

#    self.exposure_wait_variable = Tk.IntVar(value = 0)
    self.exposure_wait_variable = nbimodule.EXPOSURE_tkvar
    self.fsu_wait_variable      = Tk.IntVar()
    self.sleep_wait_variable    = Tk.IntVar()
    self.exposure_timer_id      = -1
    self.fsu_timer_id           = -1
    self.sleep_timer_id         = -1

    self.current_message        = Tk.StringVar(value='Tra-la-la-----------------------------------12344567')
#    self.error_message          = Tk.StringVar(value='')
    self.error_message		= nbimodule.ERROR_MESSAGE_tkvar

    self.ncycles   = 5                            # get from run function
    self.filter_current = Tk.StringVar(value='?') # get from run function
    self.user_stop = Tk.BooleanVar(value=False)   # set to False in run method
    self.cycle     = Tk.IntVar(value=0)
    self.dtm_last  = datetime(1,1,1)

    self.filter_list = []

    self.parent = None
    self.wlist4disable = []
    

  def estiflux(self, fullfname, ybeg, xsiz, ysiz):   # we use long stripe or full image
    printdd('flatmaker:estiflux:', fullfname, ybeg, xsiz, ysiz)
  
    hdulist = fits.open(fullfname)
    if len(hdulist) < 2:
      str_error = 'estiflux:short hdulist of ' + fullfname
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
      printe('flatmaker : '+str_error)
      return -1,-1
    head = fits.header.Header()
    for hdu in hdulist:
      head = head + hdu.header

    img = np.float64(hdulist[1].data)
    hdulist.close()
  
    ytot,xtot = img.shape
    printdd( 'ybeg,ytot,xtot =', ybeg, ytot, xtot)
    xb_q2 = int(round(xtot/2))
    if(ybeg > ytot):
      str_error = 'estiflat: ybeg is outside the image, '+str(ybeg) + ' > '+ str(ytot)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
      printe(str_error)
      return -1,-1

# Recalculate with binning later !!!
    over = 50   #px
    printdd('over_q2:',ybeg,':',ybeg+ysiz, xb_q2,':',xb_q2+over-5)
    printdd('light  :',ybeg,':',ybeg+ysiz, xb_q2+over+5,':',xb_q2+xsiz)
    over_q2      = img[ybeg:ybeg+ysiz, xb_q2+1:xb_q2+over-6]
    lightzone_q2 = img[ybeg:ybeg+ysiz, xb_q2+over+5:xb_q2+xsiz  ]
    printdd('median(over_q2) =', np.median(over_q2))
    count = np.median(lightzone_q2 - np.median(over_q2))
    flux = count / head.get('EXPTIME',99999) 
    return  head.get('DATE-AVG'),flux

  def fsu_callback(self, pvname=None, value=-1, char_value=-1, **kw):
    printdd('flatmaker:fsu_callback value=', value)
    if value == 0:
      printdd('Callback fsu: Unlock')
      self.fsu_wait_variable.set(1)      # for wait_variable unlock

  def pvcallback_filter(self, pvname=None, value=-1, char_value=-1, **kw):
    self.filter_current.set(self.fsu.filterid2name(value))

  def initExchange(self):
#    NBI_EXPOSURE.add_callback(self.exposure_callback, with_ctrlvars=False)
    FSU_CMD_RUN.add_callback(self.fsu_callback, with_ctrlvars=False)
    FSU_FILTER.add_callback (self.pvcallback_filter, with_ctrlvars=False)
    self.filter_current.set (self.fsu.filterid2name(FSU_FILTER.get()))

  def sleep_timeout(self):
      self.sleep_wait_variable.set(1)
      printdd('Timer ready')

  def smart_sleep(self, sleep_time_sec):
      self.sleep_timer_id =  self.parent.after(sleep_time_sec*1000, self.sleep_timeout)
      printdd('Sleep. Wait ', sleep_time_sec, 'sec...')
      self.parent.wait_variable(self.sleep_wait_variable)
      printdd('Sleep unlocked')

  def pre_exposure(self):
    self.twrate.set(0.0)     # twighlight rate
    flux_last,self.dtm_last = 0,datetime(1,1,1) 

# --- First flux evaluation ------------
    self.current_message.set('Waiting for the twighlight')
    printdd(' ------------- We are waiting for the Twighlight ---------------')
    printdd(' -------------- First EVALUATION: ---- exptime0 =', self.exptime.get())
  
    exptime_start = self.exptime.get()
    while True:
  
      if self.user_stop.get():
        nbimodule.ERROR_MESSAGE_tkvar.set('User stop')
        printdd('User stop 0')
        return -1
      printdd('---- New cycle step --------- Flux evaluation on the small area ------')

      self.exptime.set(exptime_start)
      if __name__ == '__main__': return 0

      userfilename = 'focusing'
      userfilename = correct_outfile(self.nbi.get_filepath(), userfilename)
      progid = 'calibration'
      xsiz = self.xsiz.get() + 60 # overscan+
      arg_list = ['--xbeg', str(self.xbeg.get()), '--ybeg', str(self.ybeg.get()), '--xsiz', str(xsiz), '--ysiz', str(self.ysiz.get()),'--fname', userfilename,\
      '--time', str(self.exptime.get()),  '--preclear', '1',  '--readout', '1','--imagetype', 'flat',\
      '--progid', progid, '--filter', self.fsu.get_filter_name(), '--filterid', self.fsu.get_filterID(), '--fres','0']

      try:
        print(arg_list)
        res = self.nbi.load(arg_list)
        print(res)
        if res < 0:
          printe('flat load pre-exposure error')
          return -1

      except Exception as e:
        printe('bad arglist', arg_list)
        self.user_stop.set(True)
        return -1

      self.nbi.write_nbi_header()
      print(1)
      self.fitsheader.write_all_keys()
      print(2)

      exposure_and_reading_time_ms = self.nbi.calc_exposure_and_reading_time_sec() * 1000
      self.exposure_timer_id =  self.parent.after(exposure_and_reading_time_ms, self.exposure_timeout) # 100 sec
      print(3)
      res = self.nbi.start_expose()
      printdd('res=',res)
      if res < 0:
        printe('flat start pre-exposure error')
        return -1

# Wait here until expose ready
      self.current_message.set('Wait exposure...')
      printdd('Wait exposure...')
# Wait here until expose ready
      self.parent.wait_variable(self.exposure_wait_variable)
      if self.exposure_wait_variable.get() < 0: return -1
      self.current_message.set('Ready')
      printdd('Unlocked')
      self.parent.after_cancel(self.exposure_timer_id)

      fullpath = os.path.join(self.nbi.get_filepath(), self.nbi.get_filename())
      self.ds9nbi.view(fullpath)
#      upload2bd.upload(self.nbi.get_filepath(), self.nbi.get_filename())
      yb = 0
      dtm_string,flux = self.estiflux(fullpath, yb, self.xsiz.get(), self.ysiz.get())
      printdd('flux =', flux)
## TMP !!!
#      flux = 97000.0 
## TMP !!!
      if flux < 0: return -1
      self.flux.set(round(flux))
      dtm = datetime.strptime(dtm_string[0:19], "%Y-%m-%dT%H:%M:%S")
      if flux_last != 0:
         twrate = (flux - flux_last) / (dtm - self.dtm_last).seconds
         self.twrate.set( round(twrate,1))
         printdd(' ------------- FLUX EVALUATION: CURRENT TWRATE =',self.twrate.get())

      flux_last, self.dtm_last = flux, dtm
#      deltasec = (datetime.utcnow() - datetime.strptime(dtm_string[0:19], "%Y-%m-%dT%H:%M:%S")).seconds
      deltasec = (datetime.utcnow() - dtm).seconds
      printdd('deltasec =', deltasec)
      flux_predicted = flux + self.twrate.get()*deltasec
      if flux_predicted <= 0: flux_predicted = 0.00001
      exptime = round2float(self.count.get()/flux_predicted,2)
      exptime_stable = round2float(self.count.get()/flux,2)
      self.exptime.set(exptime)
      printdd(' -------------- EVALUATION: exptime',self.exptime.get(),'exptime_stable =', exptime_stable, '---- flux predicted', flux_predicted, 'flux=',flux,'\n')

      printdd (exptime_stable, self.exptime.get(), self.mintime)

      if (exptime_stable < self.mintime): self.exptime.set(exptime_stable)
      if (self.exptime.get() < self.mintime):
        printdd ('exptime < mintime')

      if (self.exptime.get() > self.maxtime): self.exptime.set(exptime_stable)
      if (self.exptime.get() > self.maxtime):
        printdd ('exptime > maxtime')

      if self.user_stop.get():
        self.current_message.set('User stop')
        printdd('User stop')
        return -1

      if self.sunset.get():
        printdd('sunset')
        if self.exptime.get() < self.mintime:
          msg = 'Light. Waiting... Sleep...'
          printdd(msg)
          self.current_message.set(msg)
          self.parent.update()
          self.smart_sleep(self.wait_interval.get()) # sec
#          sleep(self.wait_interval.get()) # secs
          self.parent.update()
          continue
        else:
          self.current_message.set('Sunset!')
          if self.exptime.get() > self.maxtime:
            msg = 'Too late!'
            printdd(msg)
            self.current_message.set(msg)
            return -2
          else: return 0

      if not self.sunset.get():
        printdd('sunrise')
        if self.exptime.get() > self.maxtime:
          msg = 'Dark. Waiting... Sleep...'
          printdd(msg)
          self.current_message.set(msg)
          self.parent.update()
          self.smart_sleep(self.wait_interval.get()) # sec
#          sleep(self.wait_interval.get()) # secs
          self.parent.update()
          continue
        else:
          msg = 'Sunrise!'
          printdd(msg)
          self.current_message.set(msg)
          if self.exptime.get() < self.mintime:
            msg = 'Too late!'
            printdd(msg)
            self.current_message.set(msg)
            return -2
          else: return 0
    return 0
# ----------------- MAIN MEASURING ----------------

  def main_exposure(self):
    nshifts = 0
    exptime_stable = self.exptime.get()
    flux_last = self.flux.get()
    for cycle in range(0, self.ncycles):

      self.cycle.set(cycle)
      if self.user_stop.get():
        nbimodule.ERROR_MESSAGE_tkvar.set('User stop')
        printdd('User stop')
        break
      printdd('\n\nCCCCCCCCCCCCCCCC YYYYYYYYYYYYYYYYY CCCCCCCCCCCCCCCCC LLLLLLLLLL EEEEEEEEEE CYCLE '+str(cycle))

      if self.exptime.get() > self.maxtime: self.exptime.set(exptime_stable)
      if self.exptime.get() > self.maxtime:
        msg = 'Too dark...'
        nbimodule.ERROR_MESSAGE_tkvar.set(msg)
        printdd(msg)
        break

      if self.exptime.get() < self.mintime: self.exptime.set(exptime_stable)
      if self.exptime.get() < self.mintime:
        msg = 'Too light...'
        nbimodule.ERROR_MESSAGE_tkvar.set(msg)
        printdd(msg)
        break

      userfilename = composeSAIoutfile()
      progid = 'calibration'

      arg_list = ['--fname', userfilename,\
      '--time', str(self.exptime.get()),  '--preclear', '1',  '--readout', '1','--imagetype', 'flat',\
      '--progid', progid, '--filter', self.fsu.get_filter_name(),'--filterid', self.fsu.get_filterID(), '--fres','1']

      if __name__ == '__main__': return 0

      try:
        res = self.nbi.load(arg_list)
        if res < 0:
          printe('flat load main-exposure error')
          return -1

      except Exception as e:
        printe('bad arglist', arg_list)
        self.user_stop.set(True)
        return -1
      self.nbi.write_nbi_header()
      self.fitsheader.write_all_keys()

      exposure_and_reading_time_ms = self.nbi.calc_exposure_and_reading_time_sec() * 1000
      self.exposure_timer_id =  self.parent.after(exposure_and_reading_time_ms, self.exposure_timeout) # 100 sec
      res = self.nbi.start_expose()
      if res < 0:
        printe('flat start main exposure error')
        return res

# Wait here until expose ready
      self.current_message.set('Wait exposure...')
      printdd('Wait exposure...')
# Wait here until expose ready
      self.parent.wait_variable(self.exposure_wait_variable)
      if self.exposure_wait_variable.get() < 0: return -1
      self.current_message.set('Ready')
      printdd('Unlocked')
      self.parent.after_cancel(self.exposure_timer_id)

      dalt_deg = sec2deg(self.dalt_dsec.get())
      printdd('shiftHor')
## TMP !!!
      res = self.ics_shiftHor(0, dalt_deg)
      if res < 0:
        printe('ics_shiftHor error')
        return res
## TMP !!!
# wait --- I'm too lazy to wait PV here...
      sleep(0.1)  
      nshifts = nshifts + 1

      fullpath = os.path.join(self.nbi.get_filepath(), self.nbi.get_filename())
      print ('%%%%%%%%%%%%%%%%%%%%%% ',fullpath)
      self.ds9nbi.view(fullpath)
#      upload2bd.upload(self.nbi.get_filepath(), self.nbi.get_filename())
      upload2bd.upload(self.nbi.get_filename(), self.nbi.get_filepath())
      print('%%%%%%%%%%%%%%%%%% uploaded ? %%%%%%%%%%%%%%%%%%%%%')
      printdd(self.ybeg.get(), self.xsiz.get(), self.ysiz.get())
      dtm_string,flux = self.estiflux(fullpath, self.ybeg.get(), self.xsiz.get(), self.ysiz.get())
## TMP !!!
#      flux = 97000 + nshifts*1000
## TMP !!!
      if flux < 0: return flux
      self.flux.set(round(flux))
#      dtm_string,flux = estiflux(path2image+fname+'.fits',args.ybeg,args.xsiz.get(),args.ysiz.get())
      dtm = datetime.strptime(dtm_string[0:19], "%Y-%m-%dT%H:%M:%S")
#     twRate calculation:
      if flux_last != 0:
        twrate = (flux - flux_last) / (dtm - self.dtm_last).seconds
        printdd ('flux,flux_last,twrate,interval_sec',flux,flux_last,twrate,(dtm - self.dtm_last).seconds)
        self.twrate.set( round(twrate,1))
#        print "---- current twilight rate is", self.twrate.get(), "ADU/s2"

        printdd(' !!!!!!!!!!! CURRENT TWRATE =',self.twrate.get())

      flux_last, self.dtm_last = flux,dtm
#      deltasec = (datetime.utcnow() - datetime.strptime(dtm_string[0:19], "%Y-%m-%dT%H:%M:%S")).seconds
      deltasec = (datetime.utcnow() - dtm).seconds
      flux_predicted = flux + self.twrate.get()*deltasec
      if flux_predicted <= 0: flux_predicted = 0.00001
      exptime = round2float( self.count.get()/flux_predicted,2 )
      exptime_stable = round2float(self.count.get()/flux,2)
      self.exptime.set(exptime)

      printdd('-------------- exptime', self.exptime.get(), 'exptime_stable =', exptime_stable)
      printdd('flux =',flux,'\n')
#6      sys.stdout.flush()

# return to the start position
    dalt_deg = sec2deg(self.dalt_dsec.get())
## TMP !!!
    res = self.ics_shiftHor(0, -1.0*nshifts*dalt_deg)
    if res < 0:
        printe('ics_shiftHor error')
        return res
## TMP !!!!!!!
    self.current_message.set('Done')     
    return 0
      
  def stop(self):
    self.user_stop.set(True)
# and to unlock smart_sleep if sleeping
    self.parent.after_cancel(self.sleep_timer_id)
    self.sleep_wait_variable.set(1) 

  def run(self):

# Check AGU
    if self.agutran.get() < self.agutran_min_tkvar.get():
        str_error = 'Take off AGU paw!'
        printe(str_error)
        self.current_message.set(str_error)
        return

# Check initCCD3 first
    if not self.nbi.isCCD3init():
      str_error = 'ccd3comm has not initiated, Press InitCCD3 first'
      printe(str_error)
      self.current_message.set(str_error)
      return

    if len(self.filter_list) < 1:
      printdd('Select scenario first')
      self.current_message.set('Select scenario first')
      return

    thread_flatmaker = threading.Thread(target=self.do_it_in_the_thread)
    thread_flatmaker.start()
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Init thread start')

  def do_it_in_the_thread(self):
    printdd(' Flatmaker RUN')

    nbimodule.ERROR_MESSAGE_tkvar.set('')
    self.user_stop.set(False)
    self.disable_due_to_exposure()
    for filter_name, exp in self.filter_list:
      printdd(' NNNNNNNNNNN EEEEEEEEEEEE WWWWWWWWWWWW SSSSSSSSSSS TTTTTTT EEEEEEEEEE PPPPPP filter,exp =', filter_name, exp)

      self.exptime.set(exp)
      res = self.set_filter(filter_name)
      if res < 0:
        self.current_message.set('Stop on error')
        break

#      sleep(5)
      res = self.pre_exposure()
      if res < 0: 
        if res == -2: continue
        else: break
      sleep(1.0)
      res = self.main_exposure()
      if res < 0: break
    self.nbi.fres()
    self.enable_after_exposure()
    print(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), 'Flatmaker thread finished')


  def exposure_timeout(self):
    step_name = 'exposure'
    msg = 'Flatmaker '+ step_name + ' timeout'
    printe(msg)
    nbimodule.ERROR_MESSAGE_tkvar.set(msg)
    self.parent.after_cancel(self.exposure_timer_id)
    self.exposure_wait_variable.set(-1)

  def fsu_timeout(self):
    step_name = 'fsu'
    msg = 'Flatmaker '+ step_name + ' timeout'
    printe(msg)
    nbimodule.ERROR_MESSAGE_tkvar.set(msg)
    self.parent.after_cancel(self.fsu_timer_id)
    self.fsu_wait_variable.set(-1)

  def set_filter(self, filter_name):
    printdd('Flatmaker:set_filter', filter_name) 
    res = self.fsu.set_mirror_target('WFI')
    if res < 0: return self.stop(res)
    res = self.fsu.set_filter_target(filter_name)
    if res < 0: return self.stop(res)
    self.fsu_timer_id =  self.parent.after(self.fsu.max_move_time_sec*1000, self.fsu_timeout) # 100 sec
    res = self.fsu.move_fsu()
    if res < 0: return self.stop(res)

    self.current_message.set('Wait FSU...')
    printdd('Wait filter...')
# Wait here until expose ready
    self.parent.wait_variable(self.fsu_wait_variable)
    if self.fsu_wait_variable.get() < 0: return -1
    self.current_message.set('Ready')
    printdd('Unlocked')
    self.parent.after_cancel(self.fsu_timer_id)

    return 0


  def select_scenario(self):
    fname =  filedialog.askopenfilename(initialdir = ".",title = "Select file",filetypes = (("flat files","*.flt"),("all files","*.*")))
# Check if nothing have been selected !!!
    printdd(' scenario fname =', fname)

#    fname = 'sunrise.flt'
    if not os.path.exists(fname):
      str_error = 'Flat scenario file does not exist' 
      printe(str_error)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
      return -1

    self.ncycles = 5
    self.sunset.set(1)
    self.filter_list = []

    try:
      with open(fname,'r') as f:
        for line in f:
          if line[0] == '#': continue
          line = line.replace('\n','').replace('  ','').replace('  ','').strip()
          li = line.split(' ')
          if len(li) != 2: continue
          if li[0].upper() == 'NCYCLES' : self.ncycles = int(float(li[1]))
          elif li[0].upper() == 'OBSERVER': continue
          elif li[0].upper() == 'SUNSET': self.sunset.set(bool(float(li[1])))
          else: self.filter_list.append((li[0],float(li[1])))
    except ValueError as e:
      str_error = 'Select scenario error: ' + str(e)
      printe(str_error)
      nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
      return -1
    return 0


# -----------------------  GUI -----------------------

  def gwindow(self, parent):

    width = 11; pdx = 2; pdy = 1; pdx_frame = 0; pdy_frame = 1
    self.parent = parent
    self.wlist4disable = []

    frame01 = Tk.Frame(parent, padx=pdx_frame, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame01.pack(fill=Tk.X)

    row = 0; column = 0
    for text in ['sunset', 'twrate', 'flux', 'exptime', 'cycle','filter', 'stop']:
      Tk.Label(frame01, text=text, width=width, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
      column = column + 1

    row = row+1; column = 0
    for var in [self.sunset, self.twrate, self.flux, self.exptime, self.cycle, self.filter_current, self.user_stop]:
      Tk.Label(frame01, textvariable=var, width=width, fg='blue', anchor=Tk.CENTER).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1

    row = row+1
    column = 0   
    for text in ['dalt', 'wait [sec]', 'xbeg', 'ybeg', 'xsiz', 'ysiz', 'count[ADU]']:
      Tk.Label(frame01, text=text, width=width, anchor=Tk.CENTER ).grid(row=row, column=column, padx=pdx, pady = pdy, sticky=Tk.W)
      column = column + 1
    row = row+1
    column = 0

    for var in [self.dalt_dsec, self.wait_interval, self.xbeg, self.ybeg, self.xsiz, self.ysiz, self.count]:
      Tk.Entry(frame01, width=width, textvariable=var).grid(row=row, column=column, padx=pdx, pady=pdy, sticky=Tk.W)
      column = column + 1
# -------------------------- Messages ---------------
    frame02 = Tk.Frame(parent, padx=pdx_frame, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame02.pack(fill=Tk.X)
    row = 0; column = 0
    Tk.Label(frame02, textvariable=self.current_message, anchor=Tk.CENTER).grid(row=row, column=column, sticky=Tk.EW)
    row = row+1
    column = 0
    Tk.Label(frame02, textvariable=self.error_message, fg='black', anchor=Tk.CENTER).grid(row=row, column=column, sticky=Tk.EW)

# ------------------------- Buttons -----------------
    frame03 = Tk.Frame(parent, padx=pdx_frame, pady=pdy_frame, borderwidth=3, relief=Tk.RIDGE)
    frame03.pack(fill=Tk.X)
#    root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("jpeg files","*.jpg"),("all files","*.*")))
    select_button = Tk.Button(frame03, text='Select', highlightbackground='lightgreen', width=width, command=self.select_scenario)
    select_button.pack(side = Tk.LEFT, padx=pdx)
    play_button   = Tk.Button(frame03, text='Play'  , highlightbackground='lightgreen', width=width, command=self.run)
    play_button.pack(side = Tk.LEFT, padx=pdx)
    stop_button   = Tk.Button(frame03, text='Stop'  , highlightbackground='lightgreen', width=width, command=self.stop)
    stop_button.pack(side = Tk.LEFT, padx=pdx)
    stop_button
    self.wlist4disable.append(select_button)
    self.wlist4disable.append(play_button)

  def disable_due_to_exposure(self):
    for w in self.wlist4disable:
      w.configure(state='disabled')

  def enable_after_exposure(self):
    for w in self.wlist4disable:
      w.configure(state='normal')

if __name__ == '__main__':
  flatmaker = Flatmaker()
  flatmaker.gwindow(root)
  root.mainloop()
