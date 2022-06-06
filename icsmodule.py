#!/usr/bin/python

# Module for nbi camera
# ICSserver is started by g1.py or on its own on the binary (nbi-machine), send active
# commands to the OCS-client (telescope)
port_ICSserver = 33000 # (It's my port on binary/ I'm ICS-server
ocs_client_connection = None # no connection from OCS as client
#OCS_IP = '192.168.15.51'  # OCS client for active commands to the telescope (protoOCS by now)
#OCS_IP = '192.168.15.12'  # OCS client for active commands to the telescope (ocsd)

ocsd_host    = '192.168.15.12'
oldocs_host  = '192.168.15.51'
ocsd_default_host = ocsd_host
OCS_IP = ocsd_host
ocs_host_choices  = [ocsd_host, oldocs_host]
# TMP !!!! get  it from BD
ics_choices = ['TDS', 'NBI', 'ANC', 'SPP', 'MGL']

import tkinter as Tk

OCS_FOCUS_RUN_tkvar = Tk.StringVar(value='')      # socket-response from OCS-client
OCS_HOR_RUN_tkvar   = Tk.StringVar(value='')      # --
OCS_EQU_RUN_tkvar   = Tk.StringVar(value='')      # --

import debug
debug_ocs = debug.debug_ocs
epics = debug.epics
debug_nbi = debug.debug_nbi
import nbimodule

from jkclass import printd, printe, sec2deg
import threading
import time
import sys
import atexit

import re
import socket

OCS_RUN_CURRENT = -1   # = no tkvar yet

WHO_IS_WORKING = -1  # Who is working - focus. equ or hor
WHO_DICT = {'FOCUS':OCS_FOCUS_RUN_tkvar, 'EQU':OCS_EQU_RUN_tkvar, 'HOR':OCS_HOR_RUN_tkvar}

def closeClient():
  global ocs_client_lbl, ocs_client_connection
  printd('Send Quit to OCS Client')
  try:
    _send('quit')
  except Exception as e:
    printe('Exception in CloseClient:', e)

  ocs_client_lbl.set(False)
  print('type(ocs_client_connection)', type(ocs_client_connection), ocs_client_connection)
#  if ocs_client_connection > 0:
  if ocs_client_connection is not None:
    ocs_client_connection.close()
  ocs_client_connection = None


def _send(cmd):
  global ocs_client_connection
  printd( 'NBI-server to OCS-client _send:->',cmd)
#  if ocs_client_connection == -1:
  if ocs_client_connection is None:
    str_error = 'OCS-client connection is closed or not established yet'
    printe('icsmodule:_send '+ str_error)
    nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
    return -1
#  print '|'+cmd+'|'
  try:
    cmd = cmd+'\r'
    ocs_client_connection.send(cmd.encode())
    return 0
  except socket.error:
    str_error = 'OCS-client connection is closed'
    printe('icsmodule_send socket.error, closed connection?')
    nbimodule.ERROR_MESSAGE_tkvar.set(str_error)
    try:
      if (ocs_client_connection is not None): ocs_client_connection.close()
    except Exception as e:
      printe( 'icsmodule:_send exception', e)

    ocs_client_connection = None
    ocs_client_lbl.set(False)
    return -1



def OCS_Client_Thread(sock, addr):
# sock = ocs_client_connection here
  global ocs_client_connection, ocs_client_lbl
  global WHO_IS_WORKING

  printd( 'OCS-Client Thread has started')
  ocs_client_lbl.set(True)

  while True:

# select.select !!! to kill it

#      data = sock.recv(1024)
      try:
        data = ocs_client_connection.recv(1024).decode('utf-8')
        res = re.sub('\n', '', re.sub('\r','', data))
        printd( 'received: |'+res+'|')
        if not data:
          printe( 'icsmodule:OCS_Client_Thread: No data received from OCS, close')
          break
        ocs_response.set(res)
        WHO_tkvar = WHO_DICT.get(WHO_IS_WORKING, None)
        if WHO_tkvar is not None: WHO_tkvar.set(res)
        WHO_IS_WORKING = -1
      except Exception as e:
        printe( 'icsmodule ClientThread exception', e)

  if ocs_client_connection is not None:
    try:
      printe('icsmodule:OCS_Client_Thread tries to close the connection')
      ocs_client_connection.close()
    except Exception as e:
      printe( 'icsmodule ClientThread exception', e, 'sock=', ocs_client_connection)
  printd( 'ICS-server: OCS-client',addr,'has disconnected')
  ocs_client_connection = None
  ocs_client_lbl.set(False)

def _recv_timeout():
  global ocs_client_connection
  data = ''
  if ocs_client_connection is None:
    printe( "nbiserver: OCS-client connection is closed or not established yet")
    return
  try:
    ocs_client_connection.settimeout(5.0)
    data = ocs_client_connection.recv(1024).decode('utf-8')
    printd( '|' + data + '|')
    if data:
      printd( 'ICS-server: received from OCS' + str(len(data)) + ' bytes: |' +  res + '|')
      return res
      
    else:
      printd( 'Socket has been closed by client')
      ocs_client_connection.close()
      ocs_client_connection = None
      ocs_client_lbl.set(False)
      return -1

  except socket.error as e:
    printe( 'OCS-client connection error:', e)
    closeClient()

  return -1



# Listen to socket as server in the thread
def runICSserver():
  printd( "run ICS Server")
  global ocs_client_connection
  global ocs_client_lbl

  sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  sock.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR, 1)
  atexit.register(closeClient)
    
#  sock.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR, 1)
  try:
    sock.bind(('0.0.0.0', port_ICSserver))
    sock.listen(1)
    printd( 'runICSserver: Socket is ready')
  except Exception as e:
    printe('sock.bind or listen exception', e)
    return -1

  while True:
    printd( "ICS-server: Waiting for OCS-client connection...\n")
#    ocs_client_connection, addr = sock.accept()
    ocs_client_connection_, addr = sock.accept()
    printd( 'ICS-server has connected with OCS-client:', addr)
    if len(addr) > 1:
      ip_client = addr[0]
      printd( "IP:", ip_client)
      if ip_client != OCS_IP:
        printd( 'No, ip:', ip_client, ' is not my client. Disconnection')
        ocs_client_connection_.close()
        continue
        
      printd( 'Yes! My client!')
      ocs_client_connection = ocs_client_connection_
      ocs_client_lbl.set(True)

# for async respond:
      thread = threading.Thread(target = OCS_Client_Thread, args=(ocs_client_connection, addr))
      thread.setDaemon(True)
      thread.start()
  return 0

def initExchange():
  thread_ics = threading.Thread(target = runICSserver)
  thread_ics.daemon = True
  thread_ics.start()

# Send active command direct  to the OCS:

def shiftHor(daz, dalt):    # in degrees
  global WHO_IS_WORKING
  WHO_IS_WORKING = 'HOR'
  commstring = 'CORRECT DAZ ' + str(daz) + ' DALT ' + str(dalt)
  printd('icsmodule: shiftHor', commstring)
  return _send(commstring)

def shiftAlt(dalt):
  return shiftHor(0.0, dalt)

def shiftAz(daz):
  return shiftHor(daz, 0.0)

def shiftAzm(dazm):
  commstring = 'CORRECT DAZM ' + str(dazm)
  printd(commstring)
  return _send(commstring)

def moveCover(frac):
  commstring = 'COVER FRAC ' + str(frac)
  printd( commstring)
  return _send(commstring)

def shiftEqu(dra, ddec, ddero):
  global WHO_IS_WORKING
  WHO_IS_WORKING = 'EQU'
  commstring = 'SHIFTEQU DRA ' + str(dra) + ' DDEC ' + str(ddec) + ' DDERO ' + str(ddero)
  printd('icsmodule: shiftEqu', commstring)
  return _send(commstring)

def shiftEquM(dram, ddec, ddero):
  global WHO_IS_WORKING
  WHO_IS_WORKING = 'EQU'
  commstring = 'SHIFTEQU DRAM ' + str(dram) + ' DDEC ' + str(ddec) + ' DDERO ' + str(ddero)
  printd('icsmodule: shiftEquM', commstring)
  return _send(commstring)

def shiftDec(ddec):
  dra   = 0.0
  ddero = 0.0
  return shiftEqu(dra, ddec, ddero)

def shiftRA(dra):
  ddec  = 0.0
  ddero = 0.0
  return shiftEqu(dra, ddec, ddero)

def shiftRAM(dram):
  ddec  = 0.0
  ddero = 0.0
  return shiftEquM(dram, ddec, ddero)

def shiftDero(ddero):
  dra   = 0.0
  ddec  = 0.0
  return shiftEqu(dra, ddec, ddero)

def setFocus(m2pos):
  global WHO_IS_WORKING
  WHO_IS_WORKING = 'FOCUS'
  commstring = 'FOCUS M2 ' + str(m2pos)
  printd( commstring)
  res = _send(commstring)
#  print('icsmodule: setFocus res =', res, type(res))
  return res

def shiftFocus(d_m2pos):
  commstring = 'FOCUS DM2 ' + str(d_m2pos)
  printd(commstring)
  return _send(commstring)

def aguOn():
  commstring = 'AUTOGUIDE GUIDE 1'
  printd(commstring)
  return _send(commstring)

def aguOff():
  commstring = 'AUTOGUIDE GUIDE 0'
  printd(commstring)
  return _send(commstring)


def newICS(ics_id):
  commstring = 'NEWICS ' + ics_id
  printd( commstring)
  return _send(commstring)

def targetDone(success=1):
  commstring = 'TARGETDONE SUCCESS ' + str(success)
#  commstring = 'TARGETDONE '
  printd('icsmodule:', commstring)
  return _send(commstring)

# -----------------------

if __name__ == '__main__':
  root = Tk.Tk()
  root.wm_title("ICS-module")

m2pos       = Tk.DoubleVar()
m2pos.set(10.2)
dalt        = Tk.DoubleVar()
daz         = Tk.DoubleVar()
dazm        = Tk.DoubleVar()
cover       = Tk.DoubleVar()
ddero       = Tk.DoubleVar()
dra         = Tk.DoubleVar()
dram        = Tk.DoubleVar()
ddec        = Tk.DoubleVar()
commstring  = Tk.StringVar()
ocs_response= Tk.StringVar()

ocs_client_lbl = Tk.BooleanVar()
ocs_client_lbl.set(False)


def setFocus_():
  return setFocus(m2pos.get())

def shiftAlt_dsec():
  dalt_deg = sec2deg(dalt.get())
  return shiftAlt(dalt_deg)

#def shiftAlt_():
#  return shiftAlt(dalt.get())

def shiftAz_dsec():
  daz_deg = sec2deg(daz.get())
  return shiftAz(daz_deg)

def shiftAzm_dsec():
  dazm_deg = sec2deg(dazm.get())
  return shiftAzm(dazm_deg)
#  return shiftAzm(daz.get())

def shiftDec_dsec():
  ddec_deg = sec2deg(ddec.get())
  return shiftDec(ddec_deg)

def shiftRA_dsec():
  dra_deg = sec2deg(dra.get())
  return shiftRA(dra_deg)

def shiftRAM_dsec():
  dram_deg = sec2deg(dram.get())
  return shiftRAM(dram_deg)

def shiftDero_():
  return shiftDero(ddero.get())

def moveCover_():
  return moveCover(cover.get())

def send_command_():
  _send(commstring.get())


def pult(parent, m2pos_current=-99):
  m2pos.set(m2pos_current)

  frame1 = Tk.Frame(parent, padx=10, pady=10, borderwidth=3, relief=Tk.RIDGE)
  frame1.pack(fill=Tk.X)

  foc_frame  = Tk.Frame(frame1, padx=10, pady=5)
  foc_frame.pack(fill=Tk.X)
  Tk.Label(foc_frame,  text = 'M2', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(foc_frame,  textvariable = m2pos, width=7).pack(side = Tk.LEFT) 
  Tk.Button(foc_frame, text="Move", fg="blue", width=7, command=setFocus_).pack(side = Tk.LEFT, padx=5) 

  alt_frame  = Tk.Frame(frame1, padx=10, pady=5)
  alt_frame.pack(fill=Tk.X)
  Tk.Label(alt_frame, text='Alt shift,\"', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(alt_frame, textvariable = dalt, width=7).pack(side = Tk.LEFT) 
  Tk.Button(alt_frame, text="Move", fg="blue", width=7, command=shiftAlt_dsec).pack(side = Tk.LEFT, padx=5) 

  az_frame  = Tk.Frame(frame1, padx=10, pady=5)
  az_frame.pack(fill=Tk.X)
  Tk.Label(az_frame, text='Az shift,\"', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(az_frame, textvariable = daz, width=7).pack(side = Tk.LEFT) 
  Tk.Button(az_frame, text="Move", fg="blue", width=7, command=shiftAz_dsec).pack(side = Tk.LEFT, padx=5) 

  azm_frame  = Tk.Frame(frame1, padx=10, pady=5)
  azm_frame.pack(fill=Tk.X)
  Tk.Label(azm_frame, text='AzM shift,\"', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(azm_frame, textvariable = dazm, width=7).pack(side = Tk.LEFT) 
  Tk.Button(azm_frame, text="Move", fg="blue", width=7, command=shiftAzm_dsec).pack(side = Tk.LEFT, padx=5) 

  dec_frame  = Tk.Frame(frame1, padx=10, pady=5)
  dec_frame.pack(fill=Tk.X)
  Tk.Label(dec_frame, text='Dec shift,\"', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(dec_frame, textvariable = ddec, width=7).pack(side = Tk.LEFT) 
  Tk.Button(dec_frame, text="Move", fg="blue", width=7, command=shiftDec_dsec).pack(side = Tk.LEFT, padx=5) 

  ra_frame  = Tk.Frame(frame1, padx=10, pady=5)
  ra_frame.pack(fill=Tk.X)
  Tk.Label(ra_frame, text='RA shift,\"', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(ra_frame, textvariable = dra, width=7).pack(side = Tk.LEFT) 
  Tk.Button(ra_frame, text="Move", fg="blue", width=7, command=shiftRA_dsec).pack(side = Tk.LEFT, padx=5) 

  ram_frame  = Tk.Frame(frame1, padx=10, pady=5)
  ram_frame.pack(fill=Tk.X)
  Tk.Label(ram_frame, text='RAM shift,\"', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(ram_frame, textvariable = dram, width=7).pack(side = Tk.LEFT) 
  Tk.Button(ram_frame, text="Move", fg="blue", width=7, command=shiftRAM_dsec).pack(side = Tk.LEFT, padx=5) 

  dero_frame  = Tk.Frame(frame1, padx=10, pady=5)
  dero_frame.pack(fill=Tk.X)
  Tk.Label(dero_frame, text='DDero', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(dero_frame, textvariable = ddero, width=7).pack(side = Tk.LEFT) 
  Tk.Button(dero_frame,text="Move", fg="blue", highlightbackground="green", width=7, command=shiftDero_).pack(side = Tk.LEFT, padx=5) 

  cover_frame  = Tk.Frame(frame1, padx=10, pady=5)
  cover_frame.pack(fill=Tk.X)
  Tk.Label(cover_frame, text='Cover frac', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Entry(cover_frame, textvariable = cover, width=7).pack(side = Tk.LEFT) 
  Tk.Button(cover_frame, text="Move", fg="blue", width=7, command=moveCover_).pack(side = Tk.LEFT, padx=5) 

  comm_frame   = Tk.Frame(frame1, padx=10, pady=5)
  comm_frame.pack(fill=Tk.X)
  Tk.Entry(comm_frame, textvariable=commstring, width=20).pack(side = Tk.LEFT)
  Tk.Button(comm_frame, text="Send", fg="blue", width=7, command=send_command_).pack(side = Tk.LEFT, padx=5)
  

# ------ Status ----
  frame2 = Tk.Frame(parent, padx=10, pady=10, borderwidth=3, relief=Tk.RIDGE)
  frame2.pack(fill=Tk.X)

  frame21 = Tk.Frame(frame2, pady=2)
  frame21.pack(side = Tk.TOP, fill=Tk.BOTH)
  Tk.Label(frame21, text='OCS-client', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Label(frame21, textvariable = ocs_client_lbl, fg="blue",  anchor=Tk.W, width=7).pack(side = Tk.LEFT)
  
  frame22 = Tk.Frame(frame2, pady=2)
  frame22.pack(side=Tk.TOP, fill=Tk.BOTH)
  Tk.Label(frame22, text='respond:', width=10, anchor=Tk.W).pack(side = Tk.LEFT)
  Tk.Label(frame22, textvariable = ocs_response,  anchor=Tk.W, fg="blue", width=25).pack(side = Tk.LEFT)

  frame3 = Tk.Frame(parent, padx=10, pady=10, borderwidth=3, relief=Tk.RIDGE)
  frame3.pack(fill=Tk.X)

  targetDonebutton  = Tk.Button(frame3, text="TargetDone", fg="blue", command=targetDone).pack(side = Tk.LEFT, padx=5)
  closeClientbutton = Tk.Button(frame3, text="CloseClient", fg="red", command=closeClient).pack(side = Tk.LEFT, padx=5)

# ----------------------------------------------------------------------

  if __name__ == '__main__':
    quitbutton = Tk.Button(frame3, text="Quit", fg="red", command=root.quit).pack(side = Tk.LEFT, padx=5)

if __name__ == '__main__':
  pult(root)  
  root.mainloop()

