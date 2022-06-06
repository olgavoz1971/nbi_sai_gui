#!/usr/bin/env python
import socket
import sys

IP_image_server = '192.168.15.51'
port_image_server = 33000
from jkclass import printd, printe

def sendFileReadyIQ(fullpath):
  printd('iq.sendFileReady', fullpath)
  sock = socket.socket()
  try:
    sock.connect((IP_image_server, port_image_server))
  except socket.error as e:
    printe( 'Image server',IP_image_server,port_image_server,e)
    return -1
  printd( 'connected', IP_image_server, port_image_server)
#  data = ("focus fname "+fname+" focbeg "
#           +str(m2pos0.get())+" focstep "+str(m2step.get())+" focnum "+str(nshits.get()+1)+
#           " dero "+str(dero_real.get())+" pa "+str(parallactic_angle.get())+"\n")
  data = 'iq '+fullpath+'\n'
  sock.send(data.encode())
  printd('->',data)
  sock.close()
  printd('sendFileReady: client socket is closed')
  return 0

def sendFileReadyFocus(fullpath, m2_beg, m2_step, nsteps, dero, parallactic_angle):
  printd('focus.sendFileReady', fullpath)
  sock = socket.socket()
  try:
    sock.connect((IP_image_server, port_image_server))
  except socket.error as e:
    printe( 'Image server',IP_image_server,port_image_server,e)
    return -1
  printd( 'connected', IP_image_server, port_image_server)
#  data = ("focus fname "+fname+" focbeg "
#           +str(m2pos0.get())+" focstep "+str(m2step.get())+" focnum "+str(nshits.get()+1)+
#           " dero "+str(dero_real.get())+" pa "+str(parallactic_angle.get())+"\n")

  data = ("focus fname "+fullpath+" focbeg "
           +str(m2_beg) + " focstep " + str(m2_step) + " focnum " + str(nsteps) +
           " dero " + str(dero) + " pa " + str(parallactic_angle) + '\n')

  sock.send(data.encode())
  printd('->',data)
  sock.close()
  printd('sendFileReady: client socket is closed')
  return 0
