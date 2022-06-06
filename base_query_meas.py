#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import numpy as np 
import time
import argparse
import psycopg2
import pandas as pd
import pandas.io.sql as sqlio
import os
import tkinter as Tk
from epics import caget

class info_window():
	def __init__(self, root=None):
		
		if root is not None:
			self.myWindow = Tk.Toplevel(root) 
		else:
			self.myWindow=Tk.Tk()
		self.myWindow.title("Measurements from base") 	

		Tk.Grid.rowconfigure(self.myWindow, 0, weight=1)
		Tk.Grid.columnconfigure(self.myWindow, 0, weight=1)

		self.frame_par = Tk.Frame(self.myWindow,pady=5)		#Frame for the image
		self.frame_par.grid(row=0, column=0)
		self.frame_text = Tk.Frame(self.myWindow,pady=5)		#Frame for buttons and scales
		self.frame_text.grid(row=1, column=0)

		Tk.Label(self.frame_par, text ="Object:").pack(side='left',padx=1,pady=2)

		self.run_b = Tk.Button(self.frame_par, text="Load!", fg="black", command=self.search)
		self.run_b.pack(side='left',padx=1,pady=2)

		self.textbox = Tk.Text(self.frame_text,height=32,width=150)
		self.textbox.pack()

		self.search()
		
		icpath= os.path.dirname(os.path.abspath(__file__) )+'/ico.png'
		self.myWindow.tk.call('wm', 'iconphoto', self.myWindow._w, Tk.PhotoImage(file=icpath) )

		if __name__ == '__main__':
			self.myWindow.mainloop()

	def search(self):
			self.textbox.delete('1.0', Tk.END)
			res  = get_meas() 
			if len(res)==0: res= 'Meases not Found!'
			self.textbox.insert(1.0, str(res) )

# print(get_last_obs('RW_Aur'))


def get_meas():
	conn = psycopg2.connect("dbname='sai2p5' user='ocsuser' host='192.168.10.87' password='?ocsuser=' ")

	track_id = caget('SAI25:TRACKID')
	# track_id = '1377529'
	mysql = """SELECT measurements.id, measurements.ord, measurements.ics ,measurements.exposure,measurements.snr,measurements.active,
				measurements.nframes, measurements.repeat_num, modes.mode, bands.band
			from measurements 
			JOIN tracks ON measurements.pointing_id = tracks.pointing_id  
			JOIN bands ON bands.id = measurements.band_id
			JOIN modes ON modes.id = measurements.mode_id				
			WHERE tracks.id = {} AND measurements.active=true ORDER BY measurements.ord;""" .format(track_id)

	meases = sqlio.read_sql_query(mysql, conn)

	sql_str = """SELECT ics_params.name, ics_params.default_value
	FROM ics_params WHERE ics_params.ics = 'TDS'"""
	default_params = sqlio.read_sql_query(sql_str, conn)
	for i, row in default_params.iterrows():
		meases[row['name']] = row['default_value']

	for i, mid in enumerate(meases['id']):
		sql_str = """SELECT ics_params.name, ics_values.value
		FROM ics_params join ics_values on ics_values.param_id =ics_params.id WHERE ics_values.meas_id = {}""".format(mid)
		meas_params = sqlio.read_sql_query(sql_str, conn)
		for  i, row in meas_params.iterrows():
			# print(row)
			meases[row['name']] = row['value']

	meases['SPEED'][meases['SPEED']=='0.05'] = 'MIN'
	meases['SPEED'][meases['SPEED']=='1']    = '1'
	meases['SPEED'][meases['SPEED']=='3']    = 'MAX'
	meases['MIRROR'] = meases['MIRROR']=='true'
	meases['INBEAM'] = meases['INBEAM']=='true'

	# meases['INBEAM'][meases['INBEAM']=='true']='v'
	# meases['INBEAM'][meases['INBEAM']=='3']='MAX'
	return meases

	conn.close()


if __name__ == '__main__':
	myw = info_window()



	# string="""WITH T AS ( 
 #        WITH t1 AS (
 #                SELECT measurements.id FROM modes,bands,pointings,measurements WHERE measurements.pointing_id=pointings.id AND 
 #                pointings.id={} AND pointings.ics=\'TDS\' AND bands.id=measurements.band_id AND modes.id=measurements.mode_id 
 #            ), 
 #            t2 AS ( 
 #                SELECT measurements.id,ics_params.name,ics_values.value FROM ics_params,modes,bands,pointings,measurements,ics_values WHERE ics_values.meas_id=measurements.id AND measurements.pointing_id=pointings.id AND 
 #                pointings.id={} AND pointings.ics=\'TDS\' AND bands.id=measurements.band_id AND modes.id=measurements.mode_id AND 
 #                ics_values.param_id=ics_params.id 
 #            ) 
 #            SELECT t2.id AS T1_INDEX, MIN(CASE WHEN 
 #            name=\'SLIT\' THEN value ELSE NULL END) AS 
 #            SLIT, MIN(CASE WHEN 
 #            name=\'SPEED\' THEN value ELSE NULL END) AS SPEED, MIN(CASE 
 #            WHEN 
 #            name=\'MIRROR\' THEN value ELSE NULL END) AS MIRROR, MIN(CASE 
 #            WHEN 
 #            name=\'INBEAM\' THEN value ELSE NULL END) AS INBEAM FROM t2 GROUP BY t2.id ORDER BY t2.id 
 #            )   
 #            SELECT measurements.ord,measurements.active,LOWER(modes.mode),bands.band,measurements.exposure,measurements.nframes,t.slit,t.speed,t.mirror,t.inbeam,measurements.id FROM modes,bands,pointings,measurements LEFT JOIN t ON t.t1_index=measurements.id WHERE measurements.pointing_id=pointings.id AND 
 #            pointings.id={} AND pointings.ics=\'TDS\' AND bands.id=measurements.band_id AND modes.id=measurements.mode_id 
 #             ORDER BY measurements.ord;""".format(1377930,1377930,1377930