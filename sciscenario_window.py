#!/usr/bin/python3

import tkinter as Tk
from jkclass import printe

from collections import namedtuple
ScenarioRow = namedtuple('ScenarioRow', ['filter', 'exp', 'mode'], verbose=False)
ScenarioRow.__new__.__defaults__ = (None, None, None)


class SciScenarioTable():

  def __init__(self):
    self.filter_string_tkvar = Tk.StringVar(value='')
    self.values_string_tkvar = Tk.StringVar(value='')
    self.nrows = 5
    self.ncols = 3
    self.sci_scenario_rows = []
    self.filter_di = {}
    self.opened = False

    self.label_2dli = []
    self.sci_scenario_row_number_current = 0


  def gwindow(self, parent):

    self.parent = parent
    self.info_frame  = Tk.Frame(parent)
    self.table_frame = Tk.Frame(parent)

    self.info_frame.pack (fill=Tk.BOTH, side=Tk.TOP, padx=10, pady=3)
    self.table_frame.pack(fill=Tk.BOTH, side=Tk.TOP, padx=10, pady=3)

    pdx = 2; pdy = 2; w = 10
    row = 0; col = 0

    filter_string = Tk.Label(self.info_frame, textvariable=self.filter_string_tkvar, fg='blue')
    filter_string.pack(anchor=Tk.W, side=Tk.TOP, pady=3)
    values_string = Tk.Label(self.info_frame, textvariable=self.values_string_tkvar)
    values_string.pack(anchor=Tk.W, side=Tk.TOP, pady=3)
 
    headcolor = 'yellow'
    width=7
    Tk.Label(self.table_frame, text = 'Filter', width=width, bg=headcolor, relief='ridge', padx=20).grid(row=row, column=col); col += 1
    Tk.Label(self.table_frame, text = 'Exp, s', width=width, bg=headcolor, relief='ridge', padx=20).grid(row=row, column=col); col += 1
    Tk.Label(self.table_frame, text = 'Mode',   width=width, bg=headcolor, relief='ridge', padx=20).grid(row=row, column=col)

    for col in range(0, self.ncols):
      col_li = []
      for row in range(1, self.nrows+1):
         if row == 1: bgcolor = 'lightgreen'
         else: bgcolor = 'lightgray'
         label = Tk.Label(self.table_frame, text = '', width=width, bg=bgcolor, relief='ridge', padx=20)
         label.grid(row=row, column=col, sticky=Tk.EW, padx=pdx, pady=pdy)
         col_li.append(label)
      self.label_2dli.append(col_li)
#    Tk.Button(self.table_frame, text='Next', command=self.next).grid(row=self.nrows+1, column=0)
    self.redraw()
    self.opened = True

    print('gwindow:')
    print('self.sci_scenario_rows:', self.sci_scenario_rows)
    print('self.sci_scenario_row_number_current:', self.sci_scenario_row_number_current)


  def load(self, scenario_string):
    li = scenario_string.split('+')
    self.sci_scenario_rows = []
    filter = None
    exp    = None
    mode   = None

    for ent in li:
      ent_li = ent.split(',')
      key = ent_li[0]
      if key.upper() == 'FILTER': filter = ent_li[1]
      if key.upper() == 'EXP'   :
        exp  = ent_li[1] 
        mode = ent_li[2]
        self.sci_scenario_rows.append(ScenarioRow(filter = filter, exp = exp, mode = mode))

    self.filter_di = {}
    for ent in self.sci_scenario_rows:
      self.filter_di[ent.filter] = self.filter_di.get(ent.filter,0) + 1
      print('%7s %7s %7s' % (ent.filter, ent.exp, ent.mode))
    for filter in self.filter_di:
      self.filter_di[filter] = [self.filter_di[filter],0]
    print(self.filter_di)
    self.sci_scenario_row_number_current = 0
    if self.opened: self.redraw()
    

  def redraw(self):
    try:
      filter_string = ''
      values_string = ''
      for filter in self.filter_di:
        filter_string += '  %5s  ' % filter
#        values_string += '  %4i %4i  ' % (self.filter_di[filter][0], self.filter_di[filter][1])
        values_string += '  %5i  ' % (self.filter_di[filter][1])
      self.filter_string_tkvar.set(filter_string)
      self.values_string_tkvar.set(values_string)
      for i in range(0, self.nrows):
        sci_scenario_row_number = i + self.sci_scenario_row_number_current
        if sci_scenario_row_number >= len(self.sci_scenario_rows):
           for r in range(0, self.ncols): self.label_2dli[r][i].configure(text='')
           continue   
        row = self.sci_scenario_rows[sci_scenario_row_number]
        self.label_2dli[0][i].configure(text=row.filter)
        self.label_2dli[1][i].configure(text=row.exp, anchor=Tk.E)
        self.label_2dli[2][i].configure(text=row.mode)
    except Exception as e:
       printe(e)

  def next(self):
# TMP method
    print('self.sci_scenario_rows:', self.sci_scenario_rows)
    print('self.sci_scenario_row_number_current:', self.sci_scenario_row_number_current)
    if self.sci_scenario_row_number_current < len(self.sci_scenario_rows):
      filter = self.sci_scenario_rows[self.sci_scenario_row_number_current].filter
      if filter in self.filter_di:
        self.filter_di[filter][1] += 1
#    for filter in self.filter_di:
#      self.filter_di[filter] -= 1
        print('sciscenario_window: ', filter, self.filter_di)
    self.sci_scenario_row_number_current += 1
    if self.opened: self.redraw()

if __name__ == '__main__':


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



  import BDmeasurements

  measurements = BDmeasurements.Measurements()
  scenario_string = measurements.get_scenario_string()
  print('-------------------\n', scenario_string)
  scenario_string = unfold(scenario_string)
  print(scenario_string)

  #scenario_string = load_scenario()
  print(scenario_string)
  scenario_string = unfold(scenario_string)
  print('\n',scenario_string)

  print('\n\n\n')
  root = Tk.Tk()
  sci_scenario_table = SciScenarioTable()
  sci_scenario_table.load(scenario_string)
  sci_scenario_table.gwindow(root)

  root.mainloop()

#measurements = Measurements()
#scenario_string = measurements.get_scenario_string()
#print('-------------------\n', scenario_string)
#scenario_string = unfold(scenario_string)
#print(scenario_string)

