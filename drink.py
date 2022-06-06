import tkinter as Tk
root = Tk.Tk()
from PIL import Image, ImageTk
#on_image = Tk.PhotoImage(width=48, height=24)
#off_image = Tk.PhotoImage(width=48, height=24)
#on_image.put(("green",), to=(0, 0, 23,23))
#off_image.put(("red",), to=(24, 0, 47, 23))
scale=8
size=100
drunk_image = Tk.PhotoImage(file = 'drunk.png').subsample(scale)
sober_image = Tk.PhotoImage(file = 'sober.png').subsample(scale)
#drunk_image = Tk.PhotoImage(Image.open('drunk.png').resize(size,size))
#sober_image = Tk.PhotoImage(Image.open('sober.png').resize(size,size))
img = Image.open('drunk.png')
#img = img.resize((size,size))
#drunk_image = ImageTk.PhotoImage(img)
#sober_i = Image.open('sober.png').resize((size,size))
#drunk_image = Tk.PhotoImage(drunk_i)
#sober_image = ImageTk.PhotoImage(sober_i)
#drunk_image = Tk.PhotoImage(file = 'drunk.png').resize(size)
#sober_image = Tk.PhotoImage(file = 'sober.png').resize(size)

drink = Tk.IntVar(value=1)
#var2 = Tk.IntVar(value=0)
cb1 = Tk.Checkbutton(root, image=sober_image, selectimage=drunk_image, indicatoron=False,
                     onvalue=1, offvalue=0, variable=drink)
#cb1 = Tk.Checkbutton(root, selectimage=drunk_image, indicatoron=False,
#                     onvalue=1, offvalue=0, variable=drink)
label = Tk.Label(root, textvariable=drink)
                     
#cb2 = Tk.Checkbutton(root, image=off_image, selectimage=on_image, indicatoron=False,
#                     onvalue=1, offvalue=0, variable=var2)
#
cb1.pack(padx=20, pady=10)
label.pack(padx=20, pady=10)
#cb2.pack(padx=20, pady=10)
root.mainloop()
