#!/usr/bin/env bash 
xpaset -p KGO quit
sleep 1s
ds9 -title KGO &
sleep 8s
xprop -name "SAOImage KGO" -f WM_CLASS 8s -set WM_CLASS "ds9kgo"
