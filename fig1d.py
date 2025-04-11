# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:23:06 2024

@author: dty
"""

import matplotlib.pyplot as plt
import numpy as np

with open('NA-TP-NP-SWP-4area.txt', 'r', encoding='utf-8') as file:
    data0 = file.readlines()

# 将读取到的每一行转换为float
per = [float(line.strip()) for line in data0]


print(per)
        
        
        
per[4]=1.0-per[4]



#%%c
labels = "North Atlantic",\
          "Tropical Eastern Pacific",\
          "North Pacific",\
          "Southwest Pacific",\
          "Others"
sizes = per
# color=['salmon', 'lightsalmon','tomato','lightpink','gray' ]
explode= (0.05, 0.05, 0.05, 0.05,0.05)
# color = ['#ff9999', '#c2c2f0', '#99ff99', '#ffcc99', 'gray']
color = ['#FFB6C1', '#98C9E0', '#A3D9A5', '#F1C6AC','grey']
# fig  = plt.figure(figsize=(20, 16),dpi=300)
fig, ax = plt.subplots(figsize=(6, 5),dpi=100)
ax.pie(sizes, labels=labels,autopct='%1.1f%%',explode=explode,labeldistance=1.1,\
       startangle=160, counterclock=True,  radius=1.5, rotatelabels=False, \
       colors=color,textprops = {'fontsize':20, 'color':'black','fontname':'Arial'\
                                 })
plt.savefig( 'Figure1d.eps', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()


