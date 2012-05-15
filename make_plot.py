#!/usr/bin/python

import csv

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
#rcParams['text.usetex'] = True
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['ytick.major.pad'] = 10.
rcParams['xtick.major.pad'] = 10.
rcParams['ytick.major.size'] = 3.
rcParams['xtick.major.size'] = 3.
rcParams['ytick.labelsize'] = 'medium'
rcParams['xtick.labelsize'] = 'medium'
rcParams['legend.fontsize'] = 'medium'
rcParams['legend.frameon'] = False
rcParams['axes.titlesize'] = 'medium'
rcParams['axes.linewidth'] = 3.0
rcParams['axes.formatter.limits'] = (-3, 3)
rcParams['font.size'] = 20.0
rcParams['font.weight'] = 'bold'
rcParams['lines.linewidth'] = 3.0
#read data

data = csv.reader(open('dt.sorted.csv','rb'))

X=np.array([])
Y=np.array([])
Z=np.array([])

for row in data:
    if len(row)==3:
        X=np.append(X,float(row[0]))
        Y=np.append(Y,float(row[1])/1e-3)
        Z=np.append(Z,float(row[2])/1e-3)

XX=X.reshape((10,10))
YY=Y.reshape((10,10))
ZZ=Z.reshape((10,10))

plot_data = np.flipud(np.transpose(ZZ))
extent = (XX[0][0],XX[-1][-1],YY[0][0],YY[-1][-1])
fig = plt.figure()
ax = fig.gca()
plt.subplots_adjust(bottom=.14)
cont = ax.contour(XX, YY, ZZ,[.1,1],
                 colors='w', # negative contours will be dashed by default
                 )
ax.text(.00125,1.5,"{:.1f}".format(.1), color = 'white', fontsize = 20, weight = 'bold')
ax.text(.0025,4.5,"{:.1f}".format(1), color = 'white', fontsize = 20, weight = 'bold')
image = ax.imshow(plot_data,extent=extent,aspect = 'auto')
ax.set_ylabel(r' (mrad)',fontsize = 20, weight = 'bold')
ax.set_title("Pulse Duration (fs)",fontsize = 20, weight = 'bold')
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markeredgewidth(3)
    line.set_markersize(5)
plt.colorbar(image,ticks = [Z.min(),1,2,3,4,5,Z.max()],format='%.2f')
plt.savefig("colormap-color.svg")


