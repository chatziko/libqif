"""
This code plots in two adjacent subplots: 
subplot1: a graph with several curves for Time w.r.t number of locations
subplot2: a graph with several curves for Utility Loss w.r.t number of locations
Each curve corresponds to the Optimal method or "Kostas" method with certain Radius (R) 
"""

# !/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

locations = [64, 81, 100, 121, 144, 169, 196, 225]
xValues = np.arange(1, 9, 1)

c_timeD = np.array([2.2, 12.69, 55.73, 206.2, 809, 2324.6])
c_timeR2 = np.array([0.243, 0.6618, 1.825, 4.7725, 11.843, 31.605, 75.835, 187.02])
c_time_R2_2 = np.array([0.2467, 0.69, 1.875, 4.461, 11.844, 31.846, 88.137, 219.68])
c_timeR2_5 = np.array([0.521, 1.66, 4.28, 12.35, 36.56, 104.32, 266.5, 625.3])
c_timeR3 = np.array([0.56, 1.84, 4.92, 13.23, 40.03, 110.3, 287.1, 618.3])
c_timeR3_5 = np.array([0.79, 2.62, 7.85, 27.726, 95.6, 211.5, 516.6, 1135.6])
c_timeR4 = np.array([0.92, 3.27, 12.5,  48.4, 132.5, 301])

utilityD = np.array([2.56243, 2.78742, 2.99209, 3.17712, 3.34382, 3.49328])
utilityR2 = np.array([2.89654, 3.19917, 3.49993, 3.79664, 4.0787, 4.34521, 4.58603, 4.82057])
utility_R2_2 = np.array([2.80522, 3.10312, 3.36236, 3.61198, 3.8486, 4.07648, 4.2899, 4.48869])
utilityR2_5 = np.array([2.73575, 3.01626, 3.27404, 3.50123, 3.71121, 3.91027, 4.0971, 4.26994])
utilityR3 = np.array([2.66593, 2.92578, 3.16351, 3.37836, 3.57117, 3.74393, 3.90332, 4.05024])
utilityR3_5 = np.array([2.63476, 2.8833, 3.11085, 3.31791, 3.50393, 3.6716, 3.82279, 3.96098])
utilityR4 = np.array([2.6157, 2.85768, 3.07928, 3.28145, 3.46201, 3.62472])


fig, axs = plt.subplots(1, 2, figsize=(8.4, 4.3))
#plt.subplots_adjust(wspace=0.095)
#Subplot 1: Time wrt nb. of loc.
tRD, = axs[0].plot(xValues[:-2], c_timeD, 'kd-')
tR2, = axs[0].plot(xValues, c_timeR2, 'bs-', markeredgecolor='b')
tR2_2, = axs[0].plot(xValues, c_time_R2_2, color='#8B4513', ls='-', marker='^', markeredgecolor='#8B4513')
tR2_5, = axs[0].plot(xValues, c_timeR2_5, 'ro-', markeredgecolor='r')
tR3, = axs[0].plot(xValues, c_timeR3, 'g*-', markeredgecolor='g')
tR3_5, = axs[0].plot(xValues, c_timeR3_5, 'm+-', markeredgecolor='m')
tR4, = axs[0].plot(xValues[:-2], c_timeR4, 'cx-', markeredgecolor='c')
axs[0].set_xticklabels(locations, rotation=45)
axs[0].set_ylabel('Time (minutes)', fontsize=16)
axs[0].set_xlabel('Number of Locations', fontsize=16)
#axs[0].set_ylim([0, 350])
#axs[0].set_title('Title of plot 1')
leg = axs[0].legend([tRD, tR2, tR2_2, tR2_5, tR3, tR3_5, tR4], ['Opt.','R=2', 'R=2.2','R=2.5', 'R=3', 'R=3.5', 'R=4'], loc=2, fontsize=10.5)
#leg.get_frame().set_linewidth(0.0)
y_vec = [50, 100, 150, 200, 250, 300, 350]
for b in range(250, 2500, 250):
	axs[0].axhline(y=b, xmin=0, xmax=1, ls='dotted', color='black')

#Subplot 1: Utility Loss wrt nb. of loc.
UD,  = axs[1].plot(xValues[:-2], utilityD, 'kd-')
UR2,  = axs[1].plot(xValues, utilityR2, 'bs-', markeredgecolor='b')
UR2_2, = axs[1].plot(xValues, utility_R2_2, color='#8B4513', ls='-', marker='^', markeredgecolor='#8B4513')
UR2_5, = axs[1].plot(xValues, utilityR2_5, 'ro-', markeredgecolor='r')
UR3, = axs[1].plot(xValues, utilityR3, 'g*-', markeredgecolor='g')
UR3_5, = axs[1].plot(xValues, utilityR3_5, 'm+-',  markeredgecolor='m')
UR4, = axs[1].plot(xValues[:-2], utilityR4, 'cx-',  markeredgecolor='c')
axs[1].set_xticklabels(locations,  rotation=45)
#axs[1].yaxis.tick_right()
axs[1].set_ylabel('Utility Loss', fontsize=16)
axs[1].yaxis.set_label_position("right")
axs[1].set_xlabel('Number of Locations', fontsize=16)
axs[1].legend([UD, UR2, UR2_2, UR2_5, UR3, UR3_5, UR4], ['Opt.','R=2', 'R=2.2','R=2.5', 'R=3', 'R=3.5', 'R=4'], loc=2, fontsize=10.5)
#axs[1].set_title('Title of plot 2')
y_vec = [2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75] 
for b in y_vec:
	axs[1].axhline(y=b, xmin=0, xmax=1, ls='dotted', color='black')

plt.show()
fig.savefig('figure.pdf', format='pdf', bbox_inches='tight')
