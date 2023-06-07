"""
This code read from file "results.txt", then plots the results in two adjacent subplots: 
subplot1: a graph with several curves for Time w.r.t number of locations
subplot2: a graph with several curves for Utility Loss w.r.t number of locations
Each curve corresponds to the Optimal method or "Kostas" method with certain Radius (R) 
"""
# !/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

l_colors = ['b', 'r', 'g', 'm']
l_mark = ['o', 's', '^', 'x']

labels = []
locations = []
Utility_L = dict()
c_Time = dict()
data_file = open('results.txt', 'rU')
counter = -1
for i, line in enumerate(data_file):
    if i == 0:
        continue
        
    line = line.strip()
    line = line.replace('\t\t', '\t')
    temp = line.split('\t')

    if temp[0] == 'flag':
        labels.append(temp[1])
        counter += 1
        continue
    else:
        if counter == 0:
            locations.append(temp[0])
        try:
            Utility_L[labels[counter]].append(float(temp[2]))
        except KeyError:
            Utility_L[labels[counter]] = [float(temp[2])]
        try:
            c_Time[labels[counter]].append(float(temp[3]))
        except KeyError:
            c_Time[labels[counter]] = [float(temp[3])]

#xValues = np.array(range(1, len(locations)+1))
CT= [None for x in labels]
UL = [None for x in labels]
fig, axs = plt.subplots(1, 2, figsize=(8.4, 4.3))

#Subplot 2: Time wrt nb. of loc.
for i in range(len(labels)):
    c_t = np.array(c_Time[labels[i]])
    CT[i], = axs[0].plot(locations, c_t, color=l_colors[i], ls='-', marker=l_mark[i], markeredgecolor=l_colors[i])
#axs[0].set_xticklabels(locations, rotation=45)
axs[0].set_ylabel('Time (minutes)', fontsize=16)
axs[0].set_xlabel('Number of Locations', fontsize=16)
leg = axs[0].legend(CT, labels, loc=0, fontsize=10.5)


#Subplot 2: Utility Loss wrt nb. of loc.
for i in range(len(labels)):
    ul = np.array(Utility_L[labels[i]])
    UL[i], = axs[1].plot(locations, ul, color=l_colors[i], ls='-', marker=l_mark[i], markeredgecolor=l_colors[i])
axs[1].set_ylabel('Utility Loss', fontsize=16)
axs[1].yaxis.set_label_position("right")
#axs[1].set_xticklabels(locations, rotation=45)
axs[1].set_xlabel('Number of Locations', fontsize=16)
leg = axs[1].legend(UL, labels, loc=0, fontsize=10.5)

fig.savefig('figure.pdf', format='pdf', bbox_inches='tight')
