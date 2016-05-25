#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
import pylab as pl

def read_data(filename):
    data_org = []
    #data_org_z = []

    data = open(filename, 'r')
    x=0
    for line in data:
        if len(line)>26:
            x+=1
            y=line[28:36]
            data_org.append(float(y))
        #data_org_z.append(float(z))
         #print type(line)
    return np.array(data_org),x

def remanage_data(data_org, therm_itr, sweep_itr):
    data = []

    for i in range(therm_itr, len(data_org), sweep_itr):
        data.append(data_org[i])
    return np.array(data)

def mean_data(data):
    std = np.std(data, ddof = 1)
    return np.average(data), std*(len(data) ** (-1.0/2.0))

therm_itr = 1500;
sweep_itr = 50;


data,x = read_data('results.txt')
#print(y)
#print(z)
#plt.scatter(y,z)
#plt.show()

#data = remanage_data(data_org, therm_itr, sweep_itr)
#print data
print x

mean, std = mean_data(data)
print mean, std
