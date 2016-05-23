#!/usr/bin/env python

import numpy as np

def read_data(filename):
    data_org = []

    data = open(filename, 'r')
    x=0
    for line in data:
        x+=1
        line=line[12:]
        data_org.append(float(line))
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


data,x = read_data('wilsonloop.txt')

#data = remanage_data(data_org, therm_itr, sweep_itr)
print data
print x

mean, std = mean_data(data)
print mean, std
