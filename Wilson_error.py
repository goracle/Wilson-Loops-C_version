#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
import pylab as pl

def read_data(filename):
    data1_org = []
    data2_org = []
    #data_org_z = []

    data = open(filename, 'r')
    x=0
    for line in data:
        if len(line)>26:
            x+=1
            if(x%2==0):
                y=line[29:36]
                data2_org.append(float(y))
            else:
                y=line[28:36]
                data1_org.append(float(y))
            #data1_org.append(float(y))
        #data_org_z.append(float(z))
         #print type(line)
    return x,np.array(data1_org),np.array(data2_org)

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


x,axa,taxa = read_data('results.txt')
#print(y)
#print(z)
#plt.scatter(y,z)
#plt.show()

#data = remanage_data(data_org, therm_itr, sweep_itr)
#print data
print x

mean, std = mean_data(axa)
print "axa=(mean,std)"
print mean, std
mean, std = mean_data(taxa)
print "2axa=(mean,std)"
print mean, std
