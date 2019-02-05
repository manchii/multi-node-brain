#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

time = np.zeros(60000)
vaxon_ref = np.zeros(60000)
vaxon_test = np.zeros(60000)

with open("../InferiorOlive_Output_ref.txt") as f:
    idx=0
    f.readline()
    for line in f:
        parsed = f.readline().split()
        time[idx]=float(parsed[1])
        vaxon_ref[idx]=float(parsed[3])
        idx+=1

with open("../InferiorOlive_Output.txt") as f:
    idx=0
    f.readline()
    for line in f:
        vaxon_test[idx]=float(f.readline().split()[4])
        idx+=1


#plt.plot(time,vaxon_ref)
plt.plot(time,vaxon_ref,linewidth=5)
plt.plot(time,vaxon_test)
plt.plot(time,vaxon_ref-vaxon_test)
plt.show()
