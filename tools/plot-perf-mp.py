#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

avg_thr = {
    1:0,
    2:0,
    4:0,
    8:0,
    16:0,
    32:0
}
threads=[1,2,4,8,16,32]


with open("knl_mp.txt") as f:
    fline=f.readline()
    while(fline!=''):
        num_threads = int(fline.split()[1])
        f.readline()
        time_samples=[]
        fline=f.readline()
        while(fline!='End simulation\n'):
            time_samples.append( int( fline.split()[2] ) )
            fline=f.readline()
        avg_thr[num_threads]=np.array(time_samples).mean()
        fline=f.readline()
        print(fline)

avg = np.zeros(len(threads))
for th,idx in zip(threads,range(len(threads))):
    avg[idx] = avg_thr[th]

plt.figure(figsize=(20,10))
plt.subplot(1,2,1)
plt.bar(np.log2(threads),avg/1000)
plt.ylabel("Time (ms)")
plt.xlabel("$2^x$ threads")
plt.title("Average Xeon KNL threads execution time \n for a 10k infoli cell population")
plt.subplot(1,2,2)
plt.plot(threads,np.ones(len(avg))*avg[0]/avg)
plt.xlabel("Number of threads")
plt.ylabel("Speed-up")
plt.title("Speed-up performance")
plt.show()
