import numpy as np
from scipy.ndimage import shift

def mean(time, heat, index=None):
    if index is None:
        d_time = time[:] - shift(time[:], 1, cval=0)
        Qav = np.sum(d_time[:] * heat[:])/np.sum(d_time[:])
        return Qav, np.sum(d_time)
    
    else:
        time = time[index]
        heat = heat[index]
        d_time = time[:] - shift(time[:], 1, cval=time[0])
        Qav = np.sum(d_time[:] * heat[:])/np.sum(d_time[:])
        return Qav, np.sum(d_time)

def peaks(input):
    res = []

    j = -1
    flag = False
    for i in range(len(input)):
        if input[i] and not flag:
            j += 1
            res.append([])
            flag = True
            res[j].append(i)
        elif input[i] and flag:
            res[j].append(i)
        elif not input[i] and flag:
            flag = False

    return res