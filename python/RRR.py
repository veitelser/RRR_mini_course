import numpy as np

def proj1(x):
    return [...]

def proj2(x):
    return [...]

def norm(x):
    return [...]

def RRR(x, beta=0.5, iter_max=100):
    err = []
    epsilon = 1.e-10
    for i in range(iter_max):
        p1 = proj1(x)
        p2 = proj2(2*p1 - x)
        diff = p2 - p1
        x += beta*diff
        err.append(norm(diff))
        if (err[-1] < epsilon):
            break