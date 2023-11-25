import numpy as np
import matplotlib.pyplot as plt

# PLS(N)

# Define some speed curves

# m, pi
n05 = np.array([
    [0.4, 3.0],
    [0.46, 2.95],
    [0.5, 2.85],
    [0.555, 2.2],
    [0.56, 1.4]
])

n075 = np.array([
    [0.4, 3.0],
    [0.48, 2.95],
    [0.52, 2.75],
    [0.55, 2.2],
    [0.56, 1.4]
]) * 1.4
n075[:,1] = n075[:,1]*0.9


n1 = np.array([
    [0.45, 2.98],
    [0.55, 2.95],
    [0.6, 2.8],
    [0.61, 2.5],
    [0.62, 1.4]
])*1.6667

speed_lines = np.array([
    n05,
    n075,
    n1
])

plt.plot(speed_lines[:,:,0].T, speed_lines[:,:,1].T)
plt.show()

# Test dataset complete

# PLS-R
#---------
# q - dependant vars
# p - independent vars
# n - sample points (constituting dependant and independant vars)
# As:
# X = [x1,x2,...,xp] (n x p), Y = [y1,y2,...,yq] (n x q)
#
# Components t1, u1 from X and Y (linear combination of x1,x2,..,xp and same for Y). Two requirements for regression:
# i)  t1 and u1 should carry the most mutation information about their own datasets => t1, u1 variances are maximized
# ii) abs(correlation) of t1 and u1 is maximized
#
# after extracting t1 and u1, PLS-R conducts regression for X for t1 and Y for u1. If accuracy condition is reached
# calculation terminates; otherwise, a second round of component extraction is performed. Making use of residual info of X
# interpreted by t1 and Y by u1, and so on...
# 
# If m number of components t1,t2...,tm are finally extracted from X, PLS-R conducts the regression of yk for t1,t2,...,tm
# and then is further expressed by a regression eq. of the previous independant variables x1,x2,...,xp, k=1,2,...,q.
#  
#