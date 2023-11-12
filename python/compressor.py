import matplotlib.pyplot as plt
from scipy import optimize
import numpy as np

#########################
# Analytical compressor #
#########################

# Constants #
#-- Pressure Ratio --#
# pi_d = design point pressure ratio
# a = shape of spine
# b = position of speed lines
# k = knee sharpness
pi_d = 26.0 
a = 1.25
b = 5.
k = 0.08

#-- Efficiency --#
# n_0 = max efficiency
# m_0 = mass flow at n_0
# da = exponent delta (positive for compressors, negative for fans)
# c, d, C, D = controls decrease in efficiency as operating point moves away from m_0, p_0
n_0 = 0.887
m_0 = 0.8
da = 0.5
c = 3
d = 4
C = 15.0
D = 0.5

def m_s(N):
    return N**b
def p_s(N):
    return N**(a*b)
def p_norm(pi):
    return (pi - 1.)/(pi_d - 1.)

def pressure_ratio(m, N):
    return 1 + (pi_d-1)*(N**(a*b) + 2*N*k*np.log(1-(m-N**b)/k))
def efficiency(pi, m):
    return n_0*(1 - C*np.abs(p_norm(pi)/(m**(a+da-1))-m)**c - D*np.abs(m/m_0 - 1)**d)

m = np.linspace(0, 1.2, 10001)
N = np.linspace(0.05, 1., 11)**(1/4)
print(N)

delta = 0.01
x_m = np.arange(0.0, 1.5, delta)
y_pi = np.arange(0.0, pi_d*1.75, delta)
X, Y = np.meshgrid(x_m, y_pi)
Z = efficiency(Y, X)

eff_contrs = np.linspace(0, n_0, 15)**(1/8)
fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, eff_contrs)
ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('Simplest default with labels')

pi = []
for N_it in N:
    pi.append(pressure_ratio(m, N_it))

for pi_it in pi:
    ax.plot(m, pi_it)

ax.plot(m[::10], m[::10]**a*(pi_d-1)+1, '-k')

def pressure_ratio_inv(m_curr, pi, N):
    def optimize_func(m_var):
        return pressure_ratio(m_var, N) - pi
    return optimize.newton(optimize_func, m_curr)
reverse_m = pressure_ratio_inv(0.9, 26., 1.)
ax.plot(reverse_m, pressure_ratio(reverse_m, 1.), 'ro')

plt.show()