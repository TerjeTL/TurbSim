import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

cp_a = 1005.
cp_g = 1148.
n_m = 0.98

p_01 = 101325.0 #Pa
T_01 = 288.0 # K
gamma = 1.4

N_corr_design = 8000.0 # rpm
m_1_corr_design = 50.0 # kg/s

N_corr = N_corr_design*0.95
m_1_corr = m_1_corr_design*0.77

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
a = 1.5
b = 5.
k = 0.03

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
D = 1.0

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

pi_comp = pressure_ratio(m_1_corr/m_1_corr_design, N_corr/N_corr_design)
n_c = efficiency(pi_comp, m_1_corr/m_1_corr_design)
dT012 = T_01/n_c*(pi_comp**((gamma-1)/gamma) - 1)
T02 = T_01 + dT012
print("Inputs:\np01:", p_01, "\nT01:", T_01, "\nN:", N_corr, "\nm:", m_1_corr)
print("---------------------------")
print("p2/p1:", pi_comp)
print("nc:", n_c)
print("T02:", T02, "\t( delta:", dT012,")")

pi_turb = 2.0 #guess
m3_c = 25.0 #const
n_t = 0.87 #const

def T03_over_T01_flow_comp(pi_turb):
    return (m3_c * pi_comp / (m_1_corr * pi_turb))**2.0

#N_c_over_sq_T03 = N_corr/T03_over_T01_flow_comp

def T03_over_T01_mech_comp(pi_turb):
    dT034_over_T03 = n_t*(1 - (1/pi_turb)**((gamma-1)/gamma))
    return (dT012 * cp_a) / (dT034_over_T03 * T_01 * cp_g * n_m)

def min_f(pi_turb):
    return T03_over_T01_flow_comp(pi_turb) - T03_over_T01_mech_comp(pi_turb)

pi_turb = optimize.newton(min_f, pi_turb)

dT034_over_T03 = n_t*(1 - (1/pi_turb)**((gamma-1)/gamma))
T03_over_T01 = T03_over_T01_flow_comp(pi_turb)
T03 = T03_over_T01*T_01

print("---------------------------")
print("p3/p4", pi_turb)
print("nt:", n_t)
print("T03:", T03)
print("m3:", m3_c)

T04_over_T03 = 1-dT034_over_T03
m4_c = m3_c * pi_turb * np.sqrt(T04_over_T03)
T04 = T04_over_T03*T03

pi_comb = 0.95

p04_over_pa = pi_comp * pi_comp / pi_turb
# table would go here
m4_c_characteristic = 40.0

print("---------------------------")
print("p4/pa", p04_over_pa)
print("T04:", T04)
print("m4:", m4_c, "\t(from charact.:)", m4_c_characteristic)

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
    ax.plot(m, pi_it,'-k')

ax.plot(m[::10], m[::10]**a*(pi_d-1)+1, '-k')

ax.plot(N_corr/N_corr_design)

plt.show()