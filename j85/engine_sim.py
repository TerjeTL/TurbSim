import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# Gas Constants
cp_a = 1005.
cp_g = 1148.
R = 0.287
pa = 1.01325 # atm
T01 = 288.0 # K
gamma = 1.4
Ma = 0.0

# Design parameters
N_corr_design = 8000.0 # rpm
m1_corr_design = 200.0 # kg/s

#####################
# Engine State Vars #
#####################
class EngineState:
    def __init__(self, N1_corr, m1_corr, m3_corr, nozzle_area) -> None:
        self.n_mech = 0.98

        self.n_inlet = 0.9
        self.N1_corr = N1_corr
        self.p01 = 0.0
        self.m1_corr = m1_corr
        
        self.pi_comp = 0.0
        self.n_comp = 0.0
        self.p02 = 0.0
        self.m2_corr = 0.0
        self.T02 = 0.0

        self.pi_comb = 0.95
        self.p03 = 0.0
        self.m3_corr = m3_corr #const characteristic
        self.T03 = 0.0

        self.pi_turb = 2.0 # initial guess
        self.n_turb = 0.87
        self.p04 = 0.0
        self.m4_corr = 0.0
        self.T04 = 0.0

        self.nozzle_exit_area = nozzle_area
        self.n_nozzle = 0.9
        self.p04_over_p05 = 0.0
        self.p05_over_pa = 0.0

        self.F_over_pa = 0.0

    def disp(self):
        print("Inputs")
        print("-----------------")
        print("pa:", pa)
        print("N_corr:", self.N1_corr)
        
        print("\nStage 1")
        print("-----------------")
        print("p01:", self.p01)
        print("m1_corr:", self.m1_corr)

        print("\nStage 2")
        print("------------------")
        print("p02/p01:", self.pi_comp)
        print("comp. eff:", self.n_comp)
        print("m2_corr:", self.m2_corr)
        print("p02:", self.p02)
        print("T02:", self.T02)

        print("\nStage 3")
        print("------------------")
        print("p03/p02:", self.pi_comb)
        print("m3_corr:", self.m3_corr)
        print("p03:", self.p03)
        print("T03:", self.T03)

        print("\nStage 4")
        print("------------------")
        print("p03/p04:", self.pi_turb)
        print("turb. eff:", self.n_turb)
        print("m4_corr:", self.m4_corr)
        print("p04:", self.p04)
        print("T04:", self.T04)

        print("\nStage 5")
        print("------------------")
        print("Nozzle Area:", self.nozzle_exit_area)
        print("p04/p05:", self.p04_over_p05)
        print("p05/pa:", self.p05_over_pa)
        print("nozzle. eff:", self.n_nozzle)

        print("\nOutput")
        print("------------------")
        print("F/pa:", self.F_over_pa)
        print("Mass flow balance check:", self.m1_corr*self.p01/np.sqrt(T01), self.m2_corr*self.p02/np.sqrt(engine.T02),
              self.m3_corr*self.p03/np.sqrt(engine.T03), self.m4_corr*self.p04/np.sqrt(engine.T04))
        print("\n")

#engine = EngineState(
#    7600.0,
#    160.0,
#    20.0,
#    1.5
#)

# for the scan
#engine = EngineState(
#    8000.0,
#    200.0,
#    20.0,
#    1.3
#)

engine = EngineState(
    7360.0,
    136.0,
    20.0,
    1.5
)

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
c = 3 # width spread
d = 4 # length spread
C = 10.0
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

##########
# Nozzle #
##########
def corr_mass_flow(pi_noz):
    eff_prod = engine.n_nozzle*(1 - (1/pi_noz)**((gamma - 1)/gamma))
    c_5sq_over_T04 = 2.0*cp_g*eff_prod
    T5_over_T04 = 1 - eff_prod

    return np.sqrt(c_5sq_over_T04) * engine.nozzle_exit_area/(R * pi_noz * T5_over_T04)

critical_pr = 1/(1 - 1/engine.n_nozzle*((gamma-1)/(gamma + 1)))**(gamma/(gamma-1))
choked_flow = corr_mass_flow(critical_pr)

def nozzle_mass_flow(pi_noz):
    unchoked = corr_mass_flow(pi_noz)
    if pi_noz < critical_pr:
        if unchoked < choked_flow:
            return unchoked
    return choked_flow

#######################
# Steady-State Solver #
#######################
def iteration(m_corr):
    m_corr = m_corr[0]
    p01_over_pa = (1 + engine.n_inlet*(gamma - 1)/2*Ma**2)**(gamma/(gamma-1))
    engine.p01 = p01_over_pa*pa

    engine.pi_comp = pressure_ratio(m_corr/m1_corr_design, engine.N1_corr/N_corr_design)
    engine.n_comp = efficiency(engine.pi_comp, m_corr/m1_corr_design)
    engine.p02 = engine.p01 * engine.pi_comp
    dT012 = T01/engine.n_comp*(engine.pi_comp**((gamma-1)/gamma) - 1)
    engine.T02 = T01 + dT012
    print("Inputs:\np01:", engine.p01, "\nT01:", T01, "\nN:", engine.N1_corr, "\nm:", m_corr)
    print("---------------------------")
    print("p2/p1:", engine.pi_comp)
    print("nc:", engine.n_comp)
    print("T02:", engine.T02, "\t( delta:", dT012,")")

    T03_over_T01 = (engine.m3_corr * engine.pi_comp * engine.pi_comb / m_corr)**2.0

    #N_c_over_sq_T03 = N_corr/T03_over_T01_flow_comp

    def T03_over_T01_mech_comp(pi_turb):
        dT034_over_T03 = engine.n_turb*(1 - (1/pi_turb)**((gamma-1)/gamma))
        return (dT012 * cp_a) / (dT034_over_T03 * T01 * cp_g * engine.n_mech)

    def min_f(pi_turb):
        mech = T03_over_T01_mech_comp(pi_turb)
        return T03_over_T01 - mech

    engine.pi_turb = optimize.newton(min_f, engine.pi_turb)
    engine.p03 = engine.p02*engine.pi_comb
    dT034_over_T03 = engine.n_turb*(1 - (1/engine.pi_turb)**((gamma-1)/gamma))
    engine.T03 = T03_over_T01*T01

    print("---------------------------")
    print("p3/p4", engine.pi_turb)
    print("nt:", engine.n_turb)
    print("T03:", engine.T03)
    print("m3:", engine.m3_corr)

    T04_over_T03 = 1-dT034_over_T03
    engine.m4_corr = engine.m3_corr * engine.pi_turb * np.sqrt(T04_over_T03)
    engine.p04 = engine.p03/engine.pi_turb
    engine.T04 = T04_over_T03*engine.T03

    p04_over_pa = engine.pi_comp * engine.pi_comb / engine.pi_turb
    m4_c_characteristic = nozzle_mass_flow(p04_over_pa)

    print("---------------------------")
    print("p4/pa", p04_over_pa)
    print("T04:", engine.T04)
    print("m4:", engine.m4_corr, "\t(from charact.:)", m4_c_characteristic)

    c5_over_sqrtT04 = 0.0
    p5_over_pa = 1.0
    if p04_over_pa < choked_flow:
        eff_prod = engine.n_nozzle*(1 - (1/p04_over_pa)**((gamma - 1)/gamma))
        c5_over_sqrtT04 = np.sqrt(2.0*cp_g*eff_prod)
    else:
        p5_over_pa = 1/critical_pr * p04_over_pa
        c5_over_sqrtT04 = np.sqrt(2*gamma*R/(gamma + 1))
    Ca_over_sqrtT01 = Ma*np.sqrt(gamma*R)/np.sqrt(1+(gamma-1)/2*Ma**2)

    engine.p04_over_p05 = p04_over_pa
    engine.p05_over_pa = p5_over_pa

    engine.F_over_pa = m_corr*p01_over_pa*(c5_over_sqrtT04*np.sqrt(T04_over_T03*T03_over_T01)-Ca_over_sqrtT01)+(p5_over_pa-1)*engine.nozzle_exit_area

    print("F/pa:", engine.F_over_pa, "\n\n")
    
    return (engine.m4_corr - m4_c_characteristic)**2.

engine.m1_corr = optimize.minimize(iteration, engine.m1_corr).x[0]
engine.disp()

#running_line_N = np.linspace(8000.0, 6000.0, 41)
#running_line_m_corr = []
#running_line_pi = []

#for N in running_line_N:
#    engine.N1_corr = N
#    try:
#        engine.m1_corr = optimize.newton(iteration, engine.m1_corr)
#        engine.disp()
#        running_line_m_corr.append(engine.m1_corr/m1_corr_design)
#        running_line_pi.append(engine.pi_comp)
#    except:
#        pass

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
ax.set_title('Compressor Map')
ax.set_xlabel("Corrected Mass Flow (normalized)")
ax.set_ylabel("p2/p1")

pi = []
for N_it in N:
    pi.append(pressure_ratio(m, N_it))

for pi_it in pi:
    ax.plot(m, pi_it, '-k')

#ax.plot(m[::10], m[::10]**a*(pi_d-1)+1, '-k')

#ax.plot(running_line_m_corr, running_line_pi, '-ro', markersize=4)
ax.plot(engine.m1_corr/m1_corr_design, engine.pi_comp, '-ro')

plt.show()

########################
# Transient Simulation #
########################
# State vector takes the form 
# X = [N, p2, p4]

#proceeds as follows:
# 1. from atm, intake, Ncorr -> P01, T01 -> P02/P01 -> n_comp, m1_corr
#    -> T2, dT12, comp. torq
#
# 2. dT23 = f(ff, m1_corr, T3), m_corr -> T3, P02/P03 -> P03
# 
# 3. P3/P4 -> trq, m3_corr, T04
#
# 4. P4, T4, A -> m_corr, F
#
# 6. Derivatives:
# mass accumulation from balance at nodes:
# m_dot = m_comp + m_ff - m_turb
# P = mRT/V -> P_dot = R(m_dot T + m T_dot)/V - we assume (m T_dot) = 0
# P_dot = RT/V * m_dot
#
# 7. Integrate N, P2, P4

V_ic = 1. #default intercomponent volume
inertia = 0.05
dT23 = engine.T03 - engine.T02 # Const for now - this is the combustion heat

def pressure_ratio_inv(m_curr, pi, N):
    def optimize_func(m_var):
        return pressure_ratio(m_var, N) - pi
    return optimize.newton(optimize_func, m_curr)

def transient_step(dt, dT23):
    #State vars
    N1_corr = engine.N1_corr
    p02 = engine.p02
    p04 = engine.p04

    #Inlet
    p01_over_pa = (1 + engine.n_inlet*(gamma - 1)/2*Ma**2)**(gamma/(gamma-1))
    p01 = p01_over_pa*pa

    #Compressor
    pi_comp = p02/p01
    m1_corr = pressure_ratio_inv(engine.m1_corr/m1_corr_design, pi_comp, N1_corr/N_corr_design)*m1_corr_design
    n_comp = efficiency(pi_comp, m1_corr/m1_corr_design)
    dT012 = T01/n_comp*(pi_comp**((gamma-1)/gamma) - 1)
    T02 = T01 + dT012
    #calculate torque

    #calculate dT23 from ff (simulate comb. chamber)
    T03 = T02 + dT23
    p03 = p02 * engine.pi_comb
    pi_turb = p03/p04

    dT034_over_T03 = engine.n_turb*(1 - (1/pi_turb)**((gamma-1)/gamma))
    T04_over_T03 = 1-dT034_over_T03
    
    T04 = T04_over_T03*T03

    #Nozzle
    p04_over_pa = pi_comp * engine.pi_comb / pi_turb
    m4_corr = nozzle_mass_flow(p04_over_pa)

    c5_over_sqrtT04 = 0.0
    p5_over_pa = 1.0
    if p04_over_pa < choked_flow:
        eff_prod = engine.n_nozzle*(1 - (1/p04_over_pa)**((gamma - 1)/gamma))
        c5_over_sqrtT04 = np.sqrt(2.0*cp_g*eff_prod)
    else:
        p5_over_pa = 1/critical_pr * p04_over_pa
        c5_over_sqrtT04 = np.sqrt(2*gamma*R/(gamma + 1))
    Ca_over_sqrtT01 = Ma*np.sqrt(gamma*R)/np.sqrt(1+(gamma-1)/2*Ma**2)

    engine.p04_over_p05 = p04_over_pa
    engine.p05_over_pa = p5_over_pa

    T04_over_T03 = T04/T03
    T03_over_T01 = T03/T01

    engine.F_over_pa = m1_corr*p01_over_pa*(c5_over_sqrtT04*np.sqrt(T04_over_T03*T03_over_T01)-Ca_over_sqrtT01)+(p5_over_pa-1)*engine.nozzle_exit_area

    # Calculate derivatives and perform integration
    dm2 = m1_corr*p01/np.sqrt(T01) - engine.m3_corr*p03/np.sqrt(T03)
    dP02 = R/V_ic*(dm2*T02)

    dm4 = engine.m3_corr*p03/np.sqrt(T03) - m4_corr*p04/np.sqrt(T04)
    dP04 = R/V_ic*(dm4*T04)

    p02 += dP02*dt
    p04 += dP04*dt

    curr_N = N1_corr*np.sqrt(T01)*2.*np.pi/60.
    trq_comp = m1_corr*p01/np.sqrt(T01)*cp_a*dT012/curr_N
    trq_turb = engine.m3_corr*p03/np.sqrt(T03)*cp_g*dT034_over_T03*engine.T03/curr_N

    dN = (engine.n_mech*trq_turb - trq_comp)/inertia
    curr_N += dN*dt

    #Populate fields
    engine.N1_corr = curr_N*60.0/(np.sqrt(T01)*2.0*np.pi)
    engine.p01 = p01
    engine.p02 = p02
    engine.p03 = p03
    engine.p04 = p04
    engine.pi_comp = p02/p01
    engine.n_comp = n_comp
    engine.pi_turb = p03/p04
    engine.T02 = T02
    engine.T03 = T03
    engine.T04 = T04
    engine.m1_corr = m1_corr
    engine.m4_corr = m4_corr

##################
# Simulation Run #
##################

dt = 0.001
t_stop = 20.0
time = [0.0]
ng_transient = [engine.N1_corr]

# To plot the running line
transient_m_corr_norm = []
transient_pi = []

def throttle_func(t):
    if t<1:
        return dT23
    elif t>10:
        return dT23*1.2
    else:
        return dT23*0.2

while time[-1] < t_stop:
    t = time[-1]
    try:
        transient_step(dt, throttle_func(t))
        #print("Time:", time[-1]+dt)
        #engine.disp()
        time.append(time[-1] + dt)
        ng_transient.append(engine.N1_corr)

        transient_m_corr_norm.append(engine.m1_corr/m1_corr_design)
        transient_pi.append(engine.pi_comp)
    except:
        break

plt.plot(time, ng_transient)
plt.ylim([0.0, 8000])
plt.xlabel("Time [s]")
plt.ylabel("Ncorr")
plt.show()

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z, eff_contrs)
ax.clabel(CS, inline=True, fontsize=10)
ax.set_title('Compressor Map')
ax.set_xlabel("Corrected Mass Flow (normalized)")
ax.set_ylabel("p2/p1")

pi = []
for N_it in N:
    pi.append(pressure_ratio(m, N_it))

for pi_it in pi:
    ax.plot(m, pi_it, '-k')

ax.plot(transient_m_corr_norm, transient_pi, '-r', markersize=4)

plt.show()