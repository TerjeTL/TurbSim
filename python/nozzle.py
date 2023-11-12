import numpy as np
import matplotlib.pyplot as plt

# Just testing out equations for the nozzle

#m_c = 
gamma = 1.3
n_j = 0.85
A = 1
R = 0.287
cp = 1005

def corr_mass_flow(pi_noz):
    eff_prod = n_j*(1 - (1/pi_noz)**((gamma - 1)/gamma))
    c_5sq_over_T04 = 2.0*cp*eff_prod
    T5_over_T04 = 1 - eff_prod

    return np.sqrt(c_5sq_over_T04) * A/(R * pi_noz * T5_over_T04)

critical_pr = 1/(1 - 1/n_j*((gamma-1)/(gamma + 1)))**(gamma/(gamma-1))
choked_flow = corr_mass_flow(critical_pr)

def nozzle_mass_flow(pi_noz):
    unchoked = corr_mass_flow(pi_noz)
    if pi_noz < critical_pr:
        if unchoked < choked_flow:
            return unchoked
    return choked_flow

vfunc = np.vectorize(nozzle_mass_flow)

pr = np.linspace(1.0, critical_pr, 251)

plt.plot(pr, vfunc(pr))

plt.show()