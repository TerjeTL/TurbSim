import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# A new component map generation method for gas turbine adaptation performance simulation

R = 1

m_c0 = 1.
pi_c0 = 1.

Spi_m = 0.0
Spi_pi = 0.0
Spi_theta = 0.0
Spi_phi = np.pi/4

def a(Spi_m):
    return R*(1+Spi_m)
def b(Spi_pi):
    return R*(1+Spi_pi)

def m_c(phi_theta, Spi_m, Spi_pi, Spi_theta, Spi_phi):
    return (a(Spi_m)*np.cos(Spi_phi)*np.cos(phi_theta) - b(Spi_pi)*np.sin(Spi_phi)*np.sin(phi_theta)
        - a(Spi_m)*np.cos(Spi_theta)*np.cos(Spi_phi) + b(Spi_pi)*np.sin(Spi_theta)*np.sin(Spi_phi)
        + m_c0)

def pi_c(phi_theta, Spi_m, Spi_pi, Spi_theta, Spi_phi):
    return (a(Spi_m)*np.cos(Spi_phi)*np.cos(phi_theta) + b(Spi_pi)*np.sin(Spi_phi)*np.sin(phi_theta)
            - a(Spi_m)*np.cos(Spi_theta)*np.sin(Spi_phi) + b(Spi_pi)*np.sin(Spi_theta)*np.cos(Spi_phi)
            + pi_c0)

def test(X, Spi_m, Spi_pi, Spi_theta, Spi_phi):
    x,y = X
    theta_p = np.arctan2(x, y)
    m = m_c(theta_p, Spi_m, Spi_pi, Spi_theta, Spi_phi)
    pi = pi_c(theta_p, Spi_m, Spi_pi, Spi_theta, Spi_phi)
    return m*pi

theta = np.linspace(0, 2*np.pi, 101)

x_p = np.array([0.4, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
y_p = np.array([1.11, 1.1, 1.08, 1.05, 1.0, 0.8, 0.4])
theta_p = np.arctan2(x_p, y_p) # This will be the independent param

z=x_p*y_p
p_opt, pcov = curve_fit(test, (x_p,y_p), z, p0=[Spi_m, Spi_pi, Spi_theta, Spi_phi], bounds=([-1, -1, -np.pi, -np.pi], [np.inf, np.inf, np.pi, np.pi]))
#p_opt, pcov = curve_fit(m_c, theta_p, x_p, p0=p_opt, bounds=([-1, -1, -np.pi, -np.pi], [np.inf, np.inf, np.pi, np.pi]), absolute_sigma=True)

x = m_c(theta, p_opt[0], p_opt[1], p_opt[2], p_opt[3])
y = pi_c(theta, p_opt[0], p_opt[1], p_opt[2], p_opt[3])

#x = m_c(theta, Spi_m, Spi_pi, Spi_theta, Spi_phi)
#y = pi_c(theta, Spi_m, Spi_pi, Spi_theta, Spi_phi)

plt.plot(x_p, y_p)
plt.plot(x, y)
plt.show()