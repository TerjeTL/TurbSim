import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op

a = 20.970729
b = 4.151917
theta = 0.6#0.353522

t = np.linspace(0, np.pi/2.0, 250)

d_theta = np.arctan(-a/b*np.tan(theta))

t = t+d_theta
t_r = np.arctan(-b/a*np.tan(theta))
t_t = np.arctan(b/(a*np.tan(theta)))
x_r, y_r = (0, 0)

# Choke defined as y = c_tc * tan(t)
c_tc = 0.0

def ellipse_pos_x(t):
    return np.choose(
        t < t_r,  # if condition
        [
            a*np.cos(t)*np.cos(theta)-b*np.sin(t)*np.sin(theta), # false cond
            x_r,  # true cond
        ])

def ellipse_pos_y(t):
    return np.choose(
        t < t_r,  # if condition
        [
            a*np.cos(t)*np.sin(theta)+b*np.sin(t)*np.cos(theta), # false cond
            c_tc*np.tan(t-d_theta),  # true cond
        ])

def ellipse_pos(t):
    return (ellipse_pos_x(t), ellipse_pos_y(t))


def speedline_massflow(pi):
    # Transform the ellipse eq to a quadratic eq (for x)
    coeff = [
        np.cos(theta)**2/a**2 + np.sin(theta)**2/b**2,
        2.0*pi*np.sin(theta)*np.cos(theta)*(1./a**2 - 1./b**2),
        pi**2 * (np.sin(theta)**2/a**2 + np.cos(theta)**2/b**2) - 1
    ]

    m_corr = np.choose( pi > y_r,  # if condition
        [
            x_r, # false cond
            np.choose( pi < y_t,  # if condition
            [
                x_t, # false cond (not good)
                -0.5*(coeff[1] - np.sqrt(coeff[1]**2 - 4.*coeff[0]*coeff[2]))/coeff[0],  # true cond
            ]),  # true cond
        ])

    return m_corr

#def speedline_pi(m):
#    def min_f(t):
#        return (ellipse_pos(t*np.pi/2.0+d_theta)[0] - m)**2 
#
#    res = op.minimize_scalar(min_f, bounds=(0., 1.))
#    print(res)
#    return (ellipse_pos(res.x*np.pi/2.0+d_theta)[1], res.x)
#(np.choose(
#        t < t_r,  # if condition
#        [
#            (a*np.cos(t)*np.cos(theta)-b*np.sin(t)*np.sin(theta), a*np.cos(t)*np.sin(theta)+b*np.sin(t)*np.cos(theta)), # false cond
#            (0,0),  # true cond
#        ]))

x_r, y_r = ellipse_pos(t_r)
x_t, y_t = ellipse_pos(t_t)
c_tc = -y_r/np.tan(t_r+d_theta)

# Rotation transformed coords
x, y = ellipse_pos(t)

dx = -(a*np.cos(theta)*np.sin(t) + b*np.sin(theta)*np.cos(t))
dy = b*np.cos(theta)*np.cos(t) - a*np.sin(theta)*np.sin(t)

pi_arr = np.linspace(15, -1, 1000)
m_arr = speedline_massflow(pi_arr)

# y = 0


fig, ax = plt.subplots()

ax.plot(x, y)
#ax.scatter(x_r, y_r)
#ax.scatter(x_t, y_t)
ax.plot(m_arr, pi_arr)

plt.show()