import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op

#15.896044  0.160872  0.056402  13.100257  4.520306

#choke = 15.847315491585901
choke = 16

a = 15.896044
b = 0.160872
theta = 0.056402#0.353522

t = np.linspace(0, np.pi, 250)

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
            a*np.cos(t)*np.cos(theta)-b*np.sin(t)*np.sin(theta), # xr true cond
        ])

def ellipse_pos_y(t):
    return np.choose(
        t < t_r,  # if condition
        [
            a*np.cos(t)*np.sin(theta)+b*np.sin(t)*np.cos(theta), # false cond
            a*np.cos(t)*np.sin(theta)+b*np.sin(t)*np.cos(theta),#c_tc*np.tan(t-d_theta),  # true cond
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

def ellipse_derivative(m): # x^2 + x + c = 0 -> 2x + 1 = 0 -> x = -[1]/[0]
    coeff = [
        2.0*np.sin(theta)**2/a**2 + np.cos(theta)**2/b**2,
        2.0*m*np.sin(theta)*np.cos(theta)*(1./a**2 - 1./b**2)
    ]

    return -coeff[1]/coeff[0]    

def ellipse_pos_y_eval(m):
    # Transform the ellipse eq to a quadratic eq (for x)
    coeff = [
        np.sin(theta)**2/a**2 + np.cos(theta)**2/b**2,
        2.0*m*np.sin(theta)*np.cos(theta)*(1./a**2 - 1./b**2),
        m**2 * (np.cos(theta)**2/a**2 + np.sin(theta)**2/b**2) - 1
    ]

    return -0.5*(coeff[1] - np.sqrt(coeff[1]**2 - 4.*coeff[0]*coeff[2]))/coeff[0]

def ellipse_pos_t(m):
    def min_f(t):
        return (ellipse_pos(t)[0] - m)**2 

    res = op.minimize_scalar(min_f, bounds=(t_r, np.pi/2.))
    return res.x
#(np.choose(
#        t < t_r,  # if condition
#        [
#            (a*np.cos(t)*np.cos(theta)-b*np.sin(t)*np.sin(theta), a*np.cos(t)*np.sin(theta)+b*np.sin(t)*np.cos(theta)), # false cond
#            (0,0),  # true cond
#        ]))

def ellipse_dy_dx_t(t):
    return (b*np.cos(theta)*np.cos(t) - a*np.sin(theta)*np.sin(t))/-(a*np.cos(theta)*np.sin(t) + b*np.sin(theta)*np.cos(t))

def ellipse_pos_from_d(d):
    def min_f(t):
        return (ellipse_dy_dx_t(t) - d)**2 

    res = op.minimize_scalar(min_f, bounds=(t_r, t_t))
    pos = ellipse_pos(res.x)
    return pos

x_r, y_r = ellipse_pos(t_r)
x_t, y_t = ellipse_pos(t_t)
c_tc = -y_r/np.tan(t_r+d_theta)

# Rotation transformed coords
x, y = ellipse_pos(t)

dx = -(a*np.cos(theta)*np.sin(t) + b*np.sin(theta)*np.cos(t))
dy = b*np.cos(theta)*np.cos(t) - a*np.sin(theta)*np.sin(t)

#pi_arr = np.linspace(15, -1, 1000)
#m_arr = speedline_massflow(pi_arr)

# y = 0

# Ok let's define a bezier down to 0.. (any may do this in the future for the whole thing with a cubic bezier)
# Meaning we have x1, x2 = x_choke, x3 = point on ellipse
# and             y1 = 0, y2 = func of ellipse, y3 = point on ellipse

data = np.array([[14.7217,0.891373],
[14.8267,	0.894628],
[15.0766,	0.901138],
[15.3035,	0.906034],
[15.4899,	0.909777],
[15.5652,	0.911095],
[15.6127,	0.911253],
[15.6694,	0.910936],
[15.7208,	0.909775],
[15.7681,	0.907296],
[15.8075,	0.904089],
[15.8512,	0.899604],
[15.8882,	0.894653]])

choke = data[-1, 0]*1.01

dy_over_dx = -0.25
P1 = np.array(
    [choke, 0.0]
)
P3 = ellipse_pos_from_d(-0.25)
P2 = np.array(
    [choke, P3[1] + dy_over_dx * (P1[0] - P3[0])]
)

def bezier(t):
    px = P2[0]+(1-t)**2*(P1[0]-P2[0])+t**2*(P3[0]-P2[0])
    py = P2[1]+(1-t)**2*(P1[1]-P2[1])+t**2*(P3[1]-P2[1])

    return (px, py)

t_bez = np.linspace(0, 1, 100)
p_bez = bezier(t_bez)

fig, ax = plt.subplots()

ax.plot(x, y)
#ax.scatter(x_r, y_r)
#ax.scatter(x_t, y_t)
#ax.plot(m_arr, pi_arr)

ax.plot(np.array([P1, P2])[:,0], np.array([P1, P2])[:,1])
ax.plot(np.array([P3, P2])[:,0], np.array([P3, P2])[:,1])
ax.plot(p_bez[0], p_bez[1])
ax.plot(data[:,0], data[:,1])
#ax.axvline(choke)

plt.show()