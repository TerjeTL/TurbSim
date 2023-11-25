import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
import pandas as pd

speedline_cpr = {}
speedline_n = {}

with open("./j85/j85_cpr.csv", "r") as file:
    reader = csv.reader(file, delimiter=",")
    
    header = ""
    m_to_cpr = []

    for i, line in enumerate(reader):
        # Process row
        if len(line) == 0:
            continue
        elif line[0] == "m": # New speedline
            if header != "":
                speedline_cpr[header] = np.array(m_to_cpr)
                m_to_cpr = []
            header = line[1]
        else:
            m_to_cpr.append([float(line[0]), float(line[1])])

    speedline_cpr[header] = np.array(m_to_cpr)

with open("./j85/j85_eff.csv", "r") as file:
    reader = csv.reader(file, delimiter=",")
    
    header = ""
    m_to_n = []

    for i, line in enumerate(reader):
        # Process row
        if len(line) == 0:
            continue
        elif line[0] == "m": # New speedline
            if header != "":
                speedline_n[header] = np.array(m_to_n)
                m_to_n = []
            header = line[1]
        else:
            m_to_n.append([float(line[0]), float(line[1])])

    speedline_n[header] = np.array(m_to_n)


fig = plt.figure()
gs = fig.add_gridspec(2,2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])
fig.suptitle('Compressor Map')

N_c = []
a_pi = []
b_pi = []
theta_pi = []
surge_line = []
choke_line = []

a_n = []
b_n = []
theta_n = []

def ellipse_fit(x, y):
    # let x = mass flow, y = cpr
    X = x[None].T
    Y = y[None].T

    # Now lets fit some ellipses using lsq ||Ax - b||^2
    # Following standard form we have no translation and so we define our eq. to be
    # Ax^2 + Bxy + Cy^2 = 1 (=> F = -1)
    A = np.hstack([X**2, X*Y, Y**2]) # A, B, C in that order
    b = np.ones_like(X)
    x = np.linalg.lstsq(A, b)[0].squeeze()
    
    # Print the equation of the ellipse in standard form
    print('Ellipse (standard form): {0:.3}x^2 + {1:.3}xy+{2:.3}y^2 = 1'.format(x[0],x[1],x[2]))

    # Now convert to canonical form
    M0 = np.matrix([
        [-1,    0.0,        0.0     ],
        [0.0,   x[0],       x[1]/2.0],
        [0.0,   x[1]/2.0,   x[2]    ]
    ])
    M = np.matrix([
        [x[0],       x[1]/2.0],
        [x[1]/2.0,   x[2]    ]
    ])

    lambdas = np.sort(np.linalg.eigvals(M))
    det_M0 = np.linalg.det(M0)
    det_M = np.linalg.det(M)

    a = np.sqrt(-det_M0/(det_M*lambdas[0]))
    b = np.sqrt(-det_M0/(det_M*lambdas[1]))
    theta = np.arctan(x[1]/(x[0]-x[2]))/2.0
    print('Ellipse (canonical form): x^2/{0:.3} + y^2/{0:.3} = 1 with rotation theta = {0:.3}'.format(a,b,theta))

    return (a, b, theta)

def ellipse_eval(t, a, b, theta):
    # Coords before rotation
    m = a*np.cos(t)
    y = b*np.sin(t)
    
    # Rotation transformed coords
    m_r = m*np.cos(theta)-y*np.sin(theta)
    y_r = m*np.sin(theta)+y*np.cos(theta)
    return (m_r, y_r)

def cpr_plot(ax, t_range, a, b, theta, surge=[0, 0], choke=[0, 0]):
        t = np.linspace(t_range[0], t_range[1], 100)
        m_r, pi_r = ellipse_eval(t, a, b, theta)

        # Insert choke, surge condition
        for i in range(len(m_r)):
            if m_r[i] < surge[0]:
                pi_r[i] = np.nan
            elif pi_r[i] < choke[1]:
                m_r[i] = choke[0]

        ax.plot(m_r, pi_r, '-.k')

def eff_plot(ax, t_range, a, b, theta):
        t = np.linspace(t_range[0], t_range[1], 100)
        m_r, pi_r = ellipse_eval(t, a, b, theta)

        ax.plot(m_r, pi_r, '-.k')

for spd, cpr_arr in speedline_cpr.items():
    n_arr = speedline_n[spd]

    # Plot the data
    ax1.plot(cpr_arr[:,0], cpr_arr[:,1])
    ax1.set(xlabel='Mass flow', ylabel='CPR')

    ax2.plot(n_arr[:,0], n_arr[:,1])
    ax2.set(xlabel='Mass flow', ylabel='eff.')

    # Fit the pressure map
    a,b,theta = ellipse_fit(cpr_arr[:,0], cpr_arr[:,1])

    # Now we can determine the bounding box of the fitted ellipse (in this case, positive bounds)
    # ref: https://math.stackexchange.com/questions/91132/how-to-get-the-limits-of-rotated-ellipse
    m_r = np.sqrt(a**2*np.cos(theta)**2 + b**2*np.sin(theta)**2)
    pi_r = -(b**2-a**2)*np.sin(2.*theta)/(2.*np.sqrt(a**2*np.cos(theta)**2 + b**2*np.sin(theta)**2))

    pi_t = np.sqrt(a**2*np.sin(theta)**2 + b**2*np.cos(theta)**2)
    m_t = -(b**2-a**2)*np.sin(2.*theta)/(2.*np.sqrt(a**2*np.sin(theta)**2 + b**2*np.cos(theta)**2))

    # Save the data
    N_c.append(float(spd))
    a_pi.append(a)
    b_pi.append(b)
    theta_pi.append(theta)
    surge_line.append(m_t)
    choke_line.append(pi_r)

    print(spd, "has choke mass flow:", m_r)

    # Fit the efficiency map
    a,b,theta = ellipse_fit(n_arr[:,0], n_arr[:,1])

    # Save the data
    a_n.append(a)
    b_n.append(b)
    theta_n.append(theta)

    # Plot the ellipses
    cpr_plot(ax1, [-np.pi/8., np.pi/8.], a_pi[-1], b_pi[-1], theta_pi[-1], [cpr_arr[0,0], 0], [m_r, pi_r])
    eff_plot(ax2, [-np.pi/32., np.pi/8.], a_n[-1], b_n[-1], theta_n[-1])

N_c = np.array(N_c)
a_pi = np.array(a_pi)
b_pi = np.array(b_pi)
theta_pi = np.array(theta_pi)
a_n = np.array(a_n)
b_n = np.array(b_n)
theta_n = np.array(theta_n)

def absmaxND(a, axis=None):
    amax = a.max(axis)
    amin = a.min(axis)
    return np.where(-amin > amax, amin, amax)
def normalized_arr(a):
    return a/absmaxND(a)

ax3.plot(N_c, np.array([normalized_arr(a_pi), normalized_arr(b_pi), normalized_arr(theta_pi), normalized_arr(a_n), normalized_arr(b_n), normalized_arr(theta_n)]).T)
ax3.set(xlabel='Nc')

plt.show()

# Test interp and extr
N_arr = np.linspace(N_c[0], N_c[-1], 16)
N_arr = np.append(N_arr, N_c)

fa_pi = interpolate.interp1d(N_c, a_pi, fill_value='extrapolate', kind='quadratic')
fb_pi = interpolate.interp1d(N_c, b_pi, fill_value='extrapolate', kind='quadratic')
ftheta_pi = interpolate.interp1d(N_c, theta_pi, fill_value='extrapolate', kind='quadratic')
fa_n = interpolate.interp1d(N_c, a_n, fill_value='extrapolate', kind='quadratic')
fb_n = interpolate.interp1d(N_c, b_n, fill_value='extrapolate', kind='quadratic')
ftheta_n = interpolate.interp1d(N_c, theta_n, fill_value='extrapolate', kind='quadratic')

fig, (ax1, ax2) = plt.subplots(1, 2)

# Plot the fitting data
for spd, cpr_arr in speedline_cpr.items():
    n_arr = speedline_n[spd]

    # Plot the data
    ax1.plot(cpr_arr[:,0], cpr_arr[:,1])
    ax1.set(xlabel='Mass flow', ylabel='CPR')

    ax2.plot(n_arr[:,0], n_arr[:,1])
    ax2.set(xlabel='Mass flow', ylabel='eff.')

# interp/extrp.
m_arr = []
pi_arr = []
n_arr = []
for N in N_arr:
    cpr_plot(ax1, [-np.pi/8., np.pi/8.], fa_pi(N), fb_pi(N), ftheta_pi(N))
    eff_plot(ax2, [-np.pi/32., np.pi/8.], fa_n(N), fb_n(N), ftheta_n(N))

    tn = np.linspace(-np.pi/8., np.pi/8., 100)
    m, pi = ellipse_eval(tn, fa_pi(N), fb_pi(N), ftheta_pi(N))
    m, n = ellipse_eval(tn, fa_n(N), fb_n(N), ftheta_n(N))
    m_arr.append(m)
    pi_arr.append(pi)
    n_arr.append(n)
m_arr = np.array(m_arr).flatten()
pi_arr = np.array(pi_arr).flatten()
n_arr = np.array(n_arr).flatten()

ax1.tricontourf(m_arr, pi_arr, n_arr, 100)

plt.show()

df = pd.DataFrame({
    'a_pi': a_pi,
    'b_pi': b_pi,
    'theta_pi': theta_pi,
    'a_n': a_n,
    'b_n': b_n,
    'theta_n': theta_n,
    'surge_m': np.array(surge_line),
    'choke_pi': np.array(choke_line)
    }, index=N_c
)

print(df)
df.to_csv("j85_compressor_map.csv")