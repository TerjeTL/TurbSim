import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

P = np.array([
    [12, 0],
    [14.7217,0.891373],
[14.8267,	0.894628],
[15.0766,	0.901138],
[15.3035,	0.906034],
[15.4899,	0.909777],
[15.5652,	0.911095],
[15.6127,	0.911253],
[15.6694,	0.910936],
#[15.7208,	0.909775],
#[15.7681,	0.907296],
#[15.8075,	0.904089],
#[15.8512,	0.899604],
#[15.8882,	0.894653],
#[16.5,0]
])

P_data = np.array([[0.0, 0.0],
              [0.2, 0.1],
              [0.35, 0.2],
              [0.85, 0.8],
              [0.6, 0.28],
              [1.0, 0.0]])

t0 = np.array([0.,1.])
t1 = np.array([-1.0,0.0])
bezier = [P[0], None, None, P[-1]]

d = np.concatenate(([0.0], np.linalg.norm(P[:-1,:]-P[1:,:], axis=1)))
D = np.cumsum(d)
t = D/D[-1]

def bezier_eval_x(x, t, y):
    c1 = bezier[0][0] + t0[0]*x[0]
    c2 = bezier[3][0] + t1[0]*x[1]

    return (1.-t)**3*bezier[0][0] + 3.*(1.-t)**2*t*c1 + 3.*(1-t)*t**2*c2 + t**3*bezier[3][0] - y

def bezier_eval_y(x, t, y):
    c1 = bezier[0][1] + t0[1]*x[0]
    c2 = bezier[3][1] + t1[1]*x[1]

    return (1.-t)**3*bezier[0][1] + 3.*(1.-t)**2*t*c1 + 3.*(1-t)*t**2*c2 + t**3*bezier[3][1] - y

def bezier_eval(c_arr, t):
    return (1.-t)**3*bezier[0] + 3.*(1.-t)**2*t*bezier[1] + 3.*(1-t)*t**2*bezier[2] + t**3*bezier[3]

def bezier_deriv_eval(c_arr, t):
    return 3.*(1-t)**2*(c_arr[1] - c_arr[0]) + 6.*(1-t)*t*(c_arr[2] - c_arr[1]) + 3.*t**2*(c_arr[3] - c_arr[2])

def bezier_dderiv_eval(c_arr, t):
    return 6.*(1-t)*(c_arr[2] - 2.*c_arr[1] + c_arr[0]) + 6.*t*(c_arr[3] - 2.*c_arr[2] + c_arr[1])

def newton_rhapson_roots(c_arr, data, t):
    # t_n+1 = t_n - |b(t_n)-p * b'(t_n)| / |b'(t_n)**2 + b(t_n)-p * b''(t_n)|
    b_eval = bezier_eval(c_arr, t)
    b_deval = bezier_deriv_eval(c_arr, t)
    b_ddeval = bezier_dderiv_eval(c_arr, t)

    dp = (b_eval-data)
    prod = dp*b_deval

    num = np.sum((b_eval-data) * b_deval)
    denom = np.sum(b_deval**2 + (b_eval-data)*b_ddeval)

    dt = 0.0

    if denom != 0.0:
        dt = -num/denom
    
    return t + dt
    

def optimize_iter(bezier, data, t, x):
    t = np.array([newton_rhapson_roots(bezier, point, t_it) for (point, t_it) in zip(data, t)])
    resy_lsq = least_squares(bezier_eval_y, x, args=(t, data[:,1]), bounds=([0.0, 0.0], [np.inf, np.inf]), loss='cauchy', f_scale=0.1)
    x = resy_lsq.x
    bezier[1] = bezier[0] + t0*x[0]
    bezier[2] = bezier[3] + t1*x[1]
    return (bezier, t, x)

def fit_bezier(bezier, t, data, x0, it=100):
    bezier[1] = bezier[0] + t0*x0[0]
    bezier[2] = bezier[3] + t1*x0[1]
    result = (bezier, t, x0)
    for i in range(it):
        result = optimize_iter(result[0], data, result[1], result[2])
    
    return result

result = fit_bezier(bezier, t, P, np.array([1.,1.]), 100)

bezier = result[0]
t = result[1]

t_plot = np.linspace(0, 1, 100)
X = np.array([bezier_eval(bezier, t) for t in t_plot])
X_at_dp = np.array([bezier_eval(bezier, t) for t in t])

fig, ax = plt.subplots()

ax.scatter(P[:,0], P[:,1], c='r')
ax.scatter(X_at_dp[:,0], X_at_dp[:,1], c='b')
ax.plot(X[:,0], X[:,1], '-b')

plt.show()