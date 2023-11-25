import numpy as np
import matplotlib.pyplot as plt

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
[15.7208,	0.909775],
[15.7681,	0.907296],
[15.8075,	0.904089],
[15.8512,	0.899604],
[15.8882,	0.894653],
[18,0]
])

PX = np.array([[0.0, 0.0],
              [0.05, 0.15],
              [0.4, 0.45],
              [0.8, 0.4],
              [1.0, 0.0]])

d = np.concatenate(([0.0], np.linalg.norm(P[:-1,:]-P[1:,:], axis=1)))
D = np.concatenate(([0.0], d[0:-1] + d[1:]))
u = D/D[-1]

def bezier_eval(control_points, t):
    return (1.-t)**3*control_points[0] + 3.*(1.-t)**2*t*control_points[1] + 3.*(1.0-t)*t**2*control_points[2] + t**3*control_points[3]

def bezier_deriv_eval(control_points, t):
    return 3.*(1-t)**2*(control_points[1] - control_points[0]) + 6.*(1-t)*t*(control_points[2] - control_points[1]) + 3.*t**2*(control_points[3] - control_points[2])

def bezier_dderiv_eval(control_points, t):
    return 6.*(1-t)*(control_points[2] - 2.*control_points[1] + control_points[0]) + 6.*t*(control_points[3] - 2.*control_points[2] + control_points[1])

def reparameterize(bezier, points, parameters):
    return [newtonRaphsonRootFind(bezier, point, u) for point, u in zip(points, parameters)]

def newtonRaphsonRootFind(bezier, p, u):
    """
       Newton's root finding algorithm calculates f(x)=0 by reiterating
       x_n+1 = x_n - f(x_n)/f'(x_n)

       We are trying to find curve parameter u for some point p that minimizes
       the distance from that point to the curve. Distance point to curve is d=q(u)-p.
       At minimum distance the point is perpendicular to the curve.
       We are solving
       f = q(u)-p * q'(u) = 0
       with
       f' = q'(u) * q'(u) + q(u)-p * q''(u)

       gives
       u_n+1 = u_n - |q(u_n)-p * q'(u_n)| / |q'(u_n)**2 + q(u_n)-p * q''(u_n)|
    """
    d = bezier_eval(bezier, u)-p
    numerator = (d * bezier_deriv_eval(bezier, u)).sum()
    denominator = (bezier_deriv_eval(bezier, u)**2 + d * bezier_dderiv_eval(bezier, u)).sum()

    if denominator == 0.0:
        return u
    else:
        return u - numerator/denominator

def fit_bezier(p_arr, u_arr, l_tan, r_tan):
    bezier = [p_arr[0], None, None, p_arr[-1]]

    # compute the A's
    A = np.zeros((len(u_arr), 2, 2))
    for i, u in enumerate(u_arr):
        A[i][0] = l_tan*(3*(1-u)**2*u)
        A[i][1] = r_tan*(3*(1-u)*u**2)

    # Create the C and X matrices
    C = np.zeros((2, 2))
    X = np.zeros(2)

    for i, (point, u) in enumerate(zip(p_arr, u_arr)):
        C[0][0] += np.dot(A[i][0], A[i][0])
        C[0][1] += np.dot(A[i][0], A[i][1])
        C[1][0] += np.dot(A[i][0], A[i][1])
        C[1][1] += np.dot(A[i][1], A[i][1])

        tmp = point - bezier_eval([p_arr[0], p_arr[0], p_arr[-1], p_arr[-1]], u)

        X[0] += np.dot(A[i][0], tmp)
        X[1] += np.dot(A[i][1], tmp)

    # Compute the determinants of C and X
    det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1]
    det_C0_X  = C[0][0] * X[1] - C[1][0] * X[0]
    det_X_C1  = X[0] * C[1][1] - X[1] * C[0][1]

    # Finally, derive alpha values
    alpha_l = 0.0 if det_C0_C1 == 0 else det_X_C1 / det_C0_C1
    alpha_r = 0.0 if det_C0_C1 == 0 else det_C0_X / det_C0_C1

    # If alpha negative, use the Wu/Barsky heuristic (see text) */
    # (if alpha is 0, you get coincident control points that lead to
    # divide by zero in any subsequent NewtonRaphsonRootFind() call. */
    segment_length = np.linalg.norm(p_arr[0] - p_arr[-1])
    epsilon = 1.0e-6 * segment_length
    if alpha_l < epsilon or alpha_r < epsilon:
        # fall back on standard (probably inaccurate) formula, and subdivide further if needed.
        bezier[1] = bezier[0] + l_tan * (segment_length / 3.0)
        bezier[2] = bezier[3] + r_tan * (segment_length / 3.0)

    else:
        # First and last control points of the Bezier curve are
        # positioned exactly at the first and last data points
        # Control points 1 and 2 are positioned an alpha distance out
        # on the tangent vectors, left and right, respectively
        bezier[1] = bezier[0] + l_tan * alpha_l
        bezier[2] = bezier[3] + r_tan * alpha_r

    return np.array(bezier)

fig, ax = plt.subplots()

def fit_curve(p_arr, u, l_tan, r_tan):
    bezier_fit =  fit_bezier(p_arr, u, l_tan, r_tan)

    for i in range(100):
        u_new = reparameterize(bezier_fit, p_arr, u)
        bezier_fit = fit_bezier(p_arr, u_new, l_tan, r_tan)
        u = u_new

        print(u)

    t = np.linspace(0,1,100)
    p_reg = np.array([bezier_eval(bezier_fit, u) for u in t])
    ax.plot(p_reg[:,0], p_reg[:,1])
    return bezier_fit

t0 = (P[1]-P[0])/np.linalg.norm(P[1]-P[0])
tend = (P[-2]-P[-1])/np.linalg.norm(P[2]-P[1])

def tangent_from_range(points):
    t1 = (points[1]-points[0])/np.linalg.norm(points[1]-points[0])
    t2 = (points[-2]-points[-1])/np.linalg.norm(points[2]-points[1])
    return (t1, t2)

tangents = tangent_from_range(P[0:8])
bezier = fit_curve(P[0:8], u, np.array([0,1]), np.array([-0.1,0])) 

#t = np.linspace(0,1,100)
#p_reg = np.array([bezier_eval(beziers[-1], u) for u in t])

ax.plot(bezier[:,0], bezier[:,1])
ax.scatter(P[:,0], P[:,1], c='r')

plt.show()