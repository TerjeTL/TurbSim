import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import scipy.special

P_data = np.array([
    [13.5, 0],
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
[16.25,0]
])
P_data2 = np.array([
    [16.2, 0.0],
    [16.9325,0.811203],
    [17.0819,0.818235],
    [17.2144,0.824846],
    [17.3552,0.831293],
    [17.4876,0.83623 ],
    [17.5822,0.839664],
    [17.6451,0.840747],
    [17.6799,0.838988],
    [17.7012,0.834836],
    [17.7319,0.824273],
    [17.95, 0.0]])

P = np.array([[0.0, 0.0],
              [0.35, 1.2],
              [0.6, 1.5],
              [0.8, 1.0],
              [1.0, 0.0]])

t0 = np.array([0.,1.])
t1 = np.array([.0,1.0])

class Bezier:
    def __init__(self, n, data, tan0, tan1) -> None:
        self.n = n
        self.data = data
        self.tan0 = tan0
        self.tan1 = tan1 

        # Maybe i'll make this configurable..
        free_ctrl_pts = np.max(n-3, 0)
        x0 = np.array([1.0, 1.0])
        if (free_ctrl_pts == 1):
            x0 = np.array([1.0, 1.0, 1.0, 1.0])

        d = np.concatenate(([0.0], np.linalg.norm(data[:-1,:]-data[1:,:], axis=1)))
        D = np.cumsum(d)
        self.t_data = D/D[-1]
        
        c1 = data[0] + tan0*x0[0]
        cn1 = data[-1] + tan1*x0[-1]
        self.bezier = np.array([data[0], c1] + [[0.,0.]]*free_ctrl_pts + [cn1, data[-1]])

        self.delta_x = self.bezier[-1,0] - self.bezier[0,0]
        self.cx = self.bezier[0,0] + 0.5*self.delta_x
        self.dx_step = self.delta_x/2.0 * 0.005

        x_pos = self.cx

        if free_ctrl_pts == 1:
            x_pos = self.data[np.argmax(self.data[:,1]),0]
            y_avg = np.average(self.data[:,1])
            self.bezier[2] = np.array([x_pos, y_avg])
            x0[1] = x_pos
            x0[2] = y_avg
        
        ctrl_point_bounds_lower = [x_pos - self.dx_step, -np.inf]
        ctrl_point_bounds_upper = [x_pos + self.dx_step, np.inf]
        self.bounds = ([0.0] + ctrl_point_bounds_lower*free_ctrl_pts + [0.0], [np.inf] + ctrl_point_bounds_upper*free_ctrl_pts + [np.inf])

        self.x = x0
        pass
    
    def update_bounds(self): # Found that this works to keep stability under control. it's majorly scuffed
        curr_x = self.x[1]
        x_max = self.bezier[-1,0]
        x_min = self.bezier[0,0]
        dx = self.delta_x/2.0 * 0.01

        if (self.bounds[1][1] - curr_x < dx):
            self.bounds[1][1] = np.fmin(self.bounds[1][1] + self.dx_step, x_max)
            self.bounds[0][1] = np.fmax(curr_x - 2.*dx, x_min)
            print(self.bounds)
        elif (self.bounds[0][1] - curr_x < dx):
            self.bounds[0][1] = np.fmax(self.bounds[0][1] - self.dx_step, x_min)
            self.bounds[1][1] = np.fmin(self.x[1] + 2.*dx, x_max)
            print(self.bounds)

    def bernstein_polynomial(self, n, i, t):
        return scipy.special.binom(n, i)*t**i*(1-t)**(n-i)

    def bezier_eval(self, bezier, t):
        B = np.array([self.bernstein_polynomial(self.n, i, t) for i in range(self.n+1)])
        C = B@bezier
        return C

    def bezier_eval_d1(self, bezier, t):
        B = np.array([self.bernstein_polynomial(self.n-1, i, t) for i in range(self.n)])
        dC = self.n*(B@bezier[1:] - B@bezier[:-1])
        return dC

    def bezier_eval_d2(self, bezier, t):
        B = np.array([self.bernstein_polynomial(self.n-2, i, t) for i in range(self.n-1)])
        ddC = self.n*(self.n-1)*(B@(bezier[2:] - 2.*bezier[1:-1] + bezier[:-2]))
        return ddC

    def update_bezier(self, x):
        bezier = self.bezier
        bezier[1] = bezier[0] + self.tan0*x[0]
        bezier[-2] = bezier[-1] + self.tan1*x[-1]
        if (len(x) > 2):
            free_ctrl_pts = np.reshape(x[1:-1], (-1, 2))
            for i, pt in enumerate(free_ctrl_pts):
                bezier[i+2] = pt
        return bezier

    def newton_rhapson_roots(self, bezier, p_loc, t):
        # t_n+1 = t_n - |b(t_n)-p * b'(t_n)| / |b'(t_n)**2 + b(t_n)-p * b''(t_n)|
        c_eval = self.bezier_eval(bezier, t)
        c_eval_d1 = self.bezier_eval_d1(bezier, t)
        c_eval_d2 = self.bezier_eval_d2(bezier, t)

        dp = (c_eval-p_loc)
        prod = dp*c_eval_d1

        num = np.sum((c_eval-p_loc) * c_eval_d1)
        denom = np.sum(c_eval_d1**2 + (c_eval-p_loc)*c_eval_d2)

        dt = 0.0

        if denom != 0.0:
            dt = -num/denom
        
        return t + dt

    def bezier_eval_y(self, x, t, y):
        bezier = self.update_bezier(x)
        return np.array([self.bezier_eval(bezier, t_it)[1] for t_it in t]) - y

    def fit_bezier(self, n_it = 100):
        dx_per_iteration = 0.005
        for i in range(n_it):
            if i > 0:
                self.t_data = np.array([self.newton_rhapson_roots(self.bezier, point, t_it) for (point, t_it) in zip(self.data, self.t_data)])
            res_lsq = least_squares(self.bezier_eval_y, self.x, args=(self.t_data, self.data[:,1]), bounds=self.bounds, loss='cauchy', f_scale=0.1)
            self.x = res_lsq.x
            self.bezier = self.update_bezier(self.x)

            if (self.n == 4):
                self.update_bounds()
        
        print(self.x)
        print(self.bounds)
    
    def plot(self, ax):
        t = np.linspace(0, 1, 100)
        X = np.array([self.bezier_eval(self.bezier, t_it) for t_it in t])
        X_data = np.array([self.bezier_eval(self.bezier, t_it) for t_it in self.t_data])

        ax.plot(X[:,0], X[:,1])
        ax.plot(self.bezier[:,0], self.bezier[:,1], '-.kx')
        ax.scatter(X_data[:,0], X_data[:,1], c='b')
        ax.scatter(self.data[:,0], self.data[:,1], c='r')

bezier_curve = Bezier(4, P_data, t0, t1)
bezier_curve.fit_bezier(10)

fig, ax = plt.subplots()
bezier_curve.plot(ax)
plt.show()