# Channel cross-sections of Table 5.6.1 (Chow, 1959, Table 2.1) with Newton-Raphson
# solvers for the normal and critical depth.
#
# Symbols: y = flow depth, b0 = bottom width, s = side slope (H:1V), d0 = pipe diameter,
# n = Manning's roughness (mn), S0 = bed slope (s0), Q = discharge, T = top width.

import numpy as np

yinitial = 1.0
errorallow = 0.0001
g = 9.81


class Section:
    """Base class; each shape defines area, wetperi, topwidth, dtopwidth and gfun of y."""
    name = ''
    params = ()                    # names of the geometric parameters
    ystart = yinitial              # initial Newton iterate
    ynormal_max = np.inf           # upper bound kept by the normal-depth iterate
    ycritical_max = np.inf         # upper bound kept by the critical-depth iterate

    @property
    def geom(self):
        return {p: getattr(self, p) for p in self.params}

    def __repr__(self):
        args = ', '.join(f'{k}={v:g}' for k, v in self.geom.items())
        return f'{type(self).__name__}({args})'

    def hyrad(self, y):
        return self.area(y) / self.wetperi(y)

    def max_discharge(self, s0, mn):
        """Largest discharge with a normal depth (unbounded for open sections)"""
        return np.inf

    # ---- Normal depth: F(y) = A R^(2/3) - n Q / sqrt(S0) = 0, dF/dy = A R^(2/3) gfun(y)
    def normal_depth(self, s0, mn, Q, maxiter=100):
        """Newton-Raphson solution of Manning's equation; returns (yn, iterations).
        yn is nan when Q exceeds max_discharge(s0, mn)."""
        if Q > self.max_discharge(s0, mn):
            return np.nan, 0
        K = mn * Q / np.sqrt(s0)
        ytop = self.ynormal_max
        yn = self.ystart
        for it in range(1, maxiter + 1):
            AR = self.area(yn) * self.hyrad(yn) ** (2 / 3)
            Fyn = AR - K
            dFdyn = AR * self.gfun(yn)
            ynew = yn - Fyn / dFdyn
            if ynew <= 0:              # keep the iterate physical
                ynew = yn / 2
            if ynew >= ytop:           # stay below the upper bound
                ynew = (yn + ytop) / 2
            if abs(ynew - yn) < errorallow:
                return ynew, it
            yn = ynew
        return yn, maxiter

    # ---- Critical depth: F(y) = A^3 / T - Q^2 / g = 0, dF/dy = 3A^2 - A^3 T' / T^2
    def critical_depth(self, Q, maxiter=100):
        """Newton-Raphson solution of Q^2 T / (g A^3) = 1; returns (yc, iterations)."""
        ytop = self.ycritical_max
        yc = self.ystart
        for it in range(1, maxiter + 1):
            A = self.area(yc)
            T = self.topwidth(yc)
            Fyc = A ** 3 / T - Q * Q / g
            dFdyc = 3 * A * A - A ** 3 * self.dtopwidth(yc) / T ** 2
            ynew = yc - Fyc / dFdyc
            if ynew <= 0:              # keep the iterate physical
                ynew = yc / 2
            if ynew >= ytop:           # stay below the upper bound
                ynew = (yc + ytop) / 2
            if abs(ynew - yc) < errorallow:
                return ynew, it
            yc = ynew
        return yc, maxiter

    def froude(self, y, Q):
        """Froude number V / sqrt(g D), with hydraulic depth D = A / T."""
        A = self.area(y)
        return Q / A / np.sqrt(g * A / self.topwidth(y))

    # ---- Checks and reports
    def check_gfun(self, ys=(0.1, 0.5, 1.0, 1.5), dy=1e-6):
        """Compare gfun with a central finite difference of ln(A R^(2/3))."""
        def lnAR(y):
            return np.log(self.area(y) * self.hyrad(y) ** (2 / 3))
        for y in ys:
            fd = (lnAR(y + dy) - lnAR(y - dy)) / (2 * dy)
            print(f'y = {y:4.2f}   table = {self.gfun(y):.6f}   finite difference = {fd:.6f}')

    def report(self, s0, mn, Q):
        yn, iters = self.normal_depth(s0, mn, Q)
        print('iter = ', iters, '  Normal Depth = %.3f' % yn)
        yc, iters = self.critical_depth(Q)
        print('iter = ', iters, '  Critical Depth = %.3f' % yc)
        print('Froude number at normal depth = %.3f' % self.froude(yn, Q))
        if np.isfinite(self.max_discharge(s0, mn)):
            print('Part-full capacity Qmax = %.3f m3/s' % self.max_discharge(s0, mn))

    # ---- Drawing helpers (used by channel_plot)
    def bank_height(self, yn, yc):
        return 1.5 * max(yn, yc)

    def style_axes(self, ax):
        pass


class Rectangle(Section):
    name = 'Rectangle'
    params = ('b0',)

    def __init__(self, b0):
        self.b0 = b0

    def area(self, y):
        return self.b0 * y

    def wetperi(self, y):
        return self.b0 + 2 * y

    def topwidth(self, y):
        return self.b0 + 0 * y

    def dtopwidth(self, y):
        """dT/dy"""
        return 0 * y

    def gfun(self, y):
        """2/(3R) dR/dy + 1/A dA/dy  (last row of Table 5.6.1)"""
        b0 = self.b0
        return (5 * b0 + 6 * y) / (3 * y * (b0 + 2 * y))

    def boundary(self, hb):
        """Channel boundary (x, y) up to bank height hb"""
        return [0, 0, self.b0, self.b0], [hb, 0, 0, hb]

    def water_edges(self, y):
        """x of the left and right water-surface edges at depth y"""
        return 0.0, self.b0

    def water_polygon(self, y):
        return [(0, y), (0, 0), (self.b0, 0), (self.b0, y)]


class Trapezoid(Section):
    name = 'Trapezoid'
    params = ('b0', 's')

    def __init__(self, b0, s):
        self.b0, self.s = b0, s

    def area(self, y):
        return y * (self.b0 + self.s * y)

    def wetperi(self, y):
        return self.b0 + 2 * y * np.sqrt(1 + self.s ** 2)

    def topwidth(self, y):
        return self.b0 + 2 * self.s * y

    def dtopwidth(self, y):
        """dT/dy"""
        return 2 * self.s + 0 * y

    def gfun(self, y):
        """2/(3R) dR/dy + 1/A dA/dy  (last row of Table 5.6.1)"""
        b0, s = self.b0, self.s
        r = np.sqrt(1 + s * s)
        return ((b0 + 2 * s * y) * (5 * b0 + 6 * y * r) + 4 * s * y * y * r) \
            / (3 * y * (b0 + s * y) * (b0 + 2 * y * r))

    def boundary(self, hb):
        """Channel boundary (x, y) up to bank height hb"""
        return [-self.s * hb, 0, self.b0, self.b0 + self.s * hb], [hb, 0, 0, hb]

    def water_edges(self, y):
        """x of the left and right water-surface edges at depth y"""
        return -self.s * y, self.b0 + self.s * y

    def water_polygon(self, y):
        return [(-self.s * y, y), (0, 0), (self.b0, 0), (self.b0 + self.s * y, y)]


class Triangle(Section):
    name = 'Triangle'
    params = ('s',)

    def __init__(self, s):
        self.s = s

    def area(self, y):
        return self.s * y * y

    def wetperi(self, y):
        return 2 * y * np.sqrt(1 + self.s ** 2)

    def topwidth(self, y):
        return 2 * self.s * y

    def dtopwidth(self, y):
        """dT/dy"""
        return 2 * self.s + 0 * y

    def gfun(self, y):
        """2/(3R) dR/dy + 1/A dA/dy  (last row of Table 5.6.1)"""
        return 8 / (3 * y)

    def boundary(self, hb):
        """Channel boundary (x, y) up to bank height hb"""
        return [-self.s * hb, 0, self.s * hb], [hb, 0, hb]

    def water_edges(self, y):
        """x of the left and right water-surface edges at depth y"""
        return -self.s * y, self.s * y

    def water_polygon(self, y):
        return [(-self.s * y, y), (0, 0), (self.s * y, y)]


# A R^(2/3) of a circular section peaks at y ≈ 0.938 d0, so a part-full
# normal depth only exists for Q up to Qmax; the iterate is kept below that peak.
YMAX = 0.9382


class Circle(Section):
    name = 'Circle'
    params = ('d0',)

    def __init__(self, d0):
        self.d0 = d0
        self.ystart = 0.5 * d0
        self.ynormal_max = YMAX * d0
        self.ycritical_max = d0        # stay inside the pipe

    def theta(self, y):
        """Angle subtended at the centre by the water surface (rad)"""
        return 2 * np.arccos(1 - 2 * y / self.d0)

    def area(self, y):
        th = self.theta(y)
        return (th - np.sin(th)) * self.d0 ** 2 / 8

    def wetperi(self, y):
        return self.theta(y) * self.d0 / 2

    def topwidth(self, y):
        return 2 * np.sqrt(y * (self.d0 - y))

    def dtopwidth(self, y):
        """dT/dy"""
        return (self.d0 - 2 * y) / np.sqrt(y * (self.d0 - y))

    def gfun(self, y):
        """2/(3R) dR/dy + 1/A dA/dy  (last row of Table 5.6.1)"""
        th = self.theta(y)
        return 4 * (2 * np.sin(th) + 3 * th - 5 * th * np.cos(th)) \
            / (3 * self.d0 * th * (th - np.sin(th)) * np.sin(th / 2))

    def max_discharge(self, s0, mn):
        """Largest discharge a part-full circular section can carry at normal depth"""
        y = YMAX * self.d0
        return self.area(y) * self.hyrad(y) ** (2 / 3) * np.sqrt(s0) / mn

    def boundary(self, hb):
        """Pipe wall (x, y); centre at (0, d0/2)"""
        phi = np.linspace(0, 2 * np.pi, 361)
        return self.d0 / 2 * np.sin(phi), self.d0 / 2 * (1 - np.cos(phi))

    def water_edges(self, y):
        """x of the left and right water-surface edges at depth y"""
        h = self.topwidth(y) / 2
        return -h, h

    def water_polygon(self, y):
        th = self.theta(y)
        phi = np.linspace(-th / 2, th / 2, 181)
        return np.column_stack([self.d0 / 2 * np.sin(phi), self.d0 / 2 * (1 - np.cos(phi))])

    def bank_height(self, yn, yc):
        return self.d0

    def style_axes(self, ax):
        ax.set_aspect('equal')
        ax.set_xlim(-0.6 * self.d0, 0.6 * self.d0)
