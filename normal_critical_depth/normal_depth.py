# Hydraulics of a trapezoidal channel: geometry, normal and critical depth by Newton-Raphson

import numpy as np

yinitial = 1.0
errorallow = 0.0001
g = 9.81


def area(y, b0, s):
    return y * (b0 + s * y)


def wetperi(y, b0, s):
    return b0 + 2 * y * np.sqrt(1 + s * s)


def hyrad(y, b0, s):
    return area(y, b0, s) / wetperi(y, b0, s)


def topwidth(y, b0, s):
    return b0 + 2 * s * y


def normal_depth(b0, s, s0, mn, Q, maxiter=100):
    """Newton-Raphson solution of Manning's equation; returns (yn, iterations)."""
    yn = yinitial
    for it in range(1, maxiter + 1):
        R = hyrad(yn, b0, s)
        Fyn = area(yn, b0, s) ** (5 / 3) / wetperi(yn, b0, s) ** (2 / 3) - mn * Q / np.sqrt(s0)
        dFdyn = (5 / 3) * topwidth(yn, b0, s) * R ** (2 / 3) \
            - (2 / 3) * R ** (5 / 3) * 2 * np.sqrt(1 + s * s)
        ynew = yn - Fyn / dFdyn
        if ynew <= 0:              # keep the iterate physical
            ynew = yn / 2
        if abs(ynew - yn) < errorallow:
            return ynew, it
        yn = ynew
    return yn, maxiter


def critical_depth(b0, s, Q, maxiter=100):
    """Newton-Raphson solution of Q^2 T / (g A^3) = 1; returns (yc, iterations)."""
    yc = yinitial
    for it in range(1, maxiter + 1):
        A = area(yc, b0, s)
        T = topwidth(yc, b0, s)
        Fyc = A ** 3 / T - Q * Q / g
        dFdyc = 3 * A * A - 2 * s * A ** 3 / T ** 2
        ynew = yc - Fyc / dFdyc
        if ynew <= 0:              # keep the iterate physical
            ynew = yc / 2
        if abs(ynew - yc) < errorallow:
            return ynew, it
        yc = ynew
    return yc, maxiter


def froude(y, b0, s, Q):
    """Froude number V / sqrt(g D), with hydraulic depth D = A / T."""
    A = area(y, b0, s)
    return Q / A / np.sqrt(g * A / topwidth(y, b0, s))


if __name__ == '__main__':
    yn, iters = normal_depth(b0=20.0, s=2.0, s0=0.001, mn=0.020, Q=40)
    print('iter = ', iters, '  Normal Depth = %.3f' % yn)
    yc, iters = critical_depth(b0=20.0, s=2.0, Q=40)
    print('iter = ', iters, '  Critical Depth = %.3f' % yc)
