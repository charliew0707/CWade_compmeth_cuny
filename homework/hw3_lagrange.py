"""
HW3 - Exercise 6.16: The L1 Lagrange Point
Charles Wade
Computational Methods, CCNY
-------------------------------------------
The L1 Lagrange point is the position between the Earth and Moon at which
a satellite experiences a gravitational balance that keeps it in perfect
synchrony with the Moon's orbit. At this point, the combined gravitational
pull of the Earth (inward) and Moon (outward) provides exactly the centripetal
force required to maintain the satellite's circular orbit at the Moon's angular
velocity.

The distance r from the Earth's center to L1 satisfies the quintic equation
derived from balancing gravitational and centripetal forces. This fifth-order
polynomial has no closed-form solution, so we solve it numerically using
Newton's method.

Constants:
    M : Earth mass (from astropy.constants)
    G : Gravitational constant (from astropy.constants)
    m : Moon mass = 7.348e22 kg
    R : Earth-Moon distance = 3.844e8 m
    w : Angular velocity of Moon = 2.662e-6 s^-1

Goal: Find r accurate to at least four significant figures.

"""
"""
Part (a): Derivation
---------------------
Place Earth at the origin and the Moon at distance R along a line, with
the satellite sitting on that same line at distance r from Earth (and
so distance R - r from the Moon).

Two gravitational forces act on the satellite:
  - Earth pulls it inward (toward smaller r), with strength GM / r^2
  - The Moon pulls it outward (toward larger r), with strength Gm / (R-r)^2

For the satellite to stay on the Earth-Moon line and orbit Earth in
lockstep with the Moon, it must move in a circle of radius r with the
same angular velocity w as the Moon. That requires a net inward force
equal to the centripetal force w^2 r (satellite mass cancels out of
the whole equation, which is why it never appears).

So: (Earth's pull) - (Moon's pull) = (force needed to stay in circular orbit)

    GM / r^2 - Gm / (R - r)^2 = w^2 r

Clearing the fractions (multiplying through by r^2 (R-r)^2) and
collecting all terms turns this into a degree-5 (quintic) polynomial
equation in r, equivalent to:

    w^2 r^5 - 2 w^2 R r^4 + w^2 R^2 r^3 - (GM - Gm) r^2 + 2GMR r - GMR^2 = 0

which is the equation f(r) = 0 solved numerically below.

"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as c
from scipy.optimize import newton
import argparse

# Constants (all floats, all in SI base units: meters, kilograms, seconds)

G = c.G.si.value        # gravitational constant, m^3 / (kg s^2)
M = c.M_earth.si.value   # Earth mass, kg
m = 7.348e22             # Moon mass, kg
W = 2.662e-6             # angular velocity of Moon, s^-1
R = 3.844e8              # Earth-Moon distance, m


# Define the polynomial and deriv:
coeffs = [W**2, -2*W**2*R, W**2*R**2, -(G*M-G*m), 2*G*M*R, -G*M*R**2]
p = np.poly1d(coeffs)
p_prime = p.deriv()


# Function that defines the quintic polynomial:

def poly(r):
    return p(r)

def poly_prime(r):
    return p_prime(r)



# Function that does Newton's method, input func, initial guess, acc etc, and expected roots
# Will build 3 methods, 1 Newtown scratch, 2 Newton with np 3 Secant

def newton_scratch(func,func_prime,r=3.3e8,delta=1,eps=1e-11,steps=1000):
    for _ in range(steps):

        delta = func(r) / func_prime(r)
        r -= delta

        if abs(delta)<abs(eps):
            return r
    else:
        return None
    

def secant_method(func, x1=3.3e8, x2=3.1e8, eps=1e-11, steps=1000):
    f1 = func(x1)
    for _ in range(steps):
        f2 = func(x2)
        delta = f2 * (x2 - x1) / (f2 - f1)
        x1, x2 = x2, x2 - delta
        f1 = f2   

        if abs(delta) < eps:
            return x2
    return None


def newton_method(func, func_prime, x0, eps=1e-11, steps=1000):
    root = newton(func, x0, fprime=func_prime, tol=eps, maxiter=steps)
    return root


# need a argparse function to allow CLI inputs, set guesses to good default values, and newt for def meth:

def parse_args():
    parser = argparse.ArgumentParser(description="Solve for the L1 Lagrange point.")
    parser.add_argument(
        "--method", choices=["newton_scratch", "secant_method", "newton_method"], default="newton_method", 
        help="Enter root method: newton_scratch, secant_method, newton_method"
        )
    parser.add_argument("--x1", type=float, default=3.5e8, help="Enter starting guess in meters, between earth and moon")
    parser.add_argument("--x2", type=float, default=3.0e8, help="Enter second r guess, secant method specific")
    return parser.parse_args()

# main to run the code, use argparse call to choose method, metho returns the root closest to initial guess:
def main():
    args = parse_args()

    if args.method == "newton_scratch":
        r = newton_scratch(poly, poly_prime, args.x1)

    elif args.method == "newton_method":
        r = newton_method(poly, poly_prime, args.x1)

    elif args.method == "secant_method":
        r = secant_method(poly, args.x1, args.x2)

    # Nicer output with useful info, units and etc attached as all calcs in SI so r is always meters
    print(f"Method: {args.method}")
    print(f"Starting guess(es): x1 = {args.x1:.3e}" + (f", x2 = {args.x2:.3e}" if args.method == "secant_method" else ""))
    print(f"L1 distance from Earth's center: r = {r:.6e} m  ({r:.5g} m)")
    print(f"Fractional distance toward Moon: r/R = {r/R:.4f}")

if __name__ == "__main__":
    main()