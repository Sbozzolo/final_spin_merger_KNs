#!/usr/bin/env python3

from KerrNewman.final_spin import estimate_kn_spin
from numpy import log10
import argparse


desc = "Estimate final spin of merger of Kerr-Newman black holes."
p = argparse.ArgumentParser(description=desc)

p.add_argument('--m1', type=float, help='First mass')
p.add_argument('--m2', type=float, help='Second mass')
p.add_argument('--M', type=float, help='Total mass')
p.add_argument('--q', type=float, help='Mass ratio')
p.add_argument('--tolerance', type=float,
               help='Tolerance in the solution', default=1e-4)
p.add_argument('--lambda1', type=float, help='First charge/mass ratio',
               default=0)
p.add_argument('--lambda2', type=float, help='Second charge/mass ratio',
               default=0)

args = p.parse_args()

# Try using m1 and m2, if they are not available, use M and q
if (args.m1 is not None):
    m1 = args.m1
if (args.m2 is not None):
    m2 = args.m2
if (args.m1 is None and args.m2 is None):
    if (args.M is None or args.q is None):
        raise RuntimeError("You did not specify enough parameters!")
    else:
        m2 = args.M / (1 + args.q)
        m1 = args.M - m2

lambda1 = args.lambda1
lambda2 = args.lambda2
tol = args.tolerance

A = estimate_kn_spin(m1, m2, lambda1, lambda2, tolerance=tol)

print(f"{A:.{int(-log10(tol))}f}")
