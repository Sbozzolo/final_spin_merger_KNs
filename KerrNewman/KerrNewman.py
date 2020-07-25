#!/usr/bin/env python3

# This code implements the method decribd in 1706.06519 to estimate the final
# spin of the merger of two charged black holes in Einstein-Maxwell's theory.
# Relevant equations are (5.1) - (5.4).

# Isco computed with equations from 1907.02158 (3.13)

from sympy import cos, sin, sqrt
import sympy
import numpy as np


class KerrNewman:
    """Kerr-Newman black hole: metric and ISCO.
    """

    def __init__(self, Q=0, a=0, M=1):
        self.Q = Q
        self.M = M
        self.a = a

    # Kerr-Newman metric
    def Delta(self, r):
        return r**2 + self.a**2 - 2 * self.M * r + self.Q**2

    def rho(self, r, theta):
        return sqrt(r**2 + self.a**2 * cos(theta)**2)

    def g_tt(self, r, theta):
        return -(1 - (2 * self.M * r - self.Q**2) / self.rho(r, theta)**2)

    def g_rr(self, r, theta):
        return self.rho(r, theta) / self.Delta(r)

    def g_tphi(self, r, theta):
        return -(2 * self.M * r - self.Q**2) * self.a * sin(theta)**2 / \
            self.rho(r, theta)**2

    def g_phiphi(self, r, theta):
        return sin(theta)**2 / self.rho(r, theta)**2 * (
            (r**2 + self.a**2)**2 - self.a**2 * self.Delta(r) * sin(theta)**2)

    def A_t(self, r, theta):
        return self.Q * r / self.rho(r, theta)**2

    def A_phi(self, r, theta):
        return -self.Q * self.a * r * sin(theta)**2 / self.rho(r, theta)**2

    # Effective potential
    def V_eff(self, r, theta, L, E, e=0):
        return 1/(self.g_rr(r, theta) * self.Delta(r)) * (
            self.g_tt(r, theta) * (L + e * self.A_phi(r, theta))**2
            + 2 * self.g_tphi(r, theta) * (L + e * self.A_phi(r,
                                                              theta)) * (E - e * self.A_t(r, theta))
            + self.g_phiphi(r, theta) * (E - e * self.A_t(r, theta))**2 - self.Delta(r))

    def rl_iscos(self,
                e=0,
                guess=None):
        """Return radius and specific angular momentum
        of closest ISCO to the horizon.
        """
        r, L, E = sympy.symbols('r L E')
        V_eff_equat = self.V_eff(r, np.pi/2, L, E, e)
        V_eff_equat_d1 = sympy.diff(V_eff_equat, r)
        V_eff_equat_d2 = sympy.diff(V_eff_equat_d1, r)

        if (guess is None):
            guess = [3 * self.M, 3 * self.M**2, self.M]

        return sympy.nsolve(
            [V_eff_equat, V_eff_equat_d1, V_eff_equat_d2],
            [r, L, E],
            guess)[:2]          # return r, L, E

    # def r_iscos_analytical(self,
    #                        e=0):
    #     """This formula works only for e = 0.
    #     """
    #     if (e != 0):
    #         raise ValueError("r_iscos_analytical works only for e=0")

    #     r = sympy.symbols('r')
    #     isco_plus_eq = self.a**2 * (
    #         self.M * r**2 *
    #         (7 * self.M + 3 * r) + 8 * self.Q**4 - 2 * self.Q**2 * r *
    #         (7 * self.M + 2 * r)) + (
    #             self.M * r**2 *
    #             (6 * self.M - r) + 4 * self.Q**4 - 9 * self.M * self.Q**2 *
    #             r) * (r *
    #                   (r - 3 * self.M) + 2 * self.Q**2) + 2 * self.a * (
    #                       4 * self.Q**2 - 3 * self.M * r) * (
    #                           self.a**2 + r * (r - 2 * self.M) +
    #                           self.Q**2) * sqrt(self.M * r - self.Q**2)
    #     isco_minus_eq = self.a**2 * (
    #         self.M * r**2 *
    #         (7 * self.M + 3 * r) + 8 * self.Q**4 - 2 * self.Q**2 * r *
    #         (7 * self.M + 2 * r)) + (
    #             self.M * r**2 *
    #             (6 * self.M - r) + 4 * self.Q**4 - 9 * self.M * self.Q**2 *
    #             r) * (r *
    #                   (r - 3 * self.M) + 2 * self.Q**2) - 2 * self.a * (
    #                       4 * self.Q**2 - 3 * self.M * r) * (
    #                           self.a**2 + r * (r - 2 * self.M) +
    #                           self.Q**2) * sqrt(self.M * r - self.Q**2)
    #     isco_plus = sympy.nsolve(isco_plus_eq, r, 5 * self.M)
    #     isco_minus = sympy.nsolve(isco_minus_eq, r, 8 * self.M)
    #     return [isco_plus, isco_minus]
