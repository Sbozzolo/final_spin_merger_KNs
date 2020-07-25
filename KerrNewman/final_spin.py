#!/usr/bin/env python3

from scipy.optimize import root_scalar
from KerrNewman.KerrNewman import KerrNewman


def estimate_kn_spin(m_1=0.5,
                     m_2=0.5,
                     lambda_1=0,
                     lambda_2=0,
                     tolerance=1e-3):

    # lambda_1, lambda_2, the charge to mass ratios of the two black holes
    # q_1, q_2, the charges of the two black holes
    # m_1, m_2, the masses of the two black holes
    q_1, q_2 = lambda_1 * m_1, lambda_2 * m_2

    # Q = q_1 + q_2, the charge of the final black hole
    #
    # M = m_1 + m_2, the mass of the final black hole (it is assumed that not
    #                much energy is lost in gravitational waves)
    # a            , the kerr parameter of the final black hole
    Q = q_1 + q_2
    M = m_1 + m_2

    # nu
    nu = m_1 * m_2 / M**2

    # The charge-to-mass ratio is the reduced charge over the reduced mass
    e = ((q_1 * q_2) / Q) / ((m_1 * m_2) / M) if Q != 0 else 0

    def equation(A):
        return (KerrNewman(Q=Q, a=A)).rl_iscos(e)[1] * nu - A

    return root_scalar(equation, bracket=[0, 0.9],
                       xtol=tolerance, rtol=tolerance,
                       method='brentq').root
