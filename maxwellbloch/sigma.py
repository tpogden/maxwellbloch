# -*- coding: utf-8 -*-

""" Transition or 'flip' operators for quantum systems.

Thomas Ogden <t@ogden.eu>
"""

import qutip as qu

def sigma(n, a, b):
    """ Returns a transition or 'flip' operator.

    Returns the flip operator between states |a> and |b>, given by |a><b|, for
    a system with n states.

    For example, sigma(2,0,1) returns the quantum object represented by:
        [[ 0.  1.]
         [ 0.  0.]]

    Args:
        n: num levels
        a: level transitioned from
        b: level transitioned to

    Returns:
        |a><b|
    """

    return qu.basis(n, a) * qu.basis(n, b).dag()


def sigma_N(n, a, b, i, N):
    """ Returns the transition or 'flip' operator for a product system.

    Returns the flip operator between states |a> and |b> on the ith substate,
    given by |a><b|_i, for a system with n states and N subsystems. This
    extends sigma to a system of N interacting subsystems (for example N=2
    atoms with n=3 levels) by making a tensor product of the subsystems.

    For example, sigma(2,0,1,1,2) returns the quantum object represented by:

    [[ 0.  1.  0.  0.]
     [ 0.  0.  0.  0.]
     [ 0.  0.  0.  1.]
     [ 0.  0.  0.  0.]]

    Args:
        n: num levels
        a: level transitioned from
        b: level transitioned to
        i: which subsystem the jump is in (e.g. the ith atom)
        N: number of subsystems (e.g. atoms)

    Returns:
        |a><b|_i
    """

    # Make a list of identity matrices. Then we'll change the ith one to be our
    # subsystem flip operator.
    sigma_list = [qu.identity(n) for j in range(N)]
    sigma_list[i] = sigma(n, a, b)

    return qu.tensor(sigma_list)
