###################################################################################
# Morikawa's Sangaku Problem
#
# Numerical evaluation of μ(r) and λ(k) = μ(k^2)^2.
#
# This script:
#   (1) defines a numerical routine for μ(r), the minimal side length of a
#       square inscribed in the curvilinear triangular region determined by
#       two tangent circles of radii 1 and r ≥ 1 and their common tangent line;
#   (2) computes μ(r) by minimizing the function z_r(x) on the interval
#       [1 - 1/sqrt(2), 1] using Sage's find_local_minimum;
#   (3) defines λ(k) = μ(k^2)^2;
#   (4) prints sample numerical values of λ(1) and μ(r) for selected r.
#
###################################################################################

from sage.all import sqrt, find_local_minimum, RR

def mu(r):
    """
    Numerical approximation of μ(r) for r ≥ 1.
    """

    if r < 1:
        raise ValueError("mu(r) is defined for r ≥ 1.")

    r = RR(r)

    # z_r(x) whose minimum on [1 - 1/sqrt(2), 1] equals μ(r)
    def z(x):
        inner = 2*sqrt(r) - x - sqrt(2*x - x^2)
        return sqrt(x^2 + (r - x - sqrt(r^2 - inner^2))^2)

    mu_min,_ = find_local_minimum(z, 1 - 1/sqrt(2), 1)
    return mu_min


def Lambda(k):
    """
    Numerical approximation of λ(k) = μ(k^2)^2.
    """

    return mu(k^2)^2


RR = RealField(300)

print("lambda(1) ≈", Lambda(RR(1)),"\n")
print("mu(1.1) ≈", mu(RR(1.1)),"\n")
print("mu(1.2) ≈", mu(RR(1.2)),"\n")
print("mu(1.5) ≈", mu(RR(1.5)),"\n")
print("mu(2.0) ≈", mu(RR(2.0)),"\n")
print("mu(2.2) ≈", mu(RR(2.2)))
