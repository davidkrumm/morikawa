###################################################################################
# Morikawa's Sangaku Problem
#
# Proposition 4.1 (Analyticity at k = 1)
#
# This script:
#   (1) constructs the polynomial P(k,y) ∈ ZZ[k,y] satisfied by λ(k);
#   (2) factors P in ZZ[k,y] as P = q*w*z into three irreducible factors and
#       verifies their total degrees are (28, 42, 74);
#   (3) specializes at k = 1 and checks numerically that w(1,y) and z(1,y)
#       have no real roots in (-1,1), so s(1,a) ≠ 0 for s = w*z;
#   (4) specializes q at k = 1, factors off the linear terms 8(y-2)(y-4)(y-8),
#       and isolates the remaining degree-7 factor n(y), verifying that n is
#       irreducible over QQ.
#   (5) constructs the number field K = QQ(a) and computes the Taylor series
#       of λ(1+t) by Newton iteration in K[[t]];
#   (6) embeds the series into a high-precision real field, converts it to a
#       series for μ(1+t) via the substitution k = sqrt(r), and computes
#       rational approximations and numerical evaluations of the truncated
#       Taylor expansion.
#
###################################################################################

# Base polynomial ring over ZZ in three variables (k, x, y).
S.<k,x,y>  = PolynomialRing(ZZ)

# Convenience abbreviations.
c1 = k^2 - x
c2 = k^4
c3 = 2*k - x
c4 = 2*x - x^2

# Polynomials h(k,x,y) and f_i(k,x).
h = -((y - x^2 - c1^2 - c2 + c3^2 + c4)^2 + 4*c3^2*c4 - 4*c1^2*c2 + 4*c1^2*c3^2 + 4*c1^2*c4)^2 + c4*(8*c1^2*c3 + 4*c3*(y - x^2 - c1^2 - c2 + c3^2 + c4))^2
f1 = -2*x + 4*k
f2 = (4*k - 2)*x + k^4 - 4*k^2
f3 = (6*k - 3)*x + k^4 - 2*k^3 - 3*k^2
f4 = -x^2 + 2*x
f5 = 4*x^3 - (2*k^2 + 6*k + 7)*x^2 + (2*k^3 + 3*k^2 + 10*k)*x - 2*k^3
f6 = 8*x^3 + (-4*k^2 - 16)*x^2 + (4*k^3 - 2*k^2 + 6)*x - 4*k^3 + 8*k^2 - 4*k
f7 = (4*k^2 - 16*k)*x^3 + (-k^4 + 4*k^3 - 10*k^2 + 40*k)*x^2 + (2*k^4 - 8*k^3 + 4*k^2 - 20*k + 2)*x + 4*k^2

# Polynomial p(k,x).
p = f4*(f1*f7 + f2*f6 - 2*f3*f5)^2 - (f3^2*f4 + f5^2 - f1*f4*f6 - f2*f7)^2;

# Eliminate x between p(k,x)=0 and h(k,x,y)=0 using the resultant in x.
R = p.resultant(h,x)

# Sanity check: after taking the resultant in x, x should no longer appear.
assert R.degree(x) == 0

# Move into the 2-variable ring ZZ[k,y] by mapping x -> 0 (and keeping k,y).
T.<k,y> = PolynomialRing(ZZ, 2)
R = T(R(k=k,x=0,y=y))

# Remove the known extraneous factor 2^26*(k^2+1)^2 to obtain the primitive
# polynomial P(k,y) used in the paper.
P = R // (2^26*(k^2 + 1)^2)

# Total degree and content of P(k,y).
assert P.total_degree() == 144
assert P.content() == 1

# Factor P in ZZ[k,y]
P_fac = P.factor()

# The factorization of P(k,y) should consist of exactly three irreducible factors.
assert len(P_fac) == 3

# Extract the three factors and verify their total degrees.
q, w, z = [f for (f, e) in P_fac]
assert q.total_degree() == 28
assert w.total_degree() == 42
assert z.total_degree() == 74

# ---------------------------------------------------------------------------
# Specialize at k = 1.
#
# We now examine the factors q(1,y), w(1,y), z(1,y).
# Since λ(1) = a lies in (0,1), it suffices to check that w(1,y) and z(1,y)
# have no real roots in (-1,1), implying s(1,a) ≠ 0 for s = w*z.
# ---------------------------------------------------------------------------

T.<y> = PolynomialRing(QQ)

# Specializations in QQ[y].
q1 = T(q(k=1, y=y))
w1 = T(w(k=1, y=y))
z1 = T(z(k=1, y=y))

# Work numerically at high precision to locate real roots.
RR = RealField(300)
w1 = w1.change_ring(RR)
z1 = z1.change_ring(RR)

# Verify that w(1,y) and z(1,y) have no real roots in (-1,1).
# Since a = λ(1) ≈ 0.1485 lies in this interval, it cannot be a root of w or z.
assert [r for r in w1.roots(multiplicities=False)
        if r < RR(1) and r > RR(-1)] == []

assert [r for r in z1.roots(multiplicities=False)
        if r < RR(1) and r > RR(-1)] == []

# ---------------------------------------------------------------------------
# Factor q(1,y)
# ---------------------------------------------------------------------------

# Remove the explicit linear factors to isolate n(y).
n = q1 // (8*(y - 2)*(y - 4)*(y - 8))

# Verify that n has degree 7 and is irreducible over QQ.
assert n.degree() == 7
assert n.change_ring(QQ).is_irreducible()

# Display the real roots of n to verify that a = λ(1)) is a root.
n_roots = n.change_ring(RR).roots(multiplicities=False)
a = n_roots[0]
print("Approximate a = lambda(1) ≈", a,"\n")

# ---------------------------------------------------------------------------
# Taylor series at k = 1 (Newton iteration in a power series ring)
#
# Goal:
#   Compute the Taylor expansion of λ(1+t) as a power series y(t) ∈ K[[t]],
#   where K = QQ(a) with a = λ(1)).
#   Then embed coefficients into a high-precision real field RR, convert to a
#   series in RR[[t]], and finally obtain μ(1+t) via
#
#       μ(1+t) = sqrt( λ( sqrt(1+t) ) )
#
# ---------------------------------------------------------------------------

# Number field K = QQ(alpha) where alpha is a root of n(y).
K.<alpha> = NumberField(n)

# Real embedding K -> RR sending alpha |-> a (numerical root chosen above).
emb = K.hom([a], RR)

# Coerce q into K[k,y] so we can work over the number field.
qK = q.change_ring(K)
Qky = qK.parent()
kK, yK = Qky.gens()

def NewtonSeriesRoot(Fpoly, Fy, prec, iters=0):
    """
    Newton iteration to find the unique power series root y(t) with y(0)=alpha
    satisfying Fpoly(y(t)) = 0 to O(t^prec).

    INPUT:
        Fpoly : polynomial in Y over S = K[[t]]
        Fy    : derivative of Fpoly w.r.t. Y
        prec  : t-adic precision (compute mod t^prec)
        iters : (optional) number of Newton steps; if 0, choose enough steps
                to stabilize to O(t^prec) under quadratic convergence.

    OUTPUT:
        yser ∈ K[[t]] with yser(0) = alpha + O(t^prec) and Fpoly(yser)=0 mod t^prec.
    """
    S = Fpoly.base_ring()  # S = K[[t]]
    t = S.gen()

    if iters == 0:
        m = 1
        p = 1
        while p < prec:
            p *= 2
            m += 1
        iters = m + 2

    yser = alpha + O(t^prec)

    for _ in range(iters):
        Fval  = Fpoly(yser)
        Fyval = Fy(yser)
        yser  = yser - Fval / Fyval
        yser  = S(yser) + O(t^prec)

    return yser

def TaylorSeriesLambda(prec=40, iters=0):
    """
    Compute the Taylor series y(t) for λ(1+t) in K[[t]] to O(t^prec),
    by solving q(1+t, y(t)) = 0 with y(0)=alpha via Newton iteration.
    """
    S.<t> = PowerSeriesRing(K, prec)   # K[[t]] / (t^prec)
    RY.<Y> = PolynomialRing(S)         # (K[[t]])[Y]

    # Form F(t,Y) = q(1+t, Y)
    Fpoly = RY(qK.subs({kK: K(1) + t, yK: Y}))
    Fy = Fpoly.derivative(Y)

    return NewtonSeriesRoot(Fpoly, Fy, prec, iters=iters)

# Expansion of λ(1+t) in K[[t]].
lambda_1p = TaylorSeriesLambda(prec=40)

# Convert λ(1+t) to a real-coefficient series via the embedding K -> RR.
prec = lambda_1p.prec()
SR.<t> = PowerSeriesRing(RR, prec)
coeffs = lambda_1p.list()  # coefficients in K

lambda_1p_RR = sum(emb(coeffs[i]) * t^i for i in range(len(coeffs))) + O(t^prec)

# Substitute t = sqrt(1+t) - 1 to obtain λ(sqrt(1+t)) as a series in RR[[t]].
u = (1 + t).sqrt() - 1
lambda_sqrt1pt = lambda_1p_RR(u)

# Finally, μ(1+t) = sqrt( λ( sqrt(1+t) ) ).
mu_1p = lambda_sqrt1pt.sqrt()

# Coefficients of μ(1+t) in RR[[t]] and rational approximations with denom ≤ 10^3.
coeffRR = mu_1p.list()  # [c0, c1, ...] in RR
print("initial coefficients of μ(1+t):", coeffRR[:10],"\n")
rat_coeffs = [c.nearby_rational(max_denominator=1000) for c in coeffRR]
print("rational approximations (denom <= 10^3):")
print(rat_coeffs[:10],"\n")

# Numeric checks for Table 1. Use truncated Taylor series for μ(1+t).
mu_poly = sum(mu_1p[i] * t^i for i in range(mu_1p.prec()))

print("mu(1+0.1) ≈", mu_poly(RR(0.1)),"\n")
print("mu(1+0.2) ≈", mu_poly(RR(0.2)),"\n")
print("mu(1+0.5) ≈", mu_poly(RR(0.5)),"\n")
print("mu(1+1.0) ≈", mu_poly(RR(1.0)),"\n")
print("mu(1+1.2) ≈", mu_poly(RR(1.2)))
