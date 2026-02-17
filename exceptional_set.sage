###################################################################################
# Morikawa's Sangaku Problem
#
# Exceptional set computation for Proposition 3.3
#
# This script constructs the polynomial P(k,y) via a resultant,
# then computes the (squarefree, primitive) discriminant Δ0(k) of P with respect
# to y. The exceptional set E0 is the set of real roots of Δ0 in [1,∞).
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

# Compute the discriminant-like resultant Δ(k) = Res_y(P, ∂P/∂y).
Delta = P.resultant(P.derivative(y),y)
assert Delta.degree(y) == 0

# View Δ as an element of ZZ[k] (drop y).
W.<k> = PolynomialRing(ZZ)
Delta = W(Delta(k=k,y=0))

# Extract the squarefree part and then the primitive part (Δ0 in the paper).
Delta_sq = Delta // Delta.gcd(Delta.derivative())
Delta0 = Delta_sq // Delta_sq.content()

# Remove the explicitly known factors k(k-1)(2k+1)(k^2+1), leaving J(k).
J = Delta0 // (k*(k - 1)*(2*k + 1)*(k^2 + 1))

# Factor J(k) over ZZ and verify the number of irreducible factors.
J_fac = factor(J)
assert len(J_fac) == 19

# Record the degrees of the irreducible factors of J for reporting in the paper.
degrees = sorted(f.degree() for f, e in J_fac)
print(degrees)

# -----------------------------------------------------------------------------
# Numerical approximation of the exceptional set E0:
#   E0 = { real roots of Δ0(k) in [1,∞) }.
# Here we compute real roots of the irreducible factors of J(k) and keep those ≥ 1,
# together with k=1 (which is also a root coming from the explicit factor (k-1)).
# -----------------------------------------------------------------------------

RR = RealField(300)
E0 = [RR(1)]
for (p, e) in J_fac:
    # Work over a high-precision real field to approximate real roots reliably.
    pRR = p.change_ring(RR)
    roots = pRR.roots(multiplicities=False)
    # Filter roots >= 1
    p_roots = [ r for r in roots if r >= 1 ]
    E0.extend(p_roots)

# Remove duplicates (exact equality in RR) and sort.
E0 = sorted(set(E0))

# Expected cardinality from the computation described in the paper.
assert len(E0) == 139

# Smallest element strictly larger than 1, and the maximum element of E0.
E0[1]
max(E0)
