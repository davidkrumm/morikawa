This repository contains the complete source code used in the paper "Algebraic and analytic structure of Morikawa’s sangaku problem."

Repository Contents:

1. numerical_mu.sage
Implements high-precision numerical evaluation of the function
\mu(r) and the auxiliary function \lambda(k)=\mu(k^2)^2 using minimization of the defining function z_r(x).

2. exceptional_set.sage
Constructs the defining polynomial P(k,y), computes its discriminant
\Delta(k)=\mathrm{Res}_y(P,P_y), and determines the finite exceptional set in [1,\infty) where the real-analytic implicit function theorem does not directly apply.

3. analytic_at_k=1.sage
Analyzes the specialization at k=1, factors the relevant irreducible component, and computes the Taylor expansion of \lambda (and hence \mu) at r=1 using Newton iteration in a power series ring.

