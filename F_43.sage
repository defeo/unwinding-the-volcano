# -*- coding: utf-8 -*-

def frob(P, Q, q):
    'Computes the matrix of q-power Frobenius in basis P, Q'
    o = P.order()
    return matrix([next([i, j]
                        for i in Zmod(o) for j in Zmod(o)
                        if not (i == 0 and j == 0) and (i.lift()*P + j.lift()*Q).xy() == point )
                   for point in ((P[0]^q, P[1]^q), (Q[0]^q, Q[1]^q))]).transpose()

def change_basis(M, P, Q):
    'Compute the matrix-vector product M * (P,Q)^t'
    M = M.lift()
    return M[0,0]*P + M[0,1]*Q, M[1,0]*P + M[1,1]*Q

def elem_step(M, P, Q, q):
    "Jérôme's elementary step (I believe)"
    P, Q = change_basis(M, P, Q)
    return frob(P.division_points(2)[0], Q, q) 


K.<z> = GF(43^2)

# E1 has order 40 over GF(43), with 4 rational 2-torsion points
E1 = EllipticCurve(K, (30, 10))

# A basis of E1[4] with Frobenius matrix [[1 0] [0 3]]
P1 = E1(3, 16)
Q1 = E1(10, 14*z + 36)
#P1 = P1 + 2*Q1

# Image curve and basis by a degree 5 isogeny
psi = E1.isogeny(E1(0).division_points(5))
E2 = psi.codomain()
P2 = psi(P1)
Q2 = psi(Q1)

# The matrix of Frobenius
Phi = frob(P1, Q1, 43)

# The list of candidates for psi
Psi = [m for m in Phi.parent() if m.is_invertible() and Phi*m == m*Phi]

# Push one step along the crater
alpha = E1.isogeny(2*P1)
E1a = alpha.codomain()
P1a = alpha(P1).division_points(2)[0]
Q1a = alpha(Q1)
Phia = frob(P1a, Q1a, 43)

alpha = E2.isogeny(2*P2)
E2a = alpha.codomain()
P2a = alpha(P2).division_points(2)[0]
Q2a = alpha(Q2)

# The new list of candidates for psi
Psia = filter(lambda m : elem_step(m, 2*P2a, Q2a, 43) == Phia, Psi)

# Push one more step along the crater (now we have to use P1a + Q1a,
# because the Frobenius matrix is [[1 0] [2 3]]
alpha = E1a.isogeny(2*(P1a + Q1a))
E1b = alpha.codomain()
P1b = alpha(P1a + Q1a).division_points(2)[0]
Q1b = alpha(Q1a)
Phib = frob(P1b, Q1b, 43)

alpha = E2a.isogeny(2*(P2a + Q2a))
E2b = alpha.codomain()
P2b = alpha(P2a + Q2a).division_points(2)[0]
Q2b = alpha(Q2a)

# The new list of candidates for psi
Psib = filter(lambda m : elem_step(m, 2*P2b, Q2b, 43) == Phib, Psia)
