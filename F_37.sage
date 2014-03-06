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

def elem_step(M, h, P, i, Q, q):
    "Jérôme's elementary step (I believe)"
    P, Q = change_basis(M, h*P, i*Q)
    return frob(P.division_points(h)[0], Q.division_points(i)[0], q) 


K.<z> = GF(37^4)

# E1 has order 32 over GF(37), with 16 rational 4-torsion points
# It is the only curve on the crater of a volcano of height 2.
# It has 2 horizontal isogenies to itself. Nice :)
E1 = EllipticCurve(K, (23, 7))

# # A basis of E1[8] with Frobenius matrix [[1 0] [0 5]]
# P1 = E1(0, 9)
# #Q1 = E1(2*z + 12, 10*z + 34)
# Q1 = E1(12*z^3 + 2*z^2 + 11*z + 16, 14*z^3 + 27*z^2 + 19*z + 20)

# # # The matrix of Frobenius
# Phi = frob(P1, Q1, 37)

# # # The list of candidates for psi
# Psi = [m for m in Phi.parent() if m.is_invertible() and Phi*m == m*Phi]

# # Push one step along the crater
# alpha_a = E1.isogeny(4*P1)
# E1a = alpha_a.codomain()
# P1a = alpha_a(P1).division_points(2)[0]
# Q1a = alpha_a(Q1)
# Phia = frob(P1a, Q1a, 37)

# # The new list of candidates for psi
# Psia = filter(lambda m : elem_step(m, 2, P1a, 1, Q1a, 37) == Phia, Psi)

# # Push one more step along the crater 
# alpha_b = E1a.isogeny(4*P1a)
# E1b = alpha_b.codomain()
# P1b = alpha_b(P1a).division_points(2)[0]
# Q1b = alpha_b(Q1a)
# Phib = frob(P1b, Q1b, 37)

# # The new list of candidates for psi
# Psib = filter(lambda m : elem_step(m, 4, P1b, 1, Q1b, 37) == Phib, Psia)

# # We go back now
# beta_b = E1b.isogeny(4*Q1b)
# E1c = beta_b.codomain()
# P1c = beta_b(P1b)
# Q1c = beta_b(Q1b).division_points(2)[0]
# Phic = frob(P1c, Q1c, 37)

# # The new list of candidates for psi
# Psic = filter(lambda m : elem_step(m, 1, P1c, 2, Q1c, 37) == Phic, Psib)

# # One last step
# beta_a = E1c.isogeny(4*Q1c)
# E1d = beta_a.codomain()
# P1d = beta_a(P1c)
# Q1d = beta_a(Q1c).division_points(2)[0]
# Phid = frob(P1d, Q1d, 37)

# # The new list of candidates for psi
# Psid = filter(lambda m : elem_step(m, 1, P1d, 4, Q1d, 37) == Phid, Psic)


################### Second Phase

# We need one more torsion power
PP1 = P1.division_points(2)[0]
QQ1 = Q1.division_points(2)[0]

# We descend one level below the floor
gamma = E1.isogeny(8*(PP1 + QQ1))
sE1 = gamma.codomain()
sP1 = gamma(PP1)
sQ1 = gamma(PP1+QQ1)
Phie = frob(2*sP1, sQ1, 37)

# And we filter 
Psie = filter(lambda m : elem_step(m, 1, 2*sP1, 1, sQ1, 37) == Phie, Psid)

# And now?
PP1a = alpha_a(PP1).division_points(2)[0]
QQ1a = alpha_a(QQ1)
gamma_a = E1a.isogeny(8*(PP1a + QQ1a))
sE1a = gamma_a.codomain()
sP1a = gamma_a(PP1a)
sQ1a = gamma_a(PP1a+QQ1a)
Phif = frob(2*sP1a, sQ1a, 37)

Psif = filter(lambda m : elem_step(m, 2, 2*sP1a, 1, sQ1a, 37) == Phif, Psie)
