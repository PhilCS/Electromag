# Auteurs : Philippe Carpentier-Savard, Guillaume Vollant-Boulé

from numpy import *
import matplotlib.pyplot as plt

def capacite_2fils(N, L, s, a, d, theta1=0, theta2=0):
	"""
	:param int N: nombre d'éléments par fil
	:param float L: longueur de chaque fil [m]
	:param float s: distance entre l'extrémité gauche des 2 fils [m]
	:param float a: rayon de chaque fil [m]
	:param float d: décalage en x de l'extrémité gauche du fil 1
	:param float theta1: angle du fil 1, pivot à l'extrémité gauche [deg]
	:param float theta2: angle du fil 2, pivot à l'extrémité gauche [deg]
	"""

	offset1 = (d, s)  # décalage (x,y) de l'extrémité gauche du fil 1
	offset2 = (0.0, 0.0)  # décalage (x,y) de l'extrémité gauche du fil 2

	theta1 = deg2rad(theta1)
	theta2 = deg2rad(theta2)

	epsilon_0 = 8.8542*10**-12
	N2 = N*2
	h = L/N  # Delta

	segments = empty(N2, dtype=ndarray) # positions [x,y] de tous les segments des 2 fils

	for i in range(0,N):
		dist_i = (i+0.5)*h*L  # distance entre le centre du segment i et l'extrémité gauche du fil
		segments[i] = add([dist_i*cos(theta1), dist_i*sin(theta1)], offset1)  # position [x,y] du segment i du fil 1
		segments[i+N] = add([dist_i*cos(theta2), dist_i*sin(theta2)], offset2)  # position [x,y] du segment i du fil 2

	A = empty((N2,N2))
	A_ii = 2*log(h/a)

	for i in range(0,N2):
		for j in range(0,N2):
			if i == j:
				A[i,j] = A_ii  # élément de la trace
			else:
				A[i,j] = h / linalg.norm(segments[i] - segments[j])

	V = full((N2,1), 4*pi*epsilon_0)

	for i in range(N,N2):
		V[i] *= -1

	rho = dot(linalg.inv(A), V)
	C = h * sum(rho[0:N,:]) / 2  # divisé par 2 puisque la diff. de potentiel entre les fils est 2*V0

	return C, rho
####

q2a = capacite_2fils(N=10, L=1.0, s=0.5, a=0.001, d=0.0)[0]
q2b = capacite_2fils(N=20, L=1.0, s=0.5, a=0.001, d=0.0)[0]
q3_C, q3_rho = capacite_2fils(N=10, L=1.0, s=0.5, a=0.001, d=0.4)
q4 = capacite_2fils(N=10, L=1.0, s=0.01, a=0.001, d=0.0, theta1=45.0)[0]

print("Q2_10 =", q2a)
print("Q2_20 =",q2b)
print("Q3 =",q3_C)
print("Q4 =",q4)

# Figures

N_range = arange(1,11,1)
fig = plt.figure()

# Figure fil du haut
fil1 = fig.add_subplot(211)
fil1.plot(N_range, q3_rho[0:10,0], marker="o")
fil1.set_title("Fil 1")
fil1.set_xticks(N_range)

# Figure fil du bas
fil2 = fig.add_subplot(212)
fil2.plot(N_range, q3_rho[10:20,0], marker="o")
fil2.set_title("Fil 2")
fil2.set_xticks(N_range)
fil2.set_xlabel("N", fontsize=14)

fig.text(0.04, 0.51, "ρ", style="italic", fontsize=14)
plt.tight_layout()
plt.show()

