# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

# ------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
"Libraries"


# ------------------------------------------------------------------------------
"Parameters needed"

Pfus = 500
Q = 10
lam = 4.94
kappa = 1.7
eps = 0.3226
delta = 0.33
q = 3
gamma_r = 0.7
n_N = 0.85


Ci = 13.144
Cbet = 0.726
Cfus = 1.38e-3
Closs = 0.086
Csl = 0.0562
Cn = 3.183
M = 2.55
Dint = 1.23
k_B = 1.03
Ncoil = 20


# ------------------------------------------------------------------------------


def Pfus_condition(beta_N, Pfus1=Pfus, kappa1=kappa, eps1=eps, q1=q):
    R3_B4 = Pfus1/(Cfus*Ci**2/Cbet**2*kappa1*eps1**4*beta_N**2/q1**2)
    return R3_B4


def Q_condition(beta_N, Q1=Q, kappa1=kappa, eps1=eps, q1=q, n_N1=n_N, gamma_r1=gamma_r):
    R42_B73 = (Closs*Q1/(gamma_r1*Cfus*(1+Q1/lam)))**0.31/(Csl*Cn**0.41*Ci**0.96
                                                           * Cbet**0.38*Closs**(-0.69)*M**0.19*kappa1**0.09*eps1**0.68*q1**(-0.96)*n_N1**0.41*beta_N**(-0.38))
    return R42_B73


def R(beta_N, Pfus1=Pfus, Q1=Q, kappa1=kappa, eps1=eps, q1=q, n_N1=n_N, gamma_r1=gamma_r):
    x = Pfus_condition(beta_N, Pfus1, kappa1, eps1, q1)
    y = Q_condition(beta_N, Q1, kappa1, eps1, q1, n_N1, gamma_r1)
    return (x**0.73/y**4)**(1/(3*0.73-0.42*4))


def B(beta_N, Pfus1=Pfus, Q1=Q, kappa1=kappa, eps1=eps, q1=q, n_N1=n_N, gamma_r1=gamma_r):
    x = Pfus_condition(beta_N, Pfus1, kappa1, eps1, q1)
    y = Q_condition(beta_N, Q1, kappa1, eps1, q1, n_N1, gamma_r1)
    return (x**0.42/y**3)**(1/(4*0.42-0.73*3))


def outer_radius_1(r2, beta_N = 1.695):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    return (R0*(1+eps)/r2)**np.floor(2*np.pi/3.5*r2) + (re/(R0*(1+eps)))**np.floor(2*np.pi/3.5*r2)

# This function takes r2 AND Ncoil as parameters
def outer_radius_2(r2, Ncoil, beta_N = 1.695):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    return (R0*(1+eps)/r2)**Ncoil + (re/(R0*(1+eps)))**Ncoil

# This function computes the required surface of the cable given its current and B inside the superconductor
def Scable(Icable, beta_N = 1.695):
    # Icable = N_strand * Ic
    # Scable = N_strand * S_strand * (2/0.7) + S_channel
    
    Bmagnet = B(beta_N)*R(beta_N)/re
    Ic = 308.3 - 20.427*Bmagnet
    Nstrand = Icable//Ic + 1
    Scable = Nstrand * (np.pi*0.4**2)*2/0.7 + np.pi*5**2 # mm^2
    return Scable

def Jcable(Icable, beta_N = 1.695):
    return Icable/Scable(Icable, beta_N) # A/mm^2

def Dcable(Icable, beta_N = 1.695):
    return 2*np.sqrt(Scable(Icable, beta_N)/np.pi) # mm

# Total number of turns of the cable
def NTFS(Icable, beta_N = 1.695):
    Bmagnet = B(beta_N)*R(beta_N)/re
    NI_tot = 5*10**6*re*Bmagnet
    
    return NI_tot//Icable +1

# Number of turn of the cable in one TF coil
def NTF(Icable, beta_N = 1.695):
    return  NTFS(Icable, beta_N)//Ncoil +1 

# This functions compute the inner radius of the TF coil

def inner_radius(Icable, beta_N = 1.695):
    # Stot = NTFS * Scable
    # Stot = pi*(re^2 - ri^2)
    re = R(beta_N)*(1-eps) - Dint
    return np.sqrt(re**2 - Scable(Icable, beta_N)*NTFS(Icable, beta_N)/np.pi/1e6)



# -------------------------------------------------------------------------------
"Plots"


def plot_BR(bmin, bmax, npoints=600):
    beta_list = np.linspace(bmin, bmax, npoints)
    B_list = [B(x) for x in beta_list]
    R_list = [R(x) for x in beta_list]
    limit = npoints*[1.695]

    plt.figure()
    plt.plot(beta_list, B_list, label='B')
    plt.plot(beta_list, R_list, label='R')
    plt.plot(limit, R_list, '--', label='ITER')
    plt.xlabel(r'$\beta_N$')
    plt.legend()
    plt.show()


def plot_r2_1(beta_N):
    R0 = R(beta_N)
    x = np.linspace(R0*(1+eps), 2*R0*(1+eps), 100)
    r2_list = [outer_radius_1(r2, beta_N) for r2 in x]
    plt.plot(x, r2_list, label='ripple effect')
    plt.ylim(0, 1e-2)
    plt.show()


# Plots the domain of validity for (r2, Ncoil)
# delta stands for the proximity to 3.5m
def plot_r2_2(beta_N, delta):
    R0 = R(beta_N)

    X = np.arange(R0*(1+eps), 15, 0.01)
    Y = np.arange(10, 26, 0.1)
    x, y = np.meshgrid(X, Y)
    
    # Looking for the couples (r2, Ncoil) that satisfy delta <= 1e-2 and/or (2\pi r2/Ncoil \sim 3.5m)
    f = np.array([[((outer_radius_2(r2, Ncoil, beta_N) <= 1e-2)*2 + (abs(2*np.pi*r2/Ncoil - 3.5) <= delta)) for r2 in X] for Ncoil in Y])
    plt.figure()
    plt.pcolormesh(x, y, f, shading = 'auto')

    plt.xlabel('r2')
    plt.ylabel('N')
    plt.title(r'$(r_2, N)$ that satisfy both dimensioning conditions')
    plt.colorbar()
    plt.grid()
    plt.show()


#------------------------------------------------------------------------------------
"Values"
beta_N = 1.695
R0 = R(beta_N)
B0 = B(beta_N)
re = R0 - eps*R0 - Dint
Bmagnet = k_B * (B0*R0/re)
NI_tot = 5*10**6*re*Bmagnet

Icable = 69e3
Ncoil = 20
r2 = 10.8 #m

S_cable = Scable(Icable, beta_N) # mm^2
J_cable = Jcable(Icable, beta_N) # A/mm^2
D_cable = Dcable(Icable, beta_N)
N_TFS   = NTFS(Icable, beta_N)
N_TF    = NTF(Icable, beta_N)
ri      = inner_radius(Icable, beta_N)
# -----------------------------------------------------------------------------------
"Orders"

# plot_BR(1.5,1.9)



# print('R(beta_N = 1.695) = {}'.format(R0))
# print('B(beta_N = 1.695) = {}'.format(B0))
# print('re = {}'.format(re))
# print('NI_tot = {}'.format(NI_tot))
