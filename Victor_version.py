# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

# ------------------------------------------------------------------------------
"Libraries"

import matplotlib.pyplot as plt
import numpy as np
import scipy.special as scs


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
beta_N = 1.685

mu_0 = 4*np.pi*10**(-7)
Ci = 13.144
Cbet = 0.726
Cfus = 1.38e-3
Closs = 0.086
Csl = 0.0562
Cn = 3.183
M = 2.55
Dint = 1.23
k_B = 1.03


Icable = 68e3
#r2 = 10.6 #m censé être 10.4 m
Ncoil = 18
haa = 6.636 #m
t_insu = 1e-3 # m, isolation additional radius

sigma_max      = 600e6 # Pa



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


def outer_radius_1(r2, beta_N = beta_N):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    return (R0*(1+eps)/r2)**np.floor(2*np.pi/3.5*r2) + (re/(R0*(1+eps)))**np.floor(2*np.pi/3.5*r2)

# This function takes r2 AND Ncoil as parameters
def outer_radius_2(r2, Ncoil, beta_N = beta_N):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    return (R0*(1+eps)/r2)**Ncoil + (re/(R0*(1+eps)))**Ncoil

# This function computes the required surface of the cable given its current and B inside the superconductor
def Scable(Icable = Icable, beta_N = beta_N):
    # Icable = N_strand * Ic
    # Scable = N_strand * S_strand * (2/0.7) + S_channel
    
    Bmagnet = k_B *B(beta_N)*R(beta_N)/re
    Ic = 308.3 - 20.427*Bmagnet
    Nstrand = Icable//Ic + 1
    Scable = Nstrand * (np.pi*0.4**2)*2*1.3 + np.pi*5**2 # mm^2
    return Nstrand, Scable

def Jcable(Icable = Icable, beta_N = beta_N):
    return Icable/Scable(Icable, beta_N)[1] # A/mm^2

def Dcable(Icable = Icable, beta_N = beta_N):
    return 2*np.sqrt(Scable(Icable, beta_N)[1]/np.pi) # mm

# Total number of turns of the cable
def NTFS(Icable = Icable, beta_N = beta_N):
    Bmagnet = k_B *B(beta_N)*R(beta_N)/re
    NI_tot = 5*10**6*re*Bmagnet
    
    return NI_tot//Icable +1

# Number of turn of the cable in one TF coil
def NTF(Icable = Icable, beta_N = beta_N):
    return  NTFS(Icable, beta_N)//Ncoil +1 

# This functions compute the inner radius of the TF coil
def inner_radius(Icable = Icable, beta_N = beta_N):
    # Stot = NTFS * Scable
    # Stot = pi*(re^2 - ri^2)
    re = R(beta_N)*(1-eps) - Dint
    return np.sqrt(re**2 - Scable(Icable, beta_N)[1]*NTFS(Icable, beta_N)/np.pi/1e6)


# Inductance of the TF coils (all of them)
def L_TFS(Icable = Icable,  beta_N = beta_N):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    k = np.log(r2/re)/2
    L = mu_0*np.sqrt(re*r2)*(NTFS(Icable, beta_N)*k)**2/2*(scs.i0(k) + 2*scs.i1(k) + scs.iv(2, k))
    return L

# Total length of the cable in the TF coils (in meters)
def Total_length_cable(Icable = Icable, beta_N = beta_N):
    re = R(beta_N)*(1-eps) - Dint
    r2 = R2(beta_N = beta_N)
    k = np.log(r2/re)/2
    P_TF = 2*np.pi*np.sqrt(re*r2)*k*(scs.i0(k)+scs.i1(k))
    return P_TF * NTFS(Icable, beta_N)

def sigma_tot(thick, Icable = Icable,  beta_N = beta_N):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    ri = inner_radius(Icable, beta_N)
    Bmagnet = k_B *B(beta_N)*R(beta_N)/re
    NI_tot = 5*10**6*re*Bmagnet
    r2 = R2(beta_N = beta_N)
    
    sigma_hoop   = Bmagnet**2*(re+ri)**2/(8*mu_0*haa*ri)*np.log(2*r2/(re+ri))
    #print(r'sigma_hoop = {} MPa'.format(sigma_hoop/10**6))
    sigma_center = (NI_tot*(Bmagnet/2)/(2*np.pi*ri))*ri/thick
    return sigma_center + sigma_hoop

def inner_radius_nose(Icable = Icable,  beta_N = beta_N):
    thick = 0.2
    while sigma_tot(thick, Icable,  beta_N)/sigma_max > 1:
        thick += 1e-3
    rpi = inner_radius(Icable = Icable, beta_N = beta_N) - thick
    return thick, rpi

def R2(beta_N = beta_N, delta = 0.25):
    R0 = R(beta_N)
    r2 = 10.5
    while outer_radius_2(r2, Ncoil, beta_N) > 1e-2 or abs(2*np.pi*r2/Ncoil - 3.5) > delta:
        r2 += 1e-3
    return r2


#------------------------------------------------------------------------------
"Functions for the presentation"

# Gives the different radii (in m) that were found during this study
# WARNING : thickness of the isolating components not taken into account yet  
def camembert_radii(Icable = Icable,  beta_N = beta_N):
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    ri = inner_radius(Icable, beta_N) 
    thickness, rpi = inner_radius_nose(Icable,  beta_N)
    r2 = R2(beta_N, 0.25)
    return rpi, ri, re, R0, r2

# Gives the different surfaces (in mm^2) that were found during this study 
# WARNING : thickness of the isolating components not taken into account yet   
def camembert_surfaces(Icable = Icable,  beta_N = beta_N):
    S_trand = np.pi*0.4**2
    S_superconductor = S_trand/4
    
    N_strand, S_cable = Scable(Icable, beta_N)
    S_Cu_Nb3Sn = 2*S_trand*N_strand
    
    # Cross-section of the cables inside one TF coil
    S_TF  = NTF(Icable, beta_N)*S_cable
    
    # Cross-section of the cables inside all TF coils
    S_TFS = NTFS(Icable, beta_N)*S_cable
    
    # Cross-section at the intern leg of the TF coils
    R0 = R(beta_N)
    re = R0 - eps*R0 - Dint
    thickness, rpi = inner_radius_nose(Icable,  beta_N)
    Stot = np.pi*(re**2 - rpi**2)*10**6
    
    return S_superconductor, S_trand, S_Cu_Nb3Sn, S_cable, S_TF, S_TFS, Stot


# Gives the different critical current density Jc (in A/mm^2) that were found during this study
# WARNING : thickness of the isolating components not taken into account yet  
def camembert_Jc(Icable = Icable,  beta_N = beta_N):
    R0 = R(beta_N)
    B0 = B(beta_N)
    re = R0 - eps*R0 - Dint
    thickness, rpi = inner_radius_nose(Icable,  beta_N)
    N_strand, S_cable = Scable(Icable, beta_N)
    Stot = np.pi*(re**2 - rpi**2)*10**6 # mm^2
    
    magnet = k_B * (B0*R0/re)
    Ic_strain = 308.3 - 20.427*Bmagnet
    
    Jc_strain = Ic_strain/(np.pi*0.4**2) # For the wire, A/mm^2
    Jc_super  = Jc_strain*4 # For the superconductor itself, A/mm^2
    Jc_Cu     = Jc_strain/2 # For Cu + strain, A/mm^2
    Jc_cable  = Jc_strain/2.6 # For the cable without inox 
    Jc_inox   = Jc_cable*S_cable/(Stot/NTFS(Icable, beta_N)) # For the cable with inox
    
    return Jc_super, Jc_strain, Jc_Cu, Jc_cable, Jc_inox
    
    
    
    

# R0 = R(beta_N)
# B0 = B(beta_N)
# re = R0 - eps*R0 - Dint
# Bmagnet = k_B * (B0*R0/re)
# NI_tot = 5*10**6*re*Bmagnet

# N_strand, S_cable = Scable(Icable, beta_N) # mm^2
# J_cable = Jcable(Icable, beta_N) # A/mm^2
# D_cable = Dcable(Icable, beta_N) # mm
# N_TFS   = NTFS(Icable, beta_N)
# N_TF    = NTF(Icable, beta_N)
# ri      = inner_radius(Icable, beta_N) # m

# # The thickness of the inox nose is found to be 0.27 m approximately
# thickness, rpi = inner_radius_nose(Icable,  beta_N)

# # Values of the different critical current densities
# Ic_strain = 308.3 - 20.427*Bmagnet

# Jc_strain = Ic_strain/(np.pi*0.4**2) # For the wire, A/mm^2
# Jc_super  = Jc_strain*4 # For the superconductor itself, A/mm^2
# Jc_Cu     = Jc_strain/2 # For Cu + strain, A/mm^2
# Jc_cable  = Jc_strain/2.6 # For the cable without inox

# # Value of the mean Laplace force inside the intern leg of the TF coils

# Fl_cable = Icable*(Bmagnet/2) # N/m

# Fl_tot   = NI_tot*(Bmagnet/2) # N/m

# Fl_coil = Fl_tot/Ncoil # N/m

# # Energy contained inside the whole TF coils
# E_TFS = L_TFS(Icable, beta_N)*Icable**2/2 # in J


# # total length of the cable
# L_cable = Total_length_cable(Icable, beta_N)


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


def plot_r2_1(beta_N = beta_N):
    R0 = R(beta_N)
    x = np.linspace(R0*(1+eps), 2*R0*(1+eps), 100)
    r2_list = [outer_radius_1(r2, beta_N) for r2 in x]
    
    plt.figure()
    plt.plot(x, r2_list, label='ripple effect')
    plt.ylim(0, 1e-2)
    plt.legend()
    plt.show()


# Plots the domain of validity for (r2, Ncoil)
# delta stands for the proximity to 3.5m
def plot_r2_2(beta_N = beta_N, delta = 0.25):
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

def plot_sigma_tot(beta_N = beta_N, Icable = Icable, npoints = 500):
    thick_list = np.logspace(-3, 2, npoints)
    sigma_list = [sigma_tot(thick, beta_N, Icable)/sigma_max for thick in thick_list]
    
    plt.figure()
    plt.plot(thick_list, sigma_list, label = r'$\sigma / \sigma_{max}$')
    plt.plot(thick_list, npoints*[1])
    plt.xscale('log')
    plt.yscale('log')
    plt.title(r'Value of $\sigma / \sigma_{max}$')
    plt.legend()
    plt.show()
    


#------------------------------------------------------------------------------------
"Values"

R0 = R(beta_N)
B0 = B(beta_N)
re = R0 - eps*R0 - Dint
Bmagnet = k_B * (B0*R0/re)
NI_tot = 5*10**6*re*Bmagnet

N_strand, S_cable = Scable(Icable, beta_N) # mm^2
J_cable = Jcable(Icable, beta_N) # A/mm^2
D_cable = Dcable(Icable, beta_N) # mm
N_TFS   = NTFS(Icable, beta_N)
N_TF    = NTF(Icable, beta_N)
ri      = inner_radius(Icable, beta_N) # m
r2      = R2(beta_N, 0.25)

# The thickness of the inox nose is found to be 0.27 m approximately
thickness, rpi = inner_radius_nose(Icable,  beta_N)

# Values of the different critical current densities
Ic_strain = 308.3 - 20.427*Bmagnet

Jc_strain = Ic_strain/(np.pi*0.4**2) # For the wire, A/mm^2
Jc_super  = Jc_strain*4 # For the superconductor itself, A/mm^2
Jc_Cu     = Jc_strain/2 # For Cu + strain, A/mm^2
Jc_cable  = Jc_strain/2.6 # For the cable without inox

# Value of the mean Laplace force inside the intern leg of the TF coils

Fl_cable = Icable*(Bmagnet/2) # N/m

Fl_tot   = NI_tot*(Bmagnet/2) # N/m

Fl_coil = Fl_tot/Ncoil # N/m

# Energy contained inside the whole TF coils
E_TFS = L_TFS(Icable, beta_N)*Icable**2/2 # in J


# total length of the cable
L_cable = Total_length_cable(Icable, beta_N)

# -----------------------------------------------------------------------------------
"Orders"

# plot_BR(1.5,1.9)



# print('R(beta_N = 1.695) = {}'.format(R0))
# print('B(beta_N = 1.695) = {}'.format(B0))
# print('re = {}'.format(re))
# print('NI_tot = {}'.format(NI_tot))
