import numpy as np

## All values are in SI units
M_sun = 1.989e30
R_sun = 695700000.0
L_sun = 3.846e26
G = 6.67408e-11
gamma = 5.0/3.0
c = 299792458.0
sigma = 5.670367e-8
a = 4.0*sigma/c
m_proton = 1.6726219e-27
m_electron = 9.10938356e-31
k = 1.38064852e-23
h = 6.62607004e-34
hbar = 1.0545718e-34
X = 0.73
Y = 0.25
Z = 0.02

def epsilon(rho, T):
	epsilon_PP = 1.07e-7*(rho/10.0**5.0)*X**2.0*(T/10.0**6.0)**4.0
	X_CNO = 0.03*X
	epsilon_CNO = 8.24e-26*(rho/10.0**5.0)*X*X_CNO*(T/10.0**6.0)**19.9
	epsilon = epsilon_PP + epsilon_CNO
	return epsilon

def kappa(rho, T):
	kappa_es = 0.02*(1.0 + X)
	kappa_ff = 10.0**24.0*(Z + 0.0001)*(rho/10.0**3.0)**0.7*T**(-3.5)
	kappa_H = 2.5e-32*(Z/0.02)*(rho/10.0**3.0)**0.5*T**9.0

	if T > 10000.0:
		kappa = (1.0/kappa_H + 1.0/max(kappa_es, kappa_ff))**(-1.0)
	else:
		kappa = (1.0/kappa_H + 1.0/min(kappa_es, kappa_ff))**(-1.0)
	return kappa

def kappa_es():
	kappa_es = 0.02*(1.0 + X)
	return kappa_es

def kappa_ff(rho, T):
	kappa_ff = 10.0**24.0*(Z + 0.0001)*(rho/10.0**3.0)**0.7*T**(-3.5)
	return kappa_ff

def kappa_H(rho, T):
	kappa_H = 2.5e-32*(Z/0.02)*(rho/10.0**3.0)**0.5*T**9.0
	return kappa_H

def drho_dr(rho, T, M, r, dT_dr):
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	dP_dT = rho*k/mu/m_proton + 4.0*a*T**3.0/3.0
	dP_drho = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(2.0/3.0)/3.0/m_electron/m_proton + k*T/mu/m_proton
	drho_dr = -(G*M*rho/r**2.0 + dP_dT*dT_dr)/dP_drho
	return drho_dr

def drho_dr_3(rho, T, M, r, dT_dr, lam):
	g = G*M/r**2.0*(1.0 + lam/r)
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	dP_dT = rho*k/mu/m_proton + 4.0*a*T**3.0/3.0
	dP_drho = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(2.0/3.0)/3.0/m_electron/m_proton + k*T/mu/m_proton
	drho_dr = -(g*rho + dP_dT*dT_dr)/dP_drho
	return drho_dr

def drho_dr_1(rho, T, M, r, dT_dr, LAM):
	g = G*M/r**2.0*(1.0 + r/LAM)
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	dP_dT = rho*k/mu/m_proton + 4.0*a*T**3.0/3.0
	dP_drho = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(2.0/3.0)/3.0/m_electron/m_proton + k*T/mu/m_proton
	drho_dr = -(g*rho + dP_dT*dT_dr)/dP_drho
	return drho_dr

def dP_dr(rho, M, r):
	dP_dr = -G*M*rho/r**2.0
	return dP_dr

def dT_dr(rho, T, M, L, r, kappa):
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	P = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(5.0/3.0)/5.0/m_electron + rho*k*T/mu/m_proton + a*T**4.0/3.0
	dT_dr = -min(3.0*kappa*rho*L/16.0/np.pi/a/c/T**3.0/r**2.0, (1.0 - 1.0/gamma)*T*G*M*rho/P/r**2.0)
	return dT_dr

def dT_dr_3(rho, T, M, L, r, kappa, lam):
	g = G*M/r**2.0*(1.0 + lam/r)
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	P = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(5.0/3.0)/5.0/m_electron + rho*k*T/mu/m_proton + a*T**4.0/3.0
	dT_dr = -min(3.0*kappa*rho*L/16.0/np.pi/a/c/T**3.0/r**2.0, (1.0 - 1.0/gamma)*T*g*rho/P)
	return dT_dr

def dT_dr_1(rho, T, M, L, r, kappa, LAM):
	g = G*M/r**2.0*(1.0 + r/LAM)
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	P = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(5.0/3.0)/5.0/m_electron + rho*k*T/mu/m_proton + a*T**4.0/3.0
	dT_dr = -min(3.0*kappa*rho*L/16.0/np.pi/a/c/T**3.0/r**2.0, (1.0 - 1.0/gamma)*T*g*rho/P)
	return dT_dr

def P_deg(rho):
	P_deg = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(5.0/3.0)/5.0/m_electron
	return P_deg

def P_gas(rho, T):
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	P_gas = rho*k*T/mu/m_proton
	return P_gas

def P_phot(T):
	P_phot = a*T**4.0/3.0
	return P_phot

def P_tot(rho, T):
	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
	P_tot = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(rho/m_proton)**(5.0/3.0)/5.0/m_electron + rho*k*T/mu/m_proton + a*T**4.0/3.0
	return P_tot

def dM_dr(rho, r):
	dM_dr = 4.0*np.pi*r**2.0*rho
	return dM_dr

def dL_dr(rho, r, epsilon):
	dL_dr = 4.0*np.pi*r**2.0*rho*epsilon
	return dL_dr

def dL_dr_pp(rho, T, r):
	epsilon_pp = 1.07e-7*(rho/10.0**5.0)*X**2.0*(T/10.0**6.0)**4.0
	dL_dr_pp = 4.0*np.pi*r**2.0*rho*epsilon_pp
	return dL_dr_pp

def dL_dr_CNO(rho, T, r):
	X_CNO = 0.03*X
	epsilon_CNO = 8.24e-26*(rho/10.0**5.0)*X*X_CNO*(T/10.0**6.0)**19.9
	dL_dr_CNO = 4.0*np.pi*r**2.0*rho*epsilon_CNO
	return dL_dr_CNO

def dtau_dr(rho, kappa):
	dtau_dr = kappa*rho
	return dtau_dr

def k_array(rho, T, M, L, r, dr): #h is step size
	k_array = np.empty((5, 4))

	rho1 = drho_dr(rho, T, M, r, dT_dr(rho, T, M, L, r, kappa(rho, T)))
	T1 = dT_dr(rho, T, M, L, r, kappa(rho, T))
	M1 = dM_dr(rho, r)
	L1 = dL_dr(rho, r, epsilon(rho, T))
	tau1 = dtau_dr(rho, kappa(rho, T))
	k_array[:,0] = [rho1, T1, M1, L1, tau1]

	rho2 = drho_dr(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, r + dr/2.0, dT_dr(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, L + dr*L1/2.0, r + dr/2.0, kappa(rho + dr*rho1/2.0, T + dr*T1/2.0)))
	T2 = dT_dr(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, L + dr*L1/2.0, r + dr/2.0, kappa(rho + dr*rho1/2.0, T + dr*T/2.0))
	M2 = dM_dr(rho + dr*rho1/2.0, r + dr/2.0)
	L2 = dL_dr(rho + dr*rho1/2.0, r + dr/2.0, epsilon(rho + dr*rho1/2.0, T + dr*T1/2.0))
	tau2 = dtau_dr(rho + dr*rho1/2.0, kappa(rho + dr*rho1/2.0, T + dr*T1/2.0))
	k_array[:,1] = [rho2, T2, M2, L2, tau2]

	rho3 = drho_dr(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*M2, r + dr/2.0, dT_dr(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*T2, L + dr/2.0*L2, r + dr/2.0, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2)))
	T3 = dT_dr(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*M2, L + dr*L2/2.0, r + dr/2.0, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2))
	M3 = dM_dr(rho + dr/2.0*rho2, r + dr/2.0)           
	L3 = dL_dr(rho + dr/2.0*rho2, r + dr/2.0, epsilon(rho + dr/2.0*rho2, T + dr/2.0*T2))
	tau3 = dtau_dr(rho + dr/2.0*rho2, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2))
	k_array[:,2] = [rho3, T3, M3, L3, tau3]

	rho4 = drho_dr(rho + dr*rho3, T + dr*T3, M + dr*M3, r + dr, dT_dr(rho + dr*rho3, T + dr*T3, M + dr*M3, L + dr*L3, r + dr, kappa(rho + dr*rho3, T + dr*T3)))
	T4 = dT_dr(rho + dr*rho3, T + dr*T3, M + dr*M3, L + dr*L3, r + dr, kappa(rho + dr*rho3, T + dr*T))
	M4 = dM_dr(rho + dr*rho3, r + dr)
	L4 = dL_dr(rho + dr*rho3, r + dr, epsilon(rho + dr*rho3, T + dr*T3))
	tau4 = dtau_dr(rho + dr*rho3, kappa(rho + dr*rho3, T + dr*T3))
	k_array[:,3] = [rho4, T4, M4, L4, tau4]

	return k_array	

def k_array_3(rho, T, M, L, r, dr, lam): #h is step size
	k_array = np.empty((5, 4))

	rho1 = drho_dr_3(rho, T, M, r, dT_dr_3(rho, T, M, L, r, kappa(rho, T), lam), lam)
	T1 = dT_dr_3(rho, T, M, L, r, kappa(rho, T), lam)
	M1 = dM_dr(rho, r)
	L1 = dL_dr(rho, r, epsilon(rho, T))
	tau1 = dtau_dr(rho, kappa(rho, T))
	k_array[:,0] = [rho1, T1, M1, L1, tau1]

	rho2 = drho_dr_3(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, r + dr/2.0, dT_dr_3(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, L + dr*L1/2.0, r + dr/2.0, kappa(rho + dr*rho1/2.0, T + dr*T1/2.0), lam), lam)
	T2 = dT_dr_3(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, L + dr*L1/2.0, r + dr/2.0, kappa(rho + dr*rho1/2.0, T + dr*T/2.0), lam)
	M2 = dM_dr(rho + dr*rho1/2.0, r + dr/2.0)
	L2 = dL_dr(rho + dr*rho1/2.0, r + dr/2.0, epsilon(rho + dr*rho1/2.0, T + dr*T1/2.0))
	tau2 = dtau_dr(rho + dr*rho1/2.0, kappa(rho + dr*rho1/2.0, T + dr*T1/2.0))
	k_array[:,1] = [rho2, T2, M2, L2, tau2]

	rho3 = drho_dr_3(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*M2, r + dr/2.0, dT_dr_3(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*T2, L + dr/2.0*L2, r + dr/2.0, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2), lam), lam)
	T3 = dT_dr_3(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*M2, L + dr*L2/2.0, r + dr/2.0, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2), lam)
	M3 = dM_dr(rho + dr/2.0*rho2, r + dr/2.0)           
	L3 = dL_dr(rho + dr/2.0*rho2, r + dr/2.0, epsilon(rho + dr/2.0*rho2, T + dr/2.0*T2))
	tau3 = dtau_dr(rho + dr/2.0*rho2, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2))
	k_array[:,2] = [rho3, T3, M3, L3, tau3]

	rho4 = drho_dr_3(rho + dr*rho3, T + dr*T3, M + dr*M3, r + dr, dT_dr_3(rho + dr*rho3, T + dr*T3, M + dr*M3, L + dr*L3, r + dr, kappa(rho + dr*rho3, T + dr*T3), lam), lam)
	T4 = dT_dr_3(rho + dr*rho3, T + dr*T3, M + dr*M3, L + dr*L3, r + dr, kappa(rho + dr*rho3, T + dr*T), lam)
	M4 = dM_dr(rho + dr*rho3, r + dr)
	L4 = dL_dr(rho + dr*rho3, r + dr, epsilon(rho + dr*rho3, T + dr*T3))
	tau4 = dtau_dr(rho + dr*rho3, kappa(rho + dr*rho3, T + dr*T3))
	k_array[:,3] = [rho4, T4, M4, L4, tau4]

	return k_array

def k_array_1(rho, T, M, L, r, dr, LAM): #h is step size
	k_array = np.empty((5, 4))

	rho1 = drho_dr_1(rho, T, M, r, dT_dr_1(rho, T, M, L, r, kappa(rho, T), LAM), LAM)
	T1 = dT_dr_1(rho, T, M, L, r, kappa(rho, T), LAM)
	M1 = dM_dr(rho, r)
	L1 = dL_dr(rho, r, epsilon(rho, T))
	tau1 = dtau_dr(rho, kappa(rho, T))
	k_array[:,0] = [rho1, T1, M1, L1, tau1]

	rho2 = drho_dr_1(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, r + dr/2.0, dT_dr_1(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, L + dr*L1/2.0, r + dr/2.0, kappa(rho + dr*rho1/2.0, T + dr*T1/2.0), LAM), LAM)
	T2 = dT_dr_1(rho + dr*rho1/2.0, T + dr*T1/2.0, M + dr*M1/2.0, L + dr*L1/2.0, r + dr/2.0, kappa(rho + dr*rho1/2.0, T + dr*T/2.0), LAM)
	M2 = dM_dr(rho + dr*rho1/2.0, r + dr/2.0)
	L2 = dL_dr(rho + dr*rho1/2.0, r + dr/2.0, epsilon(rho + dr*rho1/2.0, T + dr*T1/2.0))
	tau2 = dtau_dr(rho + dr*rho1/2.0, kappa(rho + dr*rho1/2.0, T + dr*T1/2.0))
	k_array[:,1] = [rho2, T2, M2, L2, tau2]

	rho3 = drho_dr_1(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*M2, r + dr/2.0, dT_dr_1(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*T2, L + dr/2.0*L2, r + dr/2.0, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2), LAM), LAM)
	T3 = dT_dr_1(rho + dr/2.0*rho2, T + dr/2.0*T2, M + dr/2.0*M2, L + dr*L2/2.0, r + dr/2.0, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2), LAM)
	M3 = dM_dr(rho + dr/2.0*rho2, r + dr/2.0)           
	L3 = dL_dr(rho + dr/2.0*rho2, r + dr/2.0, epsilon(rho + dr/2.0*rho2, T + dr/2.0*T2))
	tau3 = dtau_dr(rho + dr/2.0*rho2, kappa(rho + dr/2.0*rho2, T + dr/2.0*T2))
	k_array[:,2] = [rho3, T3, M3, L3, tau3]

	rho4 = drho_dr_1(rho + dr*rho3, T + dr*T3, M + dr*M3, r + dr, dT_dr_1(rho + dr*rho3, T + dr*T3, M + dr*M3, L + dr*L3, r + dr, kappa(rho + dr*rho3, T + dr*T3), LAM), LAM)
	T4 = dT_dr_1(rho + dr*rho3, T + dr*T3, M + dr*M3, L + dr*L3, r + dr, kappa(rho + dr*rho3, T + dr*T), LAM)
	M4 = dM_dr(rho + dr*rho3, r + dr)
	L4 = dL_dr(rho + dr*rho3, r + dr, epsilon(rho + dr*rho3, T + dr*T3))
	tau4 = dtau_dr(rho + dr*rho3, kappa(rho + dr*rho3, T + dr*T3))
	k_array[:,3] = [rho4, T4, M4, L4, tau4]

	return k_array

def f_pc(L_star, R_star, T_star):
	SB_law = 4.0*np.pi*sigma*(R_star**2.0)*(T_star)**4.0
	function = (L_star - SB_law)/((L_star*SB_law)**(0.5))
   	return function














