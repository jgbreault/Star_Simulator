import numpy as np
import matplotlib.pyplot as plt
from equations import *
import timeit
start = timeit.default_timer()

def bisection(int_func, a, b, max_iterations, T):
    n = 1

    c = (a+b)/2.0
    a_int = rk_int(a, T)
    c_int = rk_int(c, T)
    f_pc_a = f_pc(a_int[3,-1], a_int[5,-1], a_int[1,-1])
    f_pc_c = f_pc(c_int[3,-1], c_int[5,-1], c_int[1,-1])

    while (n <= max_iterations):
    	while_start = timeit.default_timer()

        if (f_pc_a < 0.0 and f_pc_c < 0.0) or (f_pc_a > 0.0 and f_pc_c > 0.0):
        	a = c
        	c = (a+b)/2.0
        	c_int = rk_int(c, T)
        	f_pc_c = f_pc(c_int[3,-1], c_int[5,-1], c_int[1,-1])

        else:
        	b = c
        	c = (a+b)/2.0
        	c_int = rk_int(c, T)
        	f_pc_c = f_pc(c_int[3,-1], c_int[5,-1], c_int[1,-1])

        while_end = timeit.default_timer()
        print n, 'iterations,', ' Time since previous iteration', while_end - while_start
        print c
     	n += 1
    #while loop ended:
    if f_pc_c < 0.0: #yeilds most correct Rho_c that is positive
    	return rk_int(b, T)
    else:
    	return rk_int(c, T)

def rk_int(rho_init, T_init):

    ## defining initial conditions
    r0 = 1.0 # arbitrarily small but not 0
    rho0 = rho_init # a guess to be fine tuned
    T0 = T_init # to be varied to produce main sequence
    M0 = 4.0*np.pi*r0**3.0*rho0/3.0
    L0 = 4.0*np.pi*r0**3.0*rho0/3.0*epsilon(rho0, T0)
    tau0 = 0.0 # I don't know how to deal with this properly yet

    r_v = [r0]
    rho_v = [rho0]
    T_v = [T0]
    M_v = [M0]
    L_v = [L0]
    tau_v = [tau0]

    delta_tau = 1.0
    R_max = 20.0*R_sun
    M_max = 100.0*M_sun

    i = 0
    j = 1000
    while delta_tau>0.0001:
    	if T_v[i] < 50000.0:
    		dr = 500 # chosen arbitrarily
    	else:
    		dr = 5000 + r_v[i]/2000.0 # chosen arbitrarily
        
        k_v = k_array(rho_v[i], T_v[i], M_v[i], L_v[i], r_v[i], dr)

        r_v.append(r_v[i]+dr)
        rho_v.append(rho_v[i]+dr/6.0*(k_v[0,0] + 2.0*k_v[0,1] + 2.0*k_v[0,2]+ k_v[0,3]))
        T_v.append(T_v[i]+dr/6.0*(k_v[1,0] + 2.0*k_v[1,1] + 2.0*k_v[1,2]+ k_v[1,3]))
        M_v.append(M_v[i]+dr/6.0*(k_v[2,0] + 2.0*k_v[2,1] + 2.0*k_v[2,2]+ k_v[2,3]))
        L_v.append(L_v[i]+dr/6.0*(k_v[3,0] + 2.0*k_v[3,1] + 2.0*k_v[3,2]+ k_v[3,3]))
        tau_v.append(tau_v[i]+dr/6.0*(k_v[4,0] + 2.0*k_v[4,1] + 2.0*k_v[4,2]+ k_v[4,3]))

        i += 1
        delta_tau = kappa(rho_v[i], T_v[i])*rho_v[i]**2.0/abs(drho_dr(rho_v[i], T_v[i], M_v[i], r_v[i], dT_dr(rho_v[i], T_v[i], M_v[i], L_v[i], r_v[i], kappa(rho_v[i], T_v[i]))))

        #print tau_v[i] - tau_v[i-1]
        if i==j:
            j=j+1000
            #print M_v[i]/M_sun, r_v[i]/R_sun, dr, T_v[i]
            #print delta_tau

	    if r_v[i] > R_max or M_v[i] > M_max:
	    	break

    output = np.array([rho_v, T_v, M_v, L_v, tau_v, r_v])

    return output

## Bisection Method ###########################

#T0_range = np.linspace(10.0**6.6, 10.0**7.4, 1)
T0_range = [10.0**7.35]
MS_info_v = []
for i in range(len(T0_range)):
	int_result = bisection(rk_int, 300, 500000, 13, T0_range[i])
	# end = timeit.default_timer()
	# print "TOTAL TIME: ", end - start 
	MS_info = [int_result[1,-1], int_result[2,-1], int_result[3,-1], int_result[5,-1]] #T,M,L,R
	MS_info_v.append(MS_info)
	print 'Done star: ', i+1
	np.savetxt('Main_sequence.txt', MS_info_v)

kappa_es_v = []
kappa_ff_v = []
kappa_H_v = []
kappa_tot_v = []

P_deg_v = []
P_gas_v = []
P_phot_v = []
P_tot_v = []

dL_dr_pp_v = []
dL_dr_CNO_v = []
dL_dr_v = []

for i in range(len(int_result[5,:])):
	rho_v = int_result[0,:]
	T_v = int_result[1,:]
	r_v = int_result[5,:]

	kappa_es_v.append(kappa_es())
	kappa_ff_v.append(kappa_ff(rho_v[i], T_v[i]))
	kappa_H_v.append(kappa_H(rho_v[i], T_v[i]))
	kappa_tot_v.append(kappa(rho_v[i], T_v[i]))

	P_deg_v.append(P_deg(rho_v[i]))
	P_gas_v.append(P_gas(rho_v[i], T_v[i]))
	P_phot_v.append(P_phot(T_v[i]))
	P_tot_v.append(P_tot(rho_v[i], T_v[i]))

	dL_dr_pp_v.append(dL_dr_pp(rho_v[i], T_v[i], r_v[i]))
	dL_dr_CNO_v.append(dL_dr_CNO(rho_v[i], T_v[i], r_v[i]))
   	dL_dr_v.append(dL_dr(rho_v[i], r_v[i], epsilon(rho_v[i], T_v[i])))

###########################

radius = max(int_result[5,:])
plt.plot(int_result[5,:]/radius, int_result[0,:]/max(int_result[0,:]), label = r'$\rho$' + '/' + r'$\rho_c$')#rho
plt.plot(int_result[5,:]/radius, int_result[1,:]/max(int_result[1,:]), label = 'T/Tc')#T
plt.plot(int_result[5,:]/radius, int_result[2,:]/max(int_result[2,:]), label = 'M/M*')#M
plt.plot(int_result[5,:]/radius, int_result[3,:]/max(int_result[3,:]), label = 'L/L*')#L
plt.title(r'$\rho$'+'TML vs. r/R*')
plt.xlabel('r/R*')
plt.ylabel(r'$\rho$' + '/' + r'$\rho_c$' + ', T/Tc, M/M*, L/L*')
plt.legend(loc = 0)
plt.savefig('pTML.png')
plt.show()

plt.plot(int_result[5,:]/radius, P_deg_v/max(P_tot_v), 'r-.', label = 'P_deg')
plt.plot(int_result[5,:]/radius, P_gas_v/max(P_tot_v), 'g--', label = 'P_gas')
plt.plot(int_result[5,:]/radius, P_phot_v/max(P_tot_v), 'b:', label = 'P_phot')
plt.plot(int_result[5,:]/radius, P_tot_v/max(P_tot_v), 'k', label = 'P_tot')
plt.title('P vs. r/R*')
plt.xlabel('r/R*')
plt.ylabel('P/Pc')
plt.ylim(-0.05, 1.05)
plt.legend(loc = 0)
plt.savefig('Pressure.png')
plt.show()

plt.plot(int_result[5,:]/radius, np.log10(kappa_es_v), 'b:', label = r'$\kappa_{es}$')
plt.plot(int_result[5,:]/radius, np.log10(kappa_ff_v), 'g--', label = r'$\kappa_{ff}$')
plt.plot(int_result[5,:]/radius, np.log10(kappa_H_v), 'r-.',label = r'$\kappa_{H-}$')
plt.plot(int_result[5,:]/radius, np.log10(kappa_tot_v), 'k', label = r'$\kappa$')
plt.title('log(' +r'$\kappa$' + ') vs. r/R*')
plt.xlabel('r/R*')
plt.ylabel('log(' +r'$\kappa$' +')')
plt.ylim(-2.0, 10.0)
plt.legend(loc = 0)
plt.savefig('kappa.png')
plt.show()

plt.plot(int_result[5,:]/radius, dL_dr_v/(max(int_result[3,:])/radius), 'k', label = 'dL/dr')
plt.plot(int_result[5,:]/radius, dL_dr_pp_v/(max(int_result[3,:])/radius), 'r-.', label = 'dL/dr_pp')
plt.plot(int_result[5,:]/radius, dL_dr_CNO_v/(max(int_result[3,:])/radius), 'b:', label = 'dL/dr_CNO')
plt.title('dL/dr vs. r/R*')
plt.xlabel('r/R*')
plt.ylabel('dL/dr (L*/R*)')
plt.legend(loc = 0)
plt.ylim(-0.5)
plt.savefig('dLdr.png')
plt.show()

# dT_dr_v = []
# for i in range(len(int_result[5,:])):
# 	mu = (2.0*X + 0.75*Y + 0.5*Z)**(-1.0)
# 	dT_dr_v.append(dT_dr(int_result[0,i], int_result[1,i], int_result[2,i], int_result[3,i], int_result[5,i], kappa(int_result[0,i],int_result[1,i])))
# 	P = (3.0*np.pi**2.0)**(2.0/3.0)*hbar**2.0*(int_result[0,i]/m_proton)**(5.0/3.0)/5.0/m_electron + int_result[0,i]*k*int_result[1,i]/mu/m_proton + a*int_result[1,i]**4.0/3.0
	
# 	print dT_dr_v[i]
# 	print ((1.0 - 1.0/gamma)*int_result[1,i]*G*int_result[2,i]*int_result[0,i]/P/int_result[5,i]**2.0)
# 	print 3.0*kappa(int_result[0,i],int_result[1,i])*int_result[0,i]*int_result[3,i]/16.0/np.pi/a/c/int_result[1,i]**3.0/int_result[5,i]**2.0

# 	if dT_dr_v[i] == -3.0*kappa(int_result[0,i],int_result[1,i])*int_result[0,i]*int_result[3,i]/16.0/np.pi/a/c/int_result[1,i]**3.0/int_result[5,i]**2.0:
# 		print int_result[5,i]
# 		break

dlogP_dlogT = []
for i in range(len(P_tot_v)):
    dPdr_u = dP_dr(int_result[0,i], int_result[2,i], int_result[5,i])
    dTdr_u = dT_dr(int_result[0,i], int_result[1,i], int_result[2,i], int_result[3,i], int_result[5,i], kappa(int_result[0,i], int_result[1,i]))
    T_u = T_v[i]
    P_u = P_tot_v[i]
    dlogP_dlogT.append(dPdr_u*(dTdr_u)**(-1.0)*T_u/P_u)
print int_result[:,-1]
	#print dP_dr(int_result[0,i], int_result[2,i], int_result[5,i]), (dT_dr(int_result[0,i], int_result[1,i], int_result[2,i], int_result[3,i], int_result[5,i], kappa(int_result[0,i], int_result[1,i])))**(-1.0), T_v[i], P_tot_v[i]
	#print dP_dr(int_result[0,i], int_result[2,i], int_result[5,i])*(dT_dr(int_result[0,i], int_result[1,i], int_result[2,i], int_result[3,i], int_result[5,i], kappa(int_result[0,i], int_result[1,i])))**(-1.0)*T_v[i]/P_tot_v[i]

plt.plot(int_result[5,:]/radius, dlogP_dlogT)
plt.title('dlogP/dlogT vs. r/R*')
plt.xlabel('r/R*')
plt.ylabel('dlogP/dlogT')
plt.ylim(2.0, 5.0)
plt.savefig('dlogP_dlogT.png')
plt.show()















