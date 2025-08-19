#---------------------------------------------------#
#     CODE : RRKM Rate calculation using frequencies#
# Method 2 : Beyer-Swinehart Algorithm              #
#        By: Philips Kumar Rai                      #
#    Malaviya National Technology Jaipur            #
#    Email : philipsrai786@gmail.com                #
#-------------------------------------------------- #
import matplotlib.pyplot as plt

h = 6.626e-34                                                
freq_rc = [200, 400, 700] # cm-1
freq_ts = [200, 500]      # cm-1  
dE = 100                  # cm-1                                                     
E_max = 1000              # cm-1                                               

nbins = int(E_max / dE) + 1
E = [i * dE for i in range(nbins)] 
# E=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
####################################################
# A quick Example                                  #
# Suppose E_max = 300, dE = 50                     #
# nbins = int(E_max / dE) + 1                      #
# nbins = int(300 / 50) + 1 = int(6) + 1 = 7       #
# E = [i * dE for i in range(nbins)]               #
# range(nbins) = range(7) = [0, 1, 2, 3, 4, 5, 6]  #
# On Multiplying each i by dE = 50                 #
#  E = [0, 50, 100, 150, 200, 250, 300]            #
####################################################
# ------------------#
# Sum of States N(E)#
# ------------------#

N_rc = [0] * nbins  # [0,0,0,0,0,0,0,0,0,0,0]
N_ts = [0] * nbins  # [0,0,0,0,0,0,0,0,0,0,0]
N_rc[0] = 1         
N_ts[0] = 1   
# Now N_rc/N_ts list becomes [1,0,0,0,0,0,0,0,0,0,0]     
# RC
for f in freq_rc:
    f_bin = int(f / dE) 
# f_bin will be 200/100 = 2 then 400/100 = 4, and then 700/100 = 7         
    for i in range(f_bin, nbins):
        N_rc[i] += N_rc[i - f_bin]         
# Note: x+=y means x=x+y        
# Now we know N_rc = [1,0,0,0,0,0,0,0,0,0,0]   
# nbins=11 means 11 count entries will be there corresponding to E= 0 to 1000
# ------------------- #
# Frequency = 200 cm-1#
# f_bin = 200/100 = 2 #
# ------------------- #
# i=2  (corresponding to E=200): N_rc[2] = N_rc[2] + N_rc[0] = 0 + 1 = 1
# i=3  (corresponding to E=300): N_rc[3] = N_rc[3] + N_rc[1] = 0 + 0 = 0
# i=4  (corresponding to E=400): N_rc[4] = N_rc[4] + N_rc[2] = 0 + 1 = 1
# i=5  (corresponding to E=500): N_rc[5] = N_rc[5] + N_rc[3] = 0 + 0 = 0
# i=6  (corresponding to E=600): N_rc[6] = N_rc[6] + N_rc[4] = 0 + 1 = 1
# i=7  (corresponding to E=700): N_rc[7] = N_rc[7] + N_rc[5] = 0 + 0 = 0
# i=8  (corresponding to E=800): N_rc[8] = N_rc[8] + N_rc[6] = 0 + 1 = 1
# i=9  (corresponding to E=900): N_rc[9] = N_rc[9] + N_rc[7] = 0 + 0 = 0
# i=10 (corresponding to E=1000):N_rc[10]= N_rc[10]+ N_rc[8] = 0 + 1 = 1
# After 200 cm-1 mode: N_rc = [1,0,1,0,1,0,1,0,1,0,1]

# ------------------- #
# Frequency = 400 cm-1#
# f_bin = 400/100 = 4 #
# ------------------- #
# i=4  (corresponding to E=400): N_rc[4] = N_rc[4] + N_rc[0] = 1 + 1 = 2
# i=5  (corresponding to E=500): N_rc[5] = N_rc[5] + N_rc[1] = 0 + 0 = 0
# i=6  (corresponding to E=600): N_rc[6] = N_rc[6] + N_rc[2] = 1 + 1 = 2
# i=7  (corresponding to E=700): N_rc[7] = N_rc[7] + N_rc[3] = 0 + 0 = 0
# i=8  (corresponding to E=800): N_rc[8] = N_rc[8] + N_rc[4] = 1 + 2 = 3
# i=9  (corresponding to E=900): N_rc[9] = N_rc[9] + N_rc[5] = 0 + 0 = 0
# i=10 (corresponding to E=1000):N_rc[10]= N_rc[10]+ N_rc[6] = 1 + 2 = 3
# After 400 cm-1 mode: N_rc = [1,0,1,0,2,0,2,0,3,0,3]

# -------------------#
# Frequency = 700 cm-1
# f_bin = 700/100 = 7
# -------------------#
# i=7  (corresponding to E=700): N_rc[7] = N_rc[7] + N_rc[0] = 0 + 1 = 1
# i=8  (corresponding to E=800): N_rc[8] = N_rc[8] + N_rc[1] = 3 + 0 = 3
# i=9  (corresponding to E=900): N_rc[9] = N_rc[9] + N_rc[2] = 0 + 1 = 1
# i=10 (corresponding to E=1000):N_rc[10]= N_rc[10]+ N_rc[3] = 3 + 0 = 3
# After 700 cm-1 mode: N_rc = [1,0,1,0,2,0,2,1,3,1,3]
# Final N_rc array:[1, 0, 1, 0, 2, 0, 2, 1, 3, 1, 3]
# This means:
# g(0)=1    degeneracy of counts for E = 0
# g(200)=1  degeneracy of counts for E = 200
# g(400)=2  degeneracy of counts for E = 400
# g(600)=2  degeneracy of counts for E = 600
# g(700)=1  degeneracy of counts for E = 700
# g(800)=3  degeneracy of counts for E = 800
# g(900)=1  degeneracy of counts for E = 900
# g(1000)=3 degeneracy of counts for E = 1000
print(N_rc)
# TS
for f in freq_ts:
    f_bin = int(f / dE)
    for i in range(f_bin, nbins):
        N_ts[i] += N_ts[i - f_bin]
print(N_ts)
# ----------------------------------------#
# Density of States                       #
# Central Difference Approach             #
# ----------------------------------------#

rho_rc = [0] * nbins
rho_ts = [0] * nbins

# For energies inside 0 to Emax points, use central difference formula:
#   rho(E) = (N(E+dE) - N(E-dE)) / (2*dE)

for i in range(1, nbins - 1):
    rho_rc[i] = (N_rc[i+1] - N_rc[i-1]) / (2*dE)
    rho_ts[i] = (N_ts[i+1] - N_ts[i-1]) / (2*dE)

# For E = 0, use forward difference formula:
#     rho(0) = (N(100) - N(0)) / (100)

rho_rc[0]  = (N_rc[1] - N_rc[0]) / dE
rho_ts[0]  = (N_ts[1] - N_ts[0]) / dE

# For E = Emax, use backward difference formula:
#     rho(Emax) =  (N(Emax) - N(Emax - 100)) / (100)

rho_rc[-1] = (N_rc[-1] - N_rc[-2]) / dE
rho_ts[-1] = (N_ts[-1] - N_ts[-2]) / dE

#RRKM rate
c = 2.998e10   # cm/s
E0 = 0         # barrier (cm-1)

kE = [0] * nbins # we will initialize rate constant array (corresponding to Energy array)
for i, Ei in enumerate(E):
    if Ei >= E0 and rho_rc[i] > 0: 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #    
# Condition: only compute k(E) if E >= barrier and denominator rho_rc[i] > 0                                                                                                      #
# For rho_rc[i] = 0, denominator would be zero, leading to division by zero (undefined result). Physically: a zero density means                                                  #
# no states available at that energy. So the system cannot exist at that energy, i.e., rate constant doesn't make sense. On the other hand, Ei ​is the total available energy at   #
# that bin and E0 is the activation barrier height. In RRKM theory, the rate constant k(E) is only meaningful if the system has at least the barrier energy available. If Ei < E0,# 
# the system simply cannot cross the transition state, so the rate constant must be zero.                                                                                         #
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
        idx = int((Ei - E0) // dE)         
# idx just mentions the index number of TS states array    
        kE[i] = c * N_ts[idx] / rho_rc[i]
# RRKM formula:--------------------------------------------#
# k(E) = (c * N_ts(E - E0)) / rho_rc(E)                    #
# Numerator: number of accessible TS states up to energy E #
# Denominator: density of RC states at energy E            #
# Ratio: probability flux through TS channel               #
# c is for the conversion to sec-1                         #
#--------------------------------------------------------- #
E_calc = 600 # cm-1  

idx = int((E_calc - 0.0) // dE)
# Here idx will be 6 means pick up the 6th element of sum of states array which will be sum of states corresponding to Ei - E0 here 600 cm-1, i.e., N_ts(600-0)
if idx >= nbins:
    idx = nbins - 1

print(f"At E = {E[idx]} cm-1; k(E) = {kE[idx]:.3e} sec-1")
#-----------------------------------------------------------THE END------------------------------------------------------------------------------------------------------------------




















