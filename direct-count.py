#---------------------------------------------------#
#     CODE : RRKM Rate calculation using frequencies#
# Method 1 : Traditional direct counting of states  #
#        By: Philips Kumar Rai                      #
#    Malaviya National Technology Jaipur            #
#    Email : philipsrai786@gmail.com                #
#-------------------------------------------------- #

import math
import pandas as pd

c = 3.0e10 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
# Speed of light (cm/s)
# NOTE: Planck constant h is not needed because we stay in cm-1 units throughout. In this "spectroscopic unit system", RRKM formula is: k(E) = c * N_ts(E-E0) / rho_rc(E).
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
freq_rc = [200, 400, 700]            # Reactant complex modes (cm-1)                              #
freq_ts = [200, 500]                 # Transition state modes (cm-1)                              #
dE = 100                             # Bin size (cm-1)                                            #
E_max = 1000                         # Maximum energy (cm-1)                                      #
E = list(range(0, E_max + dE, dE))   # E = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000] #
nbins = len(E)                       # number of bins                                             #
g_rc = [0] * nbins                   # RC microcanonical degeneracies (states at energy E[i])     #
g_ts = [0] * nbins                   # TS microcanonical degeneracies                             #
                                     #------------------------------------------------------------#
# RC------------------------------------------------------------------------------------------#
# Loop over all quantum numbers (n1, n2, n3) for the 3 RC modes.                              #
# Energy formula: E = n1*freq1 + n2*freq2 + n3*freq3                                          #
# Example: if n1=1, n2=0, n3=0 then E = 200 cm-1                                              #
#          if n1=0, n2=1, n3=1 then E = 400+700=1100 cm-1 (too large and will be discarded)   #
#---------------------------------------------------------------------------------------------#
for n1 in range(E_max // freq_rc[0] + 1):                           #E_max // freq_rc[0] = 1000 // 200 = 5 so n1 = 0,1,2,3,4,5. n1 counts how many quanta of 200 cm-1 mode are excited.
    for n2 in range(E_max // freq_rc[1] + 1):                       #E_max // freq_rc[1] = 1000 // 400 = 2 n2 = 0,1,2. n2 counts quanta of 400 cm-1 mode.                              
        for n3 in range(E_max // freq_rc[2] + 1):                   #E_max // freq_rc[2] = 1000 // 700 = 1 → n3 = 0,1. n3 counts quanta of 700 cm-1 mode.                              
            energy = n1*freq_rc[0] + n2*freq_rc[1] + n3*freq_rc[2]  #Compute total energy. For example, let's choose n1=1, n2=1, n3=0: then E will be 600 cm-1.                        
                                                                    #Energy = 600 <= 1000 and 600 % 100 = 0, we will count this energy.
                                                                    #energy // dE = 600 // 100 = 6 increment g_rc[6] by 1. g_rc now looks like: [1,0,..,0,..,0,1,0,0,0,0].
            if energy <= E_max and energy % dE == 0:
                g_rc[energy // dE] += 1
                # Example: if energy=400 cm-1, then bin=400//100=4
                # g_rc[4] increases by 1 (one state at 400 cm-1)
#------------------------------------------------------------------------------------
# Idea:
# Each vibrational mode contributes integer multiples of its frequency:
# E = n1*200 + n2*400 + n3*700   
# I will Loop over n1, n2, n3 to generate all possible combinations.
# If total E is less than or equal to E_max and falls on a bin boundary (multiple of dE),
# increment g_rc[bin_index].
# Example--------------------------------------------------------------------------------- 
# E_max=1000, dE=100:
# Mode frequencies: 200, 400, 700 cm-1
#------------------------------------------------------------------------------------------
# 200 cm-1 mode (n1 multiples):
#   n1 = 0 gives E=0
#   n1 = 1 gives E=200
#   n1 = 2 gives E=400
#   n1 = 3 gives E=600
#   n1 = 4 gives E=800
#   n1 = 5 gives E=1000
#
# Add 400 cm-1 mode:
#   n2=1 with n1=0 gives E=400
#   n2=1 with n1=1 gives E=600
#   n2=1 with n1=2 gives E=800
#   n2=1 with n1=3 gives E=1000
#   etc.

# Add 700 cm-1 mode:
#   n3=1 with n1=0,n2=0 gives E=700
#   n3=1 with n1=1,n2=0 gives E=900
#   n3=1 with n1=0,n2=1 gives E=1100 # Not acceptable
#   etc.

# Collect counts at each bin (in multiples of 100):
#   E=0     1 state
#   E=200   1 state
#   E=400   2 states (n1=2 or n2=1)
#   E=600   2 states (n1=3 or n1=1+n2=1)
#   E=700   1 state  (n3=1)
#   E=800   3 states (n1=4, n1=2+n2=1, n2=2)
#   E=900   1 state  (n1=1+n3=1)
#   E=1000  3 states (n1=5, n1=3+n2=1, n1=1+n2=2)
# Notes:
#  Some bins like 100, 300, 500 have 0 states because no combination of n1*200 + n2*400 + n3*700 equals these energies.
#   Final g_rc (for 0-1000 cm-1 in steps of 100):
#   g_rc = [1,0,1,0,2,0,2,1,3,1,3]
#----------------------------------------------------------------------------------------------------------------------
# TS
# Same idea but only 2 vibrational modes at TS.
# Example: m1=1, m2=0 then energy = 200 cm-1
#          m1=0, m2=2 then energy = 1000 cm-1
for m1 in range(E_max // freq_ts[0] + 1):
    for m2 in range(E_max // freq_ts[1] + 1):
        energy = m1*freq_ts[0] + m2*freq_ts[1]
        if energy <= E_max and energy % dE == 0:
            g_ts[energy // dE] += 1
            # Example: energy=500 then bin=5 means g_ts[5] += 1
 
# ----------------------------------------
# N(E) = number of states with energy <= E
# Obtained by cumulative summing g(E)
#-----------------------------------------
N_rc = []
cum = 0
for val in g_rc:
    cum += val
    N_rc.append(cum)

N_ts = []
cum = 0
for val in g_ts:
    cum += val
    N_ts.append(cum)

# ------------------
# Density of states 
# ------------------
# Defined as derivative of N(E) wrt energy.
# Numerically: central difference approach.
rho_rc = [0.0] * nbins
for i in range(1, nbins-1):
    rho_rc[i] = (N_rc[i+1] - N_rc[i-1]) / (2*dE)
# At end points (forward/backward difference)
rho_rc[0]  = (N_rc[1] - N_rc[0]) / dE
rho_rc[-1] = (N_rc[-1] - N_rc[-2]) / dE
########################################
# Units                                #
# N(E): dimensionless (count of states)#
# dE: cm-1                             #
# rho_rc: states per cm-1              #
# ######################################
# RRKM rate constant                   #
########################################
def rrkm_rate(E0):
    k = []
    for i, Ei in enumerate(E):
        if Ei < E0:
            # If energy below barrier, no reaction possible
            k.append(0.0)
        else:
            # Determine which TS bin corresponds to (Ei - E0)
            idx = (Ei - E0) // dE
            if rho_rc[i] == 0:
                k.append(0.0)
            else:
                k.append(c * N_ts[idx] / rho_rc[i])
                # Example:
                # E=400 cm-1, E0=200 cm-1 then idx=(400-200)/100=2
                # k(400) = c * N_ts[2] / rho_rc[4]
    return k

k_E0_0 = rrkm_rate(0)      # barrierless case
k_E0_200 = rrkm_rate(200)  # barrier at 200 cm-1

# -------------
# Output Table
# ------------
rows = []
for i, Ei in enumerate(E):
    rows.append({
        "E (cm-1)": Ei,                        ####################################
        "g_rc": g_rc[i],                       # degeneracy of RC at Ei           #
        "N_rc": N_rc[i],                       # cumulative RC states up to Ei    #
        "rho_rc (/cm-1)": round(rho_rc[i], 3), # density of states                #
        "g_ts": g_ts[i],                       # degeneracy of TS at Ei           #
        "N_ts": N_ts[i],                       # cumulative TS states up to Ei    #
        "k(E), E0=0": k_E0_0[i],               # microcanonical k(E), barrierless #
        "k(E), E0=200": k_E0_200[i],           # microcanonical k(E), barrier=200 #
    })                                         ####################################

df = pd.DataFrame(rows)
print(df.to_string(index=False))

