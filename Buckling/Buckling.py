# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary

# Buckling Calculations

# Import relevant modules
import sys
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Logging
nowtext = datetime.now().strftime('%Y_%m_%d_%H%M%S')
# sys.stdout = open('log/' + nowtext + '.txt', 'w')
print(datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
print()

# Define all functions



# Define all parameters
D_o = 28                    # diameter
d_tether = 15.2/25.4        # tether diameter 15.2mm to [in]
L_tether = (370*1000)/25.4  # tether length overall 370m to [in]
P = 224.809*51.264          # force P of 51.264 kN [lbs]
SF = 1.5                    # safety factor [ul]
S_y = 35000                 # Yield strength 35 ksi
E_y = 2.9*10**7             # Youngs modulus 200 [GPa]
v = 0.3                     # poissons ratio [ul]
rho = 0.284                 # density of steel [lb/in^3]
SG = 1.01                   # spool gap 1% [ul]

# Calculated initial parameters
R_o = D_o/2
S_allow = S_y/SF
L = np.ceil((SG * d_tether * L_tether) / (np.pi * \
                            (2 * D_o + 3 * d_tether)))
p_req = (2*P)/(D_o * d_tether)

# Data import, x:: is x to end
ANSYS = genfromtxt('csv/ANSYS_run3.csv', delimiter=',')
t = ANSYS[::, 0]
lmbd = ANSYS[::, 1::]
mode = [1, 5, 10]

l_crt = 4.9*R_o*((R_o/t)**0.5)

q = np.zeros((len(t), 1))
PREQ = np.zeros((len(t), 1))


for i in range(0, len(t)):
    if L > l_crt[i]:
        q[i] = (1/4)*(E_y/(1-(v**2)))*((t[i]/R_o)**3)
    else:
        q[i] = 0.807 * (E_y*(t[i]**2) / (L*R_o)) * ((((1/(1-(v**2)))**3)*((t[i])/R_o)**2) ** 0.25)

    PREQ[i] = p_req

# Analytical Results
plt.figure(1)
plt.plot(t, q)
plt.xlabel('Thickness ' + r'$t\ [in]$')
plt.ylabel('Critical Pressure ' + r'$q\ [psi]$')
plt.title('Critical Pressure vs Thickness')
plt.savefig('fig/'+nowtext+'_q_cr.png')

# ANSYS Results
plt.figure(2)
for j in range(0, 3):
    plt.plot(t, lmbd[:, j], label=(r'Mode $\psi=$' + ('%i' % mode[j])))
plt.legend()
plt.xlabel('Thickness ' + r'$t\ [in]$')
plt.ylabel('Load Factor ' + r'$\lambda$')
plt.title('ANSYS Load Factor vs Thickness')
plt.savefig('fig/'+nowtext+'_ANSYS.png')

# ANSYS Results Comparison
plt.figure(3)
plt.plot(t, lmbd[:, 0], label=(r'Mode $\psi=$' + ('%i' % mode[0])))
plt.plot(t, q, label='Analytical')
plt.plot(t, PREQ, label='Required')
plt.legend()
plt.xlabel('Thickness ' + r'$t\ [in]$')
plt.ylabel('Pressure ' + r'$p\ [psi]$')
plt.title('Critical Pressure vs Thickness')
plt.savefig('fig/'+nowtext+'_comp.png')