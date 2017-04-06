# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary


# Iterative calcultion of thickness as per ASME Sec. VII Div 2

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
# help as per https://www.cyberciti.biz/faq/howto-get-current-date-time-in-python/
nowtext = datetime.now().strftime('%Y_%m_%d_%H%M%S')
sys.stdout = open('log/' + nowtext + '.txt', 'w')
print(datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
print()

# Define all functions
# Get FS as per 4.4.2
def getFS(fic, sy):

    # 4.4.1
    if fic <= 0.55*sy:
        fs = 2
    # 4.4.2
    elif 0.55*sy < fic < sy:
        fs = 2.407 - 0.741*(fic/sy)
    # 4.4.3
    elif fic == sy:
        fs = 1.667
    else:
        print('fs not found, equal 1')
        fs = 1

    return fs

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
t_0 = 0.5                 # initial thickness guess
t_step = 0.05               # step for convergence
maxit = 100                 # max iterations avoid infinite loop

# Calculated initial parameters
R_o = D_o/2
S_allow = S_y/SF
L = np.ceil((SG * d_tether * L_tether) / (np.pi * (2 * D_o + 3 * d_tether)))
p_req = (2*P)/(D_o * d_tether)

# Initialize loop
p_a = 0
t = t_0
itnum = 1

while p_a < p_req:

    # check for first interation, oteherwise add the t_step
    if itnum != 1:
        t += t_step

    print('Iteration %i for t = %.3f in' % (itnum, t))
    print('pa = %.1f psi <= preq = %.1f psi' % (p_a, p_req))

    # calc D_0/t ratio for first guess
    Dt = D_o / t
    LD = L / D_o

    F_he = (1.6*C_y*E_y*t)/D_o

    M_x = L/((R_o*t)**0.5)

    # Check if current it >= max
    if itnum >= maxit:
        # multiple blank lines
        print('\n' * 2)
        print('Did not find solution after %i iterations' % itnum)
        break

    # check if solution converged, print messages and then it will exit loop
    elif p_a >= p_req:
        print('A thickness of %.3f in will be safe' %t)
        print('pa = %.1f psi >= preq = %.1f psi' %(p_a, p_req))

    # equivalent to c++'s i++, inc itnum to avoid infinite loop
    itnum += 1

    # multiple blank lines
    print('\n' * 2)
