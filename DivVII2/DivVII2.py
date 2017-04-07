# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary


# Iterative calcultion of thickness as per ASME Sec. VII Div 2
# Outline in section 4.4.5.1 for cylindrical shells

# Import relevant modules
import sys
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
# from numpy import genfromtxt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Logging
# help as per https://www.cyberciti.biz/faq/howto-get-current-date-time-in-python/
nowtext = datetime.now().strftime('%Y_%m_%d_%H%M%S')
sys.stdout = open('log/' + nowtext + '.txt', 'w')
print(datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
print()

# Define all functions

# Get Ch as per step 2
def getCH(mx, dot):

    # 4.4.21
    if mx >= 2*(dot**0.94):
        ch = 0.55*(dot**-1)
        print('Equation 4.4.21, Mx = %.1f , Ch =%.1f' %(mx, ch))
    # 4.4.22
    elif 13 < mx < 2*(dot**0.94):
        ch = 1.12*(mx**-1.058)
        print('Equation 4.4.22, Mx = %.1f , Ch =%.1f' % (mx, ch))
    # 4.4.23
    elif 1.5 < mx <= 13:
        ch = 0.92/(mx-0.579)
        print('Equation 4.4.23, Mx = %.1f , Ch =%.1f' % (mx, ch))
    # 4.4.24 (i.e. mx <= 1.5)
    else:
        ch = 1
        print('Equation 4.4.24, Mx = %.1f , Ch =%.1f' % (mx, ch))

    return ch

# Get Fic as per step 3
def getFIC(feh, sy):

    # 4.4.25
    if (feh / sy) >= 2.439:
        fic = sy
        print('Equation 4.4.25, Feh = %.1f , Sy =%.1f' % (feh, sy))
    # 4.4.26
    elif 0.552 < feh/sy < 2.439:
        fic = 0.7*sy*((feh/sy)**0.4)
        print('Equation 4.4.26, Feh = %.1f , Sy =%.1f' % (feh, sy))
    # 4.4.27 (feh/sy <= 0.552)
    else:
        print('Equation 4.4.27, Feh = %.1f , Sy =%.1f' % (feh, sy))
        fic = feh

    return fic

# Get FS as per 4.4.2
def getFS(fic, sy):

    # 4.4.1
    if fic <= 0.55*sy:
        fs = 2
        print('Equation 4.4.1, FS = %.1f , Fic =%.1f' % (fs, fic))
    # 4.4.2
    elif 0.55*sy < fic < sy:
        fs = 2.407 - 0.741*(fic/sy)
        print('Equation 4.4.2, FS = %.1f , Fic =%.1f' % (fs, fic))
    # 4.4.3
    elif fic == sy:
        fs = 1.667
        print('Equation 4.4.3, FS = %.1f , Fic =%.1f' % (fs, fic))
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
t_step = 0.01               # step for convergence
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

# Initiliaze arrays
arr_it = []
arr_t = []
arr_pa = []

while p_a < p_req:

    # check for first iteration, otherwise add the t_step
    if itnum != 1:
        t += t_step

    # Display current iteration values
    print('Iteration %i :' % itnum)
    print('t = %.3f in' % t)

    # calc D_0/t ratio for first guess
    Dt = D_o / t
    LD = L / D_o

    # Shell parameter 4.4.20
    M_x = L/((R_o*t)**0.5)

    # Call function to get Ch parameter
    C_h = getCH(M_x, Dt)

    # elastic hoop compressive membrane failure 4.4.19
    F_he = (1.6*C_h*E_y*t)/D_o

    # Step 3 : calculate predicted buckling stress Fic:
    F_ic = getFIC(F_he, S_y)

    # Step 4 : w/ Fic get the design factor FS
    FS = getFS(F_ic, S_y)

    # Step 5 calculate allowable pressure p-a
    # Allowable hoop compressive membrane stress 4.4.29
    F_ha = F_ic/FS
    # Allowable pressure 4.4.28
    p_a = 2*F_ha*(Dt**-1)

    # Add to summary array
    arr_it.append(itnum)
    arr_pa.append(p_a)

    if p_a < p_req:
        print('pa = %.1f psi < preq = %.1f psi' % (p_a, p_req))

    # Check if current it >= max
    if itnum >= maxit:
        # multiple blank lines
        print('\n' * 2)
        print('Did not find solution after %i iterations' % itnum)
        break

    # check if solution converged, print messages and then it will exit loop
    elif p_a >= p_req:
        print('pa = %.1f psi >= preq = %.1f psi' % (p_a, p_req))
        print('A thickness of %.3f in will be safe' % t)

    # equivalent to c++'s i++, inc itnum to avoid infinite loop
    itnum += 1

    # multiple blank lines
    print('\n' * 2)

plt.figure(1)
plt.plot(arr_it, arr_pa)
plt.xlabel('Iteration ' + r'$n$')
plt.xlim([0, itnum])
plt.ylabel('Allowable pressure ' + r'$p_a [psi]$')
plt.title('Allowable pressure vs iteration number')
plt.savefig('fig/it_vs_pa.png')