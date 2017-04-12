# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in SI units
# Script as per EN-13445-3 ,
# Section 8.5, solving for cylindrical shell thickness

# Import relevant modules
import sys
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt


# Logging
nowtext = datetime.now().strftime('%Y_%m_%d_%H%M%S')
sys.stdout = open('log/' + nowtext + '.txt', 'w')
print(datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
print()

# Define all functions

# Get pr/py  using interpolation and csv array
def getPRY(arr, val):
    # Get smallest difference between array and val
    val_idx = (np.abs(arr[:, 0] - val)).argmin()

    # add check for roundup
    if val >= arr[val_idx, 0]:
        val_idx += 1

    # Get points for linear interpolation
    x1 = arr[val_idx - 1, 0]
    x2 = arr[val_idx, 0]
    y1 = arr[val_idx - 1, 1]
    y2 = arr[val_idx, 1]

    # Calculate approximate b
    pry = (val-x1)*((y2-y1)/(x2-x1))+y1

    return pry

# Get ncyl assuming L/2R =2
def getNCYL(er):

    if er <= 0.001:
        ncyl = 7
    elif 0.001 < er <= 0.0025:
        ncyl = 6
    elif 0.0025 < er <= 0.006:
        ncyl = 5
    elif 0.006 < er <= 0.02:
        ncyl = 4
    elif 0.02 < er < 0.05:
        ncyl = 3
    else:
        # i.e >= 0.05
        ncyl = 2


    print('Because e/2R = %.4f, n_cyl = %i' %(er, ncyl))
    return ncyl


# Define all parameters
D_o = 28*25.4                   # diameter
d_tether = 15.2        			# tether diameter 15.2mm
L_tether = 370*1000    			# tether length overall 370m to [mm]
P = 51264          				# force P of 51.264 kN
SF = 1.5                    	# safety factor [ul]
S_y = 241.3                 	# Yield strength 35 ksi to [MPa]
E_y = 200000            		# Youngs modulus 200 [GPa]
v = 0.3                     	# poissons ratio [ul]
SG = 1.01                   	# spool gap 1% [ul]
t_0 = 12.7                 		# initial thickness guess 0.5 in
t_step = 0.254               	# step for convergence (10 thou in)
maxit = 200                 	# max iterations avoid infinite loop

# Define EN-13445-3 specific properties
S = 1.5
R_p02 = S_y
S_e = R_p02

# Calculated initial parameters
R = D_o/2
S_allow = S_y/SF
L = np.ceil((SG * d_tether * L_tether) / (np.pi \
                            * (2 * D_o + 3 * d_tether)))
p_req = (2*P)/(D_o * d_tether)

# Get csv data of p_m vs p_r
# Note, cols as per p_m/p_y, p_r/p_y
pm_pr = genfromtxt('csv/pm_pr.csv', delimiter=',')


# Initialize loop
p_a = 0
e_a = t_0
itnum = 1

# Initiliaze arrays
arr_it = []
arr_t = []
arr_pa = []

while p_a < p_req:

    # check for first iteration, otherwise add the t_step
    if itnum != 1:
        e_a += t_step

    # Display current iteration values
    print('Iteration %i :' % itnum)
    print('t = %.2f mm' % e_a)

    # Calculate the average circumferential
    # pressure in the cylinder p_y
    # 8.5.2.4
    p_y = (S_e*e_a)/R

    # lower failure pressure p_m
    # Get epsilion
    # 8.5.2.6
    # Split in ABC for easier equation
    # print('From chart 8.5.5: e/2: 3')
    # n_cyl = sys.stdin.readline()
    n_cyl = getNCYL(e_a/(2*R))
    Z = (np.pi*R)/L
    A = 1/((n_cyl**2) - 1 + ((Z**2)/2))
    B = 1/((((n_cyl**2)/(Z**2))+1)**2)
    C = (e_a**2*((n_cyl**2) - 1 + (Z**2)))/(12*(R**2)*(1-(v**2)))
    eps = A*(B+C)

    # 8.5.2.5
    p_m = (E_y*e_a*eps)/R

    # Calculate p_m/ p_y ratio
    p_m_y = p_m /p_y

    # Use interpolation function to
    # handle worst case scenerio >= 7
    if p_m_y >= 7:
        p_r_y =7
    else:
        p_r_y = getPRY(pm_pr, p_m_y)

    print('Ratio pm/py = %.3f results in pr/py = %.3f ' %(p_m_y, p_r_y))

    # Calculate p_r
    p_r = p_r_y*p_y

    p_a = p_r/S

    # Add to summary array
    arr_it.append(itnum)
    arr_pa.append(p_a)

    if p_a < p_req:
        print('pa = %.3f MPa < preq = %.3f MPa' % (p_a, p_req))

    # Check if current it >= max
    if itnum >= maxit:
        # multiple blank lines
        print('\n' * 2)
        print('Did not find solution after %i iterations' % itnum)
        break

    # check if solution converged, print messages
    # and then it will exit loop
    elif p_a >= p_req:
        print('pa = %.3f MPa >= preq = %.3f MPa' % (p_a, p_req))
        print('A thickness of %.2f in will be safe' % e_a)

    # equivalent to c++'s i++, inc itnum to avoid infinite loop
    itnum += 1

    # multiple blank lines
    print('\n' * 2)

plt.figure(1)
plt.plot(arr_it, arr_pa)
plt.xlabel('Iteration ' + r'$n$')
plt.xlim([0, itnum])
plt.ylabel('Allowable pressure ' + r'$p_a [MPa]$')
plt.title('Allowable pressure vs iteration number')
plt.savefig('fig/it_vs_pa.png')
