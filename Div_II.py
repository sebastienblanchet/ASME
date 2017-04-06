# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary

# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

# Define all functions
# Get index
# help from MATLAB to py https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html
def findval(arr, val, col):
    # Get smallest difference between array and val
    val_idx = (np.abs(arr[:, col] - val)).argmin()

    # add check for roundup
    if val >= arr[val_idx, col]:
        val_idx += 1
    return val_idx

# 1D interpolation equation
def interpol(x, x1, x2, y1, y2):
    y = (x-x1)*((y2-y1)/(x2-x1))+y1
    return y

# bilinear interpolation equation set y_2 = array
def bi_interpol(A, x, y, x1, x2, y1, y2):

    a_11 = A[y1, x1+1]
    a_21 = A[y2, x1+1]
    a_12 = A[y1, x2+1]
    a_22 = A[y2, x2+1]

    fx_y1 = ((x_2-x)/(x_2-x_1))*a_11 + ((x-x_1)/(x_2-x_1))*a_21
    fx_y2 = ((x_2-x)/(x_2-x_1))*a_12 + ((x-x_1)/(x_2-x_1))*a_22

    z = ((y2-y)/(y2-y1))*fx_y1 +((y-y1)/(y2-y1))*fx_y2
    return z

# Define all parameters
D_o = 28                      # diameter
d_tether = 15.2/25.4        # tether diameter 15.2mm to [in]
L_tether = (370*1000)/25.4  # tether length overall 370m to [in]
P = 224.809*51.264          # force P of 51.264 kN [lbs]
SF = 1.5                    # safety factor [ul]
S_y = 35000                 # Yield strength 35 ksi
E_y = 2.9*10**7             # Youngs modulus 200 [GPa]
v = 0.3                     # poissons ratio [ul]
rho = 0.284                 # density of steel [lb/in^3]
SG = 1.01                   # spool gap 1% [ul]
t_0 = 1/2                  # initial thickness guess

# Calculated initial parameters
S_allow = S_y/SF
L = np.ceil((SG * d_tether * L_tether) / (np.pi * (2 * D_o + 3 * d_tether)))
p_req = (2*P)/(D_o * d_tether)

# Create arrays for imported csv data of ASME Div II, Part
FigG = genfromtxt('csv/IID_FigG.csv', delimiter=',')            # Figure G, to find A
CS2 = genfromtxt('csv/IID_CS2_300F.csv', delimiter=',')         # CS2 for S_y > 30 ksi, T<=300 def F

A_calc = FigG[0:,1:2]

uniqueA = set(A_calc)

# Initialize loop
p_a = 0
t = t_0
itnum = 1
maxit = 10
t_step = 1/8

while p_a < p_req:

    if itnum !=1:
        t += t_step

    # calc D_0/t ratio for first guess
    Dot = D_o/t
    LDo = L / D_o

    # Check for chart applicable aspect ratios
    if Dot >= 4:
        # Check for extreme LDo cases
        if LDo > 50:
            LDo = 50
        # else if statement
        elif LDo < 0.05:
            LDo =0.05
        else:

            # do stuff to find A from (1)
            # print('Method 1 used for A calc')
            A = 0.002

    else:
        # print('Method 2 used for A calc')
        A = 1.1/ (Dot**2)

    # Get nearest value in CS2 table
    ib = findval(CS2, A, 0)
    # print(ib)

    # Interpolate to find B
    A_i1 = CS2[ib - 1, 0]
    A_i2 = CS2[ib, 0]
    B_i1 = CS2[ib - 1, 1]
    B_i2 = CS2[ib, 1]
    B = interpol(A, A_i1, A_i2, B_i1, B_i2)
    # print(B)

    # for now assume Dot >4
    # Calculate allowable pressure
    p_a = (4*B)/(3*Dot)

    # equivalent to c++'s i++, inc itnum to avoid infinite loop
    itnum += 1
    if itnum >= maxit:
        print('Did not find solution after %i iterations' %(itnum))
        break

# print('A thickness of %.3f in will be safe' %(t))

# Array export
# Q_Exp = np.asarray([x_0, q])
# np.savetxt('Data\DistLoadCap.csv', np.transpose(Q_Exp), delimiter=",")