# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary

# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

# Define all parameters
D = 28                      # diameter
d_tether = 15.2/25.4        # tether diameter 15.2mm to [in]
L_tether = (370*1000)/25.4  # tether length overall 370m to [in]
P = 224.809*51.264          # force P of 51.264 kN [lbs]
SF = 1.5                    # safety factor [ul]
S_y = 35000                 # Yield strength 35 ksi
E_y = 2.9*10**7             # Youngs modulus 200 [GPa]
v = 0.3                     # poissons ratio [ul]
rho = 0.284                 # density of steel [lb/in^3]
SG = 1.01                   # spool gap 1% [ul]
t_0 = 0.5                   # initial thickness guess

# Calculated initial parameters
S_allow = S_y/SF
l = np.ceil((SG * d_tether * L_tether) / (np.pi * (2 * D + 3 * d_tether)))
q = (2*P)/(D*d_tether)

# Create arrays for imported csv data of ASME Div II, Part
# Figure G, to find A
FigG = genfromtxt('csv/IID_FigG.csv', delimiter=',')
# CS2 for S_y > 30 ksi, T<=300 def F
CS2 = genfromtxt('csv/IID_CS2_300F.csv', delimiter=',')

# Get index
# help from MATLAB to py https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html
def findval(arr, val, col):
    # Get smallest difference between array and val
    val_idx = (np.abs(arr[:, col] - val)).argmin()

    # add check for roundup
    if val >= arr[val_idx,col]:
        val_idx=val_idx+1

    return val_idx

def interpol_xy (x, x1, x2, y1, y2):
    y = (x-x1)*((y2-y1)/(x2-x1))+y1
    return y

A = 0.002

ib = findval(CS2, A, 0)

print(ib)

A_i1 = CS2[ib - 1, 0]
A_i2 = CS2[ib, 0]
B_i1 = CS2[ib - 1, 1]
B_i2 = CS2[ib, 1]

B = interpol_xy(A, A_i1, A_i2, B_i1, B_i2)
print(B)




