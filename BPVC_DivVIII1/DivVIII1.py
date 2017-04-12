# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary

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
sys.stdout = open('log/' + nowtext + '.txt', 'w')
print(datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
print()

# Define all functions
# Use bilinear interpolation equation get A array
def getA( mtrx, arr_dt, dt, arr_ld, ld):

    x2 = (np.abs(arr_dt - dt)).argmin()
    y2 = (np.abs(arr_ld - ld)).argmin()

    # add check for roundup
    if dt >= arr_dt[x2]:
        x2 += 1

    if ld >= arr_ld[y2]:
        y2 += 1

    x1 = x2 - 1
    y1 = y2 - 1

    a_11 = mtrx[y1, x1]
    a_21 = mtrx[y1, x2]
    a_12 = mtrx[y2, x1]
    a_22 = mtrx[y2, x2]

    x = dt
    y = ld

    fx_y1 = ((x2-x)/(x2-x1))*a_11 + ((x-x1)/(x2-x1))*a_21
    fx_y2 = ((x2-x)/(x2-x1))*a_12 + ((x-x1)/(x2-x1))*a_22

    print()
    print('A21 = %.4f   A22 = %.4f' % (a_12, a_22))
    print('A11 = %.4f   A21 = %.4f' % (a_11, a_21))

    if a_11 == 0 or a_21 == 0:
        a = a_22
    else:
        a = ((y2-y)/(y2-y1))*fx_y1 + ((y-y1)/(y2-y1))*fx_y2

    print('A = ', a)

    return a

# Get index
def getB(arr, val):
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
    b = (val-x1)*((y2-y1)/(x2-x1))+y1

    return b

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
t_0 = 1.5                		# initial thickness guess
t_step = 0.001              # step for convergence
maxit = 250                 # max iterations avoid infinite loop


# Calculated initial parameters
S_allow = S_y/SF
L = np.ceil((SG * d_tether * L_tether) / (np.pi * \
                                          (2 * D_o + 3 * d_tether)))
p_req = (2*P)/(D_o * d_tether)

# Create arrays for imported csv data of ASME Div II, Part
# Figure G, to find A
FigG = genfromtxt('csv/IID_FigG.csv', delimiter=',')
# CS2 for S_y > 30 ksi, T<=300 def F
CS2 = genfromtxt('csv/IID_CS2_300F.csv', delimiter=',')

# Create array of unique Do/t and L/Do values from FigG table
FigG_Dt = np.unique(FigG[0:, 0])

# Add line due to importing error Nan
FigG_Dt = FigG_Dt[~np.isnan(FigG_Dt)]

FigG_LD = np.unique(FigG[0:,1])

# Create grid for 3D solution
X_Dt, Y_LD = np.meshgrid(FigG_Dt, FigG_LD)

# get nums of rows and cols
cols = FigG_Dt.__len__()
rows = FigG_LD.__len__()

# Create solution matrix
FigG_A3D = np.zeros((rows, cols))

# loop for processing csv of A data
for i in range(rows):
    for j in range(cols):

        TARGET1 = np.array(np.where((FigG[:, 0] == \
                X_Dt[i, j]) & (FigG[:, 1] == Y_LD[i, j])))

        if not TARGET1:
            continue
        else:
            FigG_A3D[i, j] = FigG[TARGET1[0, 0], 2]

# Initialize loop
p_a = 0
t = t_0
itnum = 1

# Initliaze arrays for summary csv
it = []
t_arr = []
A_arr = []
B_arr = []
Dt_arr = []
LD_arr = []
pa_arr = []

# Create while loop for iteration criterion
while p_a < p_req:

    # check for first interation, oteherwise add the t_step
    if itnum != 1:
        t += t_step

    print('Iteration %i for t = %.3f in' % (itnum, t))
    print('pa = %.1f psi <= preq = %.1f psi' % (p_a, p_req))

    # calc D_0/t ratio for first guess
    Dt = D_o / t
    LD = L / D_o

    # Check for chart applicable aspect ratios
    if Dt >= 4:
        # Check for extreme LD cases
        if LD > 50:
            LD = 50
        # else if statement
        elif LD < 0.05:
            LD =0.05
        else:
            A = getA(FigG_A3D, FigG_Dt, Dt, FigG_LD, LD)
            # A = 0.0015

    else:
        A = 1.1/ (Dt ** 2)

    # Get nearest value in CS2 table
    if A <= CS2[-1, 0]:
        B = getB(CS2, A)
    else:
        # use [-1] to get last element of an array
        B = CS2[-1, 1]

    # for now assume Dot >4
    # Calculate allowable pressure
    p_a = (4*B)/(3 * Dt)

    # Add to summary array
    it.append(itnum)
    t_arr.append(t)
    Dt_arr.append(Dt)
    LD_arr.append(LD)
    A_arr.append(A)
    B_arr.append(B)
    pa_arr.append(p_a)

    # Check if current it >= max
    if itnum >= maxit:
        # multiple blank lines
        print('\n' * 2)
        print('Did not find solution after %i iterations' % itnum)
        break
    # check if solution converged, prints message
    elif p_a >= p_req:
        print('A thickness of %.3f in will be safe' %t)
        print('pa = %.1f psi >= preq = %.1f psi' %(p_a, p_req))

    # equivalent to c++'s i++, inc itnum to avoid infinite loop
    itnum += 1

    # multiple blank lines
    print('\n' * 2)

# Do final stuff

plt_bool = 1
if plt_bool == 1:

    # Plot D/t mass vs L/D vs A Fig G ASME
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X_Dt, Y_LD, FigG_A3D, cmap=cm.coolwarm,\
                           linewidth=0, antialiased=False)
    # Customize the z axis.
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlabel(r'$D_o/t$')
    ax.set_ylabel(r'$L/D_o$')
    ax.set_zlabel(r'$A$')
    plt.title('Visualization of ASME Fig-G')
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig('fig\FigG_3D.png')

    # Convergence plot
    plt.figure(2)
    plt.plot(it, pa_arr)
    plt.xlabel('Iteration ' + r'$n$')
    plt.xlim([0, itnum])
    plt.ylabel('Allowable pressure ' + r'$p_a [psi]$')
    plt.title('Allowable pressure vs iteration number')
    plt.savefig('fig/it_vs_pa.png')

# Array export
Exp_Data = np.asarray([it, t_arr, Dt_arr, LD_arr, A_arr, B_arr, pa_arr])
np.savetxt('data\it_sum_' + nowtext + '.csv', \
           np.transpose(Exp_Data), delimiter=',')