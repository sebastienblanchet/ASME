# Sebastien Blanchet 
# Altaeros Energies, Systems Intern
# Capstan pressure

# Import plotting modules
# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import sys


# Runtime
nowtext = datetime.now().strftime('%Y_%m_%d_%H%M')
print(datetime.now().strftime('%Y/%m/%d %H:%M:%S'))
print(nowtext)
print()

# Define all parameters
D = 28                      # diameter
d_tether = 15.2/25.4        # tether diameter 15.2mm to [in]
L_tether = (370*1000)/25.4  # tether length overall 370m to [in]
P = 224.809*51.264          # force P of 51.264 kN [lbs]
SF = 1.5                    # safety factor [ul]
S_y = 35000                 # Yield strength 35 ksi
E_y = 2.9*10**7             # Young's modulus 200 [GPa]
v = 0.3                     # poisson's ratio [ul]
SG = 1.01                   # spool gap 1% [ul]

# Calculated initial parameters
S_allow = S_y/SF
R = D/2

# Length calculations
D_p1 = D+d_tether
D_p2 = D_p1+d_tether
L_p1 = np.pi*D_p1
L_p2 = np.pi*D_p2
n_wraps = L_tether/(L_p1+L_p2)
l = np.ceil(SG*n_wraps*d_tether)
l_ws = SG*d_tether

# Intialize vectors for visual
# mu_low = 0.05
# mu_high = 0.4
# mu_step = 0.05
# mu = np.arange(mu_low, mu_high + mu_step, mu_step)
mu = [0.05, 0.1, 0.25, 0.3, 0.4, 0.5]
theta = np.arange(0.01, (15*2*np.pi), 0.01)
N = np.zeros((len(theta), len(mu)))
p = np.zeros((len(theta), len(mu)))
A = d_tether*R*theta
x = l_ws * (theta/(2*np.pi))

for i in range(0, len(theta)):
    for j in range(0, len(mu)):
        N[i, j] = (P/mu[j])*(1-np.exp(-mu[j]*theta[i]))
        p[i, j] = N[i, j] / A[i]

plt.close('all')
# Plot of N(mu) vs theta
f, axarr = plt.subplots(1)
# Plot for each mu
for k in range(0, len(mu)):
    axarr.plot(theta, N[:, k], label=(r'$\mu=$' + ('%0.2f' % mu[k])))
axarr.set_title('Capstan Equation: Normal Force')
axarr.set_xlabel(r'Contact Angle $\theta$ [radians]')
axarr.set_xlim(0, max(theta))
axarr.legend()
axarr.set_ylabel(r'Normal Force $N$ [lbs]')
plt.savefig('Figures/nvar/' + nowtext + '.png')


# For 2 subplots
# for k in range(0, len(mu)):
#     axarr[0].plot(theta, N[:, k], label=(r'$\mu=$' + ('%0.2f' % mu[k])))
# axarr[0].set_title('Capstan Equation: Normal Force')
# axarr[0].set_xlabel(r'Contact Angle $\theta$ [radians]')
# axarr[0].set_xlim(0, max(theta))
# axarr[0].legend()
# axarr[0].set_ylabel(r'Normal Force $N$ [lbs]')
# axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 5))
# k = 0
# for k in range(0, len(mu)):
#     axarr[1].plot(x, N[:, k], label=(r'$\mu=$' + ('%0.2f' % mu[k])))
# axarr[1].set_xlabel(r'Distance from Applied Load $x$ [in]')
# axarr[1].set_ylabel(r'Normal Force $N$ [lbs]')
# axarr[1].set_xlim(0, max(x))
# axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0, 5))
# plt.savefig('Figures\ ' + nowtext + '_nvar.png')
# plt.tight_layout()
# plt.show()

plt.figure(2)
for h in range(0, len(mu)):
    plt.plot(theta, p[:, h], label=(r'$\mu=$' + ('%0.2f' % mu[h])))
plt.legend()
plt.title('Capstan Equation: Pressure ')
plt.xlabel(r'Contact Angle $\theta$ [radians]')
plt.xlim(0, max(theta))
plt.ylabel(r'Pressure $p$ [psi]')
plt.savefig('Figures\pvar\ ' + nowtext + '.png')
# plt.show()

# Array export
# Q_Exp = np.asarray([x_0, q])
# np.savetxt('Data\DistLoadCap.csv', np.transpose(Q_Exp), delimiter=",")
