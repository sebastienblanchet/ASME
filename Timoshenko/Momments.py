# Sebastien Blanchet, Altaeros Energies, Systems Engineering Intern
# Note: all calculations to be completed in US customary
# Thin shell analysis of 28" drum

# Import plotting modules
import numpy as np
import matplotlib.pyplot as plt

# Metric parameters
d_SI_tether = 15.2
L_SI_tether = 370
P_SI = 51.264

# Define all parameters
d_tether = d_SI_tether/25.4        		# tether diameter 15.2mm to [in]
L_tether = (L_SI_tether*1000)/25.4  	# tether length overall 370m to [in]
SF = 1.5                    			# safety factor [ul]
S_y = 50763                 			# Yield strength 350 [MPa]
E_y = 2.9*10**7             			# Youngs modulus 200 [GPa]
v = 0.3                     			# poissons ratio [ul]
rho = 0.284                 			# density of steel [lb/in^3]
SG = 1.01                   			# spool gap 1% [ul]
D = 28                      			# OD 28 in drum for thin shell
t = 1+(11/16)							# safe thickness as per ASME VIII-1

# Prelim calcs
P = 224.809*2*np.pi*P_SI          		# force P of 51.264 kN [lbs]
Sig_allow = S_y / SF
R = D/2
l = np.ceil((SG * d_tether * L_tether) / (np.pi * (2 * D + 3 * d_tether)))
FlexRig = (E_y*t**3)/(12*(1-v**3))
Beta = ((3*(1-v**2))/((R**t)**2))**0.25

# Create winch class
class Drum:
    # Initialize object to be empty (__ means private)
    __OD = None
    __Length = None
    __Thickness = None

    # Define constructor
    def __init__(self, OD):
        self.__OD = OD

    def set_OD(self,OD):
        self.__OD


x_step = 0.125
x = np.arange(0, l/2 + x_step, x_step)

Bx = Beta*x
Phi_Bx = (np.exp(-Bx))*(np.cos(Bx)+np.sin(Bx))
Psi_Bx = (np.exp(-Bx))*(np.cos(Bx)-np.sin(Bx))
Theta_Bx = (np.exp(-Bx))*(np.cos(Bx))
Zeta_Bx = (np.exp(-Bx))*(np.sin(Bx))

# Results for P applied at center
w_x = (P/(8*(Beta**3)*FlexRig))*Phi_Bx
M_x = (P/(4*Beta))*Psi_Bx
Q_x = (P/2)*Theta_Bx

plt.close('all')
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(3, sharex=True)
# Plot w(x)
axarr[0].plot(x, w_x)
axarr[0].set_title('Deflection, Moment, Shear Distribution')
axarr[0].set_ylabel(r'$w(x)\ [in]$')
# Plot M(x)
axarr[1].plot(x, M_x)
axarr[1].set_ylabel(r'$M(x)\ [lbs\cdot in]$')
# Plot Q(x)
axarr[2].plot(x, Q_x)
axarr[2].set_ylabel(r'$Q(x)\ [lbs]$')
axarr[2].set_xlabel(r'$x\ [in]$')
axarr[2].set_xlim([0, l/4])
# Save fig and show
plt.savefig('Figures\wMQ_vs_x.png')
plt.show()