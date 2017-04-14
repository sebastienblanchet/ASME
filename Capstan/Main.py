# Sebastien Blanchet 
# Altaeros Energies, Systems Intern
# Capstan pressure

# Import plotting modules
import matplotlib.pyplot as plt
import numpy as np

# Define all parameters
D = 28                    # diameter
d_tether = 15.2/25.4        # tether diameter 15.2mm to [in]
L_tether = (370*1000)/25.4  # tether length overall 370m to [in]
P = 224.809*51.264          # force P of 51.264 kN [lbs]
SF = 1.5                    # safety factor [ul]
S_y = 35000                 # Yield strength 35 ksi
E_y = 2.9*10**7             # Young's modulus 200 [GPa]
v = 0.3                     # poisson's ratio [ul]
SG = 1.01                   # spool gap 1% [ul]
t_0 = 1.5                	# initial thickness guess
t_step = 0.001              # step for convergence
maxit = 250                 # max iterations avoid infinite loop

mu0 = 0.05
mu1 = 0.1
mu2 = 0.25
# mu_low = 0.05
# mu_high = 0.4
# mu_step = 0.05
# mu = np.arange(mu_low, mu_high + mu_step, mu_step)

# Calculated initial parameters
S_allow = S_y/SF
# p_req = (2*P)/(D_o * d_tether)
D_p1 = D+d_tether
D_p2 = D_p1+d_tether
L_p1 = np.pi*D_p1
L_p2 = np.pi*D_p2

n_wraps = L_tether/(L_p1+L_p2)
l = SG*n_wraps*d_tether
l_ws = SG*d_tether

wraps = np.arange(0, 30, 0.1)
x_0 = SG * d_tether * wraps




# Plot 1 : and l vs t
plt.figure(1)
plt.plot(x_0, q, label='Pressure q')
plt.legend()
plt.title('Capstan Equation')
plt.xlabel('Distance from applied force [in]')
plt.ylabel('Pressure q [psi]')
plt.savefig('Figures\p_var_x.png')
plt.show()


# f,ax=plt.subplots(1)
# x=linspace(0,3*pi,1001)
# y=sin(x)
# ax.plot(x/pi,y)
# ax.xaxis.set_major_formatter(FormatStrFormatter('%g $\pi$'))
# ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0))

Q_Exp = np.asarray([x_0, q])
np.savetxt('Data\DistLoadCap.csv', np.transpose(Q_Exp), delimiter=",")







