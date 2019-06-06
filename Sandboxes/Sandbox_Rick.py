import matplotlib.pyplot as plt
from scipy.interpolate import spline
import numpy as np

alpha = np.array([-10, -5, 0, 5, 10, 15, 19])
alphatip = np.array([-8, -3, 2, 7, 12, 17, 21])
# NACA 4415
CL_root = np.array([-0.75, -0.2, 0.35, 0.9, 1.31, 1.53, 1.5])
CD_root = np.array([0.0135, 0.008, 0.0075, 0.008, 0.014, 0.034, 0.087])

# NACA 4412
CL_tip = np.array([-0.6, -0.05, 0.5, 1.03, 1.44, 1.62, 1.58])
CD_tip = np.array([0.015, 0.0085, 0.0075, 0.008, 0.018, 0.04, 0.09])

# Plot smoothening
alphanew = np.linspace(alpha.min(), alpha.max(), 200)
alphatipnew = np.linspace(alphatip.min(), alphatip.max(), 200)
CLr = spline(alpha, CL_root, alphanew)
CDr = spline(alpha, CD_root, alphanew)
CLt = spline(alphatip, CL_tip, alphatipnew)
CDt = spline(alpha, CD_tip, alphanew)
#CDr, CLt, CDt = spline(alpha,CL_root,6), spline(alpha,CD_root), spline(alpha,CL_tip), spline(alpha,CD_tip)

plt.plot(alphanew, CLr, '-b', label = 'Airfoil at root')
plt.plot(alphanew, CLt, '--g', label = 'Airfoil at tip')
plt.plot(alphatipnew, CLt, 'sr', label = 'Airfoil at tip, with 2 degrees outwash')
plt.xlim(-7,18)
plt.xlabel(r'$\alpha$ [deg]', size = 18)
plt.ylabel(r'$C_L$ [-]', size = 18)
plt.legend(loc= 'bottom right')
plt.show()

plt.plot(CDr, CLr, '-b', label = 'Airfoil at root')
plt.plot(CDt, CLt, '--g', label = 'Airfoil at tip')
plt.xlim(0,0.04)
plt.ylim(-0.8,1.7)
plt.xlabel(r'$C_D$ [-]', size = 18)
plt.ylabel(r'$C_L$ [-]', size = 18)
plt.legend(loc= 'top left')
plt.show()

plt.plot(alphanew, CLr/CDr, '-b', label = 'Airfoil at root')
plt.plot(alphanew, CLt/CDt, '--g', label = 'Airfoil at tip')
plt.plot(alphatipnew, CLt/CDt, 'sr', label = 'Airfoil at tip, with 2 degrees outwash')
plt.xlim(-7,18)
plt.ylim(-50,140)
plt.xlabel(r'$\alpha$ [deg]', size = 18)
plt.ylabel(r'$C_L / C_D$ [-]', size = 18)
plt.legend(loc = 'top left')
plt.show()
