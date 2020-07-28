"""Duplicates Phillips' "Mechanics of Flight" Example 2.3.1"""

import pyprop
import numpy as np
from pyprop.helpers import to_rpm
import matplotlib.pyplot as plt

# Declare prop getter functions
def airfoil_CL(**kwargs):
    alpha = kwargs.get("alpha", 0.0)
    a_b = np.asarray(alpha+np.radians(2.1))
    return np.where(a_b <= 0.25, 2.0*np.pi*a_b, 0.5*np.pi*np.cos(a_b)/np.cos(0.25))

def airfoil_Cm(**kwargs):
    return 0.0

def airfoil_CD(**kwargs):
    alpha = kwargs.get("alpha", 0.0)
    a_b = np.asarray(alpha+np.radians(2.1))
    return np.where(a_b <= 0.25, 0.224*a_b**2+0.006, np.where(a_b <= 0.3, 16.6944*a_b**2-1.0234, 0.5*np.pi*np.sin(a_b)/np.cos(0.25)))

# Declare prop input
phillips_prop = {
    "airfoils" : {
        "phillips" : {
            "type" : "functional",
            "CL" : airfoil_CL,
            "Cm" : airfoil_Cm,
            "CD" : airfoil_CD
        }
    },
    "geometry" : {
        "n_blades" : 2,
        "hub_radius" : 0.05,
        "weight" : 4.0,
        "diameter" : 1.0,
        "geom_pitch" : 0.5,
        "chord" : ["elliptic", 0.075],
        "rotation" : "CCW",
        "airfoil" : "phillips",
        "grid" : 100
    }
}

# Reproduce Figures 2.3.7-9
num_pitches = 10
num_Js = 100
Js = np.linspace(0.0, 1.4, num_Js)
K_cs = np.linspace(0.3, 1.2, num_pitches)
C_T = np.zeros((num_pitches, num_Js))
C_P = np.zeros((num_pitches, num_Js))
eta = np.zeros((num_pitches, num_Js))

w = 2000 # arbitrary since we have no Reynolds dependence

# Loop through pitches
for i, K_c in enumerate(K_cs):
    phillips_prop["geometry"]["geom_pitch"] = K_c
    prop = pyprop.BladeElementProp("phillips_prop", phillips_prop)

    # Loop through advance ratios
    for j, J in enumerate(Js):
        V = prop.get_velocity(to_rpm(w), J)
        C_T[i,j] = prop.get_thrust_coef(w, V)
        C_P[i,j] = prop.get_power_coef(w, V)

    # Create figure 2.3.5
    if K_c == 0.5:
        prop.plot_angles_over_zeta(w, 0.0)
        V = prop.get_velocity(to_rpm(w), 0.25)
        prop.plot_angles_over_zeta(w, V)

# Calculate efficiency
eta = C_T*Js[np.newaxis,:]/C_P

# Plot thrust coefficient
plt.figure()
for i, K_c in enumerate(K_cs):
    plt.plot(Js, C_T[i,:], label=str(round(K_c, 1)))
plt.xlabel("J")
plt.ylabel("C_T")
plt.legend()
plt.gca().set_xlim([0.0, 1.4])
plt.gca().set_ylim([0.0, 0.11])
plt.show()

# Plot power coefficient
plt.figure()
for i, K_c in enumerate(K_cs):
    plt.plot(Js, C_P[i,:], label=str(round(K_c, 1)))
plt.xlabel("J")
plt.ylabel("C_P")
plt.legend()
plt.gca().set_xlim([0.0, 1.4])
plt.gca().set_ylim([0.0, 0.09])
plt.show()

# Plot efficiency
plt.figure()
for i, K_c in enumerate(K_cs):
    plt.plot(Js, eta[i,:], label=str(round(K_c, 1)))
plt.xlabel("J")
plt.ylabel("Propulsive Efficiency")
plt.legend()
plt.gca().set_xlim([0.0, 1.4])
plt.gca().set_ylim([0.0, 1.0])
plt.show()