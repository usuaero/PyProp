"""Tests the BladeElementProp class. Does this work?"""

import numpy as np
import pyprop
from pyprop.helpers import to_rpm, to_rads

def test_phillips_2_3_1():
    # Tests result is the same as Example 2.3.1 in Phillips' Mech. of Flight

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
    num_pitches = 5
    num_Js = 5
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

    true_eta = np.asarray([[0., 0.54478846, 1.20303466, 1.22594647, 1.25060314],
                           [0., 0.7565339 , 1.16017598, 1.14345178, 1.17712644],
                           [0., 0.67360335, 0.88773116, 1.09397494, 1.12164942],
                           [0., 0.53434467, 0.85811094, 0.69144413, 1.08282994],
                           [0., 0.37232873, 0.81087786, 0.91686828, 1.11095682]])
    assert np.allclose(eta, true_eta)

    true_C_T = np.asarray([[0.04275367, 0.0035241 , -0.05944844, -0.13747256, -0.22662599],
                           [0.07170009, 0.03766199, -0.01659641, -0.08368984, -0.16028094],
                           [0.09189391, 0.06917407,  0.02143569, -0.03709349, -0.10372785],
                           [0.0989737 , 0.0909446 ,  0.0546933 ,  0.00296234, -0.055687  ],
                           [0.09460433, 0.09929382,  0.0836305 ,  0.03738553, -0.01477035]])
    assert np.allclose(C_T, true_C_T)

    true_C_P = np.asarray([[0.00913487, 0.00226406, -0.03459078, -0.11774265, -0.2536987 ],
                           [0.02172797, 0.0174238 , -0.01001356, -0.07685005, -0.19062804],
                           [0.04248322, 0.03594241,  0.01690262, -0.03560243, -0.12946915],
                           [0.07690066, 0.05956943,  0.0446158 ,  0.00449849, -0.0719982 ],
                           [0.11960517, 0.09333912,  0.07219503,  0.04281401, -0.01861322]])
    assert np.allclose(C_P, true_C_P)