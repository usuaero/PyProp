import pyprop

# Initialize
#opter = pyprop.Optimizer()

# Perform random search
#opter.random_search(n_units=1000, airspeed=20.0, goal_type="thrust_to_weight_ratio", goal_val=0.3, airframe_weight=1.0, filename="search_results.txt")
#opter.random_search(n_units=10, airspeed=20.0, goal_type="power_to_weight_ratio", goal_val=70.0, airframe_weight=1.0)

## Jaden's propulsion unit
#motor = pyprop.Motor(name="EFLITE Power 32", Kv=770, resistance=0.02, I_no_load=2.4, weight=7.0)
#esc = pyprop.ESC(name="KDE", I_max=55, resistance=0.0, weight=2.12)
#batt = pyprop.Battery(name="Venom 16000mAh", capacity=16000.0, resistance=0.024, voltage=14.1, num_cells=4, weight=49.6, I_max=240)
##prop = pyprop.create_component_from_database(component="fit_prop", name="apc_12x6")
#prop = pyprop.DatabaseDataProp("apc_12x6")
#
#unit = pyprop.PropulsionUnit(prop, motor, batt, esc)
#unit.plot_thrust_curves([0.0, 50.0], n_thr=50)
#T = unit.calc_cruise_thrust(0.0, 1.0)
#print("{0} lbf".format(T))

## Database data prop
#prop = pyprop.DatabaseDataProp("apc_11x3")
#prop.plot_coefs()

# BET Prop
prop_file = "dev/test_prop.json"
prop = pyprop.BladeElementProp("test_prop", prop_file)
prop.plot_coefs()
C_T = prop.get_thrust_coef(1256.637, 50.0)
print(0.0023769*(1256.637/2*3.1415926535)**2*prop.diameter**4*C_T)