import pyprop

# Initialize
opter = pyprop.Optimizer()

# Perform random search
#opter.random_search(n_units=1, airspeed=20.0, goal_type="thrust_to_weight_ratio", goal_val=0.3, airframe_weight=1.0)
#opter.random_search(n_units=1, airspeed=20.0, goal_type="power_to_weight_ratio", goal_val=70.0, airframe_weight=1.0)

# Jaden's propulsion unit
motor = pyprop.Motor(type="user_defined", name="EFLITE Power 32", Kv=770, resistance=0.02, I_no_load=2.4, weight=7.0)
esc = pyprop.ESC(type="user_defined", name="KDE", I_max=55, resistance=0.0, weight=2.12)
batt = pyprop.Battery(type="user_defined", name="Venom 16000mAh", capacity=16000.0, resistance=0.006, voltage=3.7, num_cells=4, weight=12.4, I_max=240)
prop = pyprop.DatabaseFitProp(name="apc_12x6")

unit = pyprop.PropulsionUnit(prop, motor, batt, esc)
unit.plot_thrust_curves([0.0, 50.0], n_thr=50)