import pyprop

# Initialize
opter = pyprop.Optimizer()

# Perform random search
#opter.random_search(n_units=1000,
#                    airspeed=20.0,
#                    goal_type="thrust_to_weight_ratio",
#                    goal_val=0.3,
#                    airframe_weight=1.0,
#                    filename="search_results.txt",
#                    prop_constraints={"prop_type" : "data"},
#                    motor_constraints={"Kv" : [200, 1000]},
#                    battery_constraints={"num_cells" : [4, 5]})
#opter.random_search(n_units=1000,
#                    airspeed=20.0,
#                    goal_type="power_to_weight_ratio",
#                    goal_val=70.0,
#                    airframe_weight=1.0,
#                    filename="power_search_results.txt",
#                    prop_constraints={
#                        "prop_type" : "data",
#                        "diameter" : [17, 21],
#                        "pitch" : [6, 12]
#                    },
#                    motor_constraints={
#                        "Kv" : [100, 550]
#                    })

# Jaden's propulsion unit
motor = pyprop.Motor(name="EFLITE Power 32", Kv=770, resistance=0.02, I_no_load=2.4, weight=7.0)
#motor = pyprop.create_component_from_database(component="motor", Kv=770)
esc = pyprop.ESC(name="KDE", I_max=55, resistance=0.0, weight=2.12)
batt = pyprop.Battery(name="Venom 16000mAh", capacity=16000.0, resistance=0.024, voltage=14.1, num_cells=4, weight=49.6, I_max=240)
#batt = pyprop.create_component_from_database(component="battery", capacity=[2000, 5000])
#prop = pyprop.create_component_from_database(component="prop", diameter=[11, 16], pitch=[4, 9], manufacturer="APC")
#prop = pyprop.create_component_from_database(component="prop", type="data", diameter=[11, 16], pitch=[4, 9], manufacturer="APC")
prop = pyprop.create_component_from_database(component="prop", name="apc_12x6")

unit = pyprop.PropulsionUnit(prop, motor, batt, esc)
#unit.plot_thrust_curves([0.0, 50.0], n_thr=50)
T = unit.calc_cruise_thrust(0.0, 1.0)
print("Static thrust: {0} lbf".format(T))
t = unit.calc_batt_life(69.0, 4.4)
print("Cruise flight time: {0} min".format(t))