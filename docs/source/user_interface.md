# User Interface

PyProp is accessed as a Python module. Below is an example of how PyProp can be used to optimize a propulsion system using the database of COTS components:

```python
import pyprop

# Optimize for thrust to weight ratio of 0.3
opt = pyprop.Optimizer()
best_unit = opter.random_search(n_units=1000,
                                airspeed=20.0,
                                goal_type="thrust_to_weight_ratio",
                                goal_val=0.3,
                                airframe_weight=1.0,
                                filename="search_results.txt",
                                motor_constraints={"Kv" : [200, 1000]})

# Examine selected performance
best_unit.plot_thrust_curves()
print(best_unit.calc_batt_life(20.0, 4.0))

# Try with a different prop
new_prop = pyprop.create_component_from_database(component="prop", name="apc_12x6", prop_type="data")
new_unit = pyprop.PropulsionUnit(new_prop, best_unit.motor, best_unit.battery, best_unit.esc)

# Examine performance
new_unit.plot_thrust_curves()
print(new_unit.calc_batt_life(20.0, 4.0))
```

PyProp can also be used to analyze the performance of custom propulsion systems:

```python
import pyprop
import numpy as np

# Initialize airfoil coefficient functions
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
my_prop = {
    "airfoils" : {
        "my_airfoil" : {
            "type" : "functional",
            "CL" : airfoil_CL,
            "Cm" : airfoil_Cm,
            "CD" : airfoil_CD
        }
    },
    "geometry" : {
        "n_blades" : 2,
        "hub_radius" : 0.05,
        "weight" : [10.0, "oz"],
        "diameter" : 1.0,
        "geom_pitch" : 0.5,
        "chord" : ["elliptic", 0.75, "in"],
        "rotation" : "CCW",
        "airfoil" : "my_airfoil",
        "grid" : 100
    }
}

# Initialize prop
prop = pyprop.BladeElementProp("my_prop", my_prop)

# Initialize other components
motor = pyprop.Motor(name="EFLITE Power 32", Kv=770, resistance=0.02, I_no_load=2.4, weight=7.0)
esc = pyprop.ESC(name="KDE", I_max=55, resistance=0.0, weight=2.12)
batt = pyprop.Battery(name="Venom 16000mAh", capacity=16000.0, resistance=0.024, voltage=14.1, num_cells=4, weight=49.6, I_max=240)

# Initialize unit
unit = pyprop.PropulsionUnit(prop, motor, batt, esc)

# Evaluate performance
unit.plot_thrust_curves([0.0, 50.0], n_thr=50)
T = unit.calc_cruise_thrust(0.0, 1.0)
print("Max static thrust: {0} lbf".format(T))
t = unit.calc_batt_life(20.0, 0.5*T)
print("{0} min flight time at half of maximum static thrust.".format(t))
```

The user-defined battery, motor, and ESC can be stored in a database for later use. This is done using the write_to_database() method of each component. These components will then be available for use by the optimizers.

In the above example, a custom "functional"-type airfoil is used. For mor information on how to use these, see the documentation for [AirfoilDatabase](airfoildatabase.readthedocs.io/).

## Units
PyProp uses units that are standard to the RC community. In almost all cases, the following units are assumed and cannot be specified otherwise:

| Measurement         | Unit       |
| ------------------- | ---------- |
| Length              | inches     |
| Velocity            | ft/s       |
| Angular deflection  | degrees    |
| Angular rate        | rad/s      |
| Density             | slug/ft^3  |
| Weight/Force        | ounces     |
| Torque              | ft*lbf     |
| Power               | Watts      |
| Motor Kv            | rpm/Volt   |
| Current             | Amps       |
| Resistance          | Ohms       |

For convenience, the function `pyprop.to_rads()` can be used to convert a value in rpm to rad/s.

The only exception to the above defaults is when specifying a BladeElementProp. In this case, the user may specify the inputs in whichever system they like and then PyProp will convert these to English standard units (users familiar with MachUpX will recognize this syntax). For float values, this is done by making the value a list where the first element is the actual value and the second element is a string specifying the units. For example:

    "area" : [6.75, "ft^2"]

For vector inputs, such as position and velocity, the units are simply appended to the vector:

    "velocity" : [100.0, 13.0, 0.0, "ft/s"]

For array inputs, such as a chord distribution, the units are appended as another row in the array:

    "chord" : [[0.0, 1.225],
             [0.5, 1.0066],
             [1.0, 0.81935],
             ["-", "in"]]

When specifying column units in files, these are also appended as another row:

    # File: chord.csv
    0.0, 1.225
    0.5, 1.0066
    1.0, 0.81935
    "-", "in"
    
When specifying units in an array and one of the measurements is dimensionless (e.g. span fraction), "-" should be used.

The following measurements can be defined with the accompanying units:

| Measurement                   | Acceptable Units                  |
| ----------------------------- | --------------------------------- |
| Position/displacement/length: | "ft", "m", "in", "cm"             |
| Area:                         | "ft^2", "m^2"                     |
| Velocity:                     | "ft/s", "m/s", "mph", "kph", "kn" |
| Angular deflection/position:  | "deg", "rad"                      |
| Angular rate:                 | "rpm", "deg/s", "rad/s"           |
| Density:                      | "slug/ft^3", "kg/m^3"             |
| Weight/Force:                 | "lbf", "N", "oz"                  |
| Torque:                       | "ft lbf", "Nm"                    |
| Power:                        | "lbf ft/s", "W", "hp"             |
| Mass:                         | "slug", "kg", "lbm", "g"          |

*Units must be given exactly as they appear in this table. Otherwise, PyProp will not recognize the unit specification.*

**The distinction between units of weight and units of mass is strictly adhered to within PyProp. Please make sure you know when to use which.**

Except for angular measurements, the first unit listed is the default. The default for angular rates is always "rad/s". For airfoil parameters, measurements may only be expressed in radians.  These apply only to inputting a BladeElementProp. Any results generated by PyProp will be output in the units listed in the first table.