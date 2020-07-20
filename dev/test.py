import pyprop

# Initialize random search params
input_dict = {
    "computation" : {
        "units" : 100000,
        "processes" : 8,
        "outlierStdDevs" : 10
    },
    "condition" : {
        "altitude" : 0,
        "airspeed" : 10
    },
    "goal" : {
        "thrust" : 0,
        "thrustToWeightRatio" : 0.3
    },
    "aircraft" : {
        "emptyWeight" : 1,
        "components" : {
            "propeller" : {
                "name" : "",
                "manufacturer" : ""
            },
            "motor" : {
                "name" : "",
                "manufacturer" : ""
            },
            "esc" : {
                "name" : "",
                "manufacturer" : ""
            },
            "battery" : {
                "name" : "",
                "manufacturer" : ""
            }
        }
    }
}

# Initialize
pyprop.initialize_units(unit_sys="SI")
batt = pyprop.Battery(type="user_defined", cell_weight=[50, "N"])

# Load
opter = pyprop.Optimizer()
opter.random_search(input_dict)