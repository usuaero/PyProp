import pyprop

# Initialize random search params
input_dict = {
    "computation" : {
        "units" : 100,
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

# Load
opter = pyprop.Optimizer()
opter.random_search(input_dict)