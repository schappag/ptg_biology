simulation_config = {
    # Calcium
    "unit_Ca": "mg/dL",
    "lower_Ca": 4.8,
    "upper_Ca": 5.2,
    "measurement_type": "ionized",

    # Phosphate
    "unit_P": "mg/dL",
    "lower_P": 3.2,
    "upper_P": 4.0,

    # Vitamin D
    "unit_D": "ng/L",
    "lower_D": 18,
    "upper_D": 61,

    # PTH
    "unit_iPTH": "pg/mL",

    # Patient factors
    "vintage": 5,             # dialysis years
    "P_control": "poor",      # poor / good / excellent
    "D_control": "no",        # yes / no

    # Simulation horizon
    "tf": 12,                 # horizon value
    "tf_unit": "m",           # h / d / m

    # Therapy targets
    "target_P": 4,
    "target_D": 180,
    "time_to_target_P": 3,    # months
    "time_to_target_D": 3     # months
}
