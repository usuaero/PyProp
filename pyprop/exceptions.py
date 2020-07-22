"""Custom exceptions in PyProp."""

class MaxCurrentExceededError(Exception):
    """An error thrown when the max current that can be handled by one of the components is exceeded.

    Members
    -------
    component_type : str
        The type of component ("motor", "ESC", or "battery").

    component_name : str
        The name of the component.

    I_max : float
        The maximum current allowed through the component.

    I_drawn : float
        The amount of current trying to be drawn.
    """

    def __init__(self, component_type, component_name, I_max, I_drawn):
        self.component_type = component_type
        self.component_name = component_name
        self.I_max = I_max
        self.I_drawn = I_drawn
        super().__init__("{0} were drawn through the {1}, {2}. The maximum allowed is {3}.".format(self.I_drawn, self.component_type, self.component_name, self.I_max))


class ThrottleNotFoundError(Exception):
    """An error that occurs when the solver for the throttle value does not converge or the throttle
    value determined is not possible.

    Members
    -------
    error_type : str
        The type of error. Either "not_converged" or "throttle_invalid".

    v_cruise : float
        The cruise velocity in ft/s.

    T_req : float
        Required thrust in lbf.

    final_val : float
        The throttle value converged to (if "throttle_invalid").
    """

    def __init__(self, error_type, v_cruise, T_req, final_val=None):
        self.error_type = error_type
        self.v_cruise = v_cruise
        self.T_req = T_req
        self.final_val = final_val

        if self.error_type == "not_converged":
            super().__init__("The solver failed to converge to a throttle value.")
        elif self.error_type == "throttle_invalid":
            super().__init__("The solver converged to a throttle value outside the range of 0 to 1. Final value: {0}%".format(self.final_val*100.0))


class TorquesNotMatchedError(Exception):
    """An error that occurs when the propeller and motor torques cannot be matched at a speed.

    Members
    -------
    v_cruise : float
        Cruise velocity in ft/s.

    throttle : float
        Throttle value.
    """

    def __init__(self, v_cruise, throttle):
        self.v_cruise = v_cruise
        self.throttle = throttle

        super().__init__("The solver could not match the propeller and motor torques at the given conditions.")


class InvalidRuntimeError(Exception):
    """An error that occurs when a negative runtime is predicted.

    Members
    -------
    runtime : float
        The predicted runtime in minutes.
    """

    def __init__(self, runtime):
        self.runtime = runtime

        super().__init__("A runtime of {0} minutes was predicted, which is invalid.".format(self.runtime))


class DatabaseRecordNotFoundError(Exception):
    """"An error which occurs when the user tries to pull a record from the database which does not exist.

    Method
    ------
    command : str
        The sql command which failed.
    """

    def __init__(self, command):
        self.command = command

        super().__init__("The command '{0}' failed to retrieve a component from the database.".format(command))