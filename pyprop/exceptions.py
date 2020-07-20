"""Custom exceptions in PyProp."""

class UnitsNotInitializedError(Exception):
    """An error thrown when the user has not yet initialized the module units."""

    def __init__(self):
        super().__init__("Units of measurement have not yet been specified. Please initialize using pyprop.initialize_units().")