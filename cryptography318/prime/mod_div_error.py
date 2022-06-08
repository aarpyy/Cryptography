class ModularDivisionError(ArithmeticError):
    """
    Error class for errors caused by attempting to divide strictly Integral
    value by a non-divisor.
    """

    def __init__(self, dividend=None, divisor=None):
        if dividend is not None and divisor is not None:
            super().__init__(f"{divisor} does not divide {dividend}")
        else:
            super().__init__()
