"""Miscellaneous maths utilities."""


def eval_exponent_str(exp_str: str) -> int:
    """Convert a string in the form "X**Y" to an int.

    Will also accept a pure integer string ("X") as input.
    """
    try:
        return int(exp_str)
    except ValueError:
        base, exponent = exp_str.split("**", 1)
        try:
            return int(base) ** int(exponent)
        except ValueError as err:
            raise ValueError(
                f"'{exp_str}' is not a valid exponent string to be evaluated."
            ) from err
