"""Gather check functions and decorators."""
import functools


# for some reason I have concerns on using type annotations in Python
def argstype(*types):
    """
    Decorate a function to check for args types.

    @argstype(type1, type2, (type3, type4))

    For each function argument provide a type or a tuple of types.
    """
    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            for arg, type_ in zip(args, types):
                if not isinstance(arg, type_):
                    raise TypeError(
                        'expected type {} for {}'.format(type_, arg)
                        )

            return func(*args, **kwargs)
        return wrapper
    return decorator


def kwargstype(*types):
    """
    Decorate a function to check for kwargs types.

    @kwargstype(type1, type2, (type3, type4))

    For each function named argument provide a type or a tuple of types.
    """
    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            errmsg = "expected type {} for named argument '{}', got  {}"

            for kv, type_ in zip(kwargs.items(), types):
                if not isinstance(kv[1], type_):
                    raise TypeError(errmsg.format(type_, kv[0], type(kv[1])))
            return func(*args, **kwargs)
        return wrapper
    return decorator
