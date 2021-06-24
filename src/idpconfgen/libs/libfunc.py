"""
Functools-like mode.

Abstract implementations, engines, and alike.
"""


def vartial(func, *args, **kwargs):
    """
    Reimplements functools.partial saving the initial value.

    `vartial` is like `partial` except that the first argument is not
    completed in the vartialized function and is expected upon call.

    Example
    -------
    def func(value, par1, *pars)
    >>> myvartial = vartial(func, arg1, arg2, arg3)
    myvartial(value)
    """
    def newfunction(value, *fargs, **fkwargs):
        newkwargs = kwargs.copy()
        newkwargs.update(fkwargs)
        return func(value, *args, *fargs, **newkwargs)
    newfunction.func = func
    newfunction.args = args
    newfunction.kwargs = kwargs
    return newfunction


def are_all_valid(value, *actions):
    """
    Assert the validity of `value` to all actions.

    Actions must receive a single argument: the `value`.

    Return
    ------
    bool
        True if all `actions` evaluate to True. False otherwise.
    """
    return all(f(value) for f in actions)


def negate(x):
    return not bool(x)


def post_func(x, y, *args, **kwargs):
    return x(y(*args, **kwargs))
