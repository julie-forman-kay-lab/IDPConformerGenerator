"""Contain functions."""
from contextlib import contextmanager
from functools import partial, reduce
from operator import (
    lt, le, eq, ne, ge, gt,
    add, sub, truediv, floordiv, mul, matmul, mod, pow, neg, pos,
    contains, concat,
    and_, or_, xor, not_, truth,
    invert, lshift, rshift,
    is_, is_not,
    setitem, getitem, delitem,
    iadd,
    )  # noqa: F401


def vartial(func, *args, **keywords):
    """
    Prepare a function with args and kwargs except for the first arg.

    Functions like `functools.partial` except that the resulting
    preprepared function expects the first positional argument.

    Example
    -------
    >>> pow2 = vartial(math.pow, 2)
    >>> pow2(3)
    9
    >>> pow2(4)
    16

    This is different from:
    >>> pow_base_3 = partial(math.pow, 3)
    >>> pow_base_3(2)
    9
    >>> pow_base_3(4)
    81
    """
    def newfunc(value, *fargs, **fkeywords):
        # newkeywords = {**keywords, **fkeywords}
        newkeywords = keywords.copy()
        newkeywords.update(fkeywords)
        return func(value, *args, *fargs, **newkeywords)
    newfunc.func = func
    newfunc.args = args
    newfunc.keywords = keywords
    return newfunc


def give(value):
    """
    Preare a function to return a value when called.

    Ignore `*args` and `**kwargs`.

    Example
    -------
    >>> true = give(True)
    >>> true()
    True

    >>> five = give(5)
    >>> five(4, 6, 7, 8, some_args='some string')
    5
    """
    def newfunc(*ignore, **everything):
        return value
    newfunc.value = value
    return newfunc


# predefined set ups
true = give(True)
false = give(False)
none = give(None)


def pass_(value):
    """
    Do nothing, just pass the value.

    Example
    -------
    >>> pass_(1)
    1
    """
    return value


def f1f2(f1, f2, *a, **k):
    """
    Apply one function after the other.

    Call `f1` on the return value of `f2`.

    Args and kwargs apply to `f2`.

    Example
    -------
    >>> f1f2(str, int, 2)
    "2"
    """
    return f1(f2(*a, **k))


def f2f1(f1, f2, *a, **k):
    """
    Apply the second function after the first.

    Call `f2` on the return value of `f1`.

    Args and kwargs apply to `f1`.

    Example
    -------
    >>> f2f1(str, int, 2)
    2
    """
    return f2(f1(*a, **k))


def ternary_operator(iflogic, assertion, elselogic):
    """
    Apply ternary operator logic executing functions.

    Functions should be preconfigured and accept no arguments.

    Better if you see the code:

    `return iflogic() if assertion() else elselogic()`
    """
    return iflogic() if assertion() else elselogic()


def ternary_operator_v(x, iflogic, assertion, elselogic):
    """
    Apply ternary operator logic executing functions.

    Functions should receive a single value: `x`.

    Better if you see the code:

    `return iflogic(x) if assertion(x) else elselogic(x)`

    Parameters
    ----------
    x
        The value to pass to each function.
    """
    return iflogic(x) if assertion(x) else elselogic(x)


def make_iterable(value):
    """
    Transform into an iterable.

    Transforms a given `value` into an iterable if it is not.
    Else, return the value itself.

    Example
    -------
    >>> make_iterable(1)
    [1]

    >>> make_iterable([1])
    [1]
    """
    try:
        iter(value)
    except TypeError:
        return [value]
    else:
        return value


def reduce_helper(value, f, *a, **k):
    """
    Help in `reduce`.

    Helper function when applying `reduce` to a list of functions.

    Parameters
    ----------
    value : anything

    f : callable
        The function to call. This function receives `value` as first
        positional argument.

    *a, **k
        Args and kwargs passed to `f`.
    """
    return f(value, *a, **k)


def chainf(init, *funcs):
    """
    Run functions in sequence starting from an initial value.

    Example
    -------
    >>> chainf(2, [str, int, float])
    2.0
    """
    return reduce(reduce_helper, funcs, init)


def chainfs(*funcs):
    """
    Store functions be executed on a value.

    Example
    -------
    >>> do = chainfs(str, int, float)
    >>> do(2)
    2.0
    """
    def execute(value):
        return chainf(value, *funcs)
    return execute


def if_elif_else(value, condition_function_pair):
    """
    Apply logic if condition is True.

    Parameters
    ----------
    value : anything
        The initial value

    condition_function_pair : tuple
        First element is the assertion function, second element is the
        logic function to execute if assertion is true.

    Returns
    -------
    The result of the first function for which assertion is true.
    """
    for condition, func in condition_function_pair:
        if condition(value):
            return func(value)


def whileloop(cond, func, do_stopiteration=none, do_exhaust=none):
    """
    Execute while loop.

    All function accept no arguments. If state needs to be evaluated,
    `cond` and `func` need to be synchronized.

    Parameters
    ----------
    cond : callable
        The `while` loop condition.

    func : callable
        The function to call on each while loop iteration.

    do_stopiteration : callable
        The function to execute when `func` raises StopIteration error.

    do_exhaust : callable
        The function to execute when while loop exhausts.

    Return
    ------
    None
    """
    while cond():
        try:
            func()
        except StopIteration:
            do_stopiteration()
            return
    do_exhaust()
    return


def consume(gen):
    """
    Consume generator in a single statement.

    Example
    -------
    >>> consume(generator)
    """
    for _ in gen:
        pass


def mapc(f, *iterables):
    """
    Consume map function.

    Like `map()` but it is not a generator; `map` is consumed
    immediately.
    """
    return consume(map(f, *iterables))


def flatlist(list_):
    """
    Flat a list recursively.

    This is a generator.
    """
    # escape strings which would yield infinite recursion
    if isinstance(list_, str):
        yield list_
    else:
        try:
            for sublist in list_:
                yield from flatlist(sublist)
        except TypeError:  # sublist is not iterable
            yield sublist


def raise_(exception, *ignore, **everything):
    """Raise exception."""
    raise exception


@contextmanager
def context_engine(
        func,
        exceptions,
        doerror,
        doelse,
        dofinally,
        *args,
        **kwargs):
    """Make a context engine."""
    try:
        result = func(*args, **kwargs)

    except (exceptions) as err:
        doerror(err, *args, **kwargs)

    else:
        yield result
        doelse(result, *args, **kwargs)

    finally:
        dofinally(*args, **kwargs)


# If Then Else
ite = ternary_operator
itev = ternary_operator_v

# conditionals
is_none = partial(is_, None)
is_not_none = partial(is_not, None)

# To deprecate in version 1
ITE = ternary_operator
ITEX = ternary_operator_v
ternary_operator_x = ternary_operator_v
# chainfs, can be replaced with vartial(chainf, f1, f2, f3)
