"""Gather check functions and decorators."""
import functools


# for some reason I have concerns on using type annotations in Python
def argstype(*types):

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
