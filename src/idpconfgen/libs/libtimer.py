"""Function managing time."""
import functools
import time

from idpconfgen import log
from idpconfgen.logger import T, S


def record_time(process_name='', *args, **kwargs):
    """
    Record time of function execution.

    Use as decorator.
    """
    def decorator(func):
        
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            
            start = time.time()
            
            result = func(*args, **kwargs)
            
            log.info(S(f'elapsed time :{process_name}: {time.time() - start}'))
            
            return result
        return wrapper
    return decorator

