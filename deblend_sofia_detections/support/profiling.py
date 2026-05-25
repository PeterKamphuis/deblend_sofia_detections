
import os

PROFILING = False  # set to True to enable memory profiling

if PROFILING:
    from memory_profiler import profile as _mp_profile
    def profile(log_file):
        def decorator(func):
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            stream = open(log_file, 'w+')
            return _mp_profile(stream=stream)(func)
        return decorator
else:
    def profile(log_file):
        def decorator(func):
            return func
        return decorator
    


