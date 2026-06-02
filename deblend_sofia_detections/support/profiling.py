"""" This file should only be touched by developers of the code and is used for profiling the code. It contains a decorator that can be used to profile the memory usage of a function. The profiler will write the memory usage to a log file specified by the user. This can be useful for identifying memory leaks and optimizing the code. The profiler can be enabled by setting the PROFILING variable to True. The log file will be created if it does not exist and will be overwritten if it does exist. The log file will contain the memory usage of the function at each line of code. This can be useful for identifying which lines of code are using the most memory and for optimizing those lines of code. The profiler can be disabled by setting the PROFILING variable to False. In this case, the decorator will simply return the original function without any profiling. This can be useful for testing the code without profiling or for running the code in production without profiling."""
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
    


