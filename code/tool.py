### This is a header file to define the tools for assembly bias programs.
### Updated: Aug 13, 2014

import numpy as np

def radius(x, x0, y, y0):
    return np.sqrt((x - x0)**2 + (y - y0)**2)