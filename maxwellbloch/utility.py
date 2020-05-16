import numpy as np
from scipy import interpolate

def half_max_roots(x, y):
    half_max = np.max(y)/2
    spline = interpolate.UnivariateSpline(x, y-half_max, s=0)
    r1, r2 = spline.roots()
    return half_max, r1, r2
    
def full_width_at_half_max(x, y):
    half_max, r1, r2 = half_max_roots(x, y)
    return r2 - r1
