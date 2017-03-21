# -*- coding: utf-8 -*-

""" These are time functions provided for using the time-dependent solver.

Q: Why are there multiple versions of each?
A: The solver will want one list of arguments even if there are multiple
time-dependent parts to the Hamiltonian. (Say one laser is ramped on then CW
    and another is a Gaussian pulse.) To distinguish the arguments we've got
    multiple versions. Yes this is wasteful, I would like a better way to do
    it but this works.

Thomas Ogden <t@ogden.eu>
"""

import sys
from numpy import exp, log, sqrt, pi, sinc

def sech(t): return 2 / (exp(t) + exp(-t))

def square_1(t, args):

    on_1 = args['on_1']
    off_1 = args['off_1']
    ampl_1 = args['ampl_1']

    return ampl_1 * (t >= on_1) * (t <= off_1)

def square_2(t, args):

    on_2 = args['on_2']
    off_2 = args['off_2']
    ampl_2 = args['ampl_2']

    return ampl_2 * (t >= on_2) * (t <= off_2)

def square_3(t, args):

    on_3 = args['on_3']
    off_3 = args['off_3']
    ampl_3 = args['ampl_3']

    return ampl_3 * (t >= on_3) * (t <= off_3)

def gaussian_1(t, args):

    ampl_1 = args['ampl_1']
    fwhm_1 = args['fwhm_1']
    centre_1 = args['centre_1']

    return ampl_1 * exp(-4 * log(2) * ((t - centre_1) / fwhm_1)**2)

def gaussian_2(t, args):

    ampl_2 = args['ampl_2']
    fwhm_2 = args['fwhm_2']
    centre_2 = args['centre_2']

    return ampl_2 * exp(-4 * log(2) * ((t - centre_2) / fwhm_2)**2)

def gaussian_3(t, args):

    ampl_3 = args['ampl_3']
    fwhm_3 = args['fwhm_3']
    centre_3 = args['centre_3']

    return ampl_3 * exp(-4 * log(2) * ((t - centre_3) / fwhm_3)**2)

def gaussian_4(t, args):

    ampl_4 = args['ampl_4']
    fwhm_4 = args['fwhm_4']
    centre_4 = args['centre_4']

    return ampl_4 * exp(-4 * log(2) * ((t - centre_4) / fwhm_4)**2)

def gaussian_5(t, args):

    ampl_5 = args['ampl_5']
    fwhm_5 = args['fwhm_5']
    centre_5 = args['centre_5']

    return ampl_5 * exp(-4 * log(2) * ((t - centre_5) / fwhm_5)**2)

def gaussian_6(t, args):

    ampl_6 = args['ampl_6']
    fwhm_6 = args['fwhm_6']
    centre_6 = args['centre_6']

    return ampl_6 * exp(-4 * log(2) * ((t - centre_6) / fwhm_6)**2)

def gaussian_7(t, args):

    ampl_7 = args['ampl_7']
    fwhm_7 = args['fwhm_7']
    centre_7 = args['centre_7']

    return ampl_7 * exp(-4 * log(2) * ((t - centre_7) / fwhm_7)**2)

def ramp_on_1(t, args):

    ampl_1 = args['ampl_1']
    fwhm_1 = args['fwhm_1']
    on_1 = args['on_1']

    return ampl_1 * (exp(-4 * log(2) * ((t - on_1) / fwhm_1)**2) * (t <= on_1) +
                     (t > on_1))

def ramp_on_2(t, args):

    ampl_2 = args['ampl_2']
    fwhm_2 = args['fwhm_2']
    on_2 = args['on_2']

    return ampl_2 * (exp(-4 * log(2) * ((t - on_2) / fwhm_2)**2) * (t <= on_2) +
                     (t > on_2))

def ramp_on_3(t, args):

    ampl_3 = args['ampl_3']
    fwhm_3 = args['fwhm_3']
    on_3 = args['on_3']

    return ampl_3 * (exp(-4 * log(2) * ((t - on_3) / fwhm_3)**2) * (t <= on_3) +
                     (t > on_3))

def ramp_off_1(t, args):

    ampl_1 = args['ampl_1']
    fwhm_1 = args['fwhm_1']
    off_1 = args['off_1']

    return ampl_1 * (exp(-4 * log(2) * ((t - off_1) / fwhm_1)**2) * (t >= off_1) +
                     (t < off_1))

def ramp_off_2(t, args):

    ampl_2 = args['ampl_2']
    fwhm_2 = args['fwhm_2']
    off_2 = args['off_2']

    return ampl_2 * (exp(-4 * log(2) * ((t - off_2) / fwhm_2)**2) * (t >= off_2) +
                     (t < off_2))

def ramp_off_3(t, args):

    ampl_3 = args['ampl_3']
    fwhm_3 = args['fwhm_3']
    off_3 = args['off_3']

    return ampl_3 * (exp(-4 * log(2) * ((t - off_3) / fwhm_3)**2) * (t >= off_3) +
                     (t < off_3))

def ramp_onoff_1(t, args):

    ampl_1 = args['ampl_1']
    fwhm_1 = args['fwhm_1']
    on_1 = args['on_1']
    off_1 = args['off_1']

    ramp_on = (exp(-4 * log(2) * ((t - on_1) / fwhm_1)**2) *
               (t <= on_1) + (t > on_1))

    ramp_off = (exp(-4 * log(2) * ((t - off_1) / fwhm_1)**2) *
                (t >= off_1) + (t < off_1))


    ramp_onoff = ramp_on + ramp_off - 1

    return ampl_1 * ramp_onoff

def ramp_onoff_2(t, args):

    ampl_2 = args['ampl_2']
    fwhm_2 = args['fwhm_2']
    on_2 = args['on_2']
    off_2 = args['off_2']

    ramp_on = (exp(-4 * log(2) * ((t - on_2) / fwhm_2)**2) *
               (t <= on_2) + (t > on_2))

    ramp_off = (exp(-4 * log(2) * ((t - off_2) / fwhm_2)**2) *
                (t >= off_2) + (t < off_2))

    ramp_onoff = ramp_on + ramp_off - 1

    return ampl_2 * ramp_onoff

def ramp_onoff_3(t, args):

    ampl_3 = args['ampl_3']
    fwhm_3 = args['fwhm_3']
    on_3 = args['on_3']
    off_3 = args['off_3']

    ramp_on = (exp(-4 * log(2) * ((t - on_3) / fwhm_3)**2) *
               (t <= on_3) + (t > on_3))

    ramp_off = (exp(-4 * log(2) * ((t - off_3) / fwhm_3)**2) *
                (t >= off_3) + (t < off_3))

    ramp_onoff = ramp_on + ramp_off - 1

    return ampl_3 * ramp_onoff

def ramp_offon_1(t, args):

    ampl_1 = args['ampl_1']
    fwhm_1 = args['fwhm_1']
    off_1 = args['off_1']
    on_1 = args['on_1']

    ramp_on = (exp(-4 * log(2) * ((t - on_1) / fwhm_1)**2) *
               (t <= on_1) + (t > on_1))

    ramp_off = (exp(-4 * log(2) * ((t - off_1) / fwhm_1)**2) *
                (t >= off_1) + (t < off_1))

    ramp_onoff = ramp_on + ramp_off

    return ampl_1 * ramp_onoff

def ramp_offon_2(t, args):

    ampl_2 = args['ampl_2']
    fwhm_2 = args['fwhm_2']
    off_2 = args['off_2']
    on_2 = args['on_2']

    ramp_on = (exp(-4 * log(2) * ((t - on_2) / fwhm_2)**2) *
               (t <= on_2) + (t > on_2))

    ramp_off = (exp(-4 * log(2) * ((t - off_2) / fwhm_2)**2) *
                (t >= off_2) + (t < off_2))

    ramp_onoff = ramp_on + ramp_off

    return ampl_2 * ramp_onoff

def ramp_offon_3(t, args):

    ampl_3 = args['ampl_3']
    fwhm_3 = args['fwhm_3']
    off_3 = args['off_3']
    on_3 = args['on_3']

    ramp_on = (exp(-4 * log(2) * ((t - on_3) / fwhm_3)**2) *
               (t <= on_3) + (t > on_3))

    ramp_off = (exp(-4 * log(2) * ((t - off_3) / fwhm_3)**2) *
                (t >= off_3) + (t < off_3))

    ramp_onoff = ramp_on + ramp_off

    return ampl_3 * ramp_onoff

def sinc_1(t, args):

    ampl = args['ampl_1']
    width = args['width_1']

    return ampl * sinc(width * t) / sqrt(pi / 2.)

def sech_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1 * sech((t - centre_1) / width_1)

def sech_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']

    return ampl_2*sech((t - centre_2)/width_2)

def intp_1(t, args):

    from scipy.interpolate import interp1d

    tlist = args['tlist_1']
    ylist = args['ylist_1']

    yintp = interp1d(tlist, ylist, bounds_error=False, fill_value=0.0)

    return yintp(t)

def intp_2(t, args):

    from scipy.interpolate import interp1d

    tlist = args['tlist_2']
    ylist = args['ylist_2']

    yintp = interp1d(tlist, ylist, bounds_error=False, fill_value=0.0)

    return yintp(t)

def intp_3(t, args):

    from scipy.interpolate import interp1d

    tlist = args['tlist_3']
    ylist = args['ylist_3']

    yintp = interp1d(tlist, ylist, bounds_error=False, fill_value=0.0)

    return yintp(t)
