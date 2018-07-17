# -*- coding: utf-8 -*-

""" These closures provide indexed time functions provided for using the 
    time-dependent solver.

    e.g. square(1) will return a function that takes a time t and args 
        {'on_1', 'off_1', 'ampl_1'}

Thomas Ogden <t@ogden.eu>
"""

import sys

from numpy import exp, log, sqrt, pi
from numpy import sinc as npsinc
from scipy.interpolate import interp1d

def square(index):

    def func(t, args):
        on = args['on_' + str(index)]
        off = args['off_' + str(index)]
        ampl = args['ampl_' + str(index)]
        return ampl*(t >= on)*(t <= off)

    func.__name__ = 'square_' + str(index)
    return func

def gaussian(index):

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        fwhm = args['fwhm_' + str(index)]
        centre = args['centre_' + str(index)]
        return ampl * exp(-4 * log(2) * ((t - centre) / fwhm)**2)

    func.__name__ = 'gaussian_' + str(index)
    return func

def ramp_on(index):

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        fwhm = args['fwhm_' + str(index)]
        on = args['on_' + str(index)]
        return ampl * (exp(-4 * log(2) * ((t - on) / fwhm)**2) * (t <= on) +
                        (t > on))

    func.__name__ = 'ramp_on_' + str(index)
    return func

def ramp_off(index):

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        fwhm = args['fwhm_' + str(index)]
        off = args['off_' + str(index)]
        return ampl * (exp(-4 * log(2) * ((t - off) / fwhm)**2) * (t >= off) +
                     (t < off))

    func.__name__ = 'ramp_on_' + str(index)
    return func

def ramp_onoff(index):

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        fwhm = args['fwhm_' + str(index)]
        on = args['on_' + str(index)]
        off = args['off_' + str(index)]
        ramp_on = (exp(-4*log(2)*((t - on)/fwhm)**2) * (t <= on) + (t > on))
        ramp_off = (exp(-4*log(2)*((t - off)/fwhm)**2) *
                    (t >= off) + (t < off))
        ramp_onoff = ramp_on + ramp_off - 1
        return ampl * ramp_onoff

    func.__name__ = 'ramp_onoff_' + str(index)
    return func

def ramp_offon(index):

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        fwhm = args['fwhm_' + str(index)]
        off = args['off_' + str(index)]
        on = args['on_' + str(index)]
        ramp_on = (exp(-4 * log(2) * ((t - on) / fwhm)**2) *
            (t <= on) + (t > on))
        ramp_off = (exp(-4 * log(2) * ((t - off) / fwhm)**2) *
            (t >= off) + (t < off))
        ramp_offon = ramp_on + ramp_off
        return ampl * ramp_offon

    func.__name__ = 'ramp_offon_' + str(index)
    return func

def sinc(index):

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        width = args['width_' + str(index)]
        return ampl * npsinc(width * t) / sqrt(pi / 2.)

    func.__name__ = 'sinc_' + str(index)
    return func

def sech(index):

    def sech_(t): return 2 / (exp(t) + exp(-t))

    def func(t, args):
        ampl = args['ampl_' + str(index)]
        width = args['width_' + str(index)]
        centre = args['centre_' + str(index)]
        return ampl * sech_((t - centre) / width)

    func.__name__ = 'sech_' + str(index)
    return func

def intp(index):

    def func(t, args):
        tlist = args['tlist_' + str(index)]
        ylist = args['ylist_' + str(index)]
        yintp = interp1d(tlist, ylist, bounds_error=False, fill_value=0.0)
        return yintp(t)

    func.__name__ = 'intp_' + str(index)
    return func
