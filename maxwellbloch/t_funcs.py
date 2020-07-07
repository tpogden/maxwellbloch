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
    """ Return a Gaussian pulse time function.

    Notes:
        The amplitude of the Guassian pulse can be defined directly via the
        ampl_{index} argument or indirectly via n_pi_{index}, the desired pulse
        area in multiples of pi. A ValueError will be raised if both are set.
    """
    def func(t, args):
        """
        Args:
            t:      time
            args:   A dict containing fwhm_{index}, centre_{index},
                    and either ampl_{index} OR n_pi_{index}.
        """
        fwhm = args['fwhm_{}'.format(index)]
        centre = args['centre_{}'.format(index)]
        ampl_idx = 'ampl_{}'.format(index)
        n_pi_idx = 'n_pi_{}'.format(index)
        if ampl_idx in args:
            if n_pi_idx in args:
                raise KeyError('t_args can contain ampl or n_pi, not both.')
            else:
                ampl = args[ampl_idx]
        else:
            if n_pi_idx in args:
                n_pi = args[n_pi_idx]
                ampl = n_pi*sqrt(4.*pi*log(2)/fwhm**2)/(2*pi)
            else:
                raise KeyError('t_args must contain ampl or n_pi.')
        return ampl*exp(-4*log(2)*((t - centre)/fwhm)**2)

    func.__name__ = 'gaussian_{}'.format(index)
    return func

def sech(index):
    """ Return a sech pulse time function.

    Notes:
        - The amplitude of the sech pulse can be defined directly via the
          ampl_{index} argument or indirectly via n_pi_{index}, the desired 
          pulse area in multiples of pi. A ValueError will be raised if both are 
          set.
        - The width of the pulse can be defined directly via the width_{index}
          argument or indirectly via the fwhm_{index}, the full-width at half
          max. A ValueError will be raised if both are set.
    """
    def sech_(t): return 2/(exp(t) + exp(-t))
    SECH_FWHM_CONV = 1./2.6339157938

    def func(t, args):
        """
        Args:
            t:      time 
            args:   A dict containing centre_{index}, EITHER ampl_{index} OR 
                n_pi_{index} and EITHER width_{index} OR fwhm_{index}.
        """
        centre = args[f'centre_{index}']
        width_idx = f'width_{index}'
        fwhm_idx = f'fwhm_{index}'
        ampl_idx = f'ampl_{index}'
        n_pi_idx = f'n_pi_{index}'
        if width_idx in args:
            if fwhm_idx in args: 
                raise KeyError('t_args can contain width or fwhm, not both.')
            else:
                width = args[width_idx]
        else:
            if fwhm_idx in args:
               width = args[fwhm_idx]*SECH_FWHM_CONV 
            else:
                raise KeyError('t_args must contain width or fwhm.')
        if ampl_idx in args:
            if n_pi_idx in args:
                raise KeyError('t_args can contain ampl or n_pi, not both.')
            else:
                ampl = args[ampl_idx]
        else:
            if n_pi_idx in args:
                n_pi = args[n_pi_idx]
                ampl = n_pi/width/(2*pi)
            else:
                raise KeyError('t_args must contain ampl or n_pi.')
        return ampl*sech_((t - centre)/width)

    func.__name__ = 'sech_' + str(index)
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

def intp(index):

    def func(t, args):
        tlist = args['tlist_' + str(index)]
        ylist = args['ylist_' + str(index)]
        yintp = interp1d(tlist, ylist, bounds_error=False, fill_value=0.0)
        return yintp(t)

    func.__name__ = 'intp_' + str(index)
    return func
