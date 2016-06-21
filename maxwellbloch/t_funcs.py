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

sech = lambda t: 2/(exp(t) + exp(-t))

def square_1(t, args):

    on_1 = args['on_1']
    off_1 = args['off_1']
    ampl_1 = args['ampl_1']

    return ampl_1*(t >= on_1)*(t <= off_1)

def square_2(t, args):

    on_2 = args['on_2']
    off_2 = args['off_2']
    ampl_2 = args['ampl_2']

    return ampl_2*(t >= on_2)*(t <= off_2)

def square_3(t, args):

    on_3 = args['on_3']
    off_3 = args['off_3']
    ampl_3 = args['ampl_3']

    return ampl_3*(t >= on_3)*(t <= off_3)

def gaussian_fwhm_1(t, args):

    ampl_1 = args['ampl_1']
    fwhm_1 = args['fwhm_1']
    centre_1 = args['centre_1']

    return ampl_1*exp(-4*log(2)*((t - centre_1)/fwhm_1)**2)    

def gaussian_fwhm_2(t, args):

    ampl_2 = args['ampl_2']
    fwhm_2 = args['fwhm_2']
    centre_2 = args['centre_2']

    return ampl_2*exp(-4*log(2)*((t - centre_2)/fwhm_2)**2)    

def gaussian_fwhm_3(t, args):

    ampl_3 = args['ampl_3']
    fwhm_3 = args['fwhm_3']
    centre_3 = args['centre_3']

    return ampl_3*exp(-4*log(2)*((t - centre_3)/fwhm_3)**2)    

def gaussian_fwhm_4(t, args):

    ampl_4 = args['ampl_4']
    fwhm_4 = args['fwhm_4']
    centre_4 = args['centre_4']

    return ampl_4*exp(-4*log(2)*((t - centre_4)/fwhm_4)**2)    

def gaussian_fwhm_5(t, args):

    ampl_5 = args['ampl_5']
    fwhm_5 = args['fwhm_5']
    centre_5 = args['centre_5']

    return ampl_5*exp(-4*log(2)*((t - centre_5)/fwhm_5)**2)    

def gaussian_fwhm_6(t, args):

    ampl_6 = args['ampl_6']
    fwhm_6 = args['fwhm_6']
    centre_6 = args['centre_6']

    return ampl_6*exp(-4*log(2)*((t - centre_6)/fwhm_6)**2)    

def gaussian_fwhm_7(t, args):

    ampl_7 = args['ampl_7']
    fwhm_7 = args['fwhm_7']
    centre_7 = args['centre_7']

    return ampl_7*exp(-4*log(2)*((t - centre_7)/fwhm_7)**2)    

def ramp_on_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1*(exp(-2*log(2)*((t - centre_1)/width_1)**2)*(t <= centre_1) +
                   (t > centre_1))

def ramp_on_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']

    return ampl_2*(exp(-2*log(2)*((t - centre_2)/width_2)**2)*(t <= centre_2) +
                   (t > centre_2))

def ramp_on_3(t, args):

    ampl_3 = args['ampl_3']
    width_3 = args['width_3']
    centre_3 = args['centre_3']

    return ampl_3*(exp(-2*log(2)*((t - centre_3)/width_3)**2)*(t <= centre_3) +
                   (t > centre_3))

def ramp_off_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1*(exp(-2*log(2)*((t - centre_1)/width_1)**2)*(t >= centre_1) +
                   (t < centre_1))

def ramp_off_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']

    return ampl_2*(exp(-2*log(2)*((t - centre_2)/width_2)**2)*(t >= centre_2) +
                   (t < centre_2))

def ramp_off_3(t, args):

    ampl_3 = args['ampl_3']
    width_3 = args['width_3']
    centre_3 = args['centre_3']

    return ampl_3*(exp(-2*log(2)*((t - centre_3)/width_3)**2)*(t >= centre_3) +
                   (t < centre_3))

def ramp_onoff_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    on_1 = args['on_1']
    off_1 = args['off_1']

    ramp_on = (exp(-2*log(2)*((t - on_1)/width_1)**2)*
                        (t <= on_1) + (t > on_1))

    ramp_off = (exp(-2*log(2)*((t - off_1)/width_1)**2)*
                        (t >= off_1) + (t < off_1))

    return ampl_1*(ramp_on + ramp_off - 1.)

def ramp_onoff_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    on_2 = args['on_2']
    off_2 = args['off_2']

    ramp_on = (exp(-2*log(2)*((t - on_2)/width_2)**2)*
                        (t <= on_2) + (t > on_2))

    ramp_off = (exp(-2*log(2)*((t - off_2)/width_2)**2)*
                        (t >= off_2) + (t < off_2))

    return ampl_2*(ramp_on + ramp_off - 1.)

def ramp_onoff_3(t, args):

    ampl_3 = args['ampl_3']
    width_3 = args['width_3']
    on_3 = args['on_3']
    off_3 = args['off_3']

    ramp_on = (exp(-2*log(2)*((t - on_3)/width_3)**2)*
                        (t <= on_3) + (t > on_3))

    ramp_off = (exp(-2*log(2)*((t - off_3)/width_3)**2)*
                        (t >= off_3) + (t < off_3))

    return ampl_3*(ramp_on + ramp_off - 1.)

def sinc_1(t, args):

    ampl = args['ampl_1']
    width = args['width_1']

    return ampl*sinc(width*t)/sqrt(pi/2.)

def sech_1(t, args):

    ampl_1 = args['ampl_1']
    width_1 = args['width_1']
    centre_1 = args['centre_1']

    return ampl_1*sech((t - centre_1)/width_1)

def sech_2(t, args):

    ampl_2 = args['ampl_2']
    width_2 = args['width_2']
    centre_2 = args['centre_2']
    
    return ampl_2*sech((t - centre_2)/width_2)
    