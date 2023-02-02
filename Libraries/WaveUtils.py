#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 18:16:43 2022

@author: auclaije
"""
import warnings
import numpy as np
from scipy import optimize
from pars import rho_w, rho_i, g, E, v


def PM(u, f):
    alpha_s = 0.2044
    beta_s  = 1.2500

    fp = 0.877 * g / (2 * np.pi * u)
    Hs = 0.0246 * u**2
    PM=alpha_s*Hs**2*(fp**4/f**5)*np.exp(-beta_s*(fp/f)**4)

    return PM


def Jonswap(Hs, fp, f):
    sigma_s = 0.09 * np.ones_like(f)
    sigma_s[f <= fp] = 0.07

    # Spectral default parameters
    alpha_s = 0.2044
    beta_s  = 1.2500
    gamma_s = 3.3

    Gf=gamma_s**(np.exp((-(f-fp)**2)/(2*sigma_s**2*fp**2)))
    PM=alpha_s*Hs**2*(fp**4/f**5)*np.exp(-beta_s*(fp/f)**4)
    Ei=Gf*PM

    return Ei


def PowerLaw(Hs, fp, f, df, n):
    intVal = (f**n * df).sum()
    Ei = (Hs**2 / ( 16 * intVal)) * f**n

    return Ei


def SpecVars(u=10, v=False):
    Tp = 2 * np.pi * u / (0.877 * g)
    fp = 1 / Tp
    Hs = 0.0246 * u**2
    k = (2 * np.pi * fp)**2 / g
    wl = 2 * np.pi / k

    if v:
        print(f'For winds of {u:.2f}m/s, expected waves are {Hs:.2f}m high,\n'
              f'with a period of {Tp:.2f}s, corresponding to a frequency of {fp:.2f}Hz,\n'
              f'and wavenumber of {k:.2f}/m or wavelength of {wl:.2f}m.')
    else:
        return(Hs, Tp, fp, k, wl)


def parseDType(**kwargs):
    MFac = 0
    EFac = 0
    d = 5000
    for key, value in kwargs.items():
        if key == 'DispType' and value == 'Open':
            MFac = 0
            EFac = 0
        elif key == 'DispType' and value == 'ML':
            MFac = 1
            EFac = 0
        elif key == 'DispType' and value == 'El':
            MFac = 0
            EFac = 1
        elif key == 'DispType' and (value == 'MLEl' or value == 'ElML'):
            MFac = 1
            EFac = 1
        elif key == 'd':
            d = value
        else:
            print('Unknown optional input')
            return
    return MFac, EFac, d


def omega(k, h=0, **kwargs):
    MFac, EFac, d = parseDType(**kwargs)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        A = g * k + EFac * (E * h**3 * k**5) / (12 * ( 1 - v**2) * rho_w)
    B = MFac * rho_i * h * k / rho_w + 1 / np.tanh(k * d)
    return (A / B) ** 0.5


def disp_rel(k, f, h=0, **kwargs):
    return 2 * np.pi * f - omega(k, h, **kwargs)


def disp_rel_exp(k, f, h, d, DispType):
    return 2 * np.pi * f - omega(k, h, d=d, DispType=DispType)


def calc_k(f, h=0, **kwargs):
    k0 = (2 * np.pi * f)**2 / g
    if h == 0:
        DispType = 'Open'
    else:
        DispType = 'ML'
    d = 5000
    for key, value in kwargs.items():
        if key == 'DispType':
            DispType = value
        elif key == 'd':
            d = value
        else:
            print('Unknown optional input')
            return

    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            k = optimize.newton(disp_rel_exp, k0, args=[f, h, d, DispType])
    except RuntimeError as e:
        if e.args[0].find('Failed') == -1:
            raise
        else:
            k = np.nan * f
    return k


def calc_cg(k, h=0, **kwargs):
    MFac, EFac, d = parseDType(**kwargs)
    A = g * k + EFac * (E * h**3 * k**5) / (12 * (1 - v**2) * rho_w)
    B = MFac * rho_i * h * k / rho_w + 1 / np.tanh(k * d)
    cg = 0.5 * A**-0.5 * B**-0.5 \
        * ( g + EFac * (5 * E * h**3 * k**4) / (12 * (1 - v**2) * rho_w)) \
        - 0.5 * A**0.5 * B**-1.5 * (MFac * rho_i * h / rho_w + d * (1 - 1 / np.tanh(k * d)**2))
    return cg
