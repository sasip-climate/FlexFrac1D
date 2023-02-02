#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:50:26 2022

@author: auclaije
"""

import numpy as np
import matplotlib.pyplot as plt
from WaveUtils import calc_k, calc_cg


def plotDisp(f, h, **kwargs):

    ko = np.empty_like(f)
    kel = np.empty_like(f)
    kml = np.empty_like(f)
    k2 = np.empty_like(f)

    ko = calc_k(f, 0, **kwargs)
    kel = calc_k(f, h, DispType='El', **kwargs)
    kml = calc_k(f, h, DispType='ML', **kwargs)
    k2 = calc_k(f, h, DispType='ElML', **kwargs)

    ind = kml < 10
    fml = f[ind]
    kml = kml[ind]

    fig, hax = plt.subplots()
    hax.plot(ko, f, kel, f, kml, fml, '*', k2, f, '+')
    hax.legend(['Open Water', 'Elastic', 'Mass-Loading', 'Both'])
    hax.set(xlabel='Wavenumber (1/m)', ylabel='Frequency (Hz)')
    hax.set_xlim([0, 1])

    fig, hax = plt.subplots()
    hax.semilogx(ko, f, kel, f, kml, fml, '*', k2, f, '+')
    hax.legend(['Open Water', 'Elastic', 'Mass-Loading', 'Both'])
    hax.set(xlabel='Wavenumber (1/m)', ylabel='Frequency (Hz)')


def plot_cg(f, h):

    ko = np.empty_like(f)
    kel = np.empty_like(f)
    kml = np.empty_like(f)
    k2 = np.empty_like(f)

    ko = calc_k(f)
    kel = calc_k(f, h, DispType='El')
    kml = calc_k(f, h, DispType='ML')
    k2 = calc_k(f, h, DispType='ElML')

    cgo = calc_cg(ko)
    cgel = calc_cg(kel, h, DispType='El')
    cgml = calc_cg(kml, h, DispType='ML')
    cg2 = calc_cg(k2, h, DispType='ElML')

    ind = kml < 10
    fml = f[ind]
    cgml = cgml[ind]

    fig, hax = plt.subplots()
    hax.semilogx(cgo, f, cgel, f, cgml, fml, '*', cg2, f, '+')
    hax.legend(['Open Water', 'Elastic', 'Mass-Loading', 'Both'])
    hax.set(xlabel='Group velocity (m/s)', ylabel='Frequency (Hz)')
