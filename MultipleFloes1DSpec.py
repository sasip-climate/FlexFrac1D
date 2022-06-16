#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:48:40 2022

@author: auclaije
"""

import numpy as np
import matplotlib.pyplot as plt
import config
from tqdm import tqdm
from FlexUtils_obj import PlotFloes, BreakFloes, BreakFloesStrain, PlotLengths, PlotFSD
from WaveUtils import calc_k
from WaveSpecDef import WaveSpec
from WaveChecks import plotDisp, plot_cg
from IceDef import Floe

# 0: None, 1: Lengths, 2: Lengths and FSD, 3: Lengths, FSD and saved Floes, 4: Lengths, FSD and Floes
DoPlots = 3
repeats = 10
multiFrac = 2
FractureCriterion = 'Energy'

# Ice parameters
h = 1
x0 = 50
L = 100
DispType = 'ML'
EType = 'Flex'
# Initialize ice floe object
floe1 = Floe(h, x0, L, DispType=DispType)

# Wave Parameters
u = 5  # Wind speed (m/s)
# Initialize wave object
Spec = WaveSpec(u=u)

# calculate wave properties in ice
Spec.checkSpec(floe1)
ki = floe1.kw
floe1.setWPars(Spec)

xi = 1 / floe1.alpha
if L > 5 * xi[Spec.f == Spec.fp]:
    print('Warning: Floe is more than 5x the attenuation length of the peak')

plotDisp(Spec.f, h)
plot_cg(Spec.f, h)

# Initial setup
x = np.arange(2 * x0 + L + 1)
xProp = 4 / floe1.alpha
xProp[xProp > L] = L
tProp = (2 * x0 / Spec.cgw + xProp / floe1.cg)
tPropMax = max(tProp)
tSpecM = max(tProp[Spec.Ei > 0.1 * max(Spec.Ei)])

Floes = [floe1]
# Visualize energy propagation in the domain
if DoPlots > 0:
    for t in np.arange(Spec.Tp, 2 * tSpecM + 1 / Spec.f[0], tSpecM / 10):
        Spec.calcExt(x, t, Floes)
        if t < tSpecM * 1.1:
            Spec.plotEx(fname=(config.FigsDirSpec + f'Spec_{DispType}_L0_{L:04}_{t:04.0f}.png'), t=t)
        # Spec.set_phases(x, t, Floes)
        # Spec.plotWMean(x, floes=[floe1], fname='Spec/Waves_{t:04.0f}.png')

FL = [0] * repeats
t = np.arange(0, 1.2 * tSpecM + 2 / Spec.f[0], Spec.Tp / 20)

print(f'Launching {repeats} experiments:')
for iL in range(repeats):
    # Change the phases of each wave
    Spec.setWaves()

    Spec.calcExt(x, t[0], [floe1])
    # Spec.plotEx()
    Spec.set_phases(x, t[0], [floe1])

    wvf = Spec.calc_waves(floe1.xF)
    floe1.calc_Eel(wvf=wvf, EType=EType)
    Floes = [floe1]

    for it in tqdm(range(len(t))):

        nF = len(Floes)

        Spec.calcExt(x, t[it], Floes)
        # Spec.plotEx(t=t[it])
        Spec.set_phases(x, t[it], Floes)
        # Spec.plotWMean(x, floes=Floes)
        if FractureCriterion == 'Energy':
            Floes = BreakFloes(x, t[it], Floes, Spec, multiFrac, EType)
        elif FractureCriterion == 'Strain':
            Floes = BreakFloesStrain(x, t[it], Floes, Spec)
        else:
            raise ValueError('Non existing fracturation criterion')

        if DoPlots > 2 or len(Floes) > nF:
            lab = f'Exp_{iL:02}_E_{EType}_F_{FractureCriterion}'
            PlotFloes(x, t[it], Floes, Spec, lab, it)
        elif DoPlots > 3:
            PlotFloes(x, t[it], Floes, Spec)

    FL_temp = []
    for floe in Floes:
        FL_temp.append(floe.L)
    FL[iL] = FL_temp

n0 = Spec.calcHs()

if DoPlots > 0:
    fig, hax = PlotLengths(np.arange(repeats), FL, x0=x0, h=h, xunits='trials')

    root = (f'FloeLengths_Spec_E_{EType}_F_{FractureCriterion}_'
            f'{DispType}_Hs_{Spec.Hs}_wlp_{Spec.wlp:06.2f}_h_{h:3.1f}_L0_{L:04}')

    plt.savefig(config.FigsDirSumry + root + '.png')

if DoPlots > 1:
    fn = (f'_Spec_E_{EType}_F_{FractureCriterion}_'
          f'{DispType}_Hs_{Spec.Hs:05.2f}_wlp_{Spec.wlp:06.2f}_h_{h:3.1f}_L0_{L:04}')

    Lines = [[Spec.wlp / 2, '$\lambda_p/2$']]
    wvl_max = 2 * np.pi / calc_k(Spec.f[0], floe1.h, DispType=floe1.DispType)
    Lines.append([wvl_max / 2, '$\lambda_{max}/2$'])
    wvl_min = 2 * np.pi / calc_k(Spec.f[-1], floe1.h, DispType=floe1.DispType)
    Lines.append([wvl_min / 2, '$\lambda_{min}/2$'])

    PlotFSD(FL, h=h, n0=n0, FileName=fn, Lines=Lines)
