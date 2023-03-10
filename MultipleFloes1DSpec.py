#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:48:40 2022

@author: auclaije
"""

import numpy as np
import matplotlib.pyplot as plt
import config
import pars
from os import path
from tqdm import tqdm
from FlexUtils_obj import PlotFloes, PlotLengths, calc_xstar
from FlexUtils_obj import BreakFloes, BreakFloesStrain
from FSDUtils import PlotFSD
from WaveUtils import calc_k
from WaveSpecDef import WaveSpec
from WaveChecks import plotDisp, plot_cg
from IceDef import Floe
from treeForFrac import getFractureHistory, InitHistory

# 0: None, 1: Lengths, 2: Lengths and FSD, 3: Lengths, FSD and saved Floes, 4: Lengths, FSD and Floes
DoPlots = 3
repeats = 10
multiFrac = True
FractureCriterion = 'Energy'

# Ice parameters
h = 1
x0 = 50
L0 = 300
L = L0
dx = 0.5
DispType = 'ML'
EType = 'Flex'
# Initialize ice floe object
floe1 = Floe(h, x0, L, DispType=DispType, dx=dx)

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
x = np.arange(x0 + L + 2)
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
            Spec.plotEx(fname=(config.FigsDirSpec + f'Spec_{DispType}_L0_{L0:04}_{t:04.0f}.png'), t=t)
        # Spec.set_phases(x, t, Floes)
        # Spec.plotWMean(x, floes=[floe1], fname='Spec/Waves_{t:04.0f}.png')

FL = [0] * repeats

dt = dx / (Spec.fp * Spec.wlp)
t = np.arange(0, 2 * tSpecM + 2 / Spec.f[0], Spec.Tp / 20)  # min(Spec.Tp / 20, dt))

phi = 2 * np.pi * np.linspace(0, 1, num=repeats, endpoint=False)

print(f'Launching {repeats} experiments:')
print(f'Note: x* = {calc_xstar(floe1)}, wlm = {2*np.pi/Spec.k[-1]/2}, '
      f'wlp = {2*np.pi/Spec.kp/2}, wlM = {2*np.pi/Spec.k[0]/2}')
for iL in range(repeats):
    lab = f'Exp_{iL:02}_E_{EType}_{Spec.SpecType}_F_{FractureCriterion}_L0_{L0}'
    LoopName = f'{lab}_h_{h:3.1f}m_Hs_{Spec.Hs:04.1f}m.txt'
    DataPath = config.DataTempDir + LoopName
    FracHistPath = DataPath[:-4] + '_History.txt'
    if path.isfile(DataPath):
        print(f'Reading existing data for loop {iL:02}')
        FL[iL] = list(np.loadtxt(DataPath))
        if DoPlots > 2:
            DoPlots = 2
        continue

    # Change the phases of each wave
    if len(Spec.f) == 1:
        Spec.phi = np.array([phi[iL]])
    Spec.setWaves()
    # Reset the initial floe, history and domain
    if floe1.L > L0:
        L = L0
        floe1 = floe1.fracture(x0 + L0)[0]
        x = np.arange(x0 + L0 + 2)
    Floes = [floe1]
    InitHistory(floe1, t[0])

    tqdmlab = f'Time Loop {iL:02}' if repeats > 1 else 'Time Loop'
    for it in tqdm(range(len(t)), desc=tqdmlab):

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

        if DoPlots > 3:
            PlotFloes(x, t[it], Floes, Spec)
        elif DoPlots > 2 or len(Floes) > nF:
            PlotFloes(x, t[it], Floes, Spec, lab, it)

        if Floes[-1].x0 > 0.6 * L + x0:
            # Update floes in history (NOTE: It also update the actual floe's length)
            getFractureHistory().modifyLengthDomain(L / 2)
            # Update floe resolution and matrix
            nx = max(int(np.ceil(Floes[-1].L)), 100)
            Floes[-1].xF = Floes[-1].x0 + np.linspace(0, Floes[-1].L, nx)
            Floes[-1].initMatrix()
            # And update domain
            x = np.arange(0, x[-1] + L / 2 + 1, 1)
            L *= 1.5
            print('+', end='')

    FL_temp = []
    for floe in Floes:
        FL_temp.append(floe.L)
    np.savetxt(DataPath, FL_temp)
    FL[iL] = FL_temp

    if DoPlots > 2:
        DoPlots = 2
    FractureHistory = getFractureHistory()
    if DoPlots > 0:
        fGen = config.FigsDirSumry + 'Gen' + LoopName[3:-4] + '.png'
        FractureHistory.plotGeneration(filename=fGen)
        if DoPlots > 3:
            FractureHistory.plotGeneration()
        # Save the fracture history, giving (L, x0, gen, time, boolean existing, parent) informations
        # for each floe that has existed
        np.savetxt(FracHistPath, FractureHistory.asArray())

n0 = Spec.calcHs()

if DoPlots > 0:
    if len(Spec.f) == 1:
        xv = phi
        xu = 'phase'
    else:
        xv = np.arange(repeats)
        xu = 'trials'

    root = (f'FloeLengths_Spec_E_{EType}_{Spec.SpecType}_F_{FractureCriterion}_'
            f'{DispType}_Hs_{Spec.Hs:05.2f}_wlp_{Spec.wlp:06.2f}_h_{h:3.1f}_L0_{L:04}')

    fig, hax = PlotLengths(xv, FL, x0=x0, h=h, Spec=Spec, xunits=xu)
    plt.savefig(config.FigsDirSumry + root + '.png', dpi=150)

    fig, hax = PlotLengths(xv, FL, x0=x0, h=h, Spec=Spec, xunits=xu, trim=True)
    plt.savefig(config.FigsDirSumry + root + '_trim.png', dpi=150)

if DoPlots > 1:
    fn = (f'{config.FigsDirSumry}/Spec_E_{EType}_{Spec.SpecType}_F_{FractureCriterion}_'
          f'{DispType}_Hs_{Spec.Hs:05.2f}_wlp_{Spec.wlp:06.2f}_h_{h:3.1f}_L0_{L:04}')

    wvl_lab = '$\lambda_p/2$' if len(Spec.f) > 1 else '$\lambda/2$'
    Lines = [[Spec.wlp / 2, wvl_lab]]
    wvl_max = 2 * np.pi / calc_k(Spec.f[0], floe1.h, DispType=floe1.DispType)
    if wvl_max > Spec.wlp * 1.5:
        Lines.append([wvl_max / 2, '$\lambda_{max}/2$'])
    wvl_min = 2 * np.pi / calc_k(Spec.f[-1], floe1.h, DispType=floe1.DispType)
    if wvl_min < Spec.wlp / 1.5:
        Lines.append([wvl_min / 2, '$\lambda_{min}/2$'])
    Lines.append([calc_xstar(floe1), '$x^*$'])

    PlotFSD(FL, FileName=fn, Lines=Lines, Spec=Spec)
