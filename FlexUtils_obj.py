#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:23:52 2022

@author: auclaije
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import config
from treeForFrac import getFractureHistory
from pars import E, v, rho_w, g, strainCrit
from IceDef import Floe


def calc_xstar(In, **kwargs):
    if type(In) == Floe:
        h = In.h
        Eval = Floe.E
        vval = Floe.v
    elif isinstance(In, (int, float)):
        h = In
        Eval = E
        vval = v

    for key, value in kwargs.items():
        if key == 'E':
            Eval = value
        elif key == 'v':
            vval = value
        else:
            print(f'Unknown input: {value}')

    xs = (np.pi / 4) * (Eval * h**3 / (36 * (1 - vval**2) * rho_w * g))**(1 / 4)
    return xs


def PlotFloes(x, t, Floes, wave, *args):
    hfac = 0
    nFloes = len(Floes)
    ShowE = hasattr(Floes[0], 'Eel')
    Eel = []

    if wave.type == 'WaveSpec':
        fig, hax = wave.plot(x)
        n0 = wave.Hs
        n = wave.Hs
        if len(wave.f) > 1:
            wvstring = f'Hs_{wave.Hs:4.2f}_wlp_{wave.wlp:2.0f}'
            wvtstring = f'$H_s$: {wave.Hs:4.2f}m'
        else:
            wvstring = f'n0_{wave.waves[0].n0:4.2f}_wl_{wave.waves[0].wl:2.0f}'
            wvtstring = f'$\eta_0$: {wave.waves[0].n0:4.2f}m'
    else:
        fig, hax = wave.plot(x, t, floes=Floes)
        n0 = wave.n0
        n = wave.amp(t)
        wvstring = f'n0_{wave.n0:04.2f}_wl_{wave.wl:2.0f}'
        wvtstring = f'$\eta_0$: {n:4.2f}m'

    if ShowE:
        Eel = (nFloes - 1) * Floes[0].k
    for floe in Floes:
        if wave.type == 'WaveSpec':
            wvf = wave.calc_waves(floe.xF)
        else:
            wvf = wave.waves(floe.xF, t, amp=floe.a0, phi=floe.phi0, floes=[floe])
        _, _ = floe.plot(x, t, wvf, fig, hax)
        if ShowE:
            Eel += floe.Eel
        if floe.hw > hfac:
            hfac = floe.hw

    Hfac = 1.5 * hfac + n0
    _ = hax.axis([x[0], x[-1], -Hfac, Hfac])

    if Eel == []:
        E_string = ''
    elif Eel == 0:
        E_string = '0'
    else:
        oom = np.floor(np.log10(Eel))
        E_string = ' - Energy: '
        E_string += str(round(Eel * 10**(2 - oom)) / 100)
        E_string += 'x10$^{' + str(int(oom)) + '}$ Jm$^{-2}$'

    tstr = f'Time: {t:.2f}s - ' + wvtstring + E_string

    hax.set_title(tstr)
    hax.set(ylabel='Height (m)', xlabel='Distance (m)')

    if nFloes > 2:
        minL = Floes[0].L
        maxL = Floes[0].L
        # Lmean = 0
        for floe in Floes:
            # Lmean += floe.L
            minL = min(minL, floe.L)
            maxL = max(maxL, floe.L)
        # Lmean = Lmean / nFloes
        # fstr = f'N: {nFloes} - L: {Lmean:.2f}'
        plt.text(x[-1] * 0.8, Hfac * 0.8, f'N: {nFloes}')
        plt.text(x[-1] * 0.8, Hfac * 0.7, f'min: {minL:.2f}m')
        plt.text(x[-1] * 0.8, Hfac * 0.6, f'max: {maxL:.0f}m')

    if len(args) > 0:
        itlab = ''
        Explab = 'Floes'
        for arg in args:
            if type(arg) is str:
                Explab = arg
            elif isinstance(arg, (int, np.int32, np.int64)):
                itlab = f'it_{arg:04}'
        root = (f'{Explab}_{wvstring}_'
                f'h_{Floes[0].h:03}_{itlab}')

        plt.savefig(config.FigsDirFloes + root + '.png', dpi=150)
        plt.close(fig)
    else:
        plt.show()


def PlotFracE(floe, Eel_floes):
    fig, hax = plt.subplots()
    x = floe.xF[1:-1] - floe.x0

    if type(Eel_floes) is list:
        nF = len(Eel_floes[1])
        Etemp = np.zeros((nF, 2))
        for iF in range(nF):
            Etemp[iF, :] = np.array(Eel_floes[1][iF])
        Eel_floes = Etemp

    # Left floe
    left, = hax.semilogy(x, Eel_floes[:, 0], 'b')
    # Right floe
    right, = hax.semilogy(x, Eel_floes[:, 1], 'r')
    # Combined with fracture energy
    Etot = Eel_floes[:, 0] + Eel_floes[:, 1] + floe.k
    imin = np.argmin(Etot)
    tot, = hax.semilogy(x, Etot, ':m', linewidth=3)
    # Initial floe energy
    init, = hax.semilogy(x[[1, -1]], floe.Eel * np.array([1, 1]))
    # Fracture location
    minE = min(min(Eel_floes[:, 0]), min(Eel_floes[:, 1]))
    maxE = max(max(Eel_floes[:, 0]), max(Eel_floes[:, 1]))
    hax.semilogy(x[imin] * np.ones(2), [minE, maxE], 'k', linewidth=1)

    hax.set(ylabel='Elastic Energy (J/m$^2$)', xlabel='Along floe distance (m)')

    # Strain if available
    if hasattr(floe, 'strain'):
        hax2 = hax.twinx()
        hax2.plot(x, floe.strain[1:-1], 'g')
        hax2.set_ylabel('Elastic strain (m$^{-1}$)', color='green')

    hax.legend([init, left, right, tot],
               ['Initial Floe', 'Left Floe', 'Right Floe', 'Fractured Total'],
               loc='best')


def PlotSum(t, y, **kwargs):
    DoSave = False
    DoLeg = False
    multiFrac = True
    for key, value in kwargs.items():
        if key == 'pstr':
            pstr = value
            DoSave = True
        elif key == 'leg':
            leg = value
            DoLeg = True
        elif key == 'multiFrac':
            multiFrac = value

    fmt = ['x-r', 'x--m', '^:b']
    fig, hax = plt.subplots()
    for iF in np.arange(y.shape[1]):
        if y.shape[1] > 1:
            if max(y[:, 1]) < 10 * max(y[:, 0]):
                hax.plot(t, y[:, iF], fmt[iF])
            else:
                hax.semilogy(t, y[:, iF], fmt[iF])
        else:
            hax.plot(t, y[:, iF], fmt[iF])

    hax.set(ylabel='Total energy (J/m$^2$)', xlabel='Time (s)')

    if DoLeg:
        hax.legend(leg)

    if DoSave:
        plt.savefig(config.FigsDirSumry + pstr + '.png', dpi=150)
        plt.close()

    return(fig, hax)


def BreakFloes(x, t, Floes, wave, multiFrac=True, *args):
    '''
        Searches if a fracture can occur and where would it be
    ----------
    x : np.array -> mesh of the scene
    t : float -> time of simulation
    Floes : list(Floe) -> current floes list (order from left to right)
    wave : Wave class -> wave
    multiFrac: boolean -> determines if more than one fracture should be considered
    *args :
        Etype: string -> energy type among 'Flex' and 'Disp'

    Returns
    -------
    Floes : list(Floe) -> updated list of floes, still in order from left to right
    '''

    EType = args[0] if len(args) > 0 else 'Flex'
    Spec = (wave.type == 'WaveSpec')

    Broke = True
    nFrac = 0

    while Broke:

        Broke = False

        # Computes Elastic Energy and break Floes
        NewFloes = Floes.copy()
        Offset = 0
        Etot = 0
        for iF in range(len(Floes)):
            if Spec:
                wvf = wave.calc_waves(Floes[iF].xF)
            else:
                wvf = wave.waves(Floes[iF].xF, t, amp=Floes[iF].a0,
                                 phi=Floes[iF].phi0, floes=[Floes[iF]])
            Eel1 = Floes[iF].calc_Eel(EType=EType, wvf=wvf)

            # Check if it is worth looking for fractures
            maxFrac = Floes[iF].Eel / Floes[iF].k
            if maxFrac > 1:

                xFracs, floes, Etot_floes = \
                    Floes[iF].FindE_min(wave, t, EType=EType, multiFrac=(multiFrac and maxFrac > 2))

                if Etot_floes < Eel1:
                    Broke = True

                    # # Show energy and strain if it would break
                    # xf, _, _, E_lists = Floes[iF].FindE_minVerbose(1, wave, t, EType=EType, V=True)
                    # PlotFracE(Floes[iF], E_lists)

                    # Add event to fracture history
                    getFractureHistory().addChildren(Floes[iF], floes, t)

                    # Modify list of floes
                    nFrac += len(xFracs)
                    NewFloes[iF + Offset] = floes[0]
                    for k in range(1, len(xFracs) + 1):
                        Offset += 1
                        NewFloes.insert(iF + Offset, floes[k])
                    Etot += Etot_floes
                else:
                    Etot += Eel1
            else:
                Etot += Eel1

        if Broke:
            Floes = NewFloes
        else:
            # for floe in Floes:
            #     if floe.Eel > 10 * floe.k:
            #         xf, _, _, E_lists = floe.FindE_minVerbose(1, wave, t, EType=EType, V=True)
            #         PlotFracE(floe, E_lists)
            break

    return Floes


def BreakFloesStrain(x, t, Floes, wave):
    ''' Use of a breaking parametrization based on strain (cf Dumont2011)
    Inputs/Outputs: same as BreakFloes
    '''

    Spec = True if wave.type == 'WaveSpec' else False

    # Note: in the code, the nergy is computed with calc_Eel since it also computed displacement
    Broke = True
    nFrac = 0

    while Broke:

        Broke = False
        NewFloes = Floes.copy()
        Offset = 0
        for iF in range(len(Floes)):

            # Compute displacement
            if Spec:
                wvf = wave.calc_waves(Floes[iF].xF)
            else:
                # Computes the wave amplitude and information along the floe
                a_vec = wave.amp_att(Floes[iF].xF, Floes[iF].a0, [Floes[iF]])
                phi0 = Floes[iF].phi0
                kw = Floes[iF].kw
                wvf = wave.waves(Floes[iF].xF, t, amp=Floes[iF].a0,
                                 phi=Floes[iF].phi0, floes=[Floes[iF]])
            Floes[iF].calc_w(wvf)

            # Computes the strain at top and bottom edges of the floe
            Floes[iF].calc_strain()
            strain = np.abs(Floes[iF].strain)

            if np.max(strain) > strainCrit:
                Broke = True

                # Since time discretization causes intervals where strain is too high
                # we dont want to break the floe at each point of the interval
                # Thus, we only look at local maxima -> scipy.signal.find_peaks
                iFracs = find_peaks(strain, height=strainCrit)[0]

                # Computes positions of fracturation and resulting floes
                xFracs = [Floes[iF].xF[i] for i in iFracs]
                createdFloes = Floes[iF].fracture(xFracs)
                nFrac += len(iFracs)

                # Add event to fracture history
                getFractureHistory().addChildren(Floes[iF], createdFloes, t)

                # Set properties induced by wave and insert floes in list
                distanceFromLeft = 0
                nFloes = len(iFracs) + 1
                for iNF in range(nFloes):
                    # Set properties and calculate waves
                    if Spec:
                        wvf = wave.calc_waves(createdFloes[iNF].xF)
                    else:
                        createdFloes[iNF].a0 = a_vec[0] if iNF == 0 else a_vec[iFracs[iNF - 1]]
                        createdFloes[iNF].phi0 = phi0 + kw * distanceFromLeft
                        wvf = wave.waves(createdFloes[iNF].xF, t, amp=createdFloes[iNF].a0,
                                         phi=createdFloes[iNF].phi0, floes=[createdFloes[iNF]])
                    createdFloes[iNF].calc_w(wvf)

                    distanceFromLeft += createdFloes[iNF].L
                    # Insert in list
                    if iNF == 0:
                        NewFloes[iF + Offset] = createdFloes[0]
                    else:
                        Offset += 1
                        NewFloes.insert(iF + Offset, createdFloes[iNF])

        if Broke:
            Floes = NewFloes

    return Floes


def PlotLengths(t, L, **kwargs):
    nt = len(t)
    addThickness = False
    tstring = ''
    addWaves = False
    wstring = ''
    xu = 'phase'
    trim = False

    for key, value in kwargs.items():
        if key == 'waves':
            addWaves = True
            waves = value
            wstring = (f'$\lambda$={waves.wl}m, $\eta_0$={waves.n0}m, '
                       f"{'constant' if {waves.beta == 0} else 'growing'} waves")
        elif key == 'Spec':
            addWaves = True
            waves = value
            if len(waves.f) > 1:
                wstring = (f'$\lambda$={waves.wlp:4.1f}m, $H_s$={waves.Hs:4.1f}m waves')
            else:
                wstring = f'$\lambda$={waves.wlp:4.1f}m, $\eta_0$: {waves.waves[0].n0:4.1f}m waves'
        elif key == 'x0':
            x0 = value
        elif key == 'h':
            addThickness = True
            hstring = f'{value}m ice'
        elif key == 'xunits':
            xu = value
        elif key == 'trim':
            trim = value

    if addWaves and addThickness:
        hstring += ' with '

    if addWaves or addThickness:
        tstring = 'Floe lengths for ' + hstring + wstring

    if nt == 1:
        fig, hax = PlotLengths1(L, tstring)
        return(fig, hax)

    fig, hax = plt.subplots()
    frmt = ['x-c', 'x-b', 'x-m', 'x-r', 'x-y']
    for it in range(nt):
        tvec = [t[it], t[it]]
        hax.plot(tvec, [0, L[it][0]], frmt[0], linewidth=3)
        L0 = L[it][0]
        if trim:
            EndFloe = len(L[it][:-1])
        else:
            EndFloe = len(L[it])
        for iL in range(1, EndFloe):
            hax.plot(tvec, L0 + np.array([0, L[it][iL]]), frmt[iL % 5])
            L0 = L0 + L[it][iL]

    if addWaves:
        if waves.type == 'Wave':
            hax.plot(t, L0 * 1.2 + L0 * 0.1 * waves.waves(x0, t * waves.T / (2 * np.pi), amp=1))

    if addThickness or addWaves:
        hax.set_title(tstring)

    hax.set(ylabel='Floe length (m)')
    if xu == 'phase':
        hax.set(xlabel='Initial phase (rad)')
        plt.xticks(np.arange(0, 2.5, 0.5) * np.pi, ['0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'])
    elif xu == 'trials':
        hax.set(xlabel='Trial number')

    return(fig, hax)


def PlotLengths1(L, *args):
    fig, hax = plt.subplots()
    nFloes = len(L[0])
    for iF in range(nFloes):
        hax.plot([0, L[0][iF]], [iF, iF])
    hax.set(xlabel='Floe Length (m)', ylabel='Floe number')
    return(fig, hax)
