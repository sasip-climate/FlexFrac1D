#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 17:52:38 2022

@author: Jean-Pierre Auclair
"""

import numpy as np
import itertools as it
import matplotlib.pyplot as plt
from os import listdir
from scipy import stats


def PlotFSD(L, **kwargs):
    # Process input
    if isinstance(L[0], (list, np.ndarray)):
        Ll = []
        for l in L:
            Ll += l[:-1]  # Do not consider the last floe in the FSD
    else:
        Ll = L

    Ll = np.array(Ll)

    # To prevent errors in case of no fracture
    empty = False
    if Ll.size == 0:
        Ll = np.zeros((1,))
        empty = True

    # Process optional inputs
    DoSave = False
    fn = ''
    wle = False
    Lines = []

    L_min = np.floor(min(Ll))
    L_max = np.ceil(max(Ll))
    # dL = 1 if L_min < 10 else 2
    dL = np.ceil(10 * (L_max - L_min) / len(Ll)**0.5) / 10
    if dL < 1:
        dL = round(dL * 10) / 10
    else:
        dL = round(dL)
    Lmin = L_min
    Lmax = L_max
    XLims = [Lmin, Lmax]
    SetLims = False

    # Process optional inputs
    for key, value in kwargs.items():
        if key == 'FileName':
            DoSave = True
            fn = value
        elif key == 'wl':
            wl = value
            wle = True
        elif key == 'Lmin' or key == 'L_min':
            Lmin = value
        elif key == 'Lmax' or key == 'L_max':
            Lmax = value
        elif key == 'dL':
            dL = value
        elif key == 'Lines':
            Lines = value
        elif key == 'XLims':
            SetLims = True
            XLims = value

    values, edges = np.histogram(Ll, bins=np.arange(Lmin, Lmax, dL))
    tstring = (f'FSD for {len(Ll)} floes of {Lmin} to {Lmax}m in '
               f'{(Lmax-Lmin)/dL:.0f} bins of {dL}m.')
    print(tstring)

    fac = [1]
    fac.append(1 / values.sum())
    fac.append((edges[:-1] + 0.5) / values.sum())

    ylab = ['Number', 'Frequency', 'Length-fraction']

    if wle:
        Lines.append([wl / 2, '$\lambda$/2'])

    for ifac in np.arange(len(fac)):
        fig, hax = PlotHist(edges, values * fac[ifac])
        hax.set(ylabel=ylab[ifac])
        # If no fracture, put it in the title
        if empty:
            hax.set(title="No fracture for those parameters")
        # else:
        #     hax.set(title=tstring)

        if len(Lines):
            addLines(hax, Lines)

        if SetLims:
            plt.xlim(XLims)

        if DoSave:
            plt.savefig(f'{fn}_FSD_{ylab[ifac]}_{Lmin}_{Lmax}_{dL}.png', dpi=150)
            plt.close()
        else:
            plt.show()

    return(edges, values)


def PlotHist(edges, values):

    fig, hax = plt.subplots()
    plt.bar(edges[:-1], values, align='edge', width=edges[1] - edges[0])

    hax.set(xlabel='Floe length (m)')

    return(fig, hax)


def PlotFSDs(L, **kwargs):
    # Process input
    if isinstance(L[0], (list, np.ndarray)):
        Ll = []
        for l in L:
            Ll += l[:-1]  # Do not consider the last floe in the FSD
    else:
        Ll = L

    Ll = np.array(Ll)

    # To prevent errors in case of no fracture
    empty = False
    if Ll.size == 0:
        Ll = np.zeros((1,))
        empty = True

    # Process optional inputs
    DoSave = False
    fn = ''
    wle = False
    Lines = []

    L_min = np.floor(min(Ll))
    L_max = np.ceil(max(Ll))
    # dL = 1 if L_min < 10 else 2
    dL = np.ceil(10 * (L_max - L_min) / len(Ll)**0.5) / 10
    if dL < 1:
        dL = round(dL * 10) / 10
    else:
        dL = round(dL)
    Lmin = L_min
    Lmax = L_max
    XLims = [Lmin, Lmax]
    SetLims = False

    # Process optional inputs
    for key, value in kwargs.items():
        if key == 'FileName':
            DoSave = True
            fn = value
        elif key == 'wl':
            wl = value
            wle = True
        elif key == 'Lmin' or key == 'L_min':
            Lmin = value
        elif key == 'Lmax' or key == 'L_max':
            Lmax = value
        elif key == 'dL':
            dL = value
        elif key == 'Lines':
            Lines = value
        elif key == 'XLims':
            SetLims = True
            XLims = value

    values, edges = np.histogram(Ll, bins=np.arange(Lmin, Lmax, dL))
    tstring = (f'FSD for {len(Ll)} floes of {Lmin} to {Lmax}m in '
               f'{(Lmax-Lmin)/dL:.0f} bins of {dL}m.')
    print(tstring)

    fac = [1]
    fac.append(1 / values.sum())
    fac.append((edges[:-1] + 0.5) / values.sum())

    ylab = ['Number', 'Frequency', 'Length-fraction']

    if wle:
        Lines.append([wl / 2, '$\lambda$/2'])

    for ifac in np.arange(len(fac)):
        fig, hax = PlotHist(edges, values * fac[ifac])
        hax.set(ylabel=ylab[ifac])
        # If no fracture, put it in the title
        if empty:
            hax.set(title="No fracture for those parameters")
        # else:
        #     hax.set(title=tstring)

        if len(Lines):
            addLines(hax, Lines)

        if SetLims:
            plt.xlim(XLims)

        if DoSave:
            plt.savefig(f'{fn}_FSD_{ylab[ifac]}_{Lmin}_{Lmax}_{dL}.png', dpi=150)
            plt.close()
        else:
            plt.show()

    return(edges, values)


def addLines(hax, Lines):

    ylims = hax.get_ylim()
    xlims = hax.get_xlim()
    xoffset = 0.02 * (xlims[1] - xlims[0])

    colors = ['green', 'blue', 'purple', 'magenta', 'red', 'orange']
    styles = ['-', '--', '-.', ':', '-.', '--']
    yoffset = np.arange(1, 1 / len(Lines) - 1e-12, -1 / len(Lines)) * 0.9

    for iL in np.arange(len(Lines)):
        hax.plot(Lines[iL][0] * np.ones(2), ylims, color=colors[iL], linestyle=styles[iL])
        hax.text(Lines[iL][0] + xoffset, ylims[1] * yoffset[iL], Lines[iL][1], fontsize=20)


def FSD_read(dirname):
    files = listdir(f'{dirname}/database/temp')
    FS = [file for file in files if "History" not in file]

    FL = []
    for file in FS:
        temp = np.loadtxt(f'{dirname}/database/temp/{file}')
        if temp.ndim != 0:
            FL.append(list(temp))

    return(FL)


def FSD_redraw(FL, **kwargs):
    FileName = None
    Lines = []
    Lmin = None
    Lmax = None
    dL = None

    for key, value in kwargs.items():
        if key == "FileName":
            FileName = value
        elif key == "Lines":
            Lines = value
        elif key == 'Lmin':
            Lmin = value
        elif key == 'Lmax':
            Lmax = value
        elif key == 'dL':
            dL = value

    Ll = []
    for l in FL:
        Ll += l[:-1]  # Do not consider the last floe in the FSD

    nbins = np.floor((len(Ll))**0.5)
    if Lmin is None:
        Lmin = min(Ll)
    if Lmax is None:
        Lmax = max(Ll)
    if dL is None:
        dL = np.ceil(10 * (Lmax - Lmin) / nbins) / 10
    print(f'FSD for sizes {Lmin} to {Lmax}, with {len(Ll)} floes and {nbins:.0f} bins of {dL}m.')

    if FileName is None:
        PlotFSD(FL, Lines=Lines, dL=dL, Lmin=Lmin, Lmax=Lmax)
    else:
        PlotFSD(FL, FileName=FileName, Lines=Lines, dL=dL, Lmin=Lmin, Lmax=Lmax)


def FSD_BatchRead(Root, Names=False):
    folders = listdir(Root)
    FR = [folder for folder in folders if "_mF_" in folder]

    if len(FR) == 0:
        print(f'WARNING: No experiment folder in {Root}')
    FLB = [FSD_read(f'{Root}/{folder}') for folder in FR]

    if Names:
        return(FLB, FR)
    else:
        return(FLB)


def FSD_Sample(FLB, FR, keys):
    if type(keys) is str:
        keys = [keys]

    FLBS = []
    FRS = []
    for ind in np.arange(len(FR)):
        match = True
        for key in keys:
            if key not in FR[ind]:
                match = False
        if match:
            FLBS.append(FLB[ind])
            FRS.append(FR[ind])
    return(FLBS, FRS)


def FSD_Clean(FL):
    for ind in np.arange(len(FL)):
        if FL[ind][-1] < max(FL[ind]):
            print('Warning: Array is already clean')
        else:
            FL[ind] = FL[ind][:-1]
    return(FL)


def FSD_Merge(FL):
    FLM = FL[0].copy()
    for ind in np.arange(1, len(FL)):
        for floes in FL[ind]:
            FLM.append(floes)
    return(FLM)


def FSD_BatchReplot(Root, **kwargs):
    SaveDir = ''
    KeyWords = {}

    for key, value in kwargs.items():
        if key == "SaveDir":
            SaveDir = value
        elif key == "Lines":
            KeyWords["Lines"] = value
        elif key == 'Lmin':
            KeyWords["Lmin"] = value
        elif key == 'Lmax':
            KeyWords["Lmax"] = value
        elif key == 'dL':
            KeyWords["dL"] = value

    [FLB, FR] = FSD_BatchRead(Root, True)

    for ind in np.arange(len(FR)):
        KeyWords["FileName"] = f'{SaveDir}/{FR[ind]}'
        PlotFSD(FLB[ind], **KeyWords)


def PlotMFSD(FLL, **kwargs):
    # Process input
    NF = np.inf
    Labels = []
    for i in np.arange(len(FLL)):
        if isinstance(FLL[i][1], (list, np.ndarray)):
            FL = FSD_Clean(FLL[i])
            FLL[i] = np.array(FSD_Merge(FL))
        NF = min([len(FLL[i]), NF])
        Labels.append(f'FSD-{i}')

    # Process optional inputs
    DoSave = False
    fn = ''
    wle = False
    Lines = []
    colors = ['blue', 'green', 'orange', 'red', 'purple']
    alpha = 0.6

    L_min = np.floor(min(it.chain.from_iterable(FLL)))
    L_max = np.ceil(max(it.chain.from_iterable(FLL)))
    # dL = 1 if L_min < 10 else 2
    dL = np.ceil(10 * (L_max - L_min) / NF**0.5) / 10
    if dL < 1:
        dL = round(dL * 10) / 10
    else:
        dL = round(dL)
    Lmin = L_min
    Lmax = L_max
    SetXLims = False
    SetYLims = False

    # Process optional inputs
    for key, value in kwargs.items():
        if key == 'FileName':
            DoSave = True
            fn = value
        elif key == 'wl':
            wl = value
            wle = True
        elif key == 'Lmin' or key == 'L_min':
            Lmin = value
        elif key == 'Lmax' or key == 'L_max':
            Lmax = value
        elif key == 'dL':
            dL = value
        elif key == 'Lines':
            Lines = value
        elif key == 'XLims':
            SetXLims = True
            XLims = value
        elif key == 'YLims':
            SetYLims = True
            YLims = value
        elif key == 'Labels':
            Labels = value
        elif key == 'Colors':
            colors = value
        elif key =='alpha':
            alpha = value

    values = []
    edges = []
    fac = []
    for i in np.arange(len(FLL)):
        vt, et = np.histogram(FLL[i], bins=np.arange(Lmin, Lmax, dL))
        values.append(vt)
        edges.append(et)

        fact = [1]
        fact.append(1 / vt.sum())
        fact.append((et[:-1] + 0.5) / vt.sum())
        fac.append(fact)

    ylab = ['Number', 'Frequency', 'Length-fraction']

    if wle:
        Lines.append([wl / 2, '$\lambda$/2'])

    for ifac in np.arange(len(fac[0])):
        fig, hax = plt.subplots()
        for iD in np.arange(len(FLL)):
            # fig, hax = PlotHist(edges, values * fac[ifac])
            plt.bar(edges[iD][:-1], values[iD] * fac[iD][ifac], align='edge',
                    width=dL, alpha=alpha, label=Labels[iD], color=colors[iD])

        hax.set(ylabel=ylab[ifac])
        hax.set(xlabel='Floe length (m)')
        plt.legend(loc='upper right')

        if len(Lines):
            addLines(hax, Lines)

        if SetXLims:
            plt.xlim(XLims)

        if SetYLims:
            plt.ylim(YLims[ifac])

        if DoSave:
            tags = Labels[0]
            for Lab in Labels[1:]:
                tags = tags + '_' + Lab
            plt.savefig(f'{fn}_FSDS_{ylab[ifac]}_{tags}.png', dpi=150)
            plt.close()
        else:
            plt.show()

    return(edges, values)


def calc_Lx(xE, FLL):
    # Remove leading 0 if present
    if xE[0] == 0:
        xE = xE[1:]

    FNx = [[] for i in range(len(xE) + 1)]
    for FL in FLL:
        x0 = 0
        nx = 0
        for Floe in FL:
            if x0 > xE[-1]:
                FNx[-1].append(Floe)
            else:
                if x0 > xE[nx]: nx = min(nx + 1, len(xE))
                FNx[nx].append(Floe)
            x0 = x0 + Floe

    return(FNx)


def FSD_NormPlot(FLL, **kwargs):
    FileName = ''
    DType = 'normal'

    for key, value in kwargs.items():
        if key == 'FileName':
            FileName = value
        elif key == 'DType':
            DType = value

    # Calculate quantiles and least-square-fit curve
    FL = FSD_Merge(FLL)
    if DType == 'log': FL = np.log(FL)
    (quantiles, values), (slope, intercept, r) = stats.probplot(FL, dist='norm')

    # plot results
    plt.plot(values, quantiles, 'ob')
    plt.plot(quantiles * slope + intercept, quantiles, 'r')

    # define y-ticks
    ticks_perc = [1, 5, 10, 20, 50, 80, 90, 95, 99]

    # transfrom them from precentile to cumulative density
    ticks_quan = [stats.norm.ppf(i / 100.) for i in ticks_perc]

    # assign new y-ticks
    plt.yticks(ticks_quan, ticks_perc)

    # if in log, use different xticks
    if DType == 'log':
        xticks = [1, 2, 4, 8, 16, 32]
        plt.xticks(np.log(xticks), xticks)

    plt.grid()

    plt.xlabel('Floe Size (m)')
    plt.ylabel('Percentile (%)')

    if FileName == '':
        # show plot
        plt.grid()
        plt.show()
    else:
        plt.savefig(f'{FileName}_NormPlot.png', dpi=150)
        plt.close()
