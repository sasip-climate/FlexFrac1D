#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:20:57 2021

@author: auclaije
"""

# Given Young's modulus E and Poisson's ratio v (nu)
# Calculate Lamé's parameters lambda (l) and mu (u)


def Lame(E, v):
    lmbd = E*v/((1+v)*(1-2*v))
    u = E/(2*(1+v))
    return(lmbd, u)


def FracToughness(E, v, K):
    return((1-v**2)*K**2/E)
