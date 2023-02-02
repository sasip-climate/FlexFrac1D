#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:16:15 2022

@author: auclaije
"""


def prog_rep(i, n, s=10, flush=False):
    nc = (n - n % s)
    rem = i % (nc / s)
    if rem == 0 and i > 1:
        # print(i*s/nc, end='')
        print('#', end='', flush=flush)

    if i == n:
        print('', flush=flush)
