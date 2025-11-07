#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 10:14:16 2025

@author: ggslc
"""

import numpy as np


def cell_id(i, j, ni, nj):
    return j + nj * i


missing = np.nan