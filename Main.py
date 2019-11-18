# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 17:32:23 2019

@author: Adam
"""

import MethaneEOS
import thermoClass as thermo
import numpy as np
import matplotlib.pyplot as plt
mThermo = thermo.thermo()

print(mThermo.pressure(425.61,110))

mdot = 1560 #g/s
