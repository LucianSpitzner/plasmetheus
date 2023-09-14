# coding=utf-8
"""
File which stores natural constants in cgs units.
Taken from PROMETHEUS, created on 2. June 2021 by Andrea Gebek.
"""

import numpy as np

e = 4.803e-10   # Elementary charge
m_e = 9.109e-28
c = 2.998e10
G = 6.674*10**(-8)
k_B = 1.381*10**(-16)
amu = 1.661*10**(-24)
R_J = 6.99e9        #Jupiter radius
M_J = 1.898e30      #Jupiter mass
M_E = 5.974e27      #Earth mass
R_sun = 6.96e10     # Solar radius in cm
M_sun = 1.988e33
R_Io = 1.822e8      # Io radius
euler_mascheroni = 0.57721 
AU = 1.496e13   # Conversion of one astronomical unit into cm

"""
Masses for the absorbing atoms/ions (this is not a really elegant way,
but I couldn't find another solution yet)
"""

speciesInfoDict = {'NaI': ['Na', '1', 22.99 * amu], 
                   'KI': ['K', '1', 39.0983 * amu], 
                   'SiI': ['Si', '1', 28.0855 * amu], 
                   'SiII': ['Si', '2', 28.0855 * amu],
                   'SiIII': ['Si', '3', 28.0855 * amu], 
                   'SiIV': ['Si', '4', 28.0855 * amu], 
                   'MgI': ['Mg', '1', 24.305 * amu], 
                   'MgII': ['Mg', '2', 24.305 * amu],
                   'OII': ['O', '2', 15.999 * amu], 
                   'CII': ['C', '2', 12.000 * amu]}
