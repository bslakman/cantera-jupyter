
# coding: utf-8

# # Combustion and Air Pollution HW #6
# 
# In this homework we will use more detailed combustion chemistry then we previously considered in this class, along with the Zeldovich mechanism for thermal fixation of atmospheric N2, to model and predict NOx formation in heptane combustion.
# We use an open source software tool called [Cantera](http://www.cantera.org/docs/sphinx/html/index.html) which helps us solve thermodynamics and kinetics problems.
# 
# This homework was adapted from a lecture in Prof. Richard West's undergraduate chemical kinetics class in the Department of Chemical Engineering as well as existing Python notebooks in the [cantera-jupyter project](https://github.com/cantera/cantera-jupyter). 
# 
# #### Simply click inside the cells below to be able to type in and edit them, and press Shift+Enter to execute the code in a cell that is selected.

# In[1]:

import cantera as ct
import numpy as np
import scipy


# The n-heptane ($C_7H_{16}$) mechanism we will consider is from the Lawrence Livermore National Lab (LLNL), with the 3-reaction Zeldovich mechanism manually added to it. First we create a Solution object in Cantera for n-heptane gas that includes all of the thermodynamic and kinetic parameters needed to describe its combustion. (Don't worry about the warning)

# In[46]:

gas = ct.Solution('heptanesymp159.cti')


# We can print out some chemical species and reactions that are considered in the mechanism:

# In[53]:

print " ".join(gas.species_names[:20])
for i in range(20):
    print gas.reaction(i)
print "There are {0} reactions in this mechanism!".format(len(gas.reactions()))


# Let's look at the last few reactions in the mechanism, representing the Zeldovich mechanism.

# In[51]:

zeldovich = gas.reactions()[-3:]
for rxn in zeldovich:
    print rxn.equation, rxn.rate


# Just want to check at least one rate to make sure it's reasonably close to the one in the book...These came from the gri mech.

# In[102]:

T = np.linspace(300,2500,100)


# In[93]:

A = rxn.rate.pre_exponential_factor / 1000.0
n = rxn.rate.temperature_exponent
Ea = rxn.rate.activation_energy / 1000.0


# In[103]:

rate_mech = A * T**n * np.exp(-Ea/8.314/T)
rate_book = 7.1E7 * np.exp(-450.0/T)


# In[95]:

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[104]:

plt.plot(1.0/T, np.log(rate_mech))
plt.plot(1.0/T, np.log(rate_book))


# We'll specify a starting temperature and pressure of 1000 K and 10 bar, and stoichiometric combustion
# #### a. Specify the moles of O2 and N2 based on stoichiometric combustion.
# Replace the "xxxx" with the correct number of moles.

# In[118]:

gas.TPX = 1000, 10e5, 'nc7h16:1.0,o2:11.0,n2:41.58'
#gas.TPX = 1000, 10e5, 'nc7h16:1.0,o2:xxxx,n2:xxxx'


# We'll find the equilibrium flame temperature by equilibrating keeping pressure and enthalpy constant.

# In[117]:

gas.equilibrate('HP')
print "The adiabatic flame temperature is {0} K".format(int(round(gas.T)))
print "If you did part a properly, you should get a flame temperature of 2691 K."


# In[ ]:



