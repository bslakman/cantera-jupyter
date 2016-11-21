
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


# We'll specify a starting temperature and pressure of 1000 K and 3 bar, and stoichiometric combustion
# #### a. Specify the moles of O2 and N2 based on stoichiometric combustion.
# Replace the "xxxx" with the correct number of moles.

# In[219]:

gas.TPX = 1000, 3e5, 'nc7h16:1.0,o2:11.0,n2:41.58'
#gas.TPX = 1000, 3e5, 'nc7h16:1.0,o2:xxxx,n2:xxxx'


# We'll find the equilibrium flame temperature by equilibrating keeping pressure and enthalpy constant.

# In[214]:

gas.equilibrate('HP')
print "The adiabatic flame temperature is {0} K".format(int(round(gas.T)))
print "If you did part a properly, you should get a flame temperature of 2633 K."


# Let's print out the concentration of some components at equilibrium.

# In[215]:

print "Heptane: {:.2e}".format(gas['nc7h16'].X[0])
print "OH: {:.2e}".format(gas['oh'].X[0])
print "NO: {:.2e}".format(gas['no'].X[0])
print "O2: {:.2e}".format(gas['o2'].X[0])


# #### b. What is the NO concentration if the fuel is burned at an initial temperature of 600, 800, and 1200 K?
# Replace the temperature above and rerun that and the following cells.

# ### Kinetics
# We will burn the gas starting at the same temperature, but atmospheric pressure. This time we will hold volume constant using an ideal gas reactor, so pressure can change.
# The following code will print out the simulation time along with temperature, pressure, and internal energy.

# In[225]:

gas.TPX = 1000, 1e5, 'nc7h16:1.0,o2:11.0,n2:41.58'
#gas.TPX = 1000, 1e5, 'nc7h16:1.0,o2:xxxx,n2:xxxx' # fill in your calculated o2, n2 values.

reactor = ct.IdealGasReactor(gas)
reactor_network = ct.ReactorNet([reactor])

start_time = 0.0  #starting time
end_time = 0.5 # seconds
n_steps = 501
times = np.linspace(start_time, end_time, n_steps)
concentrations = np.zeros((n_steps, gas.n_species))
mole_frac = np.zeros((n_steps, gas.n_species))
pressures = np.zeros(n_steps)
temperatures = np.zeros(n_steps)

print_data = True
if print_data:
    #this just gives headings
    print('{0:>10s} {1:>10s} {2:>10s} {3:>14s}'.format(
            't [s]','T [K]','P [Pa]','u [J/kg]')) 

for n, time in enumerate(times):
    if time > 0:
        reactor_network.advance(time)
    temperatures[n] = reactor.T
    pressures[n] = reactor.thermo.P
    concentrations[n,:] = reactor.thermo.concentrations
    mole_frac[n,:] = gas.X
    if print_data:
        print('{0:10.3e} {1:10.3f} {2:10.3f} {3:14.6e}'.format(
                 reactor_network.time, reactor.T, reactor.thermo.P, reactor.thermo.u))


# #### c. What is the approximate ignition time? (hint: when do the temperature, pressure have the largest change?) Change initial temperature, pressure in the cell above, and comment on the change in ignition time.

# In[226]:

i_no = gas.species_names.index('no')
i_o2 = gas.species_names.index('o2')
i_oh = gas.species_names.index('oh')
plt.plot(times, mole_frac[:,i_no], label='NO')
#plt.plot(times, mole_frac[:,i_o2], label='O2')
plt.plot(times, mole_frac[:,i_oh], label='OH')
plt.xlim(0.01, 0.1)
plt.legend(loc='best')
print "The maximum mole fraction of NO is {:.2e}".format(np.max(concentrations[:,i_no]))


# #### d. At what conditions can we achieve a maximum NO mole fraction of 0.005?

# In[ ]:



