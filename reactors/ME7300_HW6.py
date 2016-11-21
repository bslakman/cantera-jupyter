
# coding: utf-8

# # Combustion and Air Pollution HW #6
# 
# #### Consider combustion of n-heptane ($C_7H_{16}$). Heptane and air are introduced to a constant volume chamber. First let's consider $\phi$ = 1.0 and atmospheric pressure. Upon ignition, combustion of the premixed fuel/air occurs at constant volume and constant internal energy (assuming that the chamber is adiabatic).
# 
# #### Calculate:
# 1. The adiabatic flame temperature if the fuel and air are preheated and introduced into the chamber at 500 K. 
# 2. The NOx concentration at equlibrium for $\phi$ = 0.7, 1.0, and 1.3
# 
# In this homework we will use more detailed combustion chemistry then we previously considered in this class, along with the Zeldovich mechanism for thermal fixation of atmospheric N2, to model and predict NOx formation in heptane combustion.
# We use an open source software tool called [Cantera](http://www.cantera.org/docs/sphinx/html/index.html) which helps us solve thermodynamics and kinetics problems.
# 
# This homework was adapted from a lecture in Prof. Richard West's undergraduate chemical kinetics class in the Department of Chemical Engineering as well as existing Python notebooks in the [cantera-jupyter project](https://github.com/cantera/cantera-jupyter). 
# 
# ### Simply click inside the cells below to be able to type in and edit them, and press Shift+Enter to execute the code in a cell that is selected.

# In[65]:

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import seaborn


# The n-heptane ($C_7H_{16}$) mechanism we will consider is from the Lawrence Livermore National Lab (LLNL), with the 3-reaction Zeldovich mechanism manually added to it. First we create a Solution object in Cantera for n-heptane gas that includes all of the thermodynamic and kinetic parameters needed to describe its combustion. (Don't worry about the warning)

# In[2]:

gas = ct.Solution('heptanesymp159.cti')


# We can print out some chemical species and reactions that are considered in the mechanism:

# In[38]:

print "Some of the species in the mechanism:"
print " ".join(gas.species_names[:20])
print "\n"

print "Model includes reactions that break down the fuel ('nc7h16'):"
for i in range(1052,1102,5):
    print gas.reaction(i)
print "\n"

print "And also contains small molecule chemistry:"
for i in range(20):
    print gas.reaction(i)
print "\n"
print "There are {0} total reactions in this mechanism!".format(len(gas.reactions()))


# Let's look at the last few reactions in the mechanism, representing the Zeldovich mechanism.

# In[4]:

zeldovich = gas.reactions()[-3:]
for rxn in zeldovich:
    print rxn.equation, rxn.rate


# #### 1. Calculate the adiabatic flame temperature if the fuel and air are preheated and introduced into the chamber at 500 K. 
# First, we need to initialize the fuel at the correct temperature (K), pressure (Pa) and mole fractions (these can be molar ratios, Cantera will normalize them, like in STANJAN).
# Replace the "xxxx" with the correct number of moles.

# In[68]:

gas.TPX = 500, 1e5, 'nc7h16:1.0,o2:11.0,n2:41.58'
#gas.TPX = 1000, 1e5, 'nc7h16:1.0,o2:xxxx,n2:xxxx'


# We'll find the equilibrium flame temperature by equilibrating keeping volume and internal energy
# constant.

# In[69]:

gas.equilibrate('UV')
print "The adiabatic flame temperature is {0} K".format(int(round(gas.T)))


# #### 2. Calculate the NOx concentration at equlibrium for $\phi$ = 0.7, 1.0, and 1.3
# Change the numbers of moles for these equivalence ratios,below.

# In[70]:

x_NO_stoich = gas['no'].X[0]

# For fuel lean case (0.7)
gas.TPX = 500, 1e5, 'nc7h16:1.0,o2:15.71,n2:59.4'
#gas.TPX = 1000, 1e5, 'nc7h16:1.0,o2:xxxx,n2:xxxx'
gas.equilibrate('UV')
x_NO_lean = gas['no'].X[0]

# For fuel rich case (1.3)
gas.TPX = 500, 1e5, 'nc7h16:1.0,o2:8.46,n2:31.98'
#gas.TPX = 1000, 1e5, 'nc7h16:1.0,o2:xxxx,n2:xxxx'
gas.equilibrate('UV')
x_NO_rich = gas['no'].X[0]

print "For fuel lean case, NO = {0} ppm".format(round(x_NO_lean*1e6, 1))
print "For stoichiometric case, NO = {0} ppm".format(round(x_NO_stoich*1e6, 1))
print "For fuel rich case, NO = {0} ppm".format(round(x_NO_rich*1e6, 1))


# It may look like the NO is decreasing with equivalence ratio, but there is a maximum in the fuel lean region...

# In[78]:

equivalence_ratios = np.linspace(0.1, 2, 50)
x_NO = np.zeros_like(equivalence_ratios)
for i, phi in enumerate(equivalence_ratios):
    mol_O2 = 11.0 / phi
    mol_N2 = 3.78 * mol_O2
    X_string = "nc7h16:1.0,o2:" + str(round(mol_O2, 2)) + ",n2:" + str(round(mol_N2, 2))
    gas.TPX = 500, 1e5, X_string
    gas.equilibrate('UV')
    x_NO[i] = gas['no'].X[0]
plt.plot(equivalence_ratios, x_NO*1E6)
plt.xlabel('$\phi$')
plt.ylabel('NO (ppm)')
print "The NO is maximum at equivalence ratio " + str(round(equivalence_ratios[np.argmax(x_NO)],3))


# ### Kinetics
# We will burn the gas starting at the same temperature, but atmospheric pressure. This time we will hold volume constant using an ideal gas reactor, so pressure can change.
# The following code will print out the simulation time along with temperature, pressure, and internal energy.

# In[25]:

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

# In[21]:

i_no = gas.species_names.index('no')
#i_o2 = gas.species_names.index('o2')
i_oh = gas.species_names.index('oh')
plt.plot(times, mole_frac[:,i_no], label='NO')
#plt.plot(times, mole_frac[:,i_o2], label='O2')
plt.plot(times, mole_frac[:,i_oh], label='OH')
plt.xlim(0.03, 0.04)
plt.legend(loc='best')
print "The maximum mole fraction of NO is {:.2e}".format(np.max(mole_frac[:,i_no]))


# #### d. At what conditions can we achieve a maximum NO mole fraction of 0.005?

# Just want to check at least one rate to make sure it's reasonably close to the one in the book...These came from the gri mech.

# In[5]:

T = np.linspace(300,2500,100)


# In[6]:

A = rxn.rate.pre_exponential_factor / 1000.0
n = rxn.rate.temperature_exponent
Ea = rxn.rate.activation_energy / 1000.0


# In[7]:

rate_mech = A * T**n * np.exp(-Ea/8.314/T)
rate_book = 7.1E7 * np.exp(-450.0/T)


# In[8]:

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[9]:

plt.plot(1.0/T, np.log(rate_mech))
plt.plot(1.0/T, np.log(rate_book))


# In[ ]:



