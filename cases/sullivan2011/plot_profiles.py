#
#  Purpose  : this script plots very basic profiles from output in profiles.iexpr.nc and tmser.iexpr.nc 
#               (typically averaged over the final hour)
#
#  required python dependencies are numpy, netCDF4 and matplotlib    
#
#  Author   : Steven van der Linden, Delft University of Technology (TUD)
#  Contact  : s.j.a.vanderlinden@tudelft.nl
#  Date     : 21 November 2024
#
#  This file is part of DALES.

import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.pyplot    as plt


#### Set timeframe over which to average ####
# average over final hour of simulation (9-10); specify time interval
time_begin    = 14400 - 1*3600
time_end      = 14400 + 0.0001 # << the adding of small value is required

##########
# General Figure settings
##########

mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['grid.linewidth'] = 1
# mpl.rcParams['axes.grid'] = True
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['legend.fontsize'] = 15

##########
# Plot figures
##########

#### Initialise the figures ####
fig1=plt.figure(figsize=(20,10))
ax11=fig1.add_axes([0.10,0.1,0.38,0.8])
ax12=fig1.add_axes([0.52,0.1,0.38,0.8])

fig2=plt.figure(figsize=(20,10))
ax21=fig2.add_axes([0.10,0.1,0.38,0.8])
ax22=fig2.add_axes([0.52,0.1,0.38,0.8])

fig3=plt.figure(figsize=(20,10))
ax31=fig3.add_axes([0.10,0.52,0.8,0.38])
ax32=fig3.add_axes([0.10,0.1,0.8,0.38])

#### Open files and start plotting ####
filename = 'profiles.001.nc'
data = nc.Dataset(filename,'r')

f2 = 'tmser.001.nc'
data2= nc.Dataset(f2,'r')

# Find correct indices corresponding to averaging times
# Convert to array of indices instead of array of True, False. 
indices_logic   = np.logical_and(data['time'][:]>=time_begin, data['time'][:]<=time_end)
indices_numeric = np.where(indices_logic)[0]

# Figure 1: wind, temperature
ax11.plot(((data['v'][indices_numeric]**2 + data['u'][indices_numeric]**2)**0.5).mean(axis=0),data['zt'][:],linestyle='solid', linewidth=2)

ax12.plot(data['thv'][indices_numeric].mean(axis=0),data['zt'][:],linestyle='solid', linewidth=2, label='virtual pot. temp.')
ax12.plot(data['thl'][indices_numeric].mean(axis=0),data['zt'][:],linestyle='dotted', linewidth=2, label='liquid water pot. temp.')

# Figure 2: momentum and temperature fluxes
ax21.plot(data['uwt'][indices_numeric].mean(axis=0),data['zm'][:],linestyle='solid', linewidth=2, label='total')
ax21.plot(data['uwr'][indices_numeric].mean(axis=0),data['zm'][:],linestyle='dashed', linewidth=2, label='resolved')

ax22.plot(data['wthlt'][indices_numeric].mean(axis=0),data['zm'][:],linestyle='solid', linewidth=2)
ax22.plot(data['wthlr'][indices_numeric].mean(axis=0),data['zm'][:],linestyle='dashed', linewidth=2)

ax31.plot(data2['time'][:],data2['ustar'][:],linestyle='solid', linewidth=2)
ax32.plot(data2['time'][:],data2['obukh'][:],linestyle='solid', linewidth=2)


#### Set axes, labels, etc. ####
ax11.grid()
ax12.grid()
ax11.set_ylabel(r'$z\,\mathrm{[m]}$',size=15)
ax11.set_xlabel(r'$U\,\mathrm{[m\,s^{-1}]}$',size=15)
ax12.set_xlabel(r'$\theta\,\mathrm{[K]}$',size=15)
ax11.set_ylim(0,250)
ax12.set_ylim(0,250)
lg11=ax11.legend(loc='best')

ax21.grid()
ax22.grid()
ax21.set_ylabel(r'$z\,\mathrm{[m]}$',size=15)
ax21.set_xlabel(r'$F_u\,\mathrm{[m^2\,s^{-2}]}$',size=15)
ax22.set_xlabel(r'$F_\theta\,\mathrm{[K\,m\,s^{-1}]}$',size=15)
ax21.set_ylim(0,250)
ax22.set_ylim(0,250)
lg21=ax21.legend(loc='best')

ax31.grid()
ax32.grid()
ax31.set_ylabel(r'$u_*\,\mathrm{[m\,s^{-1}]}$',size=15)
ax31.set_xticklabels([])
ax32.set_ylabel(r'$L_O\,\mathrm{[m]}$',size=15)
ax32.set_xlabel(r'$t\,\mathrm{[s]}$',size=15)
lg31=ax31.legend(loc='best')

#### Save and show ####
fig1.savefig('Wind_Temperature.png',format='png')
fig2.savefig('Fluxes_dearsmag.png',format='png')
fig3.savefig('FrictionVel_Obukhovlength.png',format='png')

plt.show()
