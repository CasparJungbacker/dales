#
#  Purpose  : this script generates all necessary input files (.inp.) for the GABLS1 case
#
#             The GABLS1 case is described in:
#             Beare, R.J. et al. An Intercomparison of Large-Eddy Simulations of the Stable Boundary Layer. Boundary-Layer Meteorol 118, 247–272 (2006). https://doi.org/10.1007/s10546-004-2820-6
#
#  Author   : Steven van der Linden, Delft University of Technology (TUD)
#  Contact  : s.j.a.vanderlinden@tudelft.nl
#  Date     : 21 November 2024
#
#  This file is part of DALES.

import numpy as np
import argparse

float_type = np.float64

#### Define the parser ####
parser = argparse.ArgumentParser(description='Give experiment number via -iexpr=XX')
parser.add_argument('-iexpr', action="store", dest='iexpr', type=int, default=001)
parser.add_argument('-expname', action="store", dest='expname', default='gabls1')

args = parser.parse_args()

# Get required parameters from namoptions.iexpr
with open('namoptions.'+str(args.iexpr)) as f:
    for line in f:
        if(str.rstrip(line.split('=')[0]) == 'kmax'):
            kmax = int(line.split('=')[1])

f.close()

print("The number of vertical levels given in namoptions.{} is equal to kmax={}".format(args.iexpr,kmax))

#### Create vertical grid ####
zsize = 400 # << Specified domain height in GABLS1
dz = zsize / kmax
z = np.linspace(0.5 * dz, zsize - 0.5 * dz, kmax)

print("An equidistant vertical grid has been constructed with vertical resolution dz={}".format(dz))

#### Initialise arrays for initial profiles ####
th = np.empty(z.shape)
u = np.empty(z.shape)
ug = np.empty(z.shape)
v = np.empty(z.shape)
vg = np.empty(z.shape)
qt = np.empty(z.shape)
tke = np.empty(z.shape)

wls = np.zeros(z.shape) # large scale subsidence velocity
emp = np.zeros(z.shape)
dthldt = np.zeros(z.shape)
dqtldt = np.zeros(z.shape)

#### Write prescribed values of GABLS1 case to these arrays ####
u[:]   = 8.
ug[:]  = 8.
v[:]   = 0.
vg[:]  = 0.
qt[:]  = 0.

dthetadz = 0.01

for k in range(kmax):
    if(z[k] <= 100.):
        th[k] = 265.
    if(z[k] > 100.):
        th[k] = 265. + dthetadz*(z[k]-100.)

for k in range(kmax):
    if(z[k] <= 250.):
        tke[k] = ( 0.4 * (1 - z[k]/250)**3 )
    if(z[k] > 250.):
        tke[k] = 0.

## Note: new str formatting in python
# In most of the cases the syntax is similar to the old %-formatting, 
# with the addition of the {} and with : used instead of %. For example, '%03.2f' can be translated to '{:03.2f}'.
# additional: > forces right alignment within availabe space

#### Create initialisation files for DALES ####
# prof.inp.iexpr is minimally needed and should contain the variables listed below.
# the first two lines are skipped during read phase.

f = open('prof.inp.'+str(args.iexpr),'w')
f.write('# EXPERIMENT: '+args.expname+'\n')
f.write('      height(m)   thl(K)     qt(kg/kg)       u(m/s)     v(m/s)     tke(m2/s2)\n')

for k in range(kmax):
    line = '{:>13.5f}{:>13.3f}{:>13.8f}{:>13.5f}{:>13.5f}{:>13.5f}\n'.format(z[k], th[k], qt[k], u[k], v[k], tke[k])
    f.write(line)

f.close()    

# lscale.inp.iexpr provides large scale forcings (if left uninitialised, variables are set to zero by default)

f = open('lscale.inp.'+str(args.iexpr),'w')

f.write('# EXPERIMENT: '+args.expname+'\n')
f.write('      height(m)   ugeo(m/s) vgeo(m/s)  wfls(m/s)    not_used   not_used     dqtdtls(kg/kg/s)    dthldt(K/s)\n')

for k in range(kmax):
    line = '{:>10.5f}{:>8.3f}{:>8.3f}{:>10.6f}{:>6.1f}{:>6.1f}{:>6.1f}{:>6.1f}\n'.format(z[k], ug[k], vg[k], wls[k], emp[k], emp[k], dqtldt[k], dthldt[k])
    f.write(line)

f.close() 

# ls_flux.inp.iexpr contains time dependent change of surface properties and/or profiles
# << GABLS1 uses a prescribed surface temperature over 9 simulation hours (isurf=2 in NAMOPTIONS)

f = open('ls_flux.inp.'+str(args.iexpr),'w')

f.write('EXPERIMENT: '+args.expname+'\n') # line not allowed to have '#'
f.write('       time      wtsurf      wqsurf    thls      qts   psurf\n') # in total have three lines.. 
f.write('       [s]     [K m/s]    [kg m/s]     [K]      [kg/kg]   [Pa]\n') # 
line = '{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(0.000, -9.000, -9.000, 265, -9.000, 100000.00) # Ts taken from ERA5, pressure taken from ERA5 but kept constant..
f.write(line)
line = '{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(32400.000, -9.000, -9.000, 262.75, -9.000, 100000.00)
f.write(line)
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('%Large scale forcings\n')
f.write('%    z [m]      ug [m/s] vg [m/s]      wmn [m/s]       dqtdx     dqtdy     dqtdtls     dthlrad\n')
f.write('#           0.00000\n')
for k in range(kmax):
     line = '{:>13.5f}{:>8.3f}{:>8.3f}{:>10.6f}{:>5.1f}{:>5.1f}{:>5.1f}{:>5.1f}\n'.format(z[k], ug[k], vg[k], wls[k], emp[k], emp[k], emp[k], emp[k])
     f.write(line)
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('%Large scale forcings\n')
f.write('%    z [m]      ug [m/s] vg [m/s]      wmn [m/s]       dqtdx     dqtdy     dqtdtls     dthlrad\n')
f.write('#           32400.000\n')
for k in range(kmax):
     line = '{:>13.5f}{:>8.3f}{:>8.3f}{:>10.6f}{:>5.1f}{:>5.1f}{:>5.1f}{:>5.1f}\n'.format(z[k], ug[k], vg[k], wls[k], emp[k], emp[k], emp[k], emp[k])
     f.write(line)

f.close()
