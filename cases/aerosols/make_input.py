import re

import numpy as np
import pandas as pd
import netCDF4 as nc

DROP = ["#z", "QC", "QR", "NC", "NR"]
DEPHY_NAMES = [""]

if __name__ == "__main__":
    # First, scalar.inp file
    scalars = pd.read_csv("scalar.inp.AERO", skiprows=1, sep='\s+')

    drop_pattern = re.compile(r'.*(CLD|RAI)')
    DROP = DROP +  list(filter(drop_pattern.search, scalars.keys()))
    scalars = scalars.drop(columns=DROP)

    with nc.Dataset("tracers.001.nc", "w") as ds:
        zt = ds.createDimension("zt", len(scalars))

        for item in scalars.keys():
            if not(item[-1] == "N"):
                name = item[:-3] + "_" + item[-3:]
            else:
                name = item
            var = ds.createVariable(name.lower(), "f", (zt,))
            var[:] = scalars[item]
            var.laero = 1

    aerosols = {
        "so4": {
            "long_name": "sulphate",
            "rho": 1841.0,
            "kappa": 0.88,
            "modes": "nus,ais,acs,cos"             
        },
        "bc": {
            "long_name": "black carbon",
            "rho": 1300.0,
            "kappa": 0.0,
            "modes": "ais,acs,cos,aii"
        },
        "pom": {
            "long_name": "organic matter",
            "rho": 1800.0,
            "kappa": 0.1,
            "modes": "ais,acs,cos,aii"
        },
        "ss": {
            "long_name": "sea salt",
            "rho": 2165.0,
            "kappa": 1.28,
            "modes": "acs,cos"
        },
        "du": {
            "long_name": "mineral dust",
            "rho": 2650.0,
            "kappa": 0.0,
            "modes": "acs,cos,aci,coi"
        }
    }

    with nc.Dataset("aerosol.001.nc", "w") as ds:
        for item in aerosols.keys():
            var = ds.createVariable(item, "i") 
            var.long_name = aerosols[item]["long_name"]
            var.rho = aerosols[item]["rho"]
            var.kappa = aerosols[item]["kappa"]
            var.modes = aerosols[item]["modes"]
