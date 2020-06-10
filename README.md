# WABBIT - RNS (reactive navier stokes physics)
## (W)avelet (A)daptive (B)lock-(B)ased solver for (I)nteractions with (T)urbulence

WABBIT-RNS is a fork from WABBIT, to implement separatly reactive navier stokes physics and other things. Due to heavy changes in some of the core routines, a merge with the WABBIT master branch would not be an easy task. But someday WABBIT-RNS will be synchronized again with the master.

New in comparison to WABBIT are:

+ reactive navier stokes physics
+ turbulence forcing

## INSTALL

1. clone from git

```
https://github.com/mario-sroka/WABBIT.git
```

2. install [MPI library](https://www.open-mpi.org/) 

3. install [HDF5 library](https://www.hdfgroup.org/downloads/hdf5/source-code/ "HDF5 Source Code")

	necessary variables for WABBIT
```
export HDF_ROOT=[...]
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HDF_ROOT/lib64
```

4. install [CANTERA library](https://www.cantera.org/), build FORTRAN module

	set up variables for WABBIT
```
export CANTERA_ROOT=[.../cantera_install]
export CANTERA_DATA=[.../cantera/data/inputs]
export PKG_CONFIG_PATH=[.../cantera_install/lib64/pkgconfig]
export PATH=$PATH:[.../cantera_install/bin]
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CANTERA_ROOT/lib64
```

## run WABBIT-RNS

Customize the .ini-files in CASES directory, run WABBIT with .ini-file name

```
wabbit [path_to_your_ini_filename.ini] 
```

## Test Cases

<table>
        <tr>
            <th>double shear layer</th>
            <th>density</th>
            <th>vorticity</th>
        </tr>
        <tr>
            <td rowspan=10>
<ul style="list-style-type:disc;">
  <li>2D inert gas</li>
  <li>Domain size: 1mx1m</li>
  <li>grid: 2048x2048</li>
  <li>periodic boundaries</li>
  <li>shock filtering <a href="https://doi.org/10.1016/j.jcp.2008.10.042">[Bogey2009]</a></li>
</ul>
</td>
            <td rowspan=10><img src="pics/rho.gif" width="60%"></td>
            <td rowspan=10><img src="pics/vort.gif" width="60%"></td>
        </tr>
</table>

## Publications

* ["An Open and Parallel Multiresolution Framework Using Block-Based Adaptive Grids"](https://link.springer.com/chapter/10.1007%2F978-3-319-98177-2_19 "Sroka2018"); Sroka, Engels, Krah, Mutzel, Schneider, Reiss; Active Flow and Combustion Control 2018
