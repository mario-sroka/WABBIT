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
