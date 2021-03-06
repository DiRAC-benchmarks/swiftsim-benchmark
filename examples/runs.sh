#!/bin/bash

# Set some global stuff
export OMP_WAIT_POLICY=PASSIVE

# Generate the initial conditions if they are not present.
if [ ! -e SodShock/sodShock.hdf5 ]
then
    echo "Generating initial conditions for the SodShock example..."
    cd SodShock
    python makeIC.py
    cd ..
fi
if [ ! -e SedovBlast/sedov.hdf5 ]
then
    echo "Generating initial conditions for the SedovBlast example..."
    cd SedovBlast/
    python makeIC_fcc.py
    cd ..
fi
if [ ! -e CosmoVolume/cosmoVolume.hdf5 ]
then
    echo "Downloading initial conditions for the CosmoVolume example..."
    cd CosmoVolume
    ./getIC.sh
    cd ..
fi


# Loop over number of cores
for cpu in {1..32}
do

    # Sod-Shock runs
    if [ ! -e SodShock_${cpu}.dump ]
    then
        ./swift -t $cpu -f SodShock/sodShock.hdf5 -m 0.01 -w 5000 -c 1. -d 1e-7 -e 0.01 > SodShock_fixed_${cpu}.dump
    fi
    
    # Sedov blast
    if [ ! -e SedovBlast_${cpu}.dump ]
    then
        ./swift -t $cpu -f SedovBlast/sedov.hdf5 -m 0.02 -w 5000 -c 1. -d 1e-7 -e 0.01 > SedovBlast_fixed_${cpu}.dump
    fi
    
    # Cosmological volume
    if [ ! -e CosmoVolume_${cpu}.dump ]
    then
        ./swift -t $cpu -f CosmoVolume/cosmoVolume.hdf5 -m 0.6 -w 5000 -c 1. -d 1e-7 -e 0.01 > CosmoVolume_fixed_${cpu}.dump
    fi

done

