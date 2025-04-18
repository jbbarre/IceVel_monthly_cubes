#!/bin/bash
# to compile cube2vel on GRICAD using cubes_py313 conda environnement
# Date: April 17th 2025
# Author: Jean Baptiste Barré 
# needs intel compiler to be installed . check nix in gricad's documentation.
# need to activate nix and intel compiler with:
#   source /applis/site/nix.sh
#   source $INTEL_ONEAPI/setvars.sh


set -e

SRC="cube2vel_feather_sp_f2py.F90"
MODULE="cube2vel_ifort"

export PATH=/home/barrejb-ext/.nix-profile/bin:$PATH

f2py --fcompiler=intelem --compiler=intelem --opt='-O2 -axAVX,SSE4.2' \
     --f90flags='-fPIC' -c $SRC -m $MODULE

echo "✅ Build complete: $MODULE.cpython-313t-x86_64-linux-gnu.so"

