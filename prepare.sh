#!/bin/bash

set -euo pipefail

# Builds the programs which are not distributed through conda (the other dependencies are in environment.yml)

# Our work - Pareto optimization, matchtig count lowerbound, sorted integers
make -C kmercamel-pareto

# Elias-Fano encoding for ex3 (uses only headers from 2i_bench)
make -C tools/EliasFano
