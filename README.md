# README.md

# Polarization-dependent 1D RSA (MATLAB)

MATLAB code to simulate a 1D Random Sequential Adsorption (RSA) process in which
the adsorption outcome depends on the boundary polarities of each available
interval and on the dipole orientation chosen at each attempt.

This repository contains:
- `interact_rsa.m`: single-run RSA simulation (returns the number of deposited unit dipoles).
- `explore_final.m` (or your script file): parameter scan over interval length `x` with many replicas
  (uses `parfor` if available) and saves results to a MAT file.

## Requirements
- MATLAB (any reasonably recent version should work).
- Parallel Computing Toolbox **optional** (needed only if you want to use `parfor`).

## Quick start

1. Clone/download the repository and open MATLAB in this folder.
2. Make sure `interact_rsa_final.m` is on the MATLAB path (placing it in the same folder is enough).
3. Run the replica driver script: explore_final.m
