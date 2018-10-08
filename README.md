# Ewald
Program to fit arrays of point charges to an Ewald potential

This is a fork of the Ewald program by Klintenberg, Derenzo and Weber. The original can be found at [Computer Physics Communications](http://cpc.cs.qub.ac.uk/summaries/ADME_v1_0.html).
The purpose of this fork is to work easily with the cryspy package.

# Dependencies

- lapack
- blas

# Installation

To install:

`cc -o Ewald Ewald.c randome.c -llapack -lblas -lm`

And add the `Ewald` binary to your `$PATH`.

To test:
```
cd test
Ewald < ewald.in.CaF2
```
Your output files in `Ewald/test` should be identical to the ones in `Ewald/test/output`.

Further instructions and usage are identical to the original program and can be accessed in [Ref. 1](https://doi.org/10.1016/S0010-4655(00)00071-0)

# Modifications
- Partial charges are allowed
- The lattice charge can be up to 10e-7 instead of a strict 0.0
- The maximum amount of atoms in the unit cell has been increased to 1,000 
- The maximum amount of atoms in the cluster has been increased to 10,000
- Outputs an additional point charge file with a convenient format with the extension .pts-cry
- Uses the dgelsy_ function from lapack instead of the deprecated dgelsx_ which is unavailable under some installations
- Beautified


# References
[1. Klintenberg, M., Derenzo, S. E. & Weber, M. J. Accurate crystal fields for embedded cluster calculations. Comput. Phys. Commun. 131, 120–128 (2000).](https://doi.org/10.1016/S0010-4655(00)00071-0)

[2. Derenzo, S. E., Klintenberg, M. K. & Weber, M. J. Determining point charge arrays that produce accurate ionic crystal fields for atomic cluster calculations. J. Chem. Phys. 112, 2074–2081 (2000).](https://doi.org/10.1063/1.480776)

