# Ewald
Program to fit arrays of point charges to an Ewald potential

This is a fork of the Ewald program by Klintenberg, Derenzo and Weber. The original can be found at [Computer Physics Communications](http://cpc.cs.qub.ac.uk/summaries/ADME_v1_0.html).
The purpose of this fork is to work easily with the cryspy package.



# Modifications
- Beautified
- Partial charges are allowed
- The lattice charge can be up to 10e-7 instead of a strict 0.0
- The maximum amount of atoms in the unit cell has been increased to 1,000 
- The maximum amount of atoms in the cluster has been increased to 10,000
- Outputs an additional point charge file with a conveninent format with the extension .pts-cry
- Uses the sgelsy_ from lapack instead of the deprecated  dgelsx_


# References
1. Derenzo, S. E., Klintenberg, M. K. & Weber, M. J. Determining point charge arrays that produce accurate ionic crystal fields for atomic cluster calculations. J. Chem. Phys. 112, 2074–2081 (2000).
2. Klintenberg, M., Derenzo, S. E. & Weber, M. J. Accurate crystal fields for embedded cluster calculations. Comput. Phys. Commun. 131, 120–128 (2000).
