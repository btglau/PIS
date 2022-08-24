# PIS
MATLAB code to calculate integrals of spherical harmonics x bessel functions, for particle-in-a-sphere models

Send the function pis_fock4.m to a cluster or your local computer with appropriate input options

The eris are listed in 8-fold reduced format,
(ij|kl), i>=j, k>=l, ij>=kl,
for use with pyscf, by replacing the _eri attribute

You may need to mex the function (pis_8f_pw2) that generates the indices
