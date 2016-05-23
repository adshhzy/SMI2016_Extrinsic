Extrinsic - Zhiyang Huang (2016)

------------------------------------

This code implements the algorithm described in

   Huang Z, Ju T. Extrinsically smooth direction fields. ComputGraph(2016), http://dx.doi.org/10.1016/j.cag.2016.05.015i

The basic function of this program is to compute extrinsically smooth direction fields “g” upon curve, surface or volume while giving a references field “f”. In particular, the default setting of “f” for curve and surface are tangent field of the curve and normal field of the surface, respectively.

The parameter-free initialization method based on least eigenvector problem and the Gauss–Seidel iterations, as described in the above paper, are implemented in this code and subjected to use.

Currently the code is only tested on Mac OS X 10.10 or above, it should work at other platforms with minor modification.

BUILDING
======================================================================================================


The code has only one dependency: Eigen/Eigen3

Eigen is a C++ template library for linear algebra, and it is a header-based library which can be use right away without binary library to link to, and no configured header file. For different platforms, it should work without modification. see http://eigen.tuxfamily.org/index.php?title=Main_Page for further instruction.

To make it easier for using, we include the Eigen3 files in the res