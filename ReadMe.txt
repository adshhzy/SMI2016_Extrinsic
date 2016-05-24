Extrinsic - Zhiyang Huang (2016)

------------------------------------

This code implements the algorithm described in

   Huang Z, Ju T. Extrinsically smooth direction fields. Comput Graph(2016), http://dx.doi.org/10.1016/j.cag.2016.05.015i

The basic function of this program is to compute extrinsically smooth direction fields “g” upon a curve, surface or volume while giving a references field “f”, such that “g” is orthogonal to “f”. In particular, the default setting of “f” for curve and surface are tangent field of the curve and normal field of the surface, respectively.

The parameter-free initialization method based on least eigenvector problem and the Gauss–Seidel iterations, as described in the above paper, are implemented in this code and subjected to use.

Currently the code is only tested on Mac OS X 10.10 or above, it should work at other platforms with minor modification.


BUILDING
======================================================================================================


The code has only one dependency: Eigen/Eigen3

Eigen is a C++ template library for linear algebra, and it is a header-based library which can be use right away without binary library to link to, and no configured header file. For different platforms, it should work without modification. see http://eigen.tuxfamily.org/index.php?title=Main_Page for further instruction.

For your convenience, we include the Eigen3 files in the repository. You can use other version of the Eigen3, but be sure to change the path (EIGEN_PATH) in Makefile if you place it in a different folder.

To build the program, just simply typing

	$make

in the root directory.  The result will be an executable called “extrinsic” (or “extrinsic.exe” on Windows if it is successfully built).


RUNNING
======================================================================================================

To run the code from the command line, type

	$./extrinsic -i input_file_name -o output_file_name [-e] [-g]

where:
1. -i: followed by the path of the input file. input_file_name is a path to the input file, currently support file format includes “.off”, “.obj”, “.curf”, “.surf”, “.volf”;
 
2. -o: followed by the path of the output file. output_file_name is a path to the output file, you don’t need to specify the extension of the output file, the program will define it to be either “.curf”, “.surf” or “.volf”, according to the input file. 

3. -e: optional argument. If -e is included in the command line, the “eigenvector initialization” will be run, otherwise the direction filed will be randomly initialized. 

4. -g: optional argument. If -g is included in the command line, Gauss–Seidel iterations will be run at the initialization field until converge or exceed max iterations.

(The detail format of “.curf”, “.surf”, “.volf” can be found below.)

*An additional output call “output_file_name.ply” will be created for a rough visualization of the generated field, which can be viewed in standard mesh viewers.

*Currently we use a build-in solver in Eigen to compute the least eigenvector (solve sparse linear equations), to further speed it up, you can use any other solvers. We recommend the SuiteSparse, whose wrapper has been already provided in Eigen3.



Generating Figures in the paper
======================================================================================================

To generate the files for related figures in the paper, type in the root directory:

	$chmod +x makefigure.sh
	$./makefigure.sh

Be sure that the program has been built in the root directory, and the “input” folder, “output” folder are in the in the root directory.

Noted that for generating the files for some figures (e.g fig 5, Biarc), we directly build the model in the program and do not require a input model.

Two output files will be generated in “output” folder as “figx_x.curf” (or “figx_x.surf”, “figx_x.volf”), and “figx_x.ply”. The “figx_x.ply” can be open at standard mesh viewers, but the rendering effect may be different from the exact corresponding figures in the paper.



FORMAT
======================================================================================================

.curf:

“.curf” file mimics the format of “.obj” (see https://en.wikipedia.org/wiki/Wavefront_.obj_file) with partial and modified key words:

v vx vy vz                  v -0.404425 -0.176438 0.554655
e v1 v2           =====>    e 1 20
vn vnx vny vnz    =====>    vn -0.00612673 -0.0156757 0.999858
vf vfx vfy vfz              vf -0.214613 0.652417 0.726838

where:
key word “v” is followed by x, y, z coordinate of a vertex;
key word “e” is followed by the index of vertices which make up the edge;
key word “vn” is followed by the x, y, z coordinate of the reference vector “f” on a vertex (curve tangent by default);
key word “vf” is followed by the x, y, z coordinate of the generated vector which is orthogonal to “vn”, on a vertex;
The order of “vn” and “vf” must follow the order of “v”, which means that the reference vector/field maps to the each vertex orderly.

————————————————————————————————————————————————————————————————————————————————————————————————————
.surf:

“.surf” file mimics the format of “.obj” with partial and modified key words:

v vx vy vz                  v -0.404425 -0.176438 0.554655
f v1 v2 v3        =====>    f 1 20 34
vn vnx vny vnz    =====>    vn -0.00612673 -0.0156757 0.999858
vf vfx vfy vfz              vf -0.214613 0.652417 0.726838

where:
key word “v” is followed by x, y, z coordinate of a vertex;
key word “f” is followed by the index of the three vertices which make up the triangle;
key word “vn” is followed by the x, y, z coordinate of the reference vector “f” on a vertex (surface normal by default);
key word “vf” is followed by the x, y, z coordinate of the generated vector which is orthogonal to “vn”, on a vertex;
The order of “vn” and “vf” must follow the order of “v”, which means that the reference vector/field maps to the each vertex orderly.

————————————————————————————————————————————————————————————————————————————————————————————————————
.volf:

“. volf” file mimics the format of “.obj” with partial and modified key words:

v vx vy vz                      v -0.404425 -0.176438 0.554655
tet v1 v2 v3 v4       =====>    tet 1 20 34 565
vn vnx vny vnz        =====>    vn -0.00612673 -0.0156757 0.999858
vf vfx vfy vfz                  vf -0.214613 0.652417 0.726838

where:
key word “v” is followed by x, y, z coordinate of a vertex;
key word “tet” is followed by the index of the four vertices which make up the tetrahedra;
key word “vn” is followed by the x, y, z coordinate of the reference vector “f” on a vertex;
key word “vf” is followed by the x, y, z coordinate of the generated vector which is orthogonal to “vn”, on a vertex;
The order of “vn” and “vf” must follow the order of “v”, which means that the reference vector/field maps to the each vertex orderly.