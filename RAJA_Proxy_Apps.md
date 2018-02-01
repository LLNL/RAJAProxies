# RAJA Proxy Applications

Each proxy application available in this project is contained in its
own subdirectory. 

These applications have been released individually and contain their own 
license and release information.

## LULESH v1.0

LULESH v1.0 (Livermore Unstructured Lagrangian Explicit Shock Hydrodynamics)
is a proxy app that models a Sedov blast wave using an explicit Lagrangian
hydrodynamics algorithm on an unstructured mesh. It was developed 
originally by:

* Jeff Keasler (keasler1@llnl.gov)
* Rich Hornung (hornung1@llnl.gov)

Details about LULESH v1.0, including versions implemented in a variety
of programming models, can be found at 
[https://codesign.llnl.gov/lulesh.php](https://codesign.llnl.gov/lulesh.php).

Depending on CMake configuration options you provide, this repository
can generate sequential and OpenMP baseline versions (i.e., non-RAJA) and
RAJA versions for sequential, OpenMP, and CUDA GPU execution. 

By default, each version of LULESH v1.0 included here runs the
Sedov problem to a pre-defined end time and reports execution timing,
a figure of merit (FOM) grind-time, and verification that the solution is 
correct. Command line arguments allow one to optionally:

 * set a fixed number of time cycles to run
 * set the number of elements in the mesh 
 * print run progress that reports information about each time cycle. 

To get information about runtime options, run an executable with the 
`-h` option.

## LULESH v2.0

LULESH v2.0 is an extension of LULESH v1.0. LULESH v2.0 also models a Sedov 
blast wave using an explicit Lagrangian hydrodynamics algorithm on an 
unstructured mesh.

It was extended from LULESH v1.0 by:

* Ian Karlin (karlin1@llnl.gov)
* Jeff Keasler (keasler1@llnl.gov)
* Rob Neely (neely4@llnl.gov)

Details about LULESH v2.0 can be found at
[https://codesign.llnl.gov/lulesh.php](https://codesign.llnl.gov/lulesh.php).
[LULESH 2.0 Updates and Changes](https://codesign.llnl.gov/pdfs/LULESH2.0_Changes.pdf) describes how LULESH v2.0 is different that LULESH v1.0.

Depending on CMake configuration options you provide, this repository
can generate sequential and OpenMP baseline versions (i.e., non-RAJA) and
RAJA versions for sequential, OpenMP, and CUDA GPU execution. 

By default, each version of LULESH v2.0 included here runs the
Sedov problem to a pre-defined end time and reports execution timing,
a figure of merit (FOM) grind-time, and verification that the solution is 
correct. Similar to LULESH v1.0, command line arguments allow one to 
optionally:

 * set a fixed number of time cycles to run
 * set the number of elements in the mesh 
 * print run progress that reports information about each time cycle. 

In addition, command line options are provided for:

 * setting the number of material regions
 * setting load imbalance and cost parameters between regions
 * generating output files for visualization with VisIt

To get information about runtime options, run an executable with the 
`-h` option.
