# chiralDFT

This package contains a massively-parallel hybrid Density Functional Theory/Monte-Carlo (DFT/MC) simulation code, designed to infer the chiral nematic behaviour of lyotropic liquid crystals from their microscopic mesogen properties.

It includes a high-performance and versatile implementation of the algorithm introduced by [Dussi et al.](http://dx.doi.org/10.1103/PhysRevE.90.020503), with additional support for [Straley’s perturbative approach](http://dx.doi.org/10.1103/PhysRevA.14.1835).

An implementation of the [RAPID](http://dx.doi.org/10.1145/237170.237244) collision detection library, developed by Manocha et al., is also provided for efficient overlap queries between sets of complex triangle meshes in hard-body simulations.

It further contains an original bounding volume hierarchy based on principal component analysis, efficiently generating and traversing binary trees of bounding structures to speed up energy calculations by several orders of magnitude for complex particle models.



## Requirements

* The [**Eigen**](https://eigen.tuxfamily.org) high-performance linear algebra library
* `gcc` >= 4.9 with an associated `mpicc` wrapper, provided by:
* a working MPI implementation, e.g. [**OpenMPI**](https://www.open-mpi.org)



## Compilation

Checkout the source code from the repository through either `svn` or direct download, then:

~~~shell
cd distruc-chiraldft-code	# enter the project folder
make librapid	  			# build RAPID mesh collision library
make -j4					# compile chiralDFT in parallel using 4 threads
~~~

This will build the `chiraldft` executable into the newly created `bin` folder.

Note that if you chose to install the **Eigen** library headers or `mpicc` linker in custom locations, you may need to change the relevant fields in the makefile to the chosen install paths.



## Executing

To run the code:

~~~shell
cd bin
mpirun -np <number_of_cores> ./chiraldft
~~~

with `<number_of_cores>` the number of available cores you wish to use for the simulation.
On most recent **Intel**® processors (Core™ i5 and later generations), setting `<number_of_cores>` to 8 will yield maximum performance, using the 4 physical cores + 4 virtual cores through [*Hyper-threading*](https://en.wikipedia.org/wiki/Hyper-threading). However, doing so will stress your CPU to 100%, so don't attempt to run this on your laptop for extended periods of time.



## Compile options

Some very limited options to be set in the `include/globals.hpp` file:

* set the `MESOGEN` symbol to whichever particle geometry you fancy (I might write a very short intro on how to modify each particle’s characteristics if anyone is interested for some reason)
* set the `FULL_RUN` switch to either 0 for a preliminary perturbative run or 1 for a full sweep of the chiral free energy landscape. Bear in mind that the full run is VERY expensive (especially for the `DNADuplex` particle type), so don’t try it on your home computer.

Don’t forget to recompile the code after modifying any of these flags.

Note that the `DNADuplex` particle type requires as input a DNA trajectory file generated by the [**oxDNA**](https://sourceforge.net/projects/oxdna/) model, formatted using the `pre_process.py` script located in the `resources/processing` folder.



## Utilities

A couple plotting utilities to be found in the `resources` folder (requires `gnuplot` >= 5):

* `utils/display_wireframe` conjures up a 3d interactive viewer for the chosen particle
* `plots/plot_landscape3d` plots the 3d free energy surface as a function of both particle density and macroscopic pitch.

Note the lists above are far from exhaustive - don’t hesitate to [contact me](mailto:maxime.tortora@chem.ox.ac.uk) for information about the more advanced options and functionalities.
