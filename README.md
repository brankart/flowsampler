## Flowsampler : Sampling of flows from a probability distribution

Flowsampler currently includes the sampling of shallow-water flows
and boundary layer flows.

Flowsampler is distributed under the terms of the CeCILL free software license agreement.
See LICENSE.txt for more information.

Flowsampler code is written in Fortran.

### Installation of Flowsampler

You need a FORTRAN-90 compiler and the NetCDF library (with f90 support) installed.

To compile the library and the examples :

- create a 'make.(MY_MACHINE)' file corresponding to your compiler in the 'macro' directory.
  This is the Makefile configurable part, which specifies
  your compiler options and where to find the NetCDF library.

```bash
cd build
ln -sf ../macro/make.(MY_MACHINE) Makefile.macro
```

- compile (library and examples) with:

```bash
cd build
make
make examples
```

- if everything goes well, the flowsampler library 'libflowsampler.a'
  has been created in the 'lib' directory, the module files (.mod)
  have been created in the 'include' directory' and the example
  exexutables (.x) have been created in the 'examples' directory.

 To check the installation :

 - run the examples in the 'example' directory, with appropriate parameter files:

```bash
cd examples
./bl_profile.x
./bl_sampler.x
./flow_sampler.x
```
