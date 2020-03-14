# Short user manual

Welcome to use NBODY6++GPU, which is a MPI parallel version of NBODY6GPU.
If you find issues during using the code but not fully answered in the follow tips and manuals in **doc/Nbody6++_manual.pdf**, 
please use the GitHub issue to ask questions.

# Install
Use the GNU configure tool to install:
	
	./configure [options]
	make
	make install

The options for configure can be found by using 
	
	./configure --help

The options include the maximum size of data array, the switchers for using parallelization methods and HDF5 support.

## some basic options

- Install path:
  
  The default install path is /usr/local, if you want to change to another path, please use:

		./configure --prefix=[YOURPATH]

  Replace *[YOURPATH]* to the full path that you want to install the code.
  There will be *bin*, *doc*, *include* and *lib* created in your install path.


- GPU:
  
  The GPU support is enabled in the default case, if you want to disable the GPU acceleration, please use:

		./configure --disable-gpu

- MPI:
  
  The MPI is enabled in the default case, if you want to disable the MPI parallelization, please use:

		./configure --disable-mpi

  Notice if no Fortran MPI compiler (mpif77) exist, the configure cannot pass, the *--disable-mpi* should be added manually

- AVX/SSE:

  The AVX is enable in default case, you can choose no AVX/SSE, SSE and AVX, the option is:

		./configure --enable-simd=[arg]

  Replace *[arg]* by *no*, *sse*, *avx*.

  Notice when the CPU don't support AVX/SSE, the code can still be compiled if the compiler support them. 
  In such case a warning will appear during the configuration. 
  Please be careful for this warning. 
  It indicates that the code cannot be used in the current machine.
  In the case for running code on a computer cluster with a job manage system, the code may be compiled in the login node without AVX/SSE support but can be used on the working nodes that support AVX/SSE. 
  Thus the configure does not switch off AVX/SSE even the CPU on the compiling node doesn't supports them.

- HDF5 output format:
  
  The HDF5 is disabled in default case, if you want to enable it, please use:

		./configure --enable-hdf5

  Notice when MPI is used, HDF5 compiler should be parallel version with Fortran enabled. 
  If you want to use non-parallel version of HDF5 compiler with no MPI, please add *--disable-mpi* during the configuration.

- Extra tools:

  There are some extra tools and libraries for reading conf.3 and do some basic analysis.
  It is disabled in the default case.
  If you want to switch on, please use:

		./configure --enable-tools

Notice all configure options should be used together. For example:

	./configure --prefix=/opt/nbody6++ --enable-tools --enable-hdf5

## Important notice:

- When large *NMAX* is used, sometimes the segmentation fault happens soon after the simulation starts. 
  This is due to the stack memory overflow. 
  In such case, you should always run 

		ulimit -s unlimited

  before start the simulation in the same shell.
  Be careful that *ulimit* only affects the current shell you enter this commander, thus it's need to be used every time when a new shell is opened.
  You can add this line in the shell initialization file, e.g., in *bash*, usually it is .bashrc or .bash_profile in your home directory. 
  After that, this commander is automatically loaded once a new shell is open.

- In some Linux systems (e.g. Ubuntu), when the large *NMAX* is used, the compiling crash with the error: 
 
		R_X86_64_PC32 against symbol

  This error is related to the memory model of the Linux system.
  The detail of memory model is shown in https://software.intel.com/en-us/node/579558 
  The solution is to add the configure option *--enable-mcmodel=large*

- The **stellar evolution package (SSE/BSE)** in the current version is based on the new update of Banerjee et al. (2019) (http://arxiv.org/abs/1902.07718)
