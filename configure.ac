#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([Nbody6++], [1.0], [longwang.astro@gmail.com])
#AC_CONFIG_MACRO_DIR([./macro])
m4_include([macro/ax_check_compile_flag.m4])
m4_include([macro/ax_fc_check_bounds.m4])
#m4_include([macro/ax_check_cuda.m4])
#AC_CONFIG_SRCDIR([helper_cuda.h])

RESULT=
EXTRASRC=
EXTRAOBJ=

#Save environment
##AC_ARG_ENABLE([local-env],
##   [AS_HELP_STRING([--enable-local-env],
##        [Use local environment listed in last part instead of configure])],
##   [],
##   [enable_local_env=no])

##AS_IF([test "x$enable_local_env" != xno ],
##     [AS_IF([test "x$FC" != x], [local_fc=$FC])
##      AS_IF([test "x$FCFLAGS" != x], [local_fcflags=$FCFLAGS])
##      AS_IF([test "x$CXXFLAGS" != x], [local_cxxflags=$CXXFLAGS])],
##    [local_fc=
##     local_fcflags=
##     local_cxxflags=])

#Features controls
AC_ARG_ENABLE([mpi],
     [AS_HELP_STRING([--disable-mpi],
         [disable MPI support])],
     [],
     [enable_mpi=yes])

AC_ARG_ENABLE([mcmodel],
     [AS_HELP_STRING([--enable-mcmodel=ARG],
         [Use compiler option -mcmodel, avaiable arguments --> small; -->medium; -->large
          (defalut is medium)])],
     [ case "${enableval}" in
       large) enable_mcmodel=large ;;
       medium) enable_mcmodel=medium ;;
       small) enable_mcmodel=small ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mcmodel) ;;
       esac],
     [enable_mcmodel=medium])

AC_ARG_ENABLE([debug],
     [AS_HELP_STRING([--enable-debug],
         [Switch on debugging mode])],
     [],
     [enable_debug=no])

AC_ARG_ENABLE([gpu],
     [AS_HELP_STRING([--disable-gpu],
         [Disable GPU (CUDA) acceleration support for regular force and potential])],
     [],
     [enable_gpu=yes])

AC_ARG_ENABLE([simd],
     [AS_HELP_STRING([--enable-simd=ARG],
         [Intel CPU SIMD instruction set acceleration for force and potential, avaiable arguments --> sse: Streaming SIMD Extensions; --> avx: Advanced Vector Extensions; --> no: switch off SIMD; (default is autodetecting)])],
     [ case "${enableval}" in
       avx) enable_simd=avx ;;
       sse) enable_simd=sse ;;
       no) enable_simd=no  ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-simd) ;;
       esac],
     [enable_simd=test])

AC_ARG_ENABLE([hdf5],
     [AS_HELP_STRING([--enable-hdf5],
         [Enable HDF5 output format support. (Ref: Cai, Maxwell, et al., 2015, submitted)])],
     [],
     [enable_hdf5=no])

AC_ARG_ENABLE([tt],
     [AS_HELP_STRING([--enable-tt],
         [Enable tidal tensor and analytical potential extension imported from NBODY6TT. (Ref: Renaud, F., et al., 2011, MNRAS, 418, 759; Renaud, F. & Gieles, M., 2015, MNRAS, 448, 3416)])],
     [],
     [enable_tt=no])

AC_ARG_ENABLE([tools],
    [AS_HELP_STRING([--enable-tools],
         [Compile extral tools for data analysis, see User_manual.pdf for detail])],
    [enable_tools=yes],
    [enable_tools=no])

# Checks for programs.
# Fortran part

AC_LANG([Fortran])

local_fc=$FC
local_fcflags=$FCFLAGS
local_fflags=$FFLAGS

AC_PROG_FC
fc_check=$FC
FC=$local_fc

# Checks for MPI and HDF5 support
AS_IF([test "x$enable_mpi" != xno],
     [AS_IF([test "x$enable_hdf5" != xno],
       [AC_CHECK_PROG([FC], [h5pfc], [h5pfc], [none])],
       [AC_CHECK_PROG([FC], [mpif77], [mpif77], [none])])],
     [AS_IF([test "x$enable_hdf5" != xno],
      [AC_CHECK_PROG([FC], [h5fc], [h5fc], [none])],
      FC=$fc_check)
      local_fcflags=$local_fcflags' -I../extra_inc/nompi'])

AS_IF([test "x$enable_debug" != xno],
     [FCFLAGS=$local_fcflags' -g'
      AX_FC_CHECK_BOUNDS
      AX_CHECK_COMPILE_FLAG([-fbacktrace],[FCFLAGS=$FCFLAGS" -fbacktrace"])
      AX_CHECK_COMPILE_FLAG([-Wall],[FCFLAGS=$FCFLAGS" -Wall"])],
     [FCFLAGS=$local_fcflags" -O3"])

AX_CHECK_COMPILE_FLAG([-fPIC],
     [FCFLAGS=$FCFLAGS" -fPIC"],
     [AX_CHECK_COMPILE_FLAG([-fpic],
         [FCFLAGS=$FCFLAGS" -fpic"],
         [AC_MSG_FAILURE(
          [compiler option -fPIC/-fpic not avaiable, please use GNU or Intel Fortran compiler])]
         )])

AX_CHECK_COMPILE_FLAG([-mcmodel=$enable_mcmodel],
   [FCFLAGS=$FCFLAGS" -mcmodel="$enable_mcmodel],
   [AC_MSG_FAILURE(
       [$FC compiler option -mcmodel=$enable_mcmodel is not avaiable])])

AC_OPENMP
FCFLAGS=$FCFLAGS' '$OPENMP_FCFLAGS

FCFLAGS=$FCFLAGS' -I../include'

# Check FC
AS_IF([test "x$FC" == xnone],
     [AS_IF([test "x$enable_mpi" != xno],
        [AS_IF([test "x$enable_hdf5" != xno],
          [AC_MSG_FAILURE(
             [h5pfc not found, please install it and check environment PATH])],
          [AC_MSG_FAILURE(
            [mpif77 not found, please intall it and check environment PATH])])],
        [AS_IF([test "x$enable_hdf5" != xno],
          [AC_MSG_FAILURE(
             [h5fc not found, please install it and check environment PATH])],
          [AC_MSG_FAILURE(
            [Fortran compiler not found, please intall it and check environment PATH])])])])


# FIXME: Replace `main' with a function in `-lmpi':
# AC_CHECK_LIB([mpi], [main])


#C++ part
AC_LANG([C++])
local_cxxflags=$CXXFLAGS
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_LIB([stdc++],[main],[LIBS=$LIBS' -lstdc++'])

#AS_IF([test "x$enable_debug" != xno],
#      [CXXFLAGS=$local_cxxflags" -O2"],
#      [CXXFLAGS=$local_cxxflags" -O3"])
CXXFLAGS=$local_cxxflags" -O3"

AC_OPENMP
CXXFLAGS=$CXXFLAGS' '$OPENMP_CXXFLAGS

CXXFLAGS=$CXXFLAGS' -I../include'

AX_CHECK_COMPILE_FLAG([-march=native],
    [CXXFLAGS=$CXXFLAGS" -march=native"],[])

	 
AX_CHECK_COMPILE_FLAG([-fPIC],
     CXXFLAGS=$CXXFLAGS" -fPIC",
     [AX_CHECK_COMPILE_FLAG([-fpic],
         CXXFLAGS=$CXXFLAGS" -fpic",
         [AC_MSG_FAILURE(
          [compiler option -fPIC/-fpic not avaiable, please use GNU or Intel C++ compiler])]
         )])

AX_CHECK_COMPILE_FLAG([-mcmodel=$enable_mcmodel],
   [CXXFLAGS=$CXXFLAGS" -mcmodel="$enable_mcmodel],
   [AC_MSG_FAILURE(
       [$CXX compiler option -mcmodel=$enable_mcmodel is not avaiable])])

# SIMD part
AS_IF([test "x$OPENMP_CXXFLAGS" != x],
   [AS_IF([test "x$enable_simd" == xtest],
       [c_simd_check="avx"]
       [AC_CHECK_FILE(["/proc/cpuinfo"],
         [AC_PROG_GREP
          SIMD_CHECK=`grep -c $c_simd_check /proc/cpuinfo`],
          SIMD_CHECK=-1)
         AC_MSG_CHECKING([whether current CPU support $c_simd_check])
         AS_IF([test "$SIMD_CHECK" -lt 0],
             [AC_MSG_RESULT([unknown])],
             [test "$SIMD_CHECK" -gt 0],
             [AC_MSG_RESULT([yes])],
             [AC_MSG_RESULT([no])])
         AX_CHECK_COMPILE_FLAG([-m$c_simd_check],
           [CXXFLAGS=$CXXFLAGS" -m"$c_simd_check
            RESULT=$RESULT"."$c_simd_check
	    enable_simd=$c_simd_check
           ],
           [c_simd_check="sse"]
           [AC_CHECK_FILE(["/proc/cpuinfo"],
              [AC_PROG_GREP
               SIMD_CHECK=`grep -c $c_simd_check /proc/cpuinfo`],
               SIMD_CHECK=-1)
            AC_MSG_CHECKING([whether current CPU support $c_simd_check])
            AS_IF([test "$SIMD_CHECK" -lt 0],
                  [AC_MSG_RESULT([unknown])],
                  [test "$SIMD_CHECK" -gt 0],
                  [AC_MSG_RESULT([yes])],
                  [AC_MSG_RESULT([no])])
            AX_CHECK_COMPILE_FLAG([-m$c_simd_check],
              [CXXFLAGS=$CXXFLAGS" -m"$c_simd_check
               RESULT=$RESULT"."$c_simd_check
	       enable_simd=$c_simd_check
              ],
              [enable_simd=no])])],
        [AS_IF([test "x$enable_simd" != xno],
               [AC_CHECK_FILE(["/proc/cpuinfo"],
                  [AC_PROG_GREP
                   SIMD_CHECK=`grep -c $enable_simd /proc/cpuinfo`],
                   SIMD_CHECK=-1)
                AC_MSG_CHECKING([whether current CPU support $enable_simd])
                AS_IF([test "$SIMD_CHECK" -lt 0],
                      [AC_MSG_RESULT([unknown])],
                      [test "$SIMD_CHECK" -gt 0],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])])
                AX_CHECK_COMPILE_FLAG([-m$enable_simd],
                  [CXXFLAGS=$CXXFLAGS" -m"$enable_simd
                   RESULT=$RESULT"."$enable_simd
                  ],
                  [AC_MSG_FAILURE(
                     [$CXX compiler option -m$enable_simd is not avaiable])])])])],
     enable_simd=no
     SIMD_CHECK=-2)

# FIXME: Replace `main' with a function in `-llvisit':
# AC_CHECK_LIB([lvisit], [main])
# FIXME: Replace `main' with a function in `-llvisit_nbody':
# AC_CHECK_LIB([lvisit_nbody], [main])
# FIXME: Replace `main' with a function in `-lvisit':
# AC_CHECK_LIB([visit], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memset sqrt strchr])

#CUDA part
# Checks for MPI support
AC_ARG_WITH([cuda], 
     [AS_HELP_STRING([--with-cuda=PREFIX],
         [Prefix of your CUDA installation])],
     [cuda_prefix=$withval],
     [cuda_prefix="/usr/local/cuda"])

AC_ARG_WITH([cuda_sdk],
     [AS_HELP_STRING([--with-cuda-sdk=PREFIX],
         [Prefix of your CUDA SDKinstallation])],
     [cuda_sdk_prefix=$withval],
     [cuda_sdk_prefix="/usr/local/cuda/samples"])

AS_IF([test "x$enable_gpu" != xno],
     [AC_CHECK_PROG([NVCC], [nvcc], [nvcc], [none])
      AS_IF([test "x$NVCC" == xnone], 
            [AC_CHECK_PROG([NVCC], [nvcc], [$cuda_prefix/bin/nvcc], [none], [$cuda_prefix/bin])
             AS_IF([test "x$NVCC" == xnone], 
                   [AC_MSG_FAILURE(
                       [Cannot find CUDA compiler nvcc, please check whether nvcc is in environment PATH or use --with-cuda=PREFIX to provide the PATH of CUDA installation])])],
            [AC_PROG_SED
             cuda_prefix=`which nvcc|sed 's:/bin/nvcc::'`])
      AC_CHECK_LIB([cudart], [main], 
          [LIBS=$LIBS' -lcudart'],
          [AC_CHECK_FILE(["$cuda_prefix/lib64/libcudart.so"],
              LIBS=$LIBS" -L$cuda_prefix/lib64 -lcudart",
              [AC_MSG_FAILURE(
                  [Cannot find CUDA library -lcudart, please provide correct cuda PREFIX by using --with-cuda=PREFIX])])])
      CUFLAGS=" -O3"
      RESULT=$RESULT".gpu"
      AC_PROG_AWK
      AC_MSG_CHECKING([for CUDA version])
      CUDA_VERSION=`$NVCC --version|awk -F ',' 'END{print $(NF-1)}'|awk '{print $2}'`
      CUDA_VERSION_I=`echo $CUDA_VERSION|awk -F '.' '{print $1}'`
      AC_MSG_RESULT([$CUDA_VERSION])
      AS_IF([test "$CUDA_VERSION_I	" -gt 4],
         CUFLAGS=$CUFLAGS" -D CUDA_5"
         [AC_CHECK_FILE(["$cuda_prefix/samples/common/inc/helper_cuda.h"],
              [CUFLAGS=$CUFLAGS" -I $cuda_prefix/samples/common/inc"],
              [AC_CHECK_FILE(["$cuda_sdk_prefix/common/inc/helper_cuda.h"],
                    [CUFLAGS=$CUFLAGS" -I $cuda_sdk_prefix/common/inc"],
                    [CUFLAGS=$CUFLAGS" -I ../extra_inc/cuda"])])],
         [AC_CHECK_FILE(["$cuda_prefix/cuda_sdk/common/inc/cutil.h"],
              [CUFLAGS=$CUFLAGS" -I $cuda_prefix/cuda_sdk/common/inc"],
              [AC_CHECK_FILE(["$cuda_sdk_prefix/common/inc/cutil.h"],
                    [CUFLAGS=$CUFLAGS" -I $cuda_prefix/common/inc"],
                    [CUFLAGS=$CUFLAGS" -I ../extra_inc/cuda"])])])
      FCFLAGS=$FCFLAGS' $(GPU_FLAGS)'
      EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
      AS_IF([test "x$OPENMP_CXXFLAGS" != x ],
           EXTRAOBJ=$EXTRAOBJ' $(CUDA_OBJECTS)', 
           EXTRAOBJ=$EXTRAOBJ' $(CUDA_OBJECTS_SINGLE)')
      LDFLAGS=$LDFLAGS' -L'$CUDA_PREFIX
      dnl AX_CHECK_COMPILE_FLAG([-Xcompiler],[],
      dnl    [AC_MSG_FAILURE(
      dnl       [Current nvcc do not support -Xcompiler])])
     ])

#FLAGS
# MPI FLAGS
AS_IF([test "x$enable_mpi" != xno],
      FCFLAGS=$FCFLAGS' $(MPI_FLAGS)'
      EXTRASRC=$EXTRASRC' $(MPI_FSOURCES)'
      RESULT=$RESULT'.mpi')

#SIMD
AS_IF([test "x$enable_simd" != xno],
      [FCFLAGS=$FCFLAGS' ${SIMD_FLAGS}']
      [CXXCOMP=`echo x$CXX|grep -c icc`]
      [AS_IF([test $CXXCOMP -gt 0],
         CXXFLAGS=$CXXFLAGS' $(INTEL_FLAGS)')])

#GPU
AS_IF([test "x$enable_simd" == xavx],
      [AS_IF([test "x$enable_gpu" != xno],
             EXTRAOBJ=$EXTRAOBJ' ${IRRAVX_OBJECTS}',
             EXTRAOBJ=$EXTRAOBJ' ${AVX_OBJECTS}'
             EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
             FCFLAGS=$FCFLAGS' $(GPU_FLAGS)')],
      [test "x$enable_simd" == xsse],
      [AS_IF([test "x$enable_gpu" != xno],
             EXTRAOBJ=$EXTRAOBJ' ${IRRSSE_OBJECTS}',
             EXTRAOBJ=$EXTRAOBJ' ${SSE_OBJECTS}'
             EXTRASRC=$EXTRASRC' $(GPU_FSOURCES)'
             FCFLAGS=$FCFLAGS' $(GPU_FLAGS)')])

#OMP
AS_IF([test "x$OPENMP_FCFLAGS" != x ],
      FCFLAGS=$FCFLAGS' ${OMP_FLAGS}',
      FCFLAGS=$FCFLAGS' -I../extra_inc/noomp')

AS_IF([test "x$OPENMP_CXXFLAGS" != x ],
      CXXFLAGS=$CXXFLAGS' ${OMP_FLAGS}',
      CXXFLAGS=$CXXFLAGS' -I ../extra_inc/noomp')

# HDF5 format
AS_IF([test "x$enable_hdf5" != xno],
       [FCFLAGS=$FCFLAGS' $(HDF5_FLAGS)'
        EXTRASRC=$EXTRASRC' $(HDF5_FSOURCES)'
        RESULT=$RESULT'.hdf5'])

# TT support
AS_IF([test "x$enable_tt" != xno],
       [FCFLAGS=$FCFLAGS' $(TT_FLAGS)'
        EXTRASRC=$EXTRASRC' $(TT_FSOURCES)'
	RESULT=$RESULT'.tt'])

# DEBUG FLAGS
AS_IF([test "x$enable_debug" != xno],
      FCFLAGS=$FCFLAGS' $(DEBUG_FLAGS)'
      RESULT=$RESULT'.db')


# params.h
AC_ARG_WITH([par], 
     [AS_HELP_STRING([--with-par=size],
         [Size of basic parameters (NMAX, KMAX, LMAX, MMAX) in params.h. Avaiable values: 1m, 512k, 256k, 128k, 64k, 32k, 16k, 8k, custom. The value of this size is NMAX, with prefix 'b' the K.S. binary size KMAX is half of NMAX and NMAX is set to 1.5*NMAX. LMAX is 600 above 64k, 500 above 16k. MMAX is 1024 above 64k, 512 above 16k. MMAX is doubled with prefix 'b'. Defaultsize is b16k])],
     [par_size=$withval],
     [par_size=b16k])

AS_IF([test "x$par_size" == x1m],
      [NMAX=1048576
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb1m],
      [NMAX=1510720
       KMAX=512000
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x512k],
      [NMAX=524288
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb512k],
      [NMAX=786432
       KMAX=262144
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x256k],
      [NMAX=262144
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb256k],
      [NMAX=393216
       KMAX=131072
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x128k],
      [NMAX=131072
       KMAX=65536
       LMAX=600
       MMAX=1024],
      [test "x$par_size" == xb128k],
      [NMAX=196608
       KMAX=65536
       LMAX=600
       MMAX=2048],
      [test "x$par_size" == x64k],
      [NMAX=65536
       KMAX=8192
       LMAX=500
       MMAX=512],
      [test "x$par_size" == xb64k],
      [NMAX=98304
       KMAX=32768
       LMAX=500
       MMAX=1024],
      [test "x$par_size" == x32k],
      [NMAX=32768
       KMAX=8192
       LMAX=500
       MMAX=512],
      [test "x$par_size" == xb32k],
      [NMAX=49152
       KMAX=16384
       LMAX=500
       MMAX=1024],
      [test "x$par_size" == x16k],
      [NMAX=16384
       KMAX=2048
       LMAX=400
       MMAX=256],
      [test "x$par_size" == xb16k],
      [NMAX=24576
       KMAX=8192
       LMAX=400
       MMAX=512],
      [test "x$par_size" == x8k],
      [NMAX=8192
       KMAX=2048
       LMAX=400
       MMAX=256],
      [test "x$par_size" == xb8k],
      [NMAX=12288
       KMAX=4096
       LMAX=400
       MMAX=512],
      [test "x$par_size" == xcustom],
      [flag_input_par=yes],
      [NMAX=24576
       KMAX=8192
       LMAX=400
       MMAX=512])

AC_ARG_WITH([nmax], 
     [AS_HELP_STRING([--with-nmax=size],
         [Set Maximum number of particles (will overwrite value set in --with-par)])],
     [NMAX=$withval],
     [flag_par_miss=yes])

AC_ARG_WITH([kmax], 
     [AS_HELP_STRING([--with-kmax=size],
         [Set Maximum number of K.S. particles (will overwrite value set in --with-par)])],
     [KMAX=$withval],
     [flag_par_miss=yes])

AC_ARG_WITH([lmax], 
     [AS_HELP_STRING([--with-lmax=size],
         [Set Maximum number of neighbor list particles (will overwrite value set in --with-par)])],
     [LMAX=$withval],
     [flag_par_miss=yes])

AC_ARG_WITH([mmax], 
     [AS_HELP_STRING([--with-mmax=size],
         [Set Maximum number of merger (stable triple) particles (will overwrite value set in --with-par)])],
     [MMAX=$withval],
     [flag_par_miss=yes])

AS_IF([test "x$flag_input_par" == xyes],
  [AS_IF([test "x$flag_par_miss" == xyes],
        [AC_MSG_FAILURE(
             [Parameters are not completely set, please make sure that all four option (--with-nmax, --with-kmax, --with-lmax, --with-mmax) are set])])])
      

# Extral tools
AS_IF([test "x$enable_tools" != xno],
   [EXTRARESULT=' $(EXTRATOOLS)'],
   [EXTRARESULT=''])

AC_SUBST([CUFLAGS])	
AC_SUBST([FCFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([RESULT])
AC_SUBST([EXTRASRC])
AC_SUBST([EXTRAOBJ])
AC_SUBST([EXTRARESULT])

AC_SUBST([NMAX])
AC_SUBST([KMAX])
AC_SUBST([LMAX])
AC_SUBST([MMAX])

AC_MSG_NOTICE([----------------------Configure Summary--------------------])
AC_MSG_NOTICE([--Features:])
AC_MSG_NOTICE([     Using MPI:      $enable_mpi])
AC_MSG_NOTICE([     Using GPU:      $enable_gpu])
AC_MSG_NOTICE([     CUDA version:   $CUDA_VERSION])
AC_MSG_NOTICE([     Using SIMD:     $enable_simd])
AS_IF([test "x$enable_simd" != xno],
  [AS_IF([test "$SIMD_CHECK" == -1],
      [AC_MSG_WARN([      No /proc/cpuinfo found, cannot comfirm whether current CPU support SIMD])],
        [test "$SIMD_CHECK" == 0],
      [AC_MSG_WARN([      SIMD option used now is not supported in current CPU])])],
        [test "$SIMD_CHECK" == -2],
      [AC_MSG_WARN([      SIMD is suppressed since OpenMP is switched off])])
AC_MSG_NOTICE([     mcmodel:        $enable_mcmodel])
AC_MSG_NOTICE([     OpenMP option:  $OPENMP_FCFLAGS])
AC_MSG_NOTICE([     Debug mode:     $enable_debug])
AC_MSG_NOTICE([     HDF5 format:    $enable_hdf5])
AC_MSG_NOTICE([     TT package:     $enable_tt])
AC_MSG_NOTICE([--Compilers:])
AC_MSG_NOTICE([     C++ compiler:      $CXX])
AC_MSG_NOTICE([     Fortran compiler:  $FC])
AC_MSG_NOTICE([     CUDA compiler:     $NVCC])
AC_MSG_NOTICE([--Targets:])
AC_MSG_NOTICE([     Program name:   nbody6++$RESULT])
AC_MSG_NOTICE([     Extra tools:    $enable_tools])
AC_MSG_NOTICE([     Install path:   $prefix])
AC_MSG_NOTICE([--Parameters:])
AC_MSG_NOTICE([     NMAX:           $NMAX])
AC_MSG_NOTICE([     KMAX:           $KMAX])
AC_MSG_NOTICE([     LMAX:           $LMAX])
AC_MSG_NOTICE([     MMAX:           $MMAX])
AC_MSG_NOTICE([-----------------------------------------------------------])
AS_IF([test "x$enable_tools" != xno],
   [AC_MSG_NOTICE([  Please add $prefix/lib into your runtime load library path:])
    AC_MSG_NOTICE([     LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$prefix/lib])
    AC_MSG_NOTICE([  in your pre-configure script of shell (such as ~/.bashrc if you use bash)])
    AC_MSG_NOTICE([-----------------------------------------------------------])])
AC_MSG_NOTICE([Reference paper: Wang, L., et al., 2015, MNRAS, 450, 4070 (http://arxiv.org/abs/1504.03687)])
AC_MSG_NOTICE([Email: Long Wang<longwang.astro@gmail.com>])
AC_MSG_NOTICE([-----------------------------------------------------------])

AC_CONFIG_FILES([Makefile build/Makefile include/params.h])
AC_OUTPUT
