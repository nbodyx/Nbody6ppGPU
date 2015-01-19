# AX_FC_CHECK_BOUNDS([ACTION-IF-SUCCESS], [ACTION-IF-FAILURE = FAILURE])
# ----------------------------------------------------------------------
# Look for a compiler flag to turn on array bounds checking for the
# Fortran (FC) compiler, and adds it to FCFLAGS.  Call
# ACTION-IF-SUCCESS (defaults to nothing) if successful (i.e. can
# compile code using new extension) and ACTION-IF-FAILURE (defaults to
# failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
#
# The known flags are:
# -fcheck=all, -fbounds-check: gfortran
#     -fbounds-check: g77, g95
# -CB, -check bounds: Intel compiler (icc, ecc, ifort)
#                 -C: Sun/Oracle compiler (f95)
#        -C, -qcheck: IBM compiler (xlf)
#           -Mbounds: Portland Group compiler
#       -C ,-Mbounds: Cray
#  -C, -check_bounds: SGI compiler
# -check_bounds, +check=all: HP Fortran
#        -C, -Rb -Rc: Absoft (-Rb: array boundaries, -Rc: array conformance)
# --chk e,s -chk (e,s): Lahey
#          -C -C=all: NAGWare
# -C, -ffortran-bounds-check: PathScale pathf90
#                 -C: f2c
#            -BOunds: Open Watcom
AC_DEFUN([AX_FC_CHECK_BOUNDS],
[AC_LANG_PUSH([Fortran])dnl
AC_CACHE_CHECK([for Fortran flag to enable array-bounds checking],
               [ac_cv_fc_check_bounds],
[ac_cv_fc_check_bounds=unknown
ac_fc_check_bounds_FCFLAGS_save=$FCFLAGS
for ac_flag in -fcheck=bounds -fbounds-check -check_bounds -Mbounds -qcheck \
               '-check bounds' +check=all --check '-Rb -Rc' -CB -C=all -C \
               -ffortran-bounds-check "--chk e,s" "-chk e -chk s" -bounds
do
  FCFLAGS="$ac_fc_check_bounds_FCFLAGS_save $ac_flag"
  # We should be able to link a correct program.
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
    [AC_LINK_IFELSE([[
      subroutine sub(a)
      integer a(:)
      a(8) = 0
      end subroutine

      program main
      integer a(1:7)
      interface
         subroutine sub(a)
         integer a(:)
         end subroutine
      end interface

      call sub(a)
      end program]],
       [# If we can run the program, require failure at run time.
       # In cross-compiling mode, we rely on the compiler not accepting
       # unknown options.
       AS_IF([test "$cross_compiling" = yes],
         [ac_cv_fc_check_bounds=$ac_flag; break],
	   [AS_IF([_AC_DO_TOKENS(./conftest$ac_exeext)],
	        [],
		     [ac_cv_fc_check_bounds=$ac_flag; break])])])])
done
rm -f conftest$ac_exeext conftest.err conftest.$ac_objext conftest.$ac_ext \
  core *.core core.conftest.*
FCFLAGS=$ac_fc_check_bounds_FCFLAGS_save
])
if test "x$ac_cv_fc_check_bounds" = xunknown; then
  m4_default([$2],
             [AC_MSG_ERROR([no Fortran flag for bounds checking found], 77)])
else
  if test "x$ac_cv_fc_check_bounds" != xnone; then
    FCFLAGS="$FCFLAGS $ac_cv_fc_check_bounds"
  fi
  $1
fi
AC_LANG_POP([Fortran])dnl
])# AC_FC_CHECK_BOUNDS

