#!/bin/ksh
set -x

# --------------------------------------
# netcdf library paths and compiler flag
# --------------------------------------
LIBNC="${NETCDF_C_ROOT}/lib -lnetcdf -L${NETCDF_FORTRAN_ROOT}/lib -lnetcdff"
INCNC="${NETCDF_C_ROOT}/include -I${NETCDF_FORTRAN_ROOT}/include"
FCFLAG="-fdefault-real-8 -O3 -ffloat-store"

# ---------------------------
# End of user-defined section
# ---------------------------
LIB=libannif.a
for f in  annif.F90  prog_testnn.F90; do
	mpif90 ${FCFLAG} -c -I$INCNC $f || exit 9
done

mpif90 ${FCFLAG} -L$LIBNC -o test_annif.x annif.o prog_testnn.o ADFirstAidKit/adStack_gcc.o

rm -f $LIB
ar -rv $LIB ADFirstAidKit/adStack_gcc.o
ar -rv $LIB annif.o
