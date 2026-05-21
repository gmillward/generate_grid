
FORTRAN=ifort

FFLAGS= -132

A_OUT=apex_prog

OBJECTS=\
	generate_apex_coordinates.o \
	divve.o \
	apxntrpb4lf.o \
	apex.o \
	ggrid.o \
	magfld.o \
        calc_apex_params_2d_2.o

$(A_OUT): $(OBJECTS)
	$(FORTRAN) -o $(A_OUT) $(OBJECTS) $(LFLAGS) 

.f.o:
	$(FORTRAN) -c $(FFLAGS) $<

calc_apex_params_2d_2.o: calc_apex_params_2d_2.f90
	$(FORTRAN) -c $(FFLAGS) calc_apex_params_2d_2.f90
