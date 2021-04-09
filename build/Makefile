FORTRAN=mpif90
F90=ifort

FLAGS=-i4 -r8 -O2 -assume byterecl

## UNCOMMENT TO RUN IN DEBUG MODE
# DEBUG_FLAGS=-g -traceback -check bounds

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff

SRC =   fft.f \
        modules.F \
        netcdf_io.f90 \
        subroutines.F \
        particleSubs.F \
        ioSubs.F

OBJS = $(addsuffix .o, $(basename $(SRC)))


lesmpi.a: $(OBJS) ../src/les.F
	$(FORTRAN) $^ -o $@  $(FLAGS) $(DEBUG_FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS)

%.o: ../src/%.f
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: ../src/%.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: ../src/%.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

