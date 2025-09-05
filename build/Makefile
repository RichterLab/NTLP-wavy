FORTRAN=mpiifort

FLAGS=-i4 -r8 -O2 -assume byterecl -march=core-avx2 -fpp -traceback -fc=ifx

## UNCOMMENT TO RUN IN DEBUG MODE
#DEBUG_FLAGS=-g -traceback -check bounds

OUTPUTINC = -I$(NETCDFBASE)/include
OUTPUTLIB = -L$(NETCDFBASE)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  = -lnetcdf -lnetcdff

SRC =   fft.f \
        modules.F \
        netcdf_io.f90 \
        subroutines.F \
        particleSubs.F \
	boundarySubs.F \
        ioSubs.F \
	les.F

OBJS = $(addsuffix .o, $(basename $(SRC)))


lesmpi.a: $(OBJS) 
	$(FORTRAN) $^ -o $@  $(FLAGS) $(DEBUG_FLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS)

%.o: ../src/%.f
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: ../src/%.f90
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)

%.o: ../src/%.F
	$(FORTRAN) $(FLAGS) $(DEBUG_FLAGS) -c $< $(OUTPUTINC) $(OUTPUTLIB)


clean:
	rm -f *.o *.mod lesmpi.a mach.file

#Dependencies:
les.o: modules.o particleSubs.o boundarySubs.o netcdf_io.o subroutines.o
particleSubs.o: modules.o
boundarySubs.o: modules.o
netcdf_io.o: modules.o
subroutines.o: modules.o
ioSubs.o: modules.o
