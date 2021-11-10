

CMD = Apinv
FC = gfortran
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
           -W -m64 -mcmodel=medium
FSRCS = kerf.f kesy.f subfunc.f
F90SRCS = lsmrDataModule.f90 lsmrblasInterface.f90\
          lsmrblas.f90 lsmrModule.f90 aprod.f90 Apinv.f90
OBJS = $(FSRCS:%.f=%.o) $(F90SRCS:%.f90=%.o) Apkernel.o
all:$(CMD)
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ -o $@

Apkernel.o:Apkernel.f90
	$(FC) -fopenmp $(FFLAGS) -c $<  -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
clean:
	rm *.o $(CMD)
