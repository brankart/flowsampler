# Makefile created by mkmf 19.3.0


MK_TEMPLATE = Makefile.template
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: ../lib/libflowsampler.a
GeneralProfile.o: ../src/BoundaryLayer/GeneralProfile.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/BoundaryLayer/GeneralProfile.F90
InternalProfile.o: ../src/BoundaryLayer/InternalProfile.F90 GeneralProfile.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/BoundaryLayer/InternalProfile.F90
ProfileSampler.o: ../src/BoundaryLayer/ProfileSampler.F90 GeneralProfile.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/BoundaryLayer/ProfileSampler.F90
dynamics.o: ../src/ShallowWater/dynamics.F90 grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/ShallowWater/dynamics.F90
grid.o: ../src/ShallowWater/grid.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/ShallowWater/grid.F90
kinematics.o: ../src/ShallowWater/kinematics.F90 grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../src/ShallowWater/kinematics.F90
./dynamics.F90: ../src/ShallowWater/dynamics.F90
	cp ../src/ShallowWater/dynamics.F90 .
./GeneralProfile.F90: ../src/BoundaryLayer/GeneralProfile.F90
	cp ../src/BoundaryLayer/GeneralProfile.F90 .
./kinematics.F90: ../src/ShallowWater/kinematics.F90
	cp ../src/ShallowWater/kinematics.F90 .
./grid.F90: ../src/ShallowWater/grid.F90
	cp ../src/ShallowWater/grid.F90 .
./ProfileSampler.F90: ../src/BoundaryLayer/ProfileSampler.F90
	cp ../src/BoundaryLayer/ProfileSampler.F90 .
./InternalProfile.F90: ../src/BoundaryLayer/InternalProfile.F90
	cp ../src/BoundaryLayer/InternalProfile.F90 .
SRC = ../src/ShallowWater/grid.F90 ../src/BoundaryLayer/InternalProfile.F90 ../src/ShallowWater/kinematics.F90 ../src/BoundaryLayer/ProfileSampler.F90 ../src/ShallowWater/dynamics.F90 ../src/BoundaryLayer/GeneralProfile.F90
OBJ = grid.o InternalProfile.o kinematics.o ProfileSampler.o dynamics.o GeneralProfile.o
OFF = ../src/ShallowWater/dynamics.F90 ../src/BoundaryLayer/GeneralProfile.F90 ../src/ShallowWater/kinematics.F90 ../src/ShallowWater/grid.F90 ../src/BoundaryLayer/ProfileSampler.F90 ../src/BoundaryLayer/InternalProfile.F90
clean: neat
	-rm -f .../lib/libflowsampler.a.cppdefs $(OBJ) ../lib/libflowsampler.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

../lib/libflowsampler.a: $(OBJ)
	$(AR) $(ARFLAGS) ../lib/libflowsampler.a $(OBJ)

