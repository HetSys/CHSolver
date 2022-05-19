PFUNIT_DIR=../../src/submodules/pFUnit/build/installed/PFUNIT-4.3
include $(PFUNIT_DIR)/include/PFUNIT.mk
FC= gfortran # can use h5fc #HDF5's compiler (that wraps gfortran) then no need for most the flags/libs
LD=$(FC)

FFLAGS= -Wall -Wextra -Wconversion-extra -std=f2008 -Wno-maybe-uninitialized

FFLAGS += -I/usr/include/hdf5/serial -I../../src/submodules/bin/jsonfortran-gnu-8.2.5/lib

FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

FLIBS = ../../src/submodules/bin/jsonfortran-gnu-8.2.5/lib/libjsonfortran.a
FLIBS += $(filter-out ../../obj/main.o ../../obj/hdf5-io.o ../../obj/solvers.o ../../obj/fd-solvers.o, $(wildcard ../../obj/*.o))
# FLIBS += $(filter-out ../../obj/main.o, $(wildcard ../../obj/*.o))

# directories
SRC_DIR=../../src
OBJ_DIR=../../obj
OUT_DIR=../../out
TEST_DIR=.
LOG_DIR=../../logs

FFLAGS += -I$(SRC_DIR) -I$(OBJ_DIR)

FACEPATH = $(SRC_DIR)/submodules/FACE/src/lib/face.F90
FLOGPATH = $(SRC_DIR)/submodules/flogging/src/logging.f90


%.o : %.F90
	$(FC) $(FFLAGS) -c -o $@ $< 

all: test_all
	$(shell mv test_all ../..)

test_all_TESTS := test_json_parser.pf test_globals.pf test_cli.pf test_setup.pf test_solver_utils.pf
test_all_OTHER_LIBRARIES := $(FLIBS)
$(eval $(call make_pfunit_test,test_all))

clean:
	$(RM) *.o *.mod *.a  *.inc *.F90


