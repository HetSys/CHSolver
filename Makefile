# compiler and linker
FC= gfortran # can use h5fc #HDF5's compiler (that wraps gfortran) then no need for most the flags/libs
LD=$(FC)
# flags and libraries
FFLAGS= -I/usr/include/hdf5/serial -Wall -Wextra -std=f2008 #h5fc-show is equiv to nf-config --fflags/flibs to find these
FLIBS= -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_fortran.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial


# executable names
EXE=chsolver2
EXE_TEST=test

# directories
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./out

FACEPATH = $(SRC_DIR)/submodules/FACE/src/lib/face.F90
FLOGPATH = $(SRC_DIR)/submodules/flogging/src/logging.f90

# differentiate between test and chSolver builds
SRC= $(SRC_DIR)/*.f90
# SRC_TEST=$(SRC_SUB) $(SRC_DIR)/$(EXE_TEST).f90


OBJ=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.f90=.o)))
# OBJ_TEST=$(addprefix $(OBJ_DIR)/, $(notdir $(SRC_TEST:.f90=.o)))

chsolver: directories $(OBJ)
	@printf "`tput bold``tput setaf 2`Linking`tput sgr0`\n"
	$(LD) $(FFLAGS) $(OBJ) $(FLIBS) -o $(EXE)
	
# Build rule for binaries (puts .mod files in SRC_DIR to simplify linting)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(FC) -J$(SRC_DIR) $(FFLAGS) -c -o $@ $< 
$(OBJ_DIR)/face.o: $(FACEPATH)
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(FC)-std=f2008 -c -o $@ $< 
$(OBJ_DIR)/logger_mod.o: $(FLOGPATH)
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(FC)-std=f2008 -c -o $@ $< 

#$(FLIBS) not needed for linking?

# # build the test file
# test: directories $(OBJ_TEST)
# 	@printf "`tput bold``tput setaf 2`Linking`tput sgr0`\n"
# 	$(LD) -o $(EXE_TEST) $(OBJ_TEST)

# create required directories
.PHONY: directories
directories:
	mkdir -p $(OBJ_DIR) $(OUT_DIR)

# removes binaries, outputs etc.
.PHONY: clean
clean:
	rm -f $(EXE) $(EXE_TEST) $(OBJ_DIR)/*.o $(SRC_DIR)/*.mod $(OUT_DIR)/**

# force rebuild of all files
.PHONY: all
all: clean chsolver

# dependencies
$(OBJ_DIR)/logger_mod.o : $(OBJ_DIR)/face.o
$(OBJ_DIR)/globals.o : $(OBJ_DIR)/logger_mod.o
$(OBJ_DIR)/main.o : $(OBJ_DIR)/json-parser.o $(OBJ_DIR)/solvers.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/progress.o $(OBJ_DIR)/hdf5-io.o $(OBJ_DIR)/setup.o $(OBJ_DIR)/globals.o
# $(OBJ_DIR)/write_netcdf.o : $(OBJ_DIR)/model_data.o
# $(OBJ_DIR)/velocity_verlet.o : $(OBJ_DIR)/model_data.o
# $(OBJ_DIR)/model_data.o : $(OBJ_DIR)/create_axis.o $(OBJ_DIR)/command_line.o
