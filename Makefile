# compiler and linker
export HDF5_FC := mpif90
export HDF5_FLINKER := mpif90
FC=h5pfc # use h5fc's alias to get correct flags
LD=$(FC)

# flags and libraries
FFLAGS= -Wall -Wextra -Wconversion-extra -std=f2008 -g -fall-intrinsics
FFLAGS+= -I./src/submodules/bin/jsonfortran-gnu-8.2.5/lib
FLIBS= -lpthread -lsz -lz -ldl -lm ./src/submodules/bin/jsonfortran-gnu-8.2.5/lib/libjsonfortran.a

# executable names
EXE=chsolver

# directories
SRC_DIR=./src
OBJ_DIR=./obj
OUT_DIR=./out
TEST_DIR=./tests/unit_tests
LOG_DIR=./logs


FACEPATH = $(SRC_DIR)/submodules/FACE/src/lib/face.F90
FLOGPATH = $(SRC_DIR)/submodules/flogging/src/logging.f90

SRC=$(filter-out $(SRC_DIR)/test.f90, $(wildcard $(SRC_DIR)/*.f90))
SRC_TEST=$(filter-out $(SRC_DIR)/main.f90, $(wildcard $(SRC_DIR)/*.f90))

OBJ= $(OBJ_DIR)/face.o $(OBJ_DIR)/logger_mod.o  $(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.f90=.o)))
OBJ_TEST=$(OBJ_DIR)/face.o $(OBJ_DIR)/logger_mod.o  $(addprefix $(OBJ_DIR)/, $(notdir $(SRC_TEST:.f90=.o)))

chsolver: directories $(OBJ)
	@printf "`tput bold``tput setaf 2`Linking`tput sgr0`\n"
	@echo $(OBJ)
	$(LD) $(FFLAGS) -o $(EXE) $(OBJ) $(FLIBS) 
	
# Build rule for binaries (puts .mod files in SRC_DIR to simplify linting)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(FC) -J$(SRC_DIR) $(FFLAGS) -c -o $@ $< 

$(OBJ_DIR)/face.o: $(FACEPATH)
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(FC) -J$(SRC_DIR) -std=f2008 -c -o $@ $<

$(OBJ_DIR)/logger_mod.o: $(FLOGPATH)
	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
	$(FC) -J$(SRC_DIR) -std=f2008 -c -o $@ $<

# test_json_parser.o: test_json_parser.mod
# %.o: %.F90
# 	@printf "`tput bold``tput setaf 6`Building %s`tput sgr0`\n" $@
# 	$(FC) -J$(SRC_DIR) $(FFLAGS) -c -o $@ $<
#$(FLIBS) not needed for linking?

$(SRC_DIR)/libchsolver.a : $(wildcard $(OBJ_DIR)/*.o)
	ar -r $@ $?


.PHONY: test
test: directories $(OBJ_TEST)
	@printf "`tput bold``tput setaf 2`Linking`tput sgr0`\n"
	@echo $(OBJ_TEST)
	$(LD) $(FFLAGS) -o test $(OBJ_TEST) $(FLIBS)

# pFUnit
# .PHONY: tests
tests: all $(SRC_DIR)/libchsolver.a
	make -C tests/unit_tests all

# create required directories
.PHONY: directories
directories:
	mkdir -p $(OBJ_DIR) $(OUT_DIR) $(LOG_DIR)

# removes binaries, outputs etc.
.PHONY: clean
clean:
	rm -f -r -d $(EXE) $(EXE_TEST) $(OBJ_DIR)/*.o $(SRC_DIR)/*.mod $(OUT_DIR)/** ./doxygen/output/* test_all
	make -C tests/unit_tests clean

# removes logfiles
.PHONY: logpurge
logpurge:
	rm -f -r -d $(LOG_DIR)/**

# force rebuild of all files
.PHONY: all
all: clean chsolver

# Generate Doxygen documentation
.PHONY: docs
docs:
	doxygen ./doxygen/doxygen-config 
	(cd ./doxygen/output/latex && make)
	cp ./doxygen/output/latex/refman.pdf CH-docs.pdf

# dependencies
$(OBJ_DIR)/logger_mod.o : $(OBJ_DIR)/face.o
$(OBJ_DIR)/comms.o : $(OBJ_DIR)/solver-utils.o
$(OBJ_DIR)/globals.o : $(OBJ_DIR)/logging.o
$(OBJ_DIR)/solver-utils.o : $(OBJ_DIR)/globals.o
$(OBJ_DIR)/solvers.o : $(OBJ_DIR)/globals.o $(OBJ_DIR)/solver-utils.o $(OBJ_DIR)/fd-solvers.o
$(OBJ_DIR)/fd-solvers.o : $(OBJ_DIR)/globals.o $(OBJ_DIR)/solver-utils.o $(OBJ_DIR)/hdf5-io.o
$(OBJ_DIR)/logging.o : $(OBJ_DIR)/logger_mod.o
$(OBJ_DIR)/json-parser.o $(OBJ_DIR)/hdf5-io.o $(OBJ_DIR)/setup.o $(OBJ_DIR)/solvers.o $(OBJ_DIR)/progress.o : $(OBJ_DIR)/globals.o 
$(OBJ_DIR)/command-line.o : $(OBJ_DIR)/globals.o
$(OBJ_DIR)/main.o : $(OBJ_DIR)/json-parser.o $(OBJ_DIR)/solvers.o $(OBJ_DIR)/globals.o $(OBJ_DIR)/progress.o $(OBJ_DIR)/hdf5-io.o $(OBJ_DIR)/setup.o $(OBJ_DIR)/globals.o

