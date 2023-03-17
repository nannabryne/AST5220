#	Makefile for project in AST5220 
#	-------------------------------


SHELL := /bin/bash
# HOME := /Users/nannabryne# /uio/hume/student-u29/nannabr

BUILD_DIR := build
# VPATH=src/
TARGETS := cmb

# Set compiler and flags
CXX = g++ -std=c++11 
CXXFLAGS = -O3 -g $(OPTIONS)

# Paths to GSL library
INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas


#	-	-	-	-	-	-	-	-	-	-	-	-	
# 	Options:
#	-	-	-	-	-	-	-	-	-	-	-	-
OPTIONS = 
OPTIONS += -D_GLIBCXX_DEBUG						# 	add bounds checking
OPTIONS += -D_SPLINE_WARNINGS_ON				#	show warnings if atempting to evaluate a spline out of bounds
# OPTIONS += -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE	# 	show info about the solution as we integrate
# OPTIONS += -D_USEOPEMP							#	add OpenMP parallelisation


# CXX += -fopenmp


#	-	-	-	-	-	-	-	-	-	-	-	-
# 	Object files and dependencies:
#	-	-	-	-	-	-	-	-	-	-	-	-

OBJS_ = main.o utils.o backgroundcosmology.o recombinationhistory.o spline.o ODESolver.o
OBJS = $(addprefix $(BUILD_DIR)/, $(OBJS_))

main.o					: backgroundcosmology.h
spline.o                : spline.h
ODESolver.o             : ODESolver.h
utils.o                 : utils.h spline.h ODESolver.h
backgroundcosmology.o   : backgroundcosmology.h utils.h spline.h ODESolver.h
recombinationhistory.o  : recombinationhistory.h backgroundcosmology.h



#	-	-	-	-	-	-	-	-	-	-	-	-
# 	Commands:
#	-	-	-	-	-	-	-	-	-	-	-	-

# Command to create build directory
dir:
	mkdir -p $(BUILD_DIR)

# Command to build all C++-programs
all: $(BUILD_DIR)/$(TARGETS)

# Command to run program
run: 
	./$(BUILD_DIR)/$(TARGETS)

# Command to clean, build and run in one go (temporary)
autorun: clean all run

# Command to link
$(BUILD_DIR)/cmb: $(OBJS)
	${CXX} -o $@ $^ $(CXXFLAGS) $(INC) $(LIBS)

# Command to compile and link
$(BUILD_DIR)/%.o: %.cpp
	${CXX} -c -o $@ $< $(CXXFLAGS) $(INC) 

# Command to delete object files and executables
clean:
	rm -rf $(BUILD_DIR)/$(TARGETS) $(BUILD_DIR)/*.o

# Command to produce all plots (needs fix)
plots: 	
	cd plotting && python3 milestone_two.py
