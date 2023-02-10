# Hans A. Winther (2020) (hans.a.winther@gmail.com)

SHELL := /bin/bash
HOME := /uio/hume/student-u29/nannabr

# Set compiler (use =c++17 if you have this availiable)
CC = g++ -std=c++11 

# Paths to GSL library
# INC  = -I/mn/stornext/u3/hansw/winther/local/include
# LIBS = -L/mn/stornext/u3/hansw/winther/local/lib -lgsl -lgslcblas
INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

#=======================================================
# Options
#=======================================================
OPTIONS = 

# Add bounds checking
OPTIONS += -D_GLIBCXX_DEBUG

# Show warnings if atempting to evaluate a spline out of bounds
OPTIONS += -D_SPLINE_WARNINGS_ON

# Show info about the solution as we integrate
# OPTIONS = -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE

# Add OpenMP parallelization
# OPTIONS += -D_USEOPEMP
# CC += -fopenmp

#=======================================================

C = -O3 -g $(OPTIONS)

#=======================================================

BUILD_DIR := build
dir:
	mkdir -p $(BUILD_DIR)

# VPATH=src/
TARGETS := cmb
all: $(BUILD_DIR)/$(TARGETS)

run: 
	./$(BUILD_DIR)/$(TARGETS)

# OBJECT FILES
# OBJS = main.o utils.o backgroundcosmology.o RecombinationHistory.o Perturbations.o PowerSpectrum.o spline.o ODESolver.o
OBJS_ = main.o utils.o backgroundcosmology.o spline.o ODESolver.o
OBJS = $(addprefix $(BUILD_DIR)/, $(OBJS_))

# DEPENDENCIES
# main.o                  : backgroundcosmology.h RecombinationHistory.h Perturbations.h PowerSpectrum.h
main.o					: backgroundcosmology.h
spline.o                : spline.h
ODESolver.o             : ODESolver.h
utils.o                 : utils.h spline.h ODESolver.h
backgroundcosmology.o   : backgroundcosmology.h utils.h spline.h ODESolver.h
# RecombinationHistory.o  : RecombinationHistory.h backgroundcosmology.h
# Perturbations.o         : Perturbations.h backgroundcosmology.h RecombinationHistory.h
# PowerSpectrum.o         : PowerSpectrum.h backgroundcosmology.h RecombinationHistory.h Perturbations.h

$(BUILD_DIR)/cmb: $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIBS)

$(BUILD_DIR)/%.o: %.cpp
	${CC} -c -o $@ $< $C $(INC) 

clean:
	rm -rf $(BUILD_DIR)/$(TARGETS) $(BUILD_DIR)/*.o
