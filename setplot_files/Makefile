# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC = gfortran
FFLAGS = -llapack -lblas -O2 -openmp

#FC = ifort
#FFLAGS = -O2 -openmp  -C -CB -CU -check all -traceback -fpe0 -ftrapuv -fp-model precise
#FFLAGS =  -C -CB -CU -fpe1 -ftrapuv  -O2 -fp-model precise -openmp

# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
VP2LIB = ../src

include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \
  $(GEOLIB)/qinit_module.f90 \

EXCLUDE_SOURCES = \
  $(GEOLIB)/qinit.f90 \
  $(GEOLIB)/b4step2.f90 \
  $(GEOLIB)/src2.f90 \
  $(GEOLIB)/setprob.f90 \
  $(GEOLIB)/bc2amr.f90 \
  $(GEOLIB)/src1d.f90 \
  $(GEOLIB)/setprob.f90 \
  $(GEOLIB)/valout.f \
  $(GEOLIB)/setaux.f90 \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
  $(VP2LIB)/qinit_module.f90 \
  $(VP2LIB)/vp2_module.f90 \

SOURCES = \
  $(VP2LIB)/qinit.f90 \
  $(VP2LIB)/src2.f90 \
  $(VP2LIB)/b4step2.f90 \
  $(VP2LIB)/bc2amr_vp2.f90 \
  $(VP2LIB)/src1d.f90 \
  $(VP2LIB)/setprob.f90 \
  $(VP2LIB)/rpn2_vp2.f90 \
  $(VP2LIB)/rpt2_vp2.f90 \
  $(VP2LIB)/vp2_riemann_utils.f \
  $(VP2LIB)/valout_vp2.f \
  $(VP2LIB)/setaux.f90 \
  $(VP2LIB)/dlap.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all
topo:
	python maketopo.py

all: 
	$(MAKE) topo
	$(MAKE) .plots
	$(MAKE) .htmls

