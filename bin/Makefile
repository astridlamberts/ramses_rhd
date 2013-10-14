#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@cea.fr
#############################################################################
# Compilation time parameters
NDIM = 3
NPRE = 8
SOLVER = rhd
#PATCH=../patch/
#PATCH=../patch/rhd/merger_EOS/
PATCH=../patch/rhd/wind/
#PATCH=../patch/rhd/test_nrj_cons_paper/
#PATCH=../patch/rhd/blast_3D/
EXEC = ramses
#############################################################################
# Fortran compiler options and directives


# --- No MPI, gfortran ----------------------------------
#F90 = gfortran
#CFLAGS = -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DWITHOUTMPI
#FFLAGS = -x f95-cpp-input $(CFLAGS)

## --- No MPI, pgf90 ----------------------------------
#F90 = pgf90 -O3
#FFLAGS = -Mpreprocess -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DWITHOUTMPI

# --- No MPI, xlf ------------------------------------
#F90 = xlf
#FFLAGS = -WF,-DNDIM=$(NDIM),-DNPRE=$(NPRE),-DSOLVER$(SOLVER),-DWITHOUTMPI -qfree=f90 -qsuffix=f=f90 -qsuffix=cpp=f90

# --- No MPI, f90 ------------------------------------
#F90 = f90
#FFLAGS = -cpp -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DWITHOUTMPI

# --- No MPI, ifort ----------------------------------
#F90 = ifort 
#FFLAGS =  -cpp -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DWITHOUTMPI

# --- MPI, pgf90 syntax ------------------------------
##F90 = mpif90 -Ktrap=fp,align,unf,inexact -Mdaz #-O3#4 -C -g
#FFLAGS =  -Mpreprocess -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER)  #-DQUADHILBERT

# --- MPI, ifort syntax ------------------------------
#F90 = mpif90
#FFLAGS = -cpp -g -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DNOSYSTEM

# --- MPI, ifort syntax, Potowatomi ------------------------------
F90 = mpif90 -f90=ifort
FFLAGS = -cpp -g -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DNOSYSTEM


# --- MPI, ifort syntax, additional checks -----------
#F90 = mpif90
#FFLAGS =  -O3 -g -fpe -traceback -ftrapuv -cpp -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DNOSYSTEM

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI = 
#LIBMPI = -lfmpi -lmpi -lelan
LIBS = $(LIBMPI)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../hydro:../pm:../poisson:../amr
#VPATH = $(PATCH):../mhd_old:../hydro:../pm:../poisson:../amr
#############################################################################
# All objects
MODOBJ = amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o  # wind_parameters.o
AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o amr_step.o update_time.o output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o refine_utils.o nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o units.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o sink_particle.o feedback.o
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o read_hydro_params.o synchro_hydro_fine.o 
# All objects
AMRLIB = $(MODOBJ) $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ)
#############################################################################
ramses:	$(AMRLIB) ramses.o
	$(F90) $(FFLAGS) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
#############################################################################
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
#############################################################################
clean :
	rm *.o *.$(MOD)
#############################################################################
