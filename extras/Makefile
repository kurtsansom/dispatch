# Set search path.  Any extras module that needs to have a special 
# version for a solver should be place in solvers/$(SOLVER)/extra,
# and if a local solver dependent version is needed, it should be
# placed in /$(SOLVER)/extra, below the experiment directory
#
VPATH := $(VPATH):$(SOLVER)/extras:$(TOP)/solvers/$(SOLVER)/extras:$(TOP)/extras
include $(TOP)/extras/spitzer_conductivity/Makefile

# Dependencies
extras_mod.o            : io_mod.o trace_mod.o patch_mod.o gpatch_mod.o connect_mod.o task_mod.o list_mod.o
selfgravity_mod.o       : patch_mod.o poisson_mod.o vector_ops.o download_mod.o
gravity_mod.o           : patch_mod.o units_mod.o scaling_mod.o link_mod.o io_mod.o
trace_particles_mod.o   : io_mod.o trace_mod.o dll_mod.o omp_timer_mod.o random_mod.o \
                          mesh_mod.o link_mod.o patch_mod.o mpi_buffer_mod.o \
                          particle_mod.o particle_list_mod.o
sinkparticles_mod.o     : io_unit_mod.o trace_mod.o download_mod.o data_hub_mod.o link_mod.o \
                          list_mod.o particle_history_mod.o gpatch_mod.o
sinkparticle_mod.o      : io_unit_mod.o trace_mod.o download_mod.o data_hub_mod.o link_mod.o \
                          list_mod.o particle_history_mod.o gpatch_mod.o
sink_patch_mod.o        : io_unit_mod.o trace_mod.o download_mod.o data_hub_mod.o link_mod.o \
                          list_mod.o gpatch_mod.o particle_solver_mod.o
connect_mod.o           :
