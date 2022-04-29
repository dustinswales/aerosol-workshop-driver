################################################################################
# Aerosol library helper function

include(model_config)

################################################################################
# Aerosol model build configuration

create_aerosol_model(NAME "ufs_aero"      # String a host model will use to find your model
                     TYPE type_ufs_aero   # Derived type describing your model
		     F90_MODULE int_ufs_aero # Fortran module where your derived type is defined
                     LIBRARY ufs_aero     # Name of the library your model will be built into
                     SOURCE UFS/ccpp-physics/physics/aerclm_def.F UFS/ccpp-physics/physics/funcphys.f90 UFS/ccpp-physics/physics/iounitdef.f UFS/ccpp-physics/physics/machine.F UFS/ccpp-physics/physics/physcons.F90 UFS/ccpp-physics/physics/physparam.f UFS/ccpp-physics/physics/radiation_aerosols.f UFS/ccpp-physics/physics/radlw_param.f UFS/ccpp-physics/physics/radsw_param.f UFS/int_ufs_aero.F90) # Space-separated list of all your aerosol model source files

################################################################################
