# UFS ccpp-physics aerosol/radiation
foreach (src ccpp-physics/physics/aerclm_def.F
	     ccpp-physics/physics/funcphys.f90
	     ccpp-physics/physics/iounitdef.f
	     ccpp-physics/physics/machine.F
	     ccpp-physics/physics/physcons.F90
	     ccpp-physics/physics/physparam.f
	     ccpp-physics/physics/radiation_aerosols.f
	     ccpp-physics/physics/radlw_param.f
	     ccpp-physics/physics/radsw_param.f)
  list(APPEND ufs_src ${CMAKE_CURRENT_SOURCE_DIR}/${src})
endforeach()
set(ufs_src ${ufs_src} PARENT_SCOPE)

