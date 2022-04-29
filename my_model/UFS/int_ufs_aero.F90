module int_ufs_aero
  type :: type_ufs_aero
     real,dimension(:),allocatable :: prsi ! pressure at model-interface (mb)
     real,dimension(:),allocatable :: prsl ! pressure at model-layer (mb)
  end type type_ufs_aero
end module int_ufs_aero
