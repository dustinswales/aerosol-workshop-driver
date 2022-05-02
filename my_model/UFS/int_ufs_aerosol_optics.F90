module int_ufs_aerosol_optics
  use aero_constants,  only : wp => real_kind
  implicit none
  type, public :: type_ufs_aerosol_optics
     logical :: &
          lsswr,          & ! Compute shortwave aerosol optical properties?
          lslwr             ! Compute longwave aerosol optical properties?
     integer :: &
          nCol,           & ! Number of horizontal columns
          nLay              ! Number of vertical layers
     real(wp), allocatable :: &
          lon(:),         & ! Longitude
          lat(:),         & ! Latitude
          lsmask(:),      & ! Land/sea mask (sea:0,land:1,sea-ice:2)
          plev(:,:),      & ! Pressure at model-interface (mb)
          play(:,:),      & ! Pressure at model-layer (mb)
          prslk(:,:),     & ! exner function = (p/p0)**rocp at model-layer (1)
          tvlay(:,:),     & ! Virtual temperature at model-layer  (k)
          rhlay(:,:),     & ! Relative humidity at model-layer (1)
          tracer(:,:,:),  & ! aerosol tracer concentration (kg/kg)
          aerfld(:,:,:)     ! GOCART aerosol climatology number concentration (kg-1)
  end type type_ufs_aerosol_optics
end module int_ufs_aerosol_optics
