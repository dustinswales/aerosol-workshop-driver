! Copyright (C) 2022 National Center for Atmospheric Research,
! National Technology & Engineering Solutions of Sandia, LLC (NTESS),
! and the U.S. Environmental Protection Agency (USEPA)
!
! SPDX-License-Identifier: Apache-2.0
!
module ufs_aerosol_optics

  use aero_constants,                  only : rk => real_kind
  use aero_grid,                       only : grid_t
  use aero_model,                      only : model_t
  use aero_state,                      only : state_t
  implicit none
  private

  public :: ufs_aerosol_optics_t

  !> UFS aerosol model parameters
  type, extends(model_t) :: ufs_aerosol_optics_t
    private
     logical :: &
          do_SW_,        & ! Compute shortwave aerosol optical properties?
          do_LW_           ! Compute longwave aerosol optical properties?
     integer :: &
          nCol_,         & ! Number of horizontal columns
          nLay_,         & ! Number of vertical layers
          nTracer_,      & ! Number of aerosol tracers
          nBandsLW_,     & ! Number of spectral bands (longwave)
          nBandsSW_        ! Number of spectral bands (shortwave) 
     real(rk), allocatable :: &
          tauLW_(:,:,:), & ! Optical depth (longwave)
          ssaLW_(:,:,:), & ! Single scattering albedo (longwave)
          gLW_(:,:,:),   & ! Asymmetry parameter (longwave)
          tauSW_(:,:,:), & ! Optical depth (shortwave)
          ssaSW_(:,:,:), & ! Single scattering albedo (shortwave)
          gSW_(:,:,:)      ! Asymmetry parameter (shortwave)
    !> Optics grid in wave number [m-1]
    type(grid_t) :: gridLW_, gridSW_
  contains
    procedure :: name => model_name
    procedure :: create_state
    procedure :: optics_grid
    procedure :: compute_optics
  end type ufs_aerosol_optics_t

  !> Aerosol state specific to this model
  type, extends(state_t) :: ufs_state_t
     private
     integer :: &
          nCol_,           & ! Number of horizontal columns
          nLay_,           & ! Number of vertical layers
          nTracer_           ! Number of aerosol tracers
     real(rk), allocatable :: &
          lon_(:),         & ! Longitude
          lat_(:),         & ! Latitude
          lsmask_(:),      & ! Land/sea mask (sea:0,land:1,sea-ice:2)
          plev_(:,:),      & ! Pressure at model-interface (mb)
          play_(:,:),      & ! Pressure at model-layer (mb)
          prslk_(:,:),     & ! exner function = (p/p0)**rocp at model-layer (1)
          tvlay_(:,:),     & ! Virtual temperature at model-layer  (k)
          rhlay_(:,:),     & ! Relative humidity at model-layer (1)
          tracer_(:,:,:),  & ! aerosol tracer concentration (kg/kg)
          aerfld_(:,:,:)     ! GOCART aerosol climatology number concentration (kg-1)
  end type ufs_state_t

  interface ufs_aerosol_optics_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates and configure the aerosol model
  function constructor( description_file ) result( model )

    use aero_array,                    only : array_t
    use aero_util,                     only : assert_msg
    use module_radsw_parameters,       only : NBDSW, wvnsw1=>wvnum1, wvnsw2=>wvnum2
    use module_radlw_parameters,       only : NBDLW, wvnlw1, wvnlw2
    use module_radiation_aerosols,     only : aer_init
#ifdef AERO_USE_NETCDF
    use netcdf,                        only : nf90_open, nf90_close,          &
                                              NF90_NOWRITE, NF90_NOERR
#endif
    type(ufs_aerosol_optics_t), pointer    :: model
    character(len=*), intent(in) :: description_file

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Set parameters/configuration data for the aerosol package here !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialize the aerosol grid with wavelength data pulled from
    ! https://acp.copernicus.org/articles/18/7815/2018/acp-18-7815-2018-f03.pdf
    real(kind=rk) :: wavelengths(4) = & ! [nm]
      (/ 440.0_rk, 675.0_rk, 870.0_rk, 1020.0_rk /)
    real(kind=rk) :: wave_numbers(4) ! [m-1]
    class(array_t), pointer :: interfaces, interfacesLW, interfacesSW
    integer :: i, netcdf_file
    integer,parameter :: nCol    = 1
    integer,parameter :: nLay    = 1
    integer,parameter :: nTracer = 1

    ! Longwave spectral grid 
    interfacesLW => array_t( (wvnlw1+wvnlw2)*0.5 )
    allocate( model )
    model%gridLW_ = grid_t( interfacesLW )
    
    ! Shortwave spectral grid
    interfacesSW => array_t( (wvnsw1+wvnsw2)*0.5 )
    allocate( model )
    model%gridSW_ = grid_t( interfacesSW )

    model%do_SW_ = .true.
    model%do_LW_ = .true.
    call aer_init(nLay,1)
    model%nCol_     = nCol
    model%nLay_     = nLay
    model%nTracer_  = nTracer
    model%nBandsLW_ = NBDLW
    model%nBandsSW_ = NBDSW
    allocate(model%tauLW_(model%nCol_, model%nLay_, model%nBandsLW_))
    allocate(model%ssaLW_(model%nCol_, model%nLay_, model%nBandsLW_))
    allocate(model%gLW_(  model%nCol_, model%nLay_, model%nBandsLW_))
    allocate(model%tauSW_(model%nCol_, model%nLay_, model%nBandsSW_))
    allocate(model%ssaSW_(model%nCol_, model%nLay_, model%nBandsSW_))
    allocate(model%gSW_(  model%nCol_, model%nLay_, model%nBandsSW_))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the aerosol model/package
  function model_name( this )

    character(len=:), allocatable :: model_name
    class(ufs_aerosol_optics_t), intent(in) :: this

    model_name = "my model"

  end function model_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a newly created aerosol state
  function create_state( this ) result( state )

    class(state_t),    pointer    :: state
    class(ufs_aerosol_optics_t), intent(in) :: this

    allocate( ufs_state_t :: state )
    select type( state )
    class is( ufs_state_t )
      allocate(state%lon_(   this%nCol_),                            &
               state%lat_(   this%nCol_),                            &
               state%lsmask_(this%nCol_),                            &
               state%plev_(  this%nCol_, this%nLay_+1),              &
               state%play_(  this%nCol_, this%nLay_),                &
               state%prslk_( this%nCol_, this%nLay_),                &
               state%tvlay_( this%nCol_, this%nLay_),                &
               state%rhlay_( this%nCol_, this%nLay_),                &
               state%tracer_(this%nCol_, this%nLay_, this%nTracer_), &
               state%aerfld_(this%nCol_, this%nLay_, this%nTracer_))
    end select

  end function create_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the aerosol optics grid, discretized in wavenumber space
  function optics_grid( this )

    !> Copy of optical property wave number grid
    type(grid_t) :: optics_grid
    !> My aerosol model
    class(ufs_aerosol_optics_t), intent(in) :: this

    optics_grid = this%gridLW_%clone( )

  end function optics_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes optical property data, given an aerosol state and destination
  !! arrays
  subroutine compute_optics(this, state, od, od_ssa, od_asym )
    use aero_array,                   only : array_t
    use module_radiation_aerosols,    only : setaer
    !
    class(ufs_aerosol_optics_t), intent(inout) :: this
    class(state_t),    intent(inout) :: state
    class(array_t),    intent(inout) :: od         ! The code breaks if I remove these?
    class(array_t),    intent(inout) :: od_ssa     ! The code breaks if I remove these?
    class(array_t),    intent(inout) :: od_asym    ! The code breaks if I remove these?
    ! Locals
    real(rk),allocatable :: aerSW_temp(:,:,:,:), aerLW_temp(:,:,:,:), aerodp_temp(:,:)

    select type( state )
    class is( ufs_state_t )
       allocate(aerSW_temp(state%nCol_, state%nLay_, this%nBandsSW_, 3))
       allocate(aerLW_temp(state%nCol_, state%nLay_, this%nBandsLW_, 3))
       call setaer(state%plev_, state%play_, state%prslk_, state%tvlay_,     &
            state%rhlay_, state%lsmask_, state%tracer_, state%aerfld_,       &
            state%lon_, state%lat_, state%nCol_, state%nLay_, state%nLay_+1, &
            this%do_SW_, this%do_LW_, aerSW_temp, aerLW_temp, aerodp_temp)
       this%tauLW_ = aerLW_temp(:,:,:,1)
       this%ssaLW_ = aerLW_temp(:,:,:,2)
       this%gLW_   = aerLW_temp(:,:,:,3)
       this%tauLW_ = aerSW_temp(:,:,:,1)
       this%ssaLW_ = aerSW_temp(:,:,:,2)
       this%gLW_   = aerSW_temp(:,:,:,3)
    end select

  end subroutine compute_optics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ufs_aerosol_optics
