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
          nBandsSW_,     & ! Number of spectral bands (shortwave) 
          nBandsLW_        ! Number of spectral bands (longwave) 
    !> Optics grid in wave number [m-1]
    type(grid_t) :: grid_
    type(grid_t) :: gridLW_
    type(grid_t) :: gridSW_
  contains
    procedure :: name => model_name
    procedure :: create_state
    procedure :: optics_grid_lw
    procedure :: optics_grid_sw
    procedure :: compute_optics_lw
    procedure :: compute_optics_sw
  end type ufs_aerosol_optics_t

  !> Aerosol state specific to this model
  type, extends(state_t) :: ufs_state_t
     private
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
    ! Dependencies
    use aero_array,                only : array_t
    use aero_util,                 only : assert_msg
    use module_radsw_parameters,   only : NBDSW, wvnsw1=>wvnum1, wvnsw2=>wvnum2
    use module_radlw_parameters,   only : NBDLW, wvnlw1, wvnlw2
    use module_radiation_aerosols, only : aer_init
    use aero_constants,            only : rk => real_kind
#ifdef AERO_USE_NETCDF
    use netcdf, only : nf90_open, nf90_close, nf90_inq_varid, nf90_inquire_dimension,&
         NF90_NOWRITE, NF90_NOERR, nf90_get_var, nf90_inquire_variable
#endif
    type(ufs_aerosol_optics_t), pointer    :: model
    character(len=*), intent(in) :: description_file
    class(array_t), pointer :: interfaces
    integer :: i, status, ncid, nlevel, ntracer
    integer,parameter :: nCol = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! COMMENT: 
    ! All of this data ingest should go into state_t, not model_t
    ! It should be that state_t and model_t are provided to compute_optics, 
    ! which uses fields in state_t to produces the desirec optical fields
    ! defined in model_t.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    integer :: varID1,varID2,varID3,varID4,varID5,varID6,varID7,varID8,varID9
    integer, dimension(3) :: dimIDs
    real(rk),allocatable :: p_lay(:), p_lev(:), tv_lay(:), exner_lay(:), tracer(:,:), aero_mr(:,:)
    real(rk) :: lon,lat,lsmask

#ifdef AERO_USE_NETCDF
    if( len_trim( description_file ) > 0 ) then
       print*, "Reading in UFS aerosol driver data from ",description_file
       ! Open file
       status = nf90_open(description_file,nf90_nowrite,ncid)
       if (status .ne. nf90_noerr) call handle_err(status)
       ! Get variable IDs
       status = nf90_inq_varid(ncid,"lon",varID1)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"lat",varID2)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"lsmask",varID3)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"p_lev",varID4)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"p_lay",varID5)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"tv_lay",varID6)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"exner_lay",varID7)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"tracer",varID8)
       if (status .ne. nf90_noerr) call handle_err(status)
       status =nf90_inq_varid(ncid,"aero_mr",varID9)
       if (status .ne. nf90_noerr) call handle_err(status)
       ! Allocate space
       status = nf90_inquire_variable(ncid,varID9,dimIDs=dimIDs)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_inquire_dimension(ncid,dimIDs(1),len = nlevel)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_inquire_dimension(ncid,dimIDs(2),len = ntracer)
       if (status .ne. nf90_noerr) call handle_err(status)
       allocate(p_lay(nlevel), p_lev(nlevel), tv_lay(nlevel), exner_lay(nlevel),&
            tracer(nlevel,ntracer), aero_mr(nlevel,ntracer))
       ! Read in data
       status = nf90_get_var(ncid,varID1,lon)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID2,lat)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID3,lsmask)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID4,p_lev)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID5,p_lay)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID6,tv_lay)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID7,exner_lay)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID8,tracer)
       if (status .ne. nf90_noerr) call handle_err(status)
       status = nf90_get_var(ncid,varID9,aero_mr)
       if (status .ne. nf90_noerr) call handle_err(status)
       ! Close file
       status = nf90_close(ncid)
       if (status .ne. nf90_noerr) call handle_err(status)
    endif
#endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! END DATA INGEST BLOCK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    ! Not sure what this statement does?
    allocate(model)

    ! Call UFS aerosol optics initialization routine.
    call aer_init(nlevel,1)

    ! Setup shortwave spectral grid
    interfaces => array_t( (wvnsw1+wvnsw2)*0.5 )
    model%gridSW_ = grid_t(interfaces)

    ! Setup longwave spectral grid
    interfaces => array_t( (wvnlw1+wvnlw2)*0.5 )
    model%gridLW_ = grid_t(interfaces)

    ! Set flags to control LW/SW schemes. Just set to true for now.
    ! *NOTE* These are used when the UFS aerosol-optics are used "inline". 
    ! It is often the case that the radiation isn't called at every 
    ! timestep, and for the SW, at every grid-point (daylit only colummns 
    ! for SW radiaiton). So the ability to flip on/off this calculation is
    ! controlled with these flags, set by the host model.
    model%do_SW_ = .true.
    model%do_LW_ = .true.

    ! Setup aerosol optical outputs.
    model%nCol_     = nCol
    model%nLay_     = nlevel
    model%nTracer_  = ntracer
    model%nBandsSW_ = NBDSW
    model%nBandsLW_ = NBDLW

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the aerosol model/package
  function model_name( this )

    character(len=:), allocatable :: model_name
    class(ufs_aerosol_optics_t), intent(in) :: this

    model_name = "ufs aerosol optics"

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

  !> Returns the longwave aerosol optics grid, discretized in wavenumber space
  function optics_grid_lw( this )

    !> Copy of optical property wave number grid
    type(grid_t) :: optics_grid_lw
    !> My aerosol model
    class(ufs_aerosol_optics_t), intent(in) :: this

    optics_grid_lw = this%gridLW_%clone( )

  end function optics_grid_lw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the shortwave aerosol optics grid, discretized in wavenumber space
  function optics_grid_sw( this )

    !> Copy of optical property wave number grid
    type(grid_t) :: optics_grid_sw
    !> My aerosol model
    class(ufs_aerosol_optics_t), intent(in) :: this

    optics_grid_sw = this%gridSW_%clone( )

  end function optics_grid_sw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes longwave aerosol optical properties
  subroutine compute_optics_lw( this, state, od, od_ssa, od_asym )
    use aero_array,                   only : array_t
    use module_radiation_aerosols,    only : setaer, NSPC1
    use aero_constants,               only : rk => real_kind
    !
    class(ufs_aerosol_optics_t), intent(inout) :: this
    class(state_t),    intent(inout) :: state
    class(array_t),    intent(inout) :: od
    class(array_t),    intent(inout) :: od_ssa
    class(array_t),    intent(inout) :: od_asym
    ! Locals
    real(rk) :: aerSW_temp( this%nCol_, this%nLay_, this%nBandsSW_, 3), &
                aerLW_temp( this%nCol_, this%nLay_, this%nBandsLW_, 3), &
                aerodp_temp(this%nCol_, NSPC1),                         &
                od_temp( this%nBandsLW_ )

    select type( state )
    class is( ufs_state_t )
       call setaer(state%plev_, state%play_, state%prslk_, state%tvlay_,   &
            state%rhlay_, state%lsmask_, state%tracer_, state%aerfld_,     &
            state%lon_, state%lat_, this%nCol_, this%nLay_, this%nLay_+1,  &
            this%do_SW_, this%do_LW_, aerSW_temp, aerLW_temp, aerodp_temp)
    end select
    ! this first hack-a-thon is just assuming a box model
    od_temp = aerLW_temp(1,1,:,1)
    call od%copy_in( od_temp )
    od_temp = od_temp(:) * aerLW_temp(1,1,:,2)
    call od%copy_in( od_temp )
    od_temp = od_temp(:) * aerLW_temp(1,1,:,3)
    call od%copy_in( od_temp )

  end subroutine compute_optics_lw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes shortwave aerosol optical properties
  subroutine compute_optics_sw( this, state, od, od_ssa, od_asym )
    use aero_array,                   only : array_t
    use module_radiation_aerosols,    only : setaer, NSPC1
    use aero_constants,               only : rk => real_kind
    !
    class(ufs_aerosol_optics_t), intent(inout) :: this
    class(state_t),    intent(inout) :: state
    class(array_t),    intent(inout) :: od
    class(array_t),    intent(inout) :: od_ssa
    class(array_t),    intent(inout) :: od_asym
    ! Locals
    real(rk) :: aerSW_temp( this%nCol_, this%nLay_, this%nBandsSW_, 3), &
                aerLW_temp( this%nCol_, this%nLay_, this%nBandsLW_, 3), &
                aerodp_temp(this%nCol_, NSPC1),                         &
                od_temp( this%nBandsSW_ )

    select type( state )
    class is( ufs_state_t )
       call setaer(state%plev_, state%play_, state%prslk_, state%tvlay_,   &
            state%rhlay_, state%lsmask_, state%tracer_, state%aerfld_,     &
            state%lon_, state%lat_, this%nCol_, this%nLay_, this%nLay_+1,  &
            this%do_SW_, this%do_LW_, aerSW_temp, aerLW_temp, aerodp_temp)
    end select
    ! this first hack-a-thon is just assuming a box model
    od_temp = aerSW_temp(1,1,:,1)
    call od%copy_in( od_temp )
    od_temp = od_temp(:) * aerSW_temp(1,1,:,2)
    call od%copy_in( od_temp )
    od_temp = od_temp(:) * aerSW_temp(1,1,:,3)
    call od%copy_in( od_temp )

  end subroutine compute_optics_sw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine handle_err(status)
#ifdef AERO_USE_NETCDF
    use netcdf, only : NF90_NOERR, nf90_strerror
#endif
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, "ERROR:", trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine handle_err

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ufs_aerosol_optics
