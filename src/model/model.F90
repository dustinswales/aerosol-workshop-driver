! Copyright (C) 2022 National Center for Atmospheric Research,
! National Technology & Engineering Solutions of Sandia, LLC (NTESS),
! and the U.S. Environmental Protection Agency (USEPA)
!
! SPDX-License-Identifier: Apache-2.0
!
module aero_model

  use aero_constants,              only : real_kind

  implicit none
  private

  public :: model_t, model_ptr

  type, abstract :: model_t
  contains
    procedure(model_name),     deferred :: name
    procedure(create_state),   deferred :: create_state
    procedure :: optics_grid
    procedure :: optics_grid_lw
    procedure :: optics_grid_sw
    procedure :: compute_optics
    procedure :: compute_optics_lw
    procedure :: compute_optics_sw
  end type model_t

  type :: model_ptr
    class(model_t), pointer :: ptr_ => null( )
  contains
    final :: model_ptr_finalize
  end type model_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the aerosol model/package
  function model_name( this )

    import :: model_t

    !> Unique name for the aerosol model
    character(len=:), allocatable :: model_name
    !> Aerosol model
    class(model_t), intent(in) :: this

  end function model_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an newly created aerosol state
  function create_state( this ) result( state )

    use aero_state,                    only : state_t
    import :: model_t

    !> New aerosol state
    class(state_t), pointer    :: state
    !> Aerosol model
    class(model_t), intent(in) :: this

  end function create_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the aerosol optics grid, discretized in wavenumber space
  function optics_grid( this )

    use aero_grid,                     only : grid_t
    use aero_util,                     only : die_msg

    !> Optical property wave number grid
    type(grid_t) :: optics_grid
    !> Aerosol model
    class(model_t), intent(in) :: this

    call die_msg( 128803475, "This aerosol model does not provide optical "// &
                             "properties on a continuous wavelength grid" )

  end function optics_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the longwave aerosol optics grid, discretized in wavenumber space
  function optics_grid_lw( this )

    use aero_grid,                     only : grid_t
    use aero_util,                     only : die_msg

    !> Optical property wave number grid
    type(grid_t) :: optics_grid_lw
    !> Aerosol model
    class(model_t), intent(in) :: this

    call die_msg( 957529823, "This aerosol model does not provide optical "// &
                             "properties on a longwave wavelength grid" )

  end function optics_grid_lw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the shortwave aerosol optics grid, discretized in wavenumber
  !! space
  function optics_grid_sw( this )

    use aero_grid,                     only : grid_t
    use aero_util,                     only : die_msg

    !> Optical property wave number grid
    type(grid_t) :: optics_grid_sw
    !> Aerosol model
    class(model_t), intent(in) :: this

    call die_msg( 173206942, "This aerosol model does not provide optical "// &
                             "properties on a shortwave wavelength grid" )

  end function optics_grid_sw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates optical property data, given an aerosol state and destination
  !! arrays
  subroutine compute_optics( this, state, od, od_ssa, od_asym )

    use aero_array,                    only : array_t
    use aero_state,                    only : state_t
    use aero_util,                     only : die_msg

    !> Aerosol model
    class(model_t),    intent(inout) :: this
    !> Aerosol state
    class(state_t),    intent(inout) :: state
    !> Aerosol optical depth [m]
    class(array_t),    intent(inout) :: od
    !> Aerosol scattering optical depth [m]
    class(array_t),    intent(inout) :: od_ssa
    !> Aerosol asymmetric scattering optical depth [m]
    class(array_t),    intent(inout) :: od_asym

    call die_msg( 324664447, "This aerosol model does not provide optical "// &
                             "properties on a continuous wavelength grid" )

  end subroutine compute_optics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates optical property data, given an aerosol state and destination
  !! arrays for the longwave
  subroutine compute_optics_lw( this, state, od, od_ssa, od_asym )

    use aero_array,                    only : array_t
    use aero_state,                    only : state_t
    use aero_util,                     only : die_msg

    !> Aerosol model
    class(model_t),    intent(inout) :: this
    !> Aerosol state
    class(state_t),    intent(inout) :: state
    !> Aerosol optical depth [m]
    class(array_t),    intent(inout) :: od
    !> Aerosol scattering optical depth [m]
    class(array_t),    intent(inout) :: od_ssa
    !> Aerosol asymmetric scattering optical depth [m]
    class(array_t),    intent(inout) :: od_asym

    call die_msg( 277354502, "This aerosol model does not provide optical "// &
                             "properties on a longwave wavelength grid" )

  end subroutine compute_optics_lw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates optical property data, given an aerosol state and destination
  !! arrays for the shortwave
  subroutine compute_optics_sw( this, state, od, od_ssa, od_asym )

    use aero_array,                    only : array_t
    use aero_state,                    only : state_t
    use aero_util,                     only : die_msg

    !> Aerosol model
    class(model_t),    intent(inout) :: this
    !> Aerosol state
    class(state_t),    intent(inout) :: state
    !> Aerosol optical depth [m]
    class(array_t),    intent(inout) :: od
    !> Aerosol scattering optical depth [m]
    class(array_t),    intent(inout) :: od_ssa
    !> Aerosol asymmetric scattering optical depth [m]
    class(array_t),    intent(inout) :: od_asym

    call die_msg( 335193061, "This aerosol model does not provide optical "// &
                             "properties on a shortwave wavelength grid" )

  end subroutine compute_optics_sw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees resources assocaited with a model pointer
  subroutine model_ptr_finalize( this )

    type(model_ptr), intent(inout) :: this

    if( associated( this%ptr_ ) ) deallocate( this%ptr_ )

  end subroutine model_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module aero_model
