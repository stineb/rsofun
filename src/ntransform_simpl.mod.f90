module md_ntransform_simpl
  !////////////////////////////////////////////////////////////////
  ! INORGANIC NITROGEN DYNAMICS MODULE, SIMPLE
  !----------------------------------------------------------------
  use md_classdefs
  use md_tile_cnmodel
  use md_params_core
  use md_rates

  implicit none

  private 
  public ntransform, getpar_modl_ntransform

  !-----------------------------------------------------------------------
  ! Module-specific model parameters
  !-----------------------------------------------------------------------
  type params_ndecay_type
    real :: kdecay_ninorg
  end type params_ndecay_type

  type( params_ndecay_type ) :: params_ndecay

contains

  subroutine ntransform( tile, tile_fluxes, landuse, aprec, doy )
    !////////////////////////////////////////////////////////////////
    ! Simple N loss as a function of soil moisture and temperature, 
    ! and of runoff. Similar as in BiomeE. Not distinguising between
    ! NO3 and NH4.
    !----------------------------------------------------------------
    use md_forcing_cnmodel, only: landuse_type

    ! arguments
    type(tile_type), dimension(nlu), intent(inout) :: tile
    type(tile_fluxes_type), dimension(nlu), intent(inout) :: tile_fluxes
    real :: aprec         ! annual total precipitation [mm/d] 
    type(landuse_type), intent(in) :: landuse
    real, optional :: aprec         ! annual total precipitation [mm/d] 
    integer, optional :: doy

    ! local variables
    real :: n_gasloss
    real :: ftemp_nloss
    integer :: lu
    
    if (present(aprec) .and. present(doy)) then
      if ( doy == 1 ) then
        !///////////////////////////////////////////////////////////////////////
        ! ANNUAL INITIALIZATION 
        !-----------------------------------------------------------------------
        ! Calculate soil PH using empirical relationship with annual precip
        ! Eq.5, Tab.5, XP08 (ntransform.cpp:65) (c++:aprec in mm/yr; F: mm/yr)
        !------------------------------------------------------------------
        ph_soil = 3810.0 / (762.0 + aprec) + 3.8

        ! Deprotonation of NH4 to NH3 depends on soil pH
        !------------------------------------------------------------------
        if (ph_soil > 6.0) then
          nh3max = 1.0
        else
          nh3max = 0.00001
        endif
        
      end if

    else

      doy = dummy
      aprec = dummy
      nh3max = 0.0

    end if
    

    ! LOOP OVER GRIDCELL LAND UNITS
    luloop: do lu=1,nlu

      !-------------------------------------------------------------------------
      ! Add N deposition to inorganic pools
      !-------------------------------------------------------------------------
      tile(lu)%soil%pno3%n14 = tile(lu)%soil%pno3%n14 + landuse%dnoy
      tile(lu)%soil%pnh4%n14 = tile(lu)%soil%pnh4%n14 + landuse%dnhx

      !------------------------------------------------------------------      
      ! NITRATE LEACHING
      !------------------------------------------------------------------      
      tile_fluxes(lu)%soil%dnleach = tile(lu)%soil%pno3%n14 * tile_fluxes(lu)%canopy%dfleach
      tile(lu)%soil%pno3%n14       = tile(lu)%soil%pno3%n14 - tile_fluxes(lu)%soil%dnleach
      tile_fluxes(lu)%soil%dnloss  = tile_fluxes(lu)%soil%dnloss + tile_fluxes(lu)%soil%dnleach

      !------------------------------------------------------------------      
      ! GASEOUS LOSS
      ! No distinction between NO3 and NH4. Only accounting for a 
      ! temperature-dependence of the gaseous loss rate.
      !------------------------------------------------------------------
      ! same temperature dependence as denitrification in Xu-Ri & Prentice, 2008
      ! unity at 22 deg C, increases to 2 at about 34 deg C
      ftemp_nloss = ftemp( tile(lu)%soil%phy%temp, "lloyd_and_taylor", ref_temp = 22.0 )

      ! NO3
      n_gasloss = (1.0 - exp(-ftemp_nloss * params_ndecay%kdecay_ninorg)) * tile(lu)%soil%pno3%n14
      tile(lu)%soil%pno3%n14 = tile(lu)%soil%pno3%n14 - n_gasloss
      tile_fluxes(lu)%soil%dnloss = tile_fluxes(lu)%soil%dnloss + n_gasloss

      ! NH4
      n_gasloss = (1.0 - exp(-ftemp_nloss * params_ndecay%kdecay_ninorg)) * tile(lu)%soil%pnh4%n14
      tile(lu)%soil%pnh4%n14 = tile(lu)%soil%pnh4%n14 - n_gasloss
      tile_fluxes(lu)%soil%dnloss = tile_fluxes(lu)%soil%dnloss + n_gasloss

    end do luloop

  end subroutine ntransform


  subroutine getpar_modl_ntransform()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads waterbalance module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use md_interface_cnmodel, only: myinterface

    ! maximum nitrification rate
    params_ndecay%kdecay_ninorg = myinterface%params_calib%maxnitr  ! re-interpret available parameter

  end subroutine getpar_modl_ntransform


end module md_ntransform_simpl
