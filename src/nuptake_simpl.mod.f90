module md_nuptake_simpl
  !////////////////////////////////////////////////////////////////
  ! SIMPLE NITROGEN UPTAKE MODULE
  !----------------------------------------------------------------
  use md_params_core, only: ndayyear, nmonth, nlu, npft
  use md_classdefs
  use md_tile_cnmodel
  use md_plant_cnmodel

  implicit none

  private
  public nuptake, getpar_modl_nuptake

  !-----------------------------------------------------------------------
  ! Module parameters
  !-----------------------------------------------------------------------
  type params_nuptake_type
    real :: kc         ! half-saturation constant with respect to root biomass (gC m-2)
    real :: kv         ! half-saturation constant with respect to inorganic soil N (gN m-2)
    real :: vmax       ! maximum rate of N uptake (at saturating root biomass and soil inorganic N) (gN d-1)

    real :: eff_nup           ! uptake efficiency for equation (gN/gC)
    real :: minimumcostfix    ! Minimum cost of N-fixation at optimal soil temperature, is 4.8 gC/gN, value from Gutschik (1981)
    real :: fixoptimum        ! Optimum soil temperature for N fixation. Taken to be equal to optimum temperature for nitrogenase activity as given by Houlton et al. (2008), Nature: 25.15+-0.66 
    real :: a_param_fix       ! shape parameter of the N fixation function taken to be equal to nitrogenase activity function given Houlton et al. (2008), Nature: -3.62+-0.52
    real :: b_param_fix       ! shape parameter of the N fixation function taken to be equal to nitrogenase activity function given Houlton et al. (2008), Nature: 0.27+-0.04 

  end type params_nuptake_type

  type( params_nuptake_type ) :: params_nuptake

  !----------------------------------------------------------------
  ! module-specific (private) variables
  !----------------------------------------------------------------
  type outtype_calc_dnup
    real :: act_nh4                         ! NH4 acquisition by active uptake [gN/m2/d]
    real :: act_no3                         ! NO3 acquisition by active uptake [gN/m2/d]
    real :: fix                             ! N acquisition by fixation [gN/m2/d]
  end type outtype_calc_dnup

contains

  subroutine nuptake( tile, tile_fluxes )
    !/////////////////////////////////////////////////////////////////
    ! SUBROUTINE NUPTAKE ASSUMING CONSTANT EXUDATION PER UNIT ROOT MASS
    !-----------------------------------------------------------------
    ! This model calculates first the passive uptake of N via the 
    ! transpiration stream.
    !-----------------------------------------------------------------
    ! arguments
    type(tile_type), dimension(nlu), intent(inout) :: tile
    type(tile_fluxes_type), dimension(nlu), intent(inout) :: tile_fluxes

    ! local variables
    integer :: lu, pft
    real :: fno3
    real :: dnup

    !-----------------------------------------------------------------
    ! Get total inorganic N. N uptake is assumed to deplete NO3 and NH4
    ! in proportion to their relative shares.
    !-----------------------------------------------------------------
    n0 = nh4 + no3
    
    ! if (n0>0.0) then
    fno3 = no3 / n0

    pftloop: do pft = 1, npft

      lu = params_pft_plant(pft)%lu_category

      !--------------------------------------------------------------------------
      ! N uptake of NO3 and NH4 in proportion to their relative pool sizes.
      !--------------------------------------------------------------------------
      fno3 = tile(lu)%soil%pno3%n14 / (tile(lu)%soil%pnh4%n14 + tile(lu)%soil%pno3%n14)
      dnup  = calc_dnup( &
        tile(lu)%plant(pft)%proot%c%c12, &
        tile(lu)%soil%pnh4%n14 + tile(lu)%soil%pno3%n14 &
        )

      if (params_pft_plant(pft)%nfixer) then
        !--------------------------------------------------------------------------
        ! Determine split of root uptake versus N fixation based on relative 
        ! respective efficiencies (=N acquisition over C cost).
        !--------------------------------------------------------------------------
        ! C cost of root N uptake
        ! Calculated as the C spent per time step (day) for sustaining the current
        ! root biomass.
        !--------------------------------------------------------------------------
        resp_rate = calc_resp_maint( 
          1.0, &
          params_plant%r_root, &
          tile(lu)%soil%phy%temp &
          )
        dcup = tile(lu)%plant(pft)%proot%c%c12 * &
          (params_pft_plant(pft)%k_decay_root / params_plant%growtheff + resp_rate)

        !--------------------------------------------------------------------------
        ! Efficiency of root uptake
        !--------------------------------------------------------------------------
        eff_nup = dnup / dcup

        !-----------------------------------------------------------------
        ! Get efficiency of BNF in gN/gC as a function of soil temperature
        !-----------------------------------------------------------------
        eff_nfix = calc_eff_nfix( tile(lu)%soil%phy%temp )


        if (eff_nfix > eff_nup) then
          !-----------------------------------------------------------------
          ! N acquisition happens at efficiency of eff_nfix, amplifying the 
          ! root N uptake by (eff_nfix / eff_nup). The difference of N acquisition
          ! is satisfied by N fixation and implies a C consumption from the 
          ! belowground allocatable pool, counted towards root respiration.
          !-----------------------------------------------------------------
          dnacq = dnup * eff_nfix / eff_nup
          dnfix = dnacq - dnup
          dcfix = dnfix / eff_nfix

          if (dcfix > tile(lu)%plant(pft)%palcb%c%c12) then
            !-----------------------------------------------------------------
            ! Withdraw all what's available from temporary belowground allocatable pool 
            ! and use if for N fixation, then update N fixation and total N acquisition.
            !-----------------------------------------------------------------
            dcfix = tile(lu)%plant(pft)%palcb%c%c12
            dnfix = eff_nfix * dcfix
            dnacq = dnfix + dnup
          end if

          !-----------------------------------------------------------------
          ! add consumed C to root respiration. Don't do anything with N.
          !-----------------------------------------------------------------
          call cmv( &
            carbon(carbon = dcfix), &
            tile(lu)%plant(pft)%palcb%c, &
            tile_fluxes(lu)%plant(pft)%drroot &
            )

        else
          !-----------------------------------------------------------------
          ! N acquisition only via root uptake (calculated above)
          !-----------------------------------------------------------------
          dnfix = 0.0
          dcfix = 0.0

        end if  ! (eff_nfix > eff_nup)

      end if  ! (params_pft_plant(pft)%nfixer)

      ! Update soil inorganic N pools
      tile(lu)%soil%pno3%n14 = tile(lu)%soil%pno3%n14 - dnup * fno3
      tile(lu)%soil%pnh4%n14 = tile(lu)%soil%pnh4%n14 - dnup * (1.0 - fno3)

      !--------------------------------------------------------------------------
      ! Update N-uptake of this PFT. N-retranslocation is not considered
      ! N-uptake.
      !--------------------------------------------------------------------------
      ! daily
      tile_fluxes(lu)%plant(pft)%dnup%n14 = dnup
      tile_fluxes(lu)%plant(pft)%dnup_fix = dnfix  

      !--------------------------------------------------------------------------
      ! Total N acquisition (root N uptake plus N fixation) to labile pool
      !--------------------------------------------------------------------------
      call ncp( tile_fluxes(lu)%plant(pft)%dnup, tile(lu)%plant(pft)%plabl%n )
      call ncp( nitrogen(nitrogen = dnfix), tile(lu)%plant(pft)%plabl%n )

    end do pftloop

  end subroutine nuptake


  function calc_dnup( croot, ninorg ) result( nup )
    !/////////////////////////////////////////////////////////////////
    ! N uptake is saturating both with respect to the root biomass and 
    ! with respect to soil inorganic N.
    !-----------------------------------------------------------------
    ! arguments
    real, intent(in) :: croot        ! root biomass (gC/m2)
    real, intent(in) :: ninorg       ! total inorganic soil N (gN/m2)

    ! function return variable
    real :: nup

    nup = calc_vn( ninorg ) * croot / (params_nuptake%kc + croot)

  end function calc_dnup


  function calc_vn( ninorg ) result( vn )
    !/////////////////////////////////////////////////////////////////
    ! N uptake capacity as a saturating function of soil inorganic N.
    !-----------------------------------------------------------------
    ! arguments
    real, intent(in) :: ninorg       ! total inorganic soil N (gN/m2)

    ! function return variable
    real :: vn

    vn = params_nuptake%vmax * ninorg / (params_nuptake%kv + ninorg)

  end function calc_vn


  function calc_eff_nfix( soiltemp ) result( out_eff_nfix )
    !////////////////////////////////////////////////////////////////
    ! Calculates N fixation efficiency (dNfix/dCex) as a function of 
    ! soil temperature. The functional form is chosen to be equal to
    ! the (fitted) relationship between soil temperature and nitro-
    ! genase activity derived by Houlton et al. (2008), Nature. The 
    ! maximum efficiency is 0.21 gN/gC, which corresponds to the 
    ! inverse of the "minimum cost of N-fixation" of 4.8 gC/gN given 
    ! Gutschik (1981). At higher and lower temperature, efficiency de-
    ! clines to zero.
    !-----------------------------------------------------------------
    ! arguments    
    real, intent(in) :: soiltemp

    ! local variables
    real, parameter :: norm = 1.25  ! used to normalise efficiency (nitrogenase activiy) function to 1 at optimal soil temperature (see Houlton et al., 2008)

    ! function return variable
    real :: out_eff_nfix             ! function return variable

    out_eff_nfix = 1.0 / params_nuptake%minimumcostfix * norm * &
      exp( params_nuptake%a_param_fix + &
           params_nuptake%b_param_fix * soiltemp * ( 1.0 - 0.5 * soiltemp / params_nuptake%fixoptimum ) )

  end function calc_eff_nfix


  subroutine getpar_modl_nuptake()
    !////////////////////////////////////////////////////////////////
    ! Subroutine reads nuptake module-specific parameters 
    ! from input file
    !----------------------------------------------------------------
    use md_interface_cnmodel, only: myinterface

    params_nuptake%kc   = myinterface%params_calib%nuptake_kc
    params_nuptake%kv   = myinterface%params_calib%nuptake_kv
    params_nuptake%vmax = myinterface%params_calib%nuptake_vmax

    ! Minimum cost of N-fixation is 4.8 gC/gN, value from Gutschik (1981)
    params_nuptake%minimumcostfix = myinterface%params_calib%minimumcostfix

    ! Optimum temperature for N fixation. Taken to be equal to optimum temperature for 
    ! nitrogenase activity as given by Houlton et al. (2008), Nature: 25.15+-0.66 
    params_nuptake%fixoptimum = 25.15  ! myinterface%params_calib%fixoptimum
 
    ! shape parameter of the N fixation function taken to be equal to nitrogenase activity 
    ! function given Houlton et al. (2008), Nature: -3.62+-0.52 
    params_nuptake%a_param_fix = -3.62  ! myinterface%params_calib%a_param_fix

    ! shape parameter of the N fixation function taken to be equal to nitrogenase activity 
    ! function given Houlton et al. (2008), Nature: 0.27+-0.04 
    params_nuptake%b_param_fix = 0.27  ! myinterface%params_calib%b_param_fix

  end subroutine getpar_modl_nuptake  


end module md_nuptake_simpl
