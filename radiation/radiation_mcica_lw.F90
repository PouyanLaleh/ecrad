! radiation_mcica_lw.F90 - Monte-Carlo Independent Column Approximation longtwave solver
!
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! Modifications
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-04-22  R. Hogan  Store surface fluxes at all g-points
!   2017-07-12  R. Hogan  Call fast adding method if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_mcica_lw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave Monte Carlo Independent Column Approximation
  ! (McICA). This implementation performs a clear-sky and a cloudy-sky
  ! calculation, and then weights the two to get the all-sky fluxes
  ! according to the total cloud cover. This method reduces noise for
  ! low cloud cover situations, and exploits the clear-sky
  ! calculations that are usually performed for diagnostic purposes
  ! simultaneously. The cloud generator has been carefully written
  ! such that the stochastic cloud field satisfies the prescribed
  ! overlap parameter accounting for this weighting.
  subroutine solver_mcica_lw(nlev,istartcol,iendcol, &
       &  config, single_level, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    use radiation_io,   only           : nulerr, radiation_abort
    use radiation_config, only         : config_type
    use radiation_single_level, only   : single_level_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type
    use radiation_two_stream, only     : calc_ref_trans_lw, &
         &                               calc_no_scattering_transmittance_lw
    use radiation_adding_ica_lw, only  : adding_ica_lw, fast_adding_ica_lw, &
         &                               calc_fluxes_no_scattering_lw
    use radiation_lw_derivatives, only : calc_lw_derivatives_ica, modify_lw_derivatives_ica
    use radiation_cloud_generator, only: cloud_generator

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(istartcol:iendcol, config%n_g_lw, nlev) :: &
         &  od
    real(jprb), intent(in), dimension(istartcol:iendcol,config%n_g_lw_if_scattering, nlev) :: &  !!! main call should be checked
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(istartcol:iendcol,config%n_bands_lw,nlev)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(istartcol:iendcol, config%n_bands_lw_if_scattering, &      !!! main call should be checked
         &  nlev) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(istartcol:iendcol, config%n_g_lw,nlev+1) :: &
         &  planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension( istartcol:iendcol, config%n_g_lw) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Diffuse reflectance and transmittance for each layer in clear
    ! and all skies
    real(jprb), dimension(istartcol:iendcol, config%n_g_lw, nlev) :: ref_clear, trans_clear, reflectance, transmittance

    ! Emission by a layer into the upwelling or downwelling diffuse
    ! streams, in clear and all skies
    real(jprb), dimension(istartcol:iendcol,config%n_g_lw, nlev) :: source_up_clear, source_dn_clear, source_up, source_dn
    ! Fluxes per g point
    real(jprb), dimension(istartcol:iendcol, config%n_g_lw, nlev+1) :: flux_up, flux_dn
    real(jprb), dimension(istartcol:iendcol, config%n_g_lw, nlev+1) :: flux_up_clear, flux_dn_clear
   ! Combined gas+aerosol+cloud optical depth, single scattering
    ! albedo and asymmetry factor
    real(jprb), dimension(istartcol:iendcol,config%n_g_lw,nlev) :: od_total, ssa_total, g_total

    ! Combined scattering optical depth
    real(jprb) :: scat_od, scat_od_total(istartcol:iendcol,config%n_g_lw)

    ! Optical depth scaling from the cloud generator, zero indicating
    ! clear skies
    real(jprb), dimension(istartcol:iendcol,config%n_g_lw,nlev) :: od_scaling

    ! Modified optical depth after McICA scaling to represent cloud
    ! inhomogeneity
    real(jprb), dimension(istartcol:iendcol,config%n_g_lw) :: od_cloud_new

    ! Total cloud cover output from the cloud generator
    real(jprb), dimension(istartcol:iendcol) :: total_cloud_cover

    ! Identify clear-sky layers
    logical, dimension(istartcol:iendcol,nlev) :: is_clear_sky_layer

    ! Index of the highest cloudy layer
    integer, dimension(istartcol:iendcol) :: i_cloud_top

    ! Number of g points
    integer :: ng, jjlev

    ! Loop indices for level, column and g point
    integer :: jlev, jcol, jg
    logical :: F_div(istartcol:iendcol)
    logical :: S_div(istartcol:iendcol,nlev)
    logical :: NF_div(istartcol:iendcol)
    logical :: tcc_div(istartcol:iendcol)


    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',0,hook_handle)

    if (.not. config%do_clear) then
      write(nulerr,'(a)') '*** Error: longwave McICA requires clear-sky calculation to be performed'
      call radiation_abort()
    end if

    ng = config%n_g_lw
   
    ! Loop through columns
    do jcol = istartcol,iendcol
      F_div(jcol) = .true.
      NF_div(jcol) = .false.
      tcc_div(jcol) = .false.
    end do

    do jcol = istartcol,iendcol
      do jlev = 1 , nlev
        S_div(jcol,jlev) = .true.
      end do
    end do
      ! Clear-sky calculation
      if (config%do_lw_aerosol_scattering) then     
        ! Scattering case: first compute clear-sky reflectance,
        ! transmittance etc at each model level
        call calc_ref_trans_lw(ng*nlev, nlev,&
             &  od, ssa, g, &
             &  planck_hl(:,:,1:nlev), planck_hl(:,:,2:nlev+1), & ! &  planck_hl(:,:,1:jlev), planck_hl(:,:,2:jlev+1), &
             &  ref_clear, trans_clear, &
             &  source_up_clear, source_dn_clear,istartcol,iendcol, F_div, S_div)
        ! Then use adding method to compute fluxes
        call adding_ica_lw(ng, nlev, &
             &  ref_clear, trans_clear, source_up_clear, source_dn_clear, &
             &  emission(:,jcol), albedo(:,jcol), &
             &  flux_up_clear, flux_dn_clear,istartcol,iendcol, F_div)
      else
        ! Non-scattering case: use simpler functions for
        ! transmission and emission
        call calc_no_scattering_transmittance_lw(ng*nlev, od, &
             &  planck_hl(:,:,1:nlev), planck_hl(:,:,2:nlev+1), &
             &  trans_clear, source_up_clear, source_dn_clear,istartcol,iendcol)
        ! Ensure that clear-sky reflectance is zero since it may be
        ! used in cloudy-sky case
        ref_clear = 0.0_jprb
       !!! i_cloud_top = nlev+1
        ! Simpler down-then-up method to compute fluxes
        call calc_fluxes_no_scattering_lw(ng, nlev, &
             &  trans_clear, source_up_clear, source_dn_clear, &
             &  emission, albedo, &
             &  flux_up_clear, flux_dn_clear,istartcol,iendcol)       
      end if
      do jcol = istartcol,iendcol
        ! Sum over g-points to compute broadband fluxes
        flux%lw_up_clear(jcol,:) = sum(flux_up_clear(jcol,:,:),1)     ! flux_up_clear is input for other subroutines
        flux%lw_dn_clear(jcol,:) = sum(flux_dn_clear(jcol,:,:),1)     ! ! flux_dn_clear is input for other subroutines
        ! Store surface spectral downwelling fluxes
        flux%lw_dn_surf_clear_g(:,jcol) = flux_dn_clear(jcol,:,nlev+1)  
      ! Do cloudy-sky calculation; add a prime number to the seed in
      ! the longwave
      call cloud_generator(ng, nlev, config%i_overlap_scheme, &
           &  single_level%iseed(jcol) + 997, &
           &  config%cloud_fraction_threshold, &
           &  cloud%fraction(jcol,:), cloud%overlap_param(jcol,:), &
           &  config%cloud_inhom_decorr_scaling, cloud%fractional_std(jcol,:), &
           &  config%pdf_sampler, od_scaling(jcol,:,:), total_cloud_cover(jcol), &
           &  use_beta_overlap=config%use_beta_overlap, &
           &  use_vectorizable_generator=config%use_vectorizable_generator)
      
      ! Store total cloud cover
      flux%cloud_cover_lw(jcol) = total_cloud_cover(jcol)
        if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then   !!! Lable1: we assumed it is always true
        ! Total-sky calculation
          i_cloud_top(jcol) = nlev+1
          do jjlev = 1,nlev
            is_clear_sky_layer(jcol,jjlev) = .true.
          end do
          F_div(jcol) = .true.
          do jlev = 1,nlev
            ! Compute combined gas+aerosol+cloud optical properties
            if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then  !!!! !!! Lable2: we assumed it is always true
              is_clear_sky_layer(jcol,jlev) = .false.
              S_div(jcol,jlev) = .true.
              ! Get index to the first cloudy layer from the top
              if (i_cloud_top(jcol) > jlev) then
                i_cloud_top(jcol) = jlev
              end if
            else
              S_div(jcol,jlev) = .false.
            end if
         end do   !!! jlev loop
        else
           i_cloud_top(jcol) = nlev+1
           F_div(jcol) = .false.
        end if
      end do     !!! jcol loop
      
      do jlev = 1,nlev
            do jg = 1,ng
              do jcol = istartcol,iendcol
                if (F_div(jcol)) then
                   if (S_div(jcol,jlev)) then
                       od_cloud_new(jcol,jg) = od_scaling(jcol,jg,jlev) &
                       &  * od_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev)
                       od_total(jcol,jg,jlev)  = od(jcol,jg,jlev) + od_cloud_new(jcol,jg)
                       ssa_total(jcol,jg,jlev) = 0.0_jprb
                       g_total(jcol,jg,jlev)   = 0.0_jprb
                   end if
                end if
              end do
            end do

            if (config%do_lw_cloud_scattering) then
              ! Scattering case: calculate reflectance and
              ! transmittance at each model level

              if (config%do_lw_aerosol_scattering) then
                ! In single precision we need to protect against the
                ! case that od_total > 0.0 and ssa_total > 0.0 but
                ! od_total*ssa_total == 0 due to underflowi
               
                do jg = 1,ng
                  do jcol = istartcol,iendcol
                      if (F_div(jcol)) then
                          if (S_div(jcol,jlev)) then
                               if (od_total(jcol,jg,jlev) > 0.0_jprb) then
                                    scat_od_total(jcol,jg) = ssa(jcol,jg,jlev)*od(jcol,jlev,jg) &
                                    &     + ssa_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev) &
                                    &     *  od_cloud_new(jcol,jg)
                                    ssa_total(jcol,jg,jlev) = scat_od_total(jcol,jg) / od_total(jcol,jg,jlev)

                                   if (scat_od_total(jcol,jg) > 0.0_jprb) then
                                      g_total(jcol,jg,jlev) = (g(jcol,jg,jlev)*ssa(jcol,jg,jlev)*od(jcol,jlev,jg) &
                                      &     +   g_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev) &
                                      &     * ssa_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev) &
                                      &     *  od_cloud_new(jcol,jg)) &
                                      &     / scat_od_total(jcol,jg)
                                  end if
                               end if
                        end if
                    end if
                  end do
                end do

              else
                do jg = 1,ng
                  do jcol = istartcol,iendcol
                        if (F_div(jcol)) then
                           if (S_div(jcol,jlev)) then
                               if (od_total(jcol,jg,jlev) > 0.0_jprb) then
                                    scat_od = ssa_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev) &
                                        &         * od_cloud_new(jcol,jg)
                                    ssa_total(jcol,jg,jlev) = scat_od / od_total(jcol,jg,jlev)
                                    if (scat_od > 0.0_jprb) then
                                        g_total(jcol,jg,jlev) = g_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev) &
                                        &     * ssa_cloud(jcol,config%i_band_from_reordered_g_lw(jg),jlev) &
                                        &     *  od_cloud_new(jcol,jg) / scat_od
                                    end if
                                end if
                           end if
                       end if
                  end do
                end do

              end if
            end if
      end do
      if (config%do_lw_cloud_scattering) then

            
              ! Compute cloudy-sky reflectance, transmittance etc at
              ! each model leve
                       call calc_ref_trans_lw(ng, nlev,&
                       &  od_total, ssa_total, g_total, &
                       &  planck_hl, planck_hl, &
                       &  reflectance, transmittance, &
                       &  source_up, source_dn, istartcol,iendcol,F_div,S_div)
      end if
            
            !!!  else
                 !!! if (total_cloud_cover(jcol) >= config%cloud_fraction_threshold) then
                  !!!  if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then       
              ! No-scattering case: use simpler functions for
              ! transmission and emission
                !!!        call calc_no_scattering_transmittance_lw(ng, od_total, &
                !!!        &  planck_hl(:,jlev,:), planck_hl(:,jlev+1, :), &
                !!!        &  transmittance(:,:,jlev), source_up(:,:,jlev), source_dn(:,:,jlev), istartcol,iendcol,F_div,S_div)
                !!!    end if
                !!! end if
      
          
      do jlev = 1,nlev
        do jg = 1,ng
          do jcol = istartcol, iendcol
            if (F_div(jcol)) then
              if (.not.(S_div(jcol,jlev))) then
            ! Clear-sky layer: copy over clear-sky values
                  reflectance(jcol,jg,jlev) = ref_clear(jcol,jg,jlev)
                  transmittance(jcol,jg,jlev) = trans_clear(jcol,jg,jlev)
                  source_up(jcol,jg,jlev) = source_up_clear(jcol,jg,jlev)
                  source_dn(jcol,jg,jlev) = source_dn_clear(jcol,jg,jlev)

              end if
            end if
          end do
        end do
      end do

     
        
        if (config%do_lw_aerosol_scattering) then
          ! Use adding method to compute fluxes for an overcast sky,
          ! allowing for scattering in all layers
          call adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission, albedo, &
               &  flux_up, flux_dn, istartcol,iendcol,F_div)
        else if (config%do_lw_cloud_scattering) then
          ! Use adding method to compute fluxes but optimize for the
          ! presence of clear-sky layers
        
          call fast_adding_ica_lw(ng, nlev, reflectance, transmittance, source_up, source_dn, &
               &  emission, albedo, &
               &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
               &  flux_up, flux_dn, istartcol,iendcol,F_div)
        else
         
          ! Simpler down-then-up method to compute fluxes
        !!!  call calc_fluxes_no_scattering_lw(ng, nlev, &
        !!!       &  transmittance, source_up, source_dn, emission, albedo, &
        !!!       &  flux_up, flux_dn,istartcol,iendcol)
        end if
        
      do jcol = istartcol,iendcol
        if (F_div(jcol)) then
        ! Store overcast broadband fluxes
           flux%lw_up(jcol,:) = sum(flux_up(jcol,:,:),1)
           flux%lw_dn(jcol,:) = sum(flux_dn(jcol,:,:),1)

        ! Cloudy flux profiles currently assume completely overcast
        ! skies; perform weighted average with clear-sky profile
           flux%lw_up(jcol,:) =  total_cloud_cover(jcol) *flux%lw_up(jcol,:) &
               &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_up_clear(jcol,:)
           flux%lw_dn(jcol,:) =  total_cloud_cover(jcol) *flux%lw_dn(jcol,:) &
               &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_clear(jcol,:)
        ! Store surface spectral downwelling fluxes
           flux%lw_dn_surf_g(:,jcol) = total_cloud_cover(jcol)*flux_dn(jcol,:,nlev+1) &
               &  + (1.0_jprb - total_cloud_cover(jcol))*flux%lw_dn_surf_clear_g(:,jcol)
         end if
      end do
        ! Compute the longwave derivatives needed by Hogan and Bozzo
        ! (2015) approximate radiation update scheme
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, transmittance, flux_up(:,:,nlev+1), &     !!! I removed jcol as an argument
               &                       flux%lw_derivatives, istartcol, iendcol, F_div)
          do jcol = istartcol,iendcol
            if (total_cloud_cover(jcol) < 1.0_jprb - config%cloud_fraction_threshold) then
               tcc_div(jcol) = .true.
            else
               tcc_div(jcol) = .false.
            end if
          end do
            ! Modify the existing derivative with the contribution from the clear sky
            call modify_lw_derivatives_ica(ng, nlev, trans_clear, flux_up_clear(:,:,nlev+1), &  !!! I removed jcol as an argument
                 &                         1.0_jprb-total_cloud_cover, flux%lw_derivatives, istartcol,iendcol, F_div,tcc_div)
        end if

       do jcol = istartcol, iendcol
         if (.not.(F_div(jcol))) then
        ! No cloud in profile and clear-sky fluxes already
        ! calculated: copy them over
            flux%lw_up(jcol,:) = flux%lw_up_clear(jcol,:)
            flux%lw_dn(jcol,:) = flux%lw_dn_clear(jcol,:)
            flux%lw_dn_surf_g(:,jcol) = flux%lw_dn_surf_clear_g(:,jcol)
            NF_div(jcol) = .true.
          end if
        end do
        if (config%do_lw_derivatives) then
          call calc_lw_derivatives_ica(ng, nlev, trans_clear, flux_up_clear(:,:,nlev+1), &
               &                       flux%lw_derivatives, istartcol, iendcol, NF_div)
 
        end if

    if (lhook) call dr_hook('radiation_mcica_lw:solver_mcica_lw',1,hook_handle)
    
  end subroutine solver_mcica_lw

end module radiation_mcica_lw
