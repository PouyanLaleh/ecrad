! radiation_adding_ica_sw.F90 - Shortwave adding method in independent column approximation
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
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_sw

  public

contains

  subroutine adding_ica_sw(ng, nlev, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct, istartcol, iendcol, SAH_mask, F_div)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in) :: ng ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels
    integer, intent(in) :: istartcol, iendcol

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! Cosine of the solar zenith angle
    real(jprb), intent(in),  dimension(istartcol:iendcol)         :: cos_sza

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: reflectance, transmittance

    ! Fraction of direct-beam solar radiation entering the top of a
    ! layer that is reflected back up or scattered forward into the
    ! diffuse stream at the base of the layer
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: ref_dir, trans_dir_diff

    ! Direct transmittance, i.e. fraction of direct beam that
    ! penetrates a layer without being scattered or absorbed
    real(jprb), intent(in),  dimension(istartcol:iendcol, ng, nlev)   :: trans_dir_dir
    
    logical, intent(in) :: SAH_mask(istartcol:iendcol)
    logical, intent(in) :: F_div(istartcol:iendcol)

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(istartcol:iendcol, ng, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(istartcol:iendcol, ng, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), dimension(istartcol:iendcol, ng, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(istartcol:iendcol, ng, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jg,  jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',0,hook_handle)

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    do jg = 1,ng
      do jcol = istartcol,iendcol
         if (SAH_mask(jcol)) then
           if(F_div(jcol)) then
           !call ftrace_region_begin(" ---->  radiation_adding_ica_sw:adding_ica_sw  <----")
             flux_dn_direct(jcol,jg,1) = incoming_toa(jcol,jg)
           end if
         end if
      end do
    end do
    do jlev = 1,nlev
      do jg = 1, ng
        do jcol = istartcol,iendcol
          if (SAH_mask(jcol)) then
            if(F_div(jcol)) then
               flux_dn_direct(jcol,jg,jlev+1) = flux_dn_direct(jcol,jg,jlev)*trans_dir_dir(jcol,jg,jlev)
            end if
          end if
        end do
      end do
    end do
    do jcol = istartcol,iendcol
      do jg = 1, ng
        if (SAH_mask(jcol)) then
          if(F_div(jcol)) then
            albedo(jcol,jg,nlev+1) = albedo_surf_diffuse(jcol,jg)
           ! At the surface, the direct solar beam is reflected back into the
           ! diffuse stream
            source(jcol,jg,nlev+1) = albedo_surf_direct(jcol,jg) * flux_dn_direct(jcol,jg,nlev+1) * cos_sza(jcol)
          ! call ftrace_region_end(" ---->  radiation_adding_ica_sw:adding_ica_sw  <----")
          end if
        end if
       end do
    end do

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
      do jg = 1,ng
        do jcol = istartcol,iendcol
          if (SAH_mask(jcol)) then
            if(F_div(jcol)) then
              ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
              inv_denominator(jcol,jg,jlev) = 1.0_jprb / (1.0_jprb-albedo(jcol,jg,jlev+1)*reflectance(jcol,jg,jlev))
              ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
              albedo(jcol,jg,jlev) = reflectance(jcol,jg,jlev) + transmittance(jcol,jg,jlev) * transmittance(jcol,jg,jlev) &
              &                                     * albedo(jcol,jg,jlev+1) * inv_denominator(jcol,jg,jlev)
              ! Shonk & Hogan (2008) Eq 11:
              source(jcol,jg,jlev) = ref_dir(jcol,jg,jlev)*flux_dn_direct(jcol,jg,jlev) &
                &  + transmittance(jcol,jg,jlev)*(source(jcol,jg,jlev+1) &
                &        + albedo(jcol,jg,jlev+1)*trans_dir_diff(jcol,jg,jlev)*flux_dn_direct(jcol,jg,jlev)) &
                &  * inv_denominator(jcol,jg,jlev)
            end if
          end if
       end do
      end do
    end do
    do jg = 1,ng
      do jcol = istartcol,iendcol
        if (SAH_mask(jcol)) then
          if(F_div(jcol)) then
          ! At top-of-atmosphere there is no diffuse downwelling radiation
            flux_dn_diffuse(jcol,jg,1) = 0.0_jprb

          ! At top-of-atmosphere, all upwelling radiation is due to
          ! scattering by the direct beam below that level
            flux_up(jcol,jg,1) = source(jcol,jg,1)
          end if
        end if
      end do
    end do

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
      do jg = 1,ng
        do jcol = istartcol,iendcol
          if (SAH_mask(jcol)) then
            if(F_div(jcol)) then
              ! Shonk & Hogan (2008) Eq 14 (after simplification):
               flux_dn_diffuse(jcol,jg,jlev+1) &
               &  = (transmittance(jcol,jg,jlev)*flux_dn_diffuse(jcol,jg,jlev) &
               &     + reflectance(jcol,jg,jlev)*source(jcol,jg,jlev+1) &
               &     + trans_dir_diff(jcol,jg,jlev)*flux_dn_direct(jcol,jg,jlev)) * inv_denominator(jcol,jg,jlev)
               ! Shonk & Hogan (2008) Eq 12:
               flux_up(jcol,jg,jlev+1) = albedo(jcol,jg,jlev+1)*flux_dn_diffuse(jcol,jg,jlev+1) &
               &            + source(jcol,jg,jlev+1)
               flux_dn_direct(jcol,jg,jlev) = flux_dn_direct(jcol,jg,jlev)*cos_sza(jcol)
            end if
          end if
        end do
      end do
    end do
    
    do jg = 1,ng
      do jcol = istartcol,iendcol
        if (SAH_mask(jcol)) then
          if(F_div(jcol)) then
            flux_dn_direct(jcol,jg,nlev+1) = flux_dn_direct(jcol,jg,nlev+1)*cos_sza(jcol)
          end if
        end if
      end do
    end do

    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',1,hook_handle)

  end subroutine adding_ica_sw

end module radiation_adding_ica_sw
