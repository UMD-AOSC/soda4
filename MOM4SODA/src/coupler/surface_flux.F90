!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module surface_flux_mod
!
! <CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov">GFDL </CONTACT>
!
!
! <OVERVIEW>
!  Driver program for the calculation of fluxes on the exchange grids. 
! </OVERVIEW>
!
! <DESCRIPTION>
!
! </DESCRIPTION>
!
! ============================================================================

use             fms_mod, only: FATAL, close_file, mpp_pe, mpp_root_pe, write_version_number
use             fms_mod, only: file_exist, check_nml_error, open_namelist_file, stdlog
use   monin_obukhov_mod, only: mo_drag, mo_profile
use  sat_vapor_pres_mod, only: escomp, descomp
use       constants_mod, only: cp_air, hlv, stefan, rdgas, rvgas, grav, vonkarm
use             mpp_mod, only: input_nml_file

implicit none
private

! ==== public interface ======================================================
public  surface_flux
! ==== end of public interface ===============================================

! <INTERFACE NAME="surface_flux">
!   <OVERVIEW>
!   For the calculation of fluxes on the exchange grids. 
!   </OVERVIEW>
!   <DESCRIPTION>
!   For the calculation of fluxes on the exchange grids. 
!   </DESCRIPTION>
!
!  <IN NAME="t_atm" TYPE="real, dimension(:)" UNITS="Kelvin">
!  Air temp lowest atmospheric level.  
!  </IN>
!  <IN NAME="q_atm" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Mixing ratio at lowest atmospheric level (kg/kg).  
!  </IN>
!  <IN NAME="u_atm" TYPE="real, dimension(:)" UNITS="m/s">
!  Zonal wind velocity at lowest atmospheric level.       
!  </IN>
!  <IN NAME="v_atm" TYPE="real, dimension(:)" UNITS="m/s">
!  Meridional wind velocity at lowest atmospheric level.    
!  </IN>
!  <IN NAME="p_atm" TYPE="real, dimension(:)" UNITS="Pascal">
!  Pressure lowest atmospheric level.    
!  </IN>
!  <IN NAME="z_atm" TYPE="real, dimension(:)" UNITS="m" >
!  Height lowest atmospheric level. 
!  </IN>
!  <IN NAME="p_surf" TYPE="real, dimension(:)" UNITS="Pascal">
!   Pressure at the earth's surface
!  </IN>
!  <IN NAME="t_surf" TYPE="real, dimension(:)" UNITS="Kelvin">
!   Temp at the earth's surface
!  </IN>
!  <IN NAME="t_ca" TYPE="real, dimension(:)" UNITS="Kelvin">
!   Air temp at the canopy 
!  </IN>
!  <OUT NAME="q_surf" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Mixing ratio at earth surface (kg/kg). 
!  </OUT>
!  <IN NAME="u_surf" TYPE="real, dimension(:)" UNITS="m/s">
!  Zonal wind velocity at earth surface.   
!  </IN>
!  <IN NAME="v_surf" TYPE="real, dimension(:)" UNITS="m/s">
!  Meridional wind velocity at earth surface. 
!  </IN>
!  <IN NAME="rough_mom" TYPE="real, dimension(:)" UNITS="m">
!  Momentum roughness length
!  </IN>
!  <IN NAME="rough_heat" TYPE="real, dimension(:)" UNITS="m">
!  Heat roughness length
!  </IN>
!  <IN NAME="rough_moist" TYPE="real, dimension(:)" UNITS="m">
! <Moisture roughness length
!  </IN>
!  <IN NAME="rough_scale" TYPE="real, dimension(:)" UNITS="dimensionless" >
!  Scale factor used to topographic roughness calculation
!  </IN>
!  <IN NAME="gust" TYPE="real, dimension(:)"  UNITS="m/s">
!   Gustiness factor 
!  </IN>
!  <OUT NAME="flux_t" TYPE="real, dimension(:)" UNITS="W/m^2">
!  Sensible heat flux 
!  </OUT>
!  <OUT NAME="flux_q" TYPE="real, dimension(:)" UNITS="kg/(m^2 s)">
!  Evaporative water flux 
!  </OUT>
!  <OUT NAME="flux_r" TYPE="real, dimension(:)" UNITS="W/m^2">
!  Radiative energy flux 
!  </OUT>
!  <OUT NAME="flux_u" TYPE="real, dimension(:)" UNITS="Pa">
!  Zonal momentum flux 
!  </OUT>
!  <OUT NAME="flux_v" TYPE="real, dimension(:)" UNITS="Pa">
! Meridional momentum flux 
!  </OUT>
!  <OUT NAME="cd_m" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Momentum exchange coefficient 
!  </OUT>
!  <OUT NAME="cd_t" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Heat exchange coefficient 
!  </OUT>
!  <OUT NAME="cd_q" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Moisture exchange coefficient 
!  </OUT>
!  <OUT NAME="w_atm" TYPE="real, dimension(:)" UNITS="m/s">
!  Absolute wind at the lowest atmospheric level
!  </OUT>
!  <OUT NAME="u_star" TYPE="real, dimension(:)" UNITS="m/s">
!  Turbulent velocity scale 
!  </OUT>
!  <OUT NAME="b_star" TYPE="real, dimension(:)" UNITS="m/s^2">
!  Turbulent buoyant scale
!  </OUT>
!  <OUT NAME="q_star" TYPE="real, dimension(:)" UNITS="dimensionless">
!  Turbulent moisture scale
!  </OUT>
!  <OUT NAME="dhdt_surf" TYPE="real, dimension(:)" UNITS="(W/m^2)/K">
!  Sensible heat flux temperature sensitivity
!  </OUT>
!  <OUT NAME="dedt_surf" TYPE="real, dimension(:)" UNITS="1/K">
!   Moisture flux temperature sensitivity 
!  </OUT>
!  <OUT NAME="dedq_surf" TYPE="real, dimension(:)" UNITS="(kg/m^2)/s">
!  Moisture flux humidity sensitivity  
!  </OUT>
!  <OUT NAME="drdt_surf" TYPE="real, dimension(:)" UNITS="(W/m^2)/K">
!  Radiative energy flux temperature sensitivity 
!  </OUT>
!  <OUT NAME="dhdt_atm" TYPE="real, dimension(:)" UNITS="(W/m^2)/K">
!  Derivative of sensible heat flux over temp at the lowest atmos level.
!  </OUT>
!  <OUT NAME="dedq_atm" TYPE="real, dimension(:)" UNITS="(kg/m^2/sec)/K">
!  Derivative of water vapor flux over temp at the lowest atmos level.
!  </OUT>
!  <OUT NAME="dtaudu_atm" TYPE="real, dimension(:)" UNITS="Pa/(m/s)">
!  Derivative of zonal wind stress w.r.t the lowest level zonal 
!  wind speed of the atmos
!  </OUT>
!  <OUT NAME="dtaudv_atm" TYPE="real, dimension(:)" UNITS="Pa/(m/s)">
!  Derivative of meridional wind stress w.r.t the lowest level meridional 
!  wind speed of the atmos
!  </OUT>
!  <OUT NAME="dt" TYPE="real">
!  Time step (it is not used presently)
!  </OUT>
!  <IN NAME="land" TYPE="logical, dimension(:)">
!  Indicates where land exists (true if exchange cell is on land). 
!  </IN>
!  <IN NAME="seawater" TYPE="logical, dimension(:)">
!  Indicates where liquid ocean water exists 
!  (true if exchange cell is on liquid ocean water). 
!  </IN>
!  <IN NAME="avail" TYPE="logical, dimension(:)">
!  True where the exchange cell is active.  
!  </IN>


interface surface_flux
!    module procedure surface_flux_0d
    module procedure surface_flux_1d
    module procedure surface_flux_2d  
end interface
! </INTERFACE>

!-----------------------------------------------------------------------

character(len=*), parameter :: version = '$Id: surface_flux.F90,v 20.0 2013/12/13 23:27:45 fms Exp $'
character(len=*), parameter :: tagname = '$Name: tikal $'
   
logical :: do_init = .true.

real, parameter :: d622   = rdgas/rvgas
real, parameter :: d378   = 1.-d622
real, parameter :: hlars  = hlv/rvgas
real, parameter :: gcp    = grav/cp_air
real, parameter :: kappa  = rdgas/cp_air
real            :: d608   = d378/d622
      ! d608 set to zero at initialization if the use of 
      ! virtual temperatures is turned off in namelist
      
      
! ---- namelist with default values ------------------------------------------
! <NAMELIST NAME="surface_flux_nml">
!   <DATA NAME="no_neg_q"  TYPE="logical"  DEFAULT=".false.">
!    If q_atm_in (specific humidity) is negative (because of numerical truncation),  
!    then override with 0. 
!   </DATA>
!   <DATA NAME="use_virtual_temp"  TYPE="logical"  DEFAULT=".true.">
!    If true, use virtual potential temp to calculate the stability of the surface layer.
!    if false, use potential temp.
!   </DATA>
!   <DATA NAME="alt_gustiness"  TYPE="logical"  DEFAULT=".false.">
!   An alternative formulation for gustiness calculation. 
!   A minimum bound on the wind speed used influx calculations, with the bound 
!   equal to gust_const 
!   </DATA>
!   <DATA NAME="old_dtaudv"  TYPE="logical"  DEFAULT=".false.">
!   The derivative of surface wind stress w.r.t. the zonal wind and
!   meridional wind are approximated by the same tendency.
!   </DATA>
!   <DATA NAME="use_mixing_ratio"  TYPE="logical"  DEFAULT=".false.">
!   An option to provide capability to run the Manabe Climate form of the 
!   surface flux (coded for legacy purposes). 
!   </DATA>
!   <DATA NAME="gust_const"  TYPE=""  DEFAULT="1.0">
!    Constant for alternative gustiness calculation.
!   </DATA>
!   <DATA NAME="gust_min"  TYPE=""  DEFAULT="0.0">
!    Minimum gustiness used when alt_gustiness = false.
!   </DATA>
!   <DATA NAME="ncar_ocean_flux"  TYPE="logical"  DEFAULT=".false.">
!    Use NCAR climate model turbulent flux calculation described by
!    Large and Yeager, NCAR Technical Document, 2004
!   </DATA>
!   <DATA NAME="coare4_ocean_flux" TYPE="logical"  DEFAULT=".false.">
!    Use COARE4 ocean fluxes calculations which is vectorized version 
!    of COARE3 code (Fairall et al, 2003) with modification based on the 
!    CLIMODE, MBL and CBLAST experiments (Edson et al., 2011). 
!    The cool skin option is retained but warm layer and surface wave options removed. 
!   </DATA>
!   <DATA NAME="ncar_ocean_flux_orig"  TYPE="logical"  DEFAULT=".false.">
!    Use NCAR climate model turbulent flux calculation described by
!    Large and Yeager, NCAR Technical Document, 2004, using the original
!    GFDL implementation, which contains a bug in the specification of 
!    the exchange coefficient for the sensible heat.  This option is available
!    for legacy purposes, and is not recommended for new experiments.   
!   </DATA>
!   <DATA NAME="raoult_sat_vap"  TYPE="logical"  DEFAULT=".false.">
!    Reduce saturation vapor pressures to account for seawater salinity.
!   </DATA>
! </NAMELIST>

logical :: no_neg_q              = .false.  ! for backwards compatibility
logical :: use_virtual_temp      = .true. 
logical :: alt_gustiness         = .false.
logical :: old_dtaudv            = .false.
logical :: use_mixing_ratio      = .false.
real    :: gust_const            =  1.0
real    :: gust_min              =  0.0
logical :: ncar_ocean_flux       = .false.
logical :: ncar_ocean_flux_orig  = .false. ! for backwards compatibility 
logical :: coare4_ocean_flux     = .false.
logical :: raoult_sat_vap        = .false.
logical :: do_simple             = .false.


namelist /surface_flux_nml/ no_neg_q,             &
                            use_virtual_temp,     &
                            alt_gustiness,        &
                            gust_const,           &
                            gust_min,             &
                            old_dtaudv,           &
                            use_mixing_ratio,     &
                            ncar_ocean_flux,      &
                            ncar_ocean_flux_orig, &
                            coare4_ocean_flux,    &
                            raoult_sat_vap,       &
                            do_simple       
   


contains


! ============================================================================
! <SUBROUTINE NAME="surface_flux_1d" INTERFACE="surface_flux">
!  <IN NAME="t_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="q_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="u_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="v_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="p_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="z_atm" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="p_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="t_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="t_ca" TYPE="real, dimension(:)"> </IN>
!  <OUT NAME="q_surf" TYPE="real, dimension(:)"> </OUT>
!  <IN NAME="u_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="v_surf" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_mom" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_heat" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_moist" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="rough_scale" TYPE="real, dimension(:)"> </IN>
!  <IN NAME="gust" TYPE="real, dimension(:)"> </IN>
!  <OUT NAME="flux_t" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="flux_q" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="flux_r" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="flux_u" TYPE="real, dimension(:)"></OUT>
!  <OUT NAME="flux_v" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="cd_m" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="cd_t" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="cd_q" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="w_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="u_star" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="b_star" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="q_star" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dhdt_surf" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dedt_surf" TYPE="real, dimension(:)"></OUT>
!  <OUT NAME="dedq_surf" TYPE="real, dimension(:)"></OUT>
!  <OUT NAME="drdt_surf" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dhdt_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dedq_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dtaudu_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dtaudv_atm" TYPE="real, dimension(:)"> </OUT>
!  <OUT NAME="dt" TYPE="real"> </OUT>
!  <IN NAME="land" TYPE="logical, dimension(:)"> </IN>
!  <IN NAME="seawater" TYPE="logical, dimension(:)"> </IN>
!  <IN NAME="avail" TYPE="logical, dimension(:)"> </IN>


!<PUBLICROUTINE INTERFACE="surface_flux">
subroutine surface_flux_1d (                                           &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t, flux_q, flux_r, flux_u, flux_v,                           &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dt,        land,      seawater,     avail  )
!</PUBLICROUTINE>
!  slm Mar 28 2002 -- remove agument drag_q since it is just cd_q*wind
! ============================================================================
  ! ---- arguments -----------------------------------------------------------
  logical, intent(in), dimension(:) :: land,  seawater, avail
  real, intent(in),  dimension(:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,  &
       rough_mom, rough_heat, rough_moist,  rough_scale, gust
  real, intent(out), dimension(:) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q
  real, intent(inout), dimension(:) :: q_surf
  real, intent(in) :: dt

  ! ---- local constants -----------------------------------------------------
  ! temperature increment and its reciprocal value for comp. of derivatives
  real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp

  ! ---- local vars ----------------------------------------------------------
  real, dimension(size(t_atm(:))) ::                          &
       thv_atm,  th_atm,   tv_atm,    thv_surf,            &
       e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
       t_surf0,  t_surf1,  u_dif,     v_dif,               &
       rho_drag, drag_t,    drag_m,   drag_q,    rho,      &
       q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust

  integer :: i, nbad


  if (do_init) call surface_flux_init

  !---- use local value of surf temp ----

  t_surf0 = 200.   !  avoids out-of-bounds in es lookup 
  where (avail)
     where (land)
        t_surf0 = t_ca
     elsewhere
        t_surf0 = t_surf
     endwhere
  endwhere

  t_surf1 = t_surf0 + del_temp

  call escomp ( t_surf0, e_sat  )  ! saturation vapor pressure
  call escomp ( t_surf1, e_sat1 )  ! perturbed  vapor pressure

  if(use_mixing_ratio) then
    ! surface mixing ratio at saturation
    q_sat   = d622*e_sat /(p_surf-e_sat )  
    q_sat1  = d622*e_sat1/(p_surf-e_sat1)  
  elseif(do_simple) then                  !rif:(09/02/09)
    q_sat   = d622*e_sat / p_surf
    q_sat1  = d622*e_sat1/ p_surf   
  else
    ! surface specific humidity at saturation
    q_sat   = d622*e_sat /(p_surf-d378*e_sat )  
    q_sat1  = d622*e_sat1/(p_surf-d378*e_sat1)     
  endif

  ! initilaize surface air humidity according to surface type
  where (land)
     q_surf0 = q_surf ! land calculates it
  elsewhere
     q_surf0 = q_sat  ! everything else assumes saturated sfc humidity
  endwhere

  if (raoult_sat_vap) where (seawater) q_surf0 = 0.98 * q_surf0

  ! check for negative atmospheric humidities
  where(avail) q_atm = q_atm_in
  if(no_neg_q) then
     where(avail .and. q_atm_in < 0.0) q_atm = 0.0
  endif

  ! generate information needed by monin_obukhov
  where (avail)
     p_ratio = (p_surf/p_atm)**kappa

     tv_atm  = t_atm  * (1.0 + d608*q_atm)     ! virtual temperature
     th_atm  = t_atm  * p_ratio                ! potential T, using p_surf as refernce
     thv_atm = tv_atm * p_ratio                ! virt. potential T, using p_surf as reference 
     thv_surf= t_surf0 * (1.0 + d608*q_surf0 ) ! surface virtual (potential) T
!     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

     u_dif = u_surf - u_atm                    ! velocity components relative to surface
     v_dif = v_surf - v_atm
  endwhere

  if(alt_gustiness) then
     do i = 1, size(avail)
        if (.not.avail(i)) cycle
        w_atm(i) = max(sqrt(u_dif(i)**2 + v_dif(i)**2), gust_const)
        ! derivatives of surface wind w.r.t. atm. wind components
        if(w_atm(i) > gust_const) then
           dw_atmdu(i) = u_dif(i)/w_atm(i)
           dw_atmdv(i) = v_dif(i)/w_atm(i)
        else
           dw_atmdu(i) = 0.0
           dw_atmdv(i) = 0.0
        endif
     enddo
  else
     if (gust_min > 0.0) then 
       where(avail)
         w_gust = max(gust,gust_min) ! minimum gustiness
       end where
     else
       where(avail)
         w_gust = gust
       end where
     endif  
           
     where(avail) 
        w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + w_gust*w_gust)
        ! derivatives of surface wind w.r.t. atm. wind components
        dw_atmdu = u_dif/w_atm
        dw_atmdv = v_dif/w_atm
     endwhere
  endif

  !  monin-obukhov similarity theory 
  call mo_drag (thv_atm, thv_surf, z_atm,                  &
       rough_mom, rough_heat, rough_moist, w_atm,          &
       cd_m, cd_t, cd_q, u_star, b_star, avail             )

  ! override with ocean fluxes from NCAR calculation
  if (ncar_ocean_flux .or. ncar_ocean_flux_orig) then
    call  ncar_ocean_fluxes (w_atm, th_atm, t_surf0, q_atm, q_surf0, z_atm, &
                             seawater, cd_m, cd_t, cd_q, u_star, b_star     )
  else
  ! override with ocean fluxes from COARE4 calculation
    if (coare4_ocean_flux) then 
      call  coare40vn_ocean_fluxes (w_atm, th_atm, t_surf0, q_atm, q_surf0, z_atm, &
                             seawater, cd_m, cd_t, cd_q, u_star, b_star     )		     
    endif
  end if

  where (avail)
     ! scale momentum drag coefficient on orographic roughness
     cd_m = cd_m*(log(z_atm/rough_mom+1)/log(z_atm/rough_scale+1))**2
     ! surface layer drag coefficients
     drag_t = cd_t * w_atm
     drag_q = cd_q * w_atm
     drag_m = cd_m * w_atm

     ! density
     rho = p_atm / (rdgas * tv_atm)  

     ! sensible heat flux
     rho_drag = cp_air * drag_t * rho
     flux_t = rho_drag * (t_surf0 - th_atm)  ! flux of sensible heat (W/m**2)
     dhdt_surf =  rho_drag                   ! d(sensible heat flux)/d(surface temperature)
     dhdt_atm  = -rho_drag*p_ratio           ! d(sensible heat flux)/d(atmos temperature)

     ! evaporation
     rho_drag  =  drag_q * rho
     flux_q    =  rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))

     where (land)
        dedq_surf = rho_drag
        dedt_surf = 0
     elsewhere
        dedq_surf = 0
        dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
     endwhere
        
     dedq_atm  = -rho_drag   ! d(latent heat flux)/d(atmospheric mixing ratio)

     q_star = flux_q / (u_star * rho)             ! moisture scale
     ! ask Chris and Steve K if we still want to keep this for diagnostics
     q_surf = q_atm + flux_q / (rho*cd_q*w_atm)   ! surface specific humidity

     ! upward long wave radiation
     flux_r    =   stefan*t_surf**4               ! (W/m**2)
     drdt_surf = 4*stefan*t_surf**3               ! d(upward longwave)/d(surface temperature)

     ! stresses
     rho_drag   = drag_m * rho
     flux_u     = rho_drag * u_dif   ! zonal      component of stress (Nt/m**2)
     flux_v     = rho_drag * v_dif   ! meridional component of stress 

  elsewhere
     ! zero-out un-available data in output only fields
     flux_t     = 0.0
     flux_q     = 0.0
     flux_r     = 0.0
     flux_u     = 0.0
     flux_v     = 0.0
     dhdt_surf  = 0.0
     dedt_surf  = 0.0
     dedq_surf  = 0.0
     drdt_surf  = 0.0
     dhdt_atm   = 0.0
     dedq_atm   = 0.0
     u_star     = 0.0
     b_star     = 0.0
     q_star     = 0.0
     q_surf     = 0.0
     w_atm      = 0.0
  endwhere

  ! calculate d(stress component)/d(atmos wind component)
  dtaudu_atm = 0.0
  dtaudv_atm = 0.0
  if (old_dtaudv) then
     where(avail)
        dtaudv_atm = -rho_drag
        dtaudu_atm = -rho_drag
     endwhere
  else
     where(avail)
        dtaudu_atm = -cd_m*rho*(dw_atmdu*u_dif + w_atm)
        dtaudv_atm = -cd_m*rho*(dw_atmdv*v_dif + w_atm)
     endwhere
  endif

end subroutine surface_flux_1d
! </SUBROUTINE>

!#######################################################################

subroutine surface_flux_0d (                                                 &
     t_atm_0,     q_atm_0,      u_atm_0,     v_atm_0,   p_atm_0, z_atm_0,    &
     p_surf_0,    t_surf_0,     t_ca_0,      q_surf_0,                       &
     u_surf_0,    v_surf_0,                                                  &
     rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0,        &
     flux_t_0,    flux_q_0,     flux_r_0,    flux_u_0,  flux_v_0,            &
     cd_m_0,      cd_t_0,       cd_q_0,                                      &
     w_atm_0,     u_star_0,     b_star_0,     q_star_0,                      &
     dhdt_surf_0, dedt_surf_0,  dedq_surf_0,  drdt_surf_0,                   &
     dhdt_atm_0,  dedq_atm_0,   dtaudu_atm_0, dtaudv_atm_0,                  &
     dt,          land_0,       seawater_0,  avail_0  )

  ! ---- arguments -----------------------------------------------------------
  logical, intent(in) :: land_0,  seawater_0, avail_0
  real, intent(in) ::                                                  &
       t_atm_0,     q_atm_0,      u_atm_0,     v_atm_0,                &
       p_atm_0,     z_atm_0,      t_ca_0,                              &
       p_surf_0,    t_surf_0,     u_surf_0,    v_surf_0,               &
       rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0
  real, intent(out) ::                                                 &
       flux_t_0,    flux_q_0,     flux_r_0,    flux_u_0,  flux_v_0,    &
       dhdt_surf_0, dedt_surf_0,  dedq_surf_0, drdt_surf_0,            &
       dhdt_atm_0,  dedq_atm_0,   dtaudu_atm_0,dtaudv_atm_0,           &
       w_atm_0,     u_star_0,     b_star_0,    q_star_0,               &
       cd_m_0,      cd_t_0,       cd_q_0
  real, intent(inout) :: q_surf_0
  real, intent(in)    :: dt

  ! ---- local vars ----------------------------------------------------------
  logical, dimension(1) :: land,  seawater, avail
  real, dimension(1) :: &
       t_atm,     q_atm,      u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       rough_mom, rough_heat, rough_moist,  rough_scale, gust
  real, dimension(1) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q
  real, dimension(1) :: q_surf


  avail = .true.

  t_atm(1)       = t_atm_0
  q_atm(1)       = q_atm_0
  u_atm(1)       = u_atm_0
  v_atm(1)       = v_atm_0
  p_atm(1)       = p_atm_0
  z_atm(1)       = z_atm_0
  t_ca(1)        = t_ca_0
  p_surf(1)      = p_surf_0
  t_surf(1)      = t_surf_0
  u_surf(1)      = u_surf_0
  v_surf(1)      = v_surf_0
  rough_mom(1)   = rough_mom_0
  rough_heat(1)  = rough_heat_0
  rough_moist(1) = rough_moist_0
  rough_scale(1) = rough_scale_0
  gust(1)        = gust_0
  q_surf(1)      = q_surf_0
  land(1)        = land_0
  seawater(1)    = seawater_0
  avail(1)       = avail_0

  call surface_flux_1d (                                                 &
       t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
       p_surf,    t_surf,     t_ca,      q_surf,                         &
       u_surf,    v_surf,                                                &
       rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
       flux_t, flux_q, flux_r, flux_u, flux_v,                           &
       cd_m,      cd_t,       cd_q,                                      &
       w_atm,     u_star,     b_star,     q_star,                        &
       dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
       dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
       dt,        land,      seawater, avail  )

  flux_t_0     = flux_t(1)
  flux_q_0     = flux_q(1)
  flux_r_0     = flux_r(1)
  flux_u_0     = flux_u(1)
  flux_v_0     = flux_v(1)
  dhdt_surf_0  = dhdt_surf(1)
  dedt_surf_0  = dedt_surf(1)
  dedq_surf_0  = dedq_surf(1)
  drdt_surf_0  = drdt_surf(1)
  dhdt_atm_0   = dhdt_atm(1)
  dedq_atm_0   = dedq_atm(1)
  dtaudu_atm_0 = dtaudu_atm(1)
  dtaudv_atm_0 = dtaudv_atm(1)
  w_atm_0      = w_atm(1)
  u_star_0     = u_star(1)
  b_star_0     = b_star(1)
  q_star_0     = q_star(1)
  q_surf_0     = q_surf(1)
  cd_m_0       = cd_m(1)
  cd_t_0       = cd_t(1)
  cd_q_0       = cd_q(1)

end subroutine surface_flux_0d

subroutine surface_flux_2d (                                           &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t,    flux_q,     flux_r,    flux_u,    flux_v,              &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dt,        land,       seawater,  avail  )

  ! ---- arguments -----------------------------------------------------------
  logical, intent(in), dimension(:,:) :: land,  seawater, avail
  real, intent(in),  dimension(:,:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       rough_mom, rough_heat, rough_moist, rough_scale, gust
  real, intent(out), dimension(:,:) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q
  real, intent(inout), dimension(:,:) :: q_surf
  real, intent(in) :: dt

  ! ---- local vars -----------------------------------------------------------
  integer :: j

  do j = 1, size(t_atm,2)
     call surface_flux_1d (                                           &
          t_atm(:,j),     q_atm_in(:,j),   u_atm(:,j),     v_atm(:,j),     p_atm(:,j),     z_atm(:,j),    &
          p_surf(:,j),    t_surf(:,j),     t_ca(:,j),      q_surf(:,j),                                   &
          u_surf(:,j),    v_surf(:,j),                                                                    &
          rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), rough_scale(:,j), gust(:,j),                 &
          flux_t(:,j),    flux_q(:,j),     flux_r(:,j),    flux_u(:,j),    flux_v(:,j),                   &
          cd_m(:,j),      cd_t(:,j),       cd_q(:,j),                                                     &
          w_atm(:,j),     u_star(:,j),     b_star(:,j),     q_star(:,j),                                  &
          dhdt_surf(:,j), dedt_surf(:,j),  dedq_surf(:,j),  drdt_surf(:,j),                               &
          dhdt_atm(:,j),  dedq_atm(:,j),   dtaudu_atm(:,j), dtaudv_atm(:,j),                              &
          dt,             land(:,j),       seawater(:,j),  avail(:,j)  )
  end do
end subroutine surface_flux_2d


! ============================================================================
!  Initialization of the surface flux module--reads the nml.     
!
subroutine surface_flux_init

! ---- local vars ----------------------------------------------------------
  integer :: unit, ierr, io

  ! read namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, surface_flux_nml, iostat=io)
      ierr = check_nml_error(io,'surface_flux_nml')
#else
  if ( file_exist('input.nml')) then
     unit = open_namelist_file ()
     ierr=1; 
     do while (ierr /= 0)
        read  (unit, nml=surface_flux_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'surface_flux_nml')
     enddo
10   call close_file (unit)
  endif
#endif

  ! write version number
  call write_version_number(version, tagname)

  unit = stdlog()
  if ( mpp_pe() == mpp_root_pe() )  write (unit, nml=surface_flux_nml)
  
  if(.not. use_virtual_temp) d608 = 0.0
  
  do_init = .false.
  
end subroutine surface_flux_init



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
! Original  code: GFDL.Climate.Model.Info
! Update Jul2007: GFDL.Climate.Model.Info (ch and ce exchange coeff bugfix)  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
                              cd, ch, ce, ustar, bstar       )
real   , intent(in)   , dimension(:) :: u_del, t, ts, q, qs, z
logical, intent(in)   , dimension(:) :: avail
real   , intent(inout), dimension(:) :: cd, ch, ce, ustar, bstar

  real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real :: cd_rt                                ! full drag coefficients @ z
  real :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real :: u, u10, tv, tstar, qstar, z0, xx, stab
  integer, parameter :: n_itts = 2
  integer               i, j

  if(ncar_ocean_flux_orig) then

      do i=1,size(u_del(:))
         if (avail(i)) then
             tv = t(i)*(1+0.608*q(i));
             u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
             u10 = u;                                                ! first guess 10m wind

             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
             cd_n10_rt = sqrt(cd_n10);
             ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
             stab = 0.5 + sign(0.5,t(i)-ts(i))
             ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c

             cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
             ch(i) = ch_n10;
             ce(i) = ce_n10;
             do j=1,n_itts                                           ! Monin-Obukhov iteration
                cd_rt = sqrt(cd(i));
                ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
                tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
                qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
                bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
                zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
                zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
                x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
                x2 = max(x2, 1.0);                                    ! undocumented NCAR
                x = sqrt(x2);

                if (zeta > 0) then
                    psi_m = -5*zeta;                                    ! L-Y eqn. 8c
                    psi_h = -5*zeta;                                    ! L-Y eqn. 8c
                else
                    psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
                    psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
                end if

                u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
                cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
                cd_n10_rt = sqrt(cd_n10);
                ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
                stab = 0.5 + sign(0.5,zeta)
                ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
                z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic

                xx = (log(z(i)/10)-psi_m)/vonkarm;
                cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
                xx = (log(z(i)/10)-psi_h)/vonkarm;
                ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !     10b (this code is wrong)  
                ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !     10c (this code is wrong)
             end do
         end if
      end do

  else

      do i=1,size(u_del(:))
         if (avail(i)) then
             tv = t(i)*(1+0.608*q(i));
             u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
             u10 = u;                                                ! first guess 10m wind

             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
             cd_n10_rt = sqrt(cd_n10);
             ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
             stab = 0.5 + sign(0.5,t(i)-ts(i))
             ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c

             cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
             ch(i) = ch_n10;
             ce(i) = ce_n10;
             do j=1,n_itts                                           ! Monin-Obukhov iteration
                cd_rt = sqrt(cd(i));
                ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
                tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
                qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
                bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
                zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
                zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
                x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
                x2 = max(x2, 1.0);                                    ! undocumented NCAR
                x = sqrt(x2);

                if (zeta > 0) then
                    psi_m = -5*zeta;                                    ! L-Y eqn. 8c
                    psi_h = -5*zeta;                                    ! L-Y eqn. 8c
                else
                    psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
                    psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
                end if

                u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
                cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
                cd_n10_rt = sqrt(cd_n10);
                ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
                stab = 0.5 + sign(0.5,zeta)
                ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
                z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic

                xx = (log(z(i)/10)-psi_m)/vonkarm;
                cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
                xx = (log(z(i)/10)-psi_h)/vonkarm;
                ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10b (corrected code)
                ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10c (corrected code)
             end do
         end if
      end do

  endif

end subroutine ncar_ocean_fluxes


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following COARE3 code (Fairall et al, 2003) with 
! modification based on the CLIMODE, MBL and CBLAST experiments 
! (Edson et al., 2011). The cool skin option is retained but warm layer 
! and surface wave options removed.           !
! Code : Senya Grodsky  (converted from MATLAB version)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
	subroutine coare40vn_ocean_fluxes (u, t, ts, Q, Qs, zu, avail, &
!     &                               Cd, Ch, Ce, ustar, bstar,n      )  ! Senya's test interface
                                    Cd, Ch, Ce, ustar, bstar       )     
!*********************************************************************

!function A=coare40vn(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi)
!
! Vectorized version of COARE3 code (Fairall et al, 2003) with 
! modification based on the CLIMODE, MBL and CBLAST experiments 
! (Edson et al., 2011). The cool skin option is retained but warm layer 
! and surface wave options removed. 
!
!********************************************************************
! An important component of this code is whether the inputed ts 
! represents the skin temperature of a near surface temperature.  
! How this variable is treated is determined by the jcool parameter:
! set jcool=1 if Ts is bulk ocean temperature (default),
!     jcool=0 if Ts is true ocean skin temperature. 
!********************************************************************

!jcool=0; !ts is treated as the true SST

! The code assumes u,t,rh,ts are vectors; 
! sensor heights zu,zt,zl, latitude lat, and PBL height zi are constants;
! air pressure P and radiation Rs,Rl may be vectors or constants. 
! Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
! available.  Input NaNs to indicate no data. Defaults should be set to 
! representative regional values if possible.
!
! Input:  
!
!     u = relative wind speed (m/s) at height zu(m)
!     t = bulk air temperature (degC) at height zt(m)
!    rh = relative humidity (!) at height zq(m)
!     P = surface air pressure (mb) (default = 1013)
!    ts = water temperature (degC) see jcool below
!    Rs = downward shortwave radiation (W/m^2) (default = 150) 
!    Rl = downward longwave radiation (W/m^2) (default = 370)
!   lat = latitude (default = +45 N)
!    zi = PBL height (m) (default = 600m)
!
! The user controls the output.  This is currently set as:
!
! Output:  A=[usr tau hsb hlb hbb hsbb tsr qsr zot zoq Cd Ch Ce  L zet dter
! dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10];
!
!  where
!
!   usr = friction velocity that includes gustiness (m/s)
!   tau = wind stress (N/m^2)
!   hsb = sensible heat flux into ocean (W/m^2)
!   hlb = latent heat flux into ocean (W/m^2)
!   hbb = buoyany flux into ocean (W/m^2)
!   hsbb = "sonic" buoyancy flux measured directly by sonic anemometer 
!   tsr = temperature scaling parameter (K)
!   qsr = specific humidity scaling parameter (g/Kg)
!   zot = thermal roughness length (m)
!   zoq = moisture roughness length (m)
!   Cd = wind stress transfer (drag) coefficient at height zu   
!   Ch = sensible heat transfer coefficient (Stanton number) at height zu   
!   Ce = latent heat transfer coefficient (Dalton number) at height zu
!    L = Obukhov length scale (m) 
!  zet = Monin-Obukhov stability parameter zu/L 
! dter = cool-skin temperature depression (degC)
! dqer = cool-skin humidity depression (degC)
!  tkt = cool-skin thickness (m)
!  Urf = wind speed at reference height (user can select height below)
!  Tfr = temperature at reference height
!  Qfr = specific humidity at reference height
! RHfr = relative humidity at reference height
! UrfN = neutral value of wind speed at reference height
!  Rnl = Upwelling IR radiation computed by COARE
!   Le = latent heat of vaporization
! rhoa = density of air
!   UN = neutral value of wind speed at zu
!  U10 = wind speed adjusted to 10 m
! UN10 = neutral value of wind speed at 10m
!Cdn_10 = neutral value of drag coefficient at 10m    
!Chn_10 = neutral value of Stanton number at 10m    
!Cen_10 = neutral value of Dalton number at 10m    
!

! Notes: 1) u is the relative wind speed, i.e., the magnitude of the
!           difference between the wind (at zu) and ocean surface current 
!           vectors.
!        2) Set jcool=0 in code if ts is true surface skin temperature,
!           otherwise ts is assumed the bulk temperature and jcool=1.
!        3) Set P=NaN to assign default value if no air pressure data 
!           available. 
!        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns 
!           default values to Rs, Rl so that cool skin option can be applied. 
!        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
!           and/or PBL height not given. 
!        6) The code to compute the heat flux caused by precipitation is 
!           included if rain data is available (default is no rain).
!        7) Code updates the cool-skin temperature depression dter and thickness
!           tkt during iteration loop for consistency.
!        8) Number of iterations set to nits = 6.

! Reference:
!
!  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!  Bulk parameterization of air sea fluxes: updates and verification for the 
!  COARE algorithm, J. Climate, 16, 571-590.

! Code history:
! 
! 1. 12/14/05 - created based on scalar version coare26sn.m with input
!    on vectorization from C. Moffat.  
! 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
!    values from the first pass through the iteration loop for the stable case
!    with very thin M-O length relative to zu (zetu>50) (as is done in the 
!    scalar coare26sn and COARE3 codes).
! 3. 7/26/11 - S = dt was corrected to read S = ut.
! 4. 7/28/11 - modification to roughness length parameterizations based 
!    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
! 5. 02/14/2017 matlab-->fortran senya@umd.edu mimics ncar_ocean_fluxes I/O parameters
!
!-----------------------------------------------------------------------

!	real   , intent(in)   , dimension(1) :: u, t, ts, Q, Qs, zu          ! Senya's test interface
!	logical, intent(in)   , dimension(1) :: avail                        ! Senya's test interface
!	real   , intent(inout), dimension(n) :: cd, ch, ce, ustar, bstar     ! Senya's test interface

	real   , intent(in)   , dimension(:) :: u, t, ts, Q, Qs, zu
        logical, intent(in)   , dimension(:) :: avail
        real   , intent(inout), dimension(:) :: cd, ch, ce, ustar, bstar
!ta, ts [K]
!Q, Qs [kg/kg]

!  	real    :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
!  	real    :: cd_rt                                ! full drag coefficients @ z
!  	real    :: zeta, x2, x, psi_m, psi_h            ! stability parameters
!  	real    :: u10, tv, tstar, qstar, z0, xx, stab
	real    :: Le, cpv, rhoa, visa, Al, bigc, wetc
	real    :: P, Rs, Rl, lat, zi, Rns , L, L50, L10
  	integer :: i, j, n
	
	real    :: gs, Rnl, us, rain, du, dt, dq, ta, ug, dter, ut
	real    :: Ch10, Ct10, Cd10, Cd1, Ct1, CC, Ribcu, Ribu, zetu
	real    :: usr, zo10, k50, gf, tsr, qsr, tkt, u10, zot10,charn
	real    :: umax, a1, a2, zet, zo, rr, zot, zoq, cdhf, cqhf, cthf
	real    :: tvsr, tssr, Bf, hsb, hlb, qout, dels, qcol, alq, xlamx
	real    :: dqer, usr50, tsr50, zet50, dter50, dqer50, tkt50
	real    :: qsr50, tau, hbb, hsbb
	
  	real,    parameter :: zt   = 2, zq = 2;         ! T and Q measurement level (m)
	integer, parameter :: nits = 6;
	

!***********  set constants **********************************************
	real, parameter :: Beta = 1.2;
	real, parameter :: von  = 0.4;
	real, parameter :: fdg  = 1.00; ! Turbulent Prandtl number
	real, parameter :: tdk  = 0; !273.16; !ta, ts [K]
	real, parameter :: tdc  = 273.16; !ta, ts [K]->[C]
	real, parameter :: grav = 9.8062; !grv(45deg);

!***********  air constants **********************************************
	real, parameter :: Rgas = 287.1;
	real, parameter :: cpa  = 1004.67;

!***********  cool skin constants  ***************************************
	integer, parameter :: jcool= 0; !ts is treated as the true SST
	real,    parameter :: be   = 0.026;
	real,    parameter :: cpw  = 4000;
	real,    parameter :: rhow = 1022;
	real,    parameter :: visw = 1e-6;
	real,    parameter :: tcw  = 0.6;


!************************************************************

!The following parameter determines which version of the moisture 
!roughness length.
!  0: Step this to 0 to use the form of the moisture rough derived by 
!  Fairall et al. (2003).
!
!  1: Step this to 1 to use the form of the moisture rough determined by the
!  CLIMODE, GASEX and CBLAST data.
!
!  Note that the thermal roughness length gives a Stanton number that is
!  very similar to COARE 3.0.
!
	integer, parameter :: climodeversion=0;
!******************************************************************

! set  variables that are not inputed to fortran codes to default values
!***********************************************************************
	P=1013;     ! pressure (mb)
	Rs=150;     ! incident shortwave radiation (W/m^2)
	Rl=370;     ! incident longwave radiation (W/m^2)
	lat=45;     ! latitude
	zi=600;     ! PBL height (m)
	
!***********  net radiation fluxes ***************************************
	Rns = 0.945*Rs; ! albedo correction
! IRup = eps*sigma*T^4 + (1-eps)*IR
! Rnl = IRup - IR
! Rll = eps*sigma*T^4 - eps*IR  as below


!****	do i=1,n
	do i=1,size(u(:))
	   if (avail(i)) then
	Rnl = 0.97*(5.67e-8*(ts(i)-0.3*jcool+tdk)**4-Rl); ! initial value

! IRup = Rnl + IR
!


! input variable u is assumed relative wind speed (magnitude of difference
! between wind and surface current vectors). to follow orginal Fairall code, set
! surface current speed us=0. if us data are available, construct u prior to
! using this code.
	us = 0;

! convert rh to specific humidity !! not needed Qs and Q are inputed in (kg/kg)
!Qs = qsat26sea(ts,P)./1000;    ! surface water specific humidity (kg/kg)
!Q  = qsat26air(t,P,rh)./1000;  ! specific humidity of air (kg/kg)

! set rain to zero
	rain = 0; ! rain rate (mm/hr) - keep as option

!***********  air variables **********************************************
	Le   = (2.501-.00237*(ts(i)-tdc) )*1e6;
	cpv  = cpa*(1+0.84*Q(i));
	rhoa = P*100./(Rgas*(t(i)+tdk)*(1+0.61*Q(i)));
	visa = 1.326e-5*(1+6.542e-3*(t(i)-tdc)+8.301e-6*(t(i)-tdc)**2 &
               -4.84e-9*(t(i)-tdc)**3);
	
!***********  cool skin variables  ***************************************
	Al   = 2.1e-5*(ts(i)-tdc+3.2)**0.79;
	bigc = 16*grav*cpw*(rhow*visw)**3/(tcw**2*rhoa**2);
	wetc = 0.622*Le*Qs(i)/(Rgas*(ts(i)+tdk)**2);

	Rnl = 0.97*(5.67e-8*(ts(i)-0.3*jcool+tdk)**4-Rl); ! net LW initial value

! IRup = Rnl + IR

!****************  begin bulk loop ********************************************

!***********  first guess ************************************************
	du = u(i)-us;
	dt = ts(i)-t(i)-.0098*zt;
	dq = Qs(i)-Q(i);
	ta = t(i)+tdk;
	ug=0.5;
	dter=0.3;
	ut    = sqrt(du**2+ug**2); 
	u10   = ut*log(10./1e-4)/log(zu(i)/1e-4);
	usr   = 0.035*u10;
	zo10  = 0.011*usr**2/grav + 0.11*visa/usr;
	Cd10  = (von/log(10./zo10))**2;
	Ch10  = 0.00115;
	Ct10  = Ch10/sqrt(Cd10);
	zot10 = 10/exp(von/Ct10);
	Cd1    = (von/log(zu(i)/zo10))**2;
	Ct1    = von/log(zt/zot10);
	CC    = von*Ct1/Cd1;
	Ribcu = -zu(i)/zi/0.004/Beta**3;
	Ribu  = -grav*zu(i)/ta*((dt-dter*jcool)+0.61*ta*dq)/ut**2;
	zetu = CC*Ribu*(1.+27./9.*Ribu/CC);
!!!!!!!!k50=find(zetu>50); ! stable with very thin M-O length relative to zu
	k50 = 0; if (zetu.gt.50.) k50 = 1;
!k=find(Ribu<0); zetu(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu(k)); clear k;
        if (Ribu.lt.0.) zetu=CC*Ribu/(1+Ribu/Ribcu);
	L10 = zu(i)/zetu;
	gf=ut/du; 
	usr = ut*von/(log(zu(i)/zo10)-psiu_26(zu(i)/L10));
	tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26(zt/L10));
	qsr = -(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10) &
              -psit_26(zq/L10));
	tkt = 0.001;

!**********************************************************
!  The following gives the new formulation for the Charnock parameter
!**********************************************************

	charn = 0.011;
	umax=22;
	a1=0.0016;
	a2=-0.0035;
	charn=a1*u10+a2;
!k=find(u10>umax);
!charn(k)=a1*umax+a2;
	if (u10.gt.umax) charn=a1*umax+a2;


!**************  bulk loop **************************************************

	do j=1,nits
    	zet=von*grav*zu(i)/ta*(tsr +0.61*ta*qsr)/(usr**2); !in the original coare40vn.m
!	zet=von*grav*zu(i)/ta*(tsr +0.61*ta*qsr/(1+0.61*Q(i)))/(usr**2);
    	zo=charn*usr**2/grav+0.11*visa/usr; ! surface roughness
    	rr=zo*usr/visa;
    	L=zu(i)/zet;
    	zot=min(1.0e-4/rr**0.55,2.4e-4/rr**1.2); ! temp roughness
    	if (climodeversion.eq.1) then
        	zoq=min(2.0e-5/rr**0.22,1.1e-4/rr**0.9);  ! moisture roughness
    	else
        	zoq=min(1.15e-4,5.5e-5/rr**0.60);         ! moisture roughness
    	end if

    	cdhf=von/(log(zu(i)/zo)-psiu_26(zu(i)/L));
    	cqhf=von/(log(zq/zoq)-psit_26(zq/L));
    	cthf=von/(log(zt/zot)-psit_26(zt/L));
    	usr=ut*cdhf;
    	qsr=-(dq-wetc*dter*jcool)*cqhf;
    	tsr=-(dt-dter*jcool)*cthf;
    	tvsr=tsr+0.61*ta*qsr;
    	tssr=tsr+0.51*ta*qsr;
    	Bf=-grav/ta*usr*tvsr;
    	ug=0.2;
!    k=find(Bf>0); ug(k)=max(.2,Beta*(Bf(k).*zi).^.333); clear k;
	if (Bf.gt.0.) ug=max(.2,Beta*(Bf*zi)**.333); 
    	ut=sqrt(du**2+ug**2);
    	gf=ut/du;
    	hsb=-rhoa*cpa*usr*tsr;
    	hlb=-rhoa*Le*usr*qsr;
    	qout=Rnl+hsb+hlb;
	dels=Rns*(0.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4)));
    	qcol=qout-dels;
    	alq=Al*qcol+be*hlb*cpw/Le;
    	xlamx=6.0;
    	tkt=min(0.01, xlamx*visw/(sqrt(rhoa/rhow)*usr));
!    k=find(alq>0); xlamx(k)=6./(1+(bigc(k).*alq(k)./usr(k).^4).^0.75).^0.333;
!    tkt(k)=xlamx(k).*visw./(sqrt(rhoa(k)./rhow).*usr(k)); clear k;
	if (alq.gt.0.) then
	xlamx=6./(1+(bigc*alq/usr**4)**0.75)**0.333;
	tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr);
	end if
    	dter=qcol*tkt/tcw;
    	dqer=wetc*dter;
    	Rnl=0.97*(5.67e-8*(ts(i)-dter*jcool+tdk)**4-Rl); ! update dter
!    if j==1; ! save first iteration solution for case of zetu>50;
!        usr50=usr(k50);tsr50=tsr(k50);qsr50=qsr(k50);L50=L(k50);
!        zet50=zet(k50);dter50=dter(k50);dqer50=dqer(k50);tkt50=tkt(k50);
!    end
	if (j.eq.1.and.k50.eq.1) then
	usr50=usr;tsr50=tsr;qsr50=qsr;L50=L;
        zet50=zet;dter50=dter;dqer50=dqer;tkt50=tkt;
	end if

    	u10 = ut + usr/von*(log(10./zu(i))-psiu_26(10./L)  &
             +psiu_26(zu(i)/L));
    	charn=a1*u10+a2;
    	if (u10.gt.umax) charn=a1*umax+a2;
	end do ! iteration

! insert first iteration solution for case with zetu>50
!usr(k50)=usr50;tsr(k50)=tsr50;qsr(k50)=qsr50;L(k50)=L50;
!zet(k50)=zet50;dter(k50)=dter50;dqer(k50)=dqer50;tkt(k50)=tkt50;
	if (k50.eq.1) then
	usr=usr50;tsr=tsr50;qsr=qsr50;L=L50;
	zet=zet50;dter=dter50;dqer=dqer50;tkt=tkt50;
	end if

!****************  compute fluxes  ********************************************
	tau=rhoa*usr*usr/gf;      ! wind stress
	hsb=rhoa*cpa*usr*tsr;     ! sensible heat flux
	hlb=rhoa*Le*usr*qsr;      ! latent heat flux
	hbb=rhoa*cpa*usr*tvsr;    ! buoyancy flux
	hsbb=rhoa*cpa*usr*tssr;   ! sonic heat flux

!*****  compute transfer coeffs relative to ut @ meas. ht  ********************
	Cd(i)=tau/rhoa/ut/max(.1,du);
	Ch(i)=-usr*tsr/ut/(dt-dter*jcool);
	Ce(i)=-usr*qsr/(dq-dqer*jcool)/ut;
	ustar(i)=usr;
!	bstar(i)=zet*usr**2/(von*zu(i));
!	bstar(i) = grav*(tsr/ta+0.608*qsr)/(0.608*q(i)+1);
	bstar(i) = -Bf/usr;
!	print *, 'tstar=',tsr,' qstar=',qsr
!***  compute 10-m neutral coeff relative to ut (output if needed x1000) ************
!	Cdn_10=1000*von**2/log(10./zo)**2;
!	Chn_10=1000*von**2*fdg/log(10./zo)/log(10./zot);
!	Cen_10=1000*von**2*fdg/log(10./zo)/log(10./zoq);

!***  compute 10-m neutral coeff relative to ut (output if needed) ************
!  Find the stability functions
!*********************************
!	zrf_u=10;             !User defined reference heights
!	zrf_t=10;
!	zrf_q=10;
!	psi=psiu_26(zu(i)/L);
!	psi10=psiu_26(10./L);
!	psirf=psiu_26(zrf_u/L);
!	psiT=psit_26(zt/L);
!	psi10T=psit_26(10./L);
!	psirfT=psit_26(zrf_t/L);
!	psirfQ=psit_26(zrf_q/L);
!	gf=ut/du;

!*********************************************************
!  Determine the wind speeds relative to ocean surface
!  Note that usr is the friction velocity that includes 
!  gustiness usr = sqrt(Cd) S, which is equation (18) in
!  Fairall et al. (1996)
!*********************************************************
!********************still in Matlab**********************
!using UPPER CASE may interfere with FORTRAN
!*********************************************************
!*
!*

!S = ut;
!U = du;
!S10 = S + usr./von.*(log(10./zu)-psi10+psi);
!U10 = S10./gf;
! or U10 = U + usr./von./gf.*(log(10/zu)-psi10+psi);
!Urf = U + usr./von./gf.*(log(zrf_u./zu)-psirf+psi);
!UN = U + psi.*usr/von./gf;
!U10N = U10 + psi10.*usr/von./gf;
!UrfN = Urf + psirf.*usr/von./gf;

!UN2 = usr/von./gf.*log(zu./zo);
!U10N2 = usr./von./gf.*log(10./zo);
!UrfN2  = usr./von./gf.*log(zrf_u./zo);

!lapse=grav/cpa;
!SST=ts-dter*jcool;

!T = t;
!T10 = T + tsr./von.*(log(10./zt)-psi10T+psiT) + lapse*(zt-10);
!Trf = T + tsr./von.*(log(zrf_t./zt)-psirfT+psiT) + lapse*(zt-zrf_t);
!TN = T + psiT.*tsr/von;
!T10N = T10 + psi10T.*tsr/von;
!TrfN = Trf + psirfT.*tsr/von;

!TN2 = SST + tsr/von.*log(zt./zot)-lapse*zt;
!T10N2 = SST + tsr/von.*log(10./zot)-lapse*10;
!TrfN2 = SST + tsr/von.*log(zrf_t./zot)-lapse*zrf_t;

!dqer=wetc.*dter*jcool;
!SSQ=Qs-dqer;
!SSQ=SSQ*1000;
!Q=Q*1000;
!qsr=qsr*1000;
!Q10 = Q + qsr./von.*(log(10./zq)-psi10T+psiT);
!Qrf = Q + qsr./von.*(log(zrf_q./zq)-psirfQ+psiT);
!QN = Q + psiT.*qsr/von./sqrt(gf);
!Q10N = Q10 + psi10T.*qsr/von;
!QrfN = Qrf + psirfQ.*qsr/von;

!QN2 = SSQ + qsr/von.*log(zq./zoq);
!Q10N2 = SSQ + qsr/von.*log(10./zoq);
!QrfN2 = SSQ + qsr/von.*log(zrf_q./zoq);
!RHrf=RHcalc(Trf,P,Qrf/1000);

!****************  output  ****************************************************

!A=[usr tau hsb hlb hbb hsbb tsr qsr zot zoq Cd Ch Ce  L zet dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10];
!   1   2   3   4   5   6    7   8   9  10  11 12 13 14  15  16   17   18  19  20  21  22   23  24  25  26  27  28  29     30     31    32
	end if !! avail
	end do !! U
	end subroutine coare40vn_ocean_fluxes
!
!------------------------------------------------------------------------------
	function psit_26(zet) result(psi)
	real, intent(in)  :: zet
	real            :: psi, dzet, x, psik, f, psic  
	
! computes temperature structure function
	dzet=min(50.,0.35*zet); ! stable
	psi=-((1+0.6667*zet)**1.5+0.6667*(zet-14.28)*exp(-dzet)+8.525);
	if (zet.lt.0.) then ! unstable
	x=(1-15*zet)**0.5;
	psik=2*log((1+x)/2);
	x=(1-34.15*zet)**0.3333;
	psic=1.5*log((1+x+x**2)/3)-sqrt(3.)*atan((1+2*x)/sqrt(3.))  &
            +4*atan(1.)/sqrt(3.);
	f=zet**2/(1+zet**2);
	psi=(1-f)*psik+f*psic;
	end if
	end function psit_26
!------------------------------------------------------------------------------
	function psiu_26(zet) result(psi)
	real, intent(in)  :: zet
	real            :: psi, x, psik, f, dzet, psic
	
! computes velocity structure function
	dzet=min(50.,0.35*zet); ! stable
	psi=-((1+zet)+0.6667*(zet-14.28)*exp(-dzet)+8.525);
	if (zet.lt.0.) then ! unstable
	x=(1-15*zet)**0.25;
	psik=2*log((1+x)/2)+log((1+x*x)/2)-2*atan(x)+2*atan(1.);
	x=(1-10.15*zet)**0.3333;

	psic=1.5*log((1+x+x**2)/3)-sqrt(3.)*atan((1+2*x)/sqrt(3.))  &
            +4*atan(1.)/sqrt(3.);
	f=zet**2./(1+zet**2);
	psi=(1-f)*psik+f*psic;
	end if
	end function psiu_26
!------------------------------------------------------------------------------
	function bucksat(T,P) result(exx)
	real, intent(in)  :: T,P
	real            :: exx, ex, es, em
	
! computes saturation vapor pressure [mb]
! given T [degC] and P [mb]
	exx=6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P);
	end function bucksat
!------------------------------------------------------------------------------
	function qsat26sea(T,P) result(qs)
	real, intent(in)  :: T,P
	real            :: qs, ex, es, em
	
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
	ex=bucksat(T,P);
	es=0.98*ex; ! reduction at sea surface
	qs=622*es/(P-0.378*es);
	end function qsat26sea
!------------------------------------------------------------------------------
	function qsat26air(T,P,rh) result(qs)
	real, intent(in)  :: T,P,rh
	real            :: qs, es, em
	
! computes saturation specific humidity [g/kg]
! given T [degC],rh [%], and P [mb]
	es=bucksat(T,P);
	em=0.01*rh*es;
	qs=622*em/(P-0.378*em);
	end function qsat26air
!------------------------------------------------------------------------------
!function g=grv(lat)
! computes g [m/sec^2] given lat in deg
!gamma=9.7803267715;
!c1=0.0052790414;
!c2=0.0000232718;
!c3=0.0000001262;
!c4=0.0000000007;
!phi=lat*pi/180;
!x=sin(phi);
!g=gamma*(1+c1*x.^2+c2*x.^4+c3*x.^6+c4*x.^8);
!end

!------------------------------------------------------------------------------
!function RHrf=RHcalc(T,P,Q)
! computes relative humidity given T,P, & Q

!es=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
!em=Q.*P./(0.378.*Q+0.622);
!RHrf=100*em./es;
!end


end module surface_flux_mod

