module godas_mod
! STEVE: barebones rewrite to simply (a) read in ocean_cor.res.nc computed by godas 3dvar_solo,
!        and (b) apply correction as increments at each model integration timestep.
!
!<CONTACT EMAIL="david.behringer@noaa.gov"> David Behringer
!<CONTACT EMAIL="Steve.Penny@noaa.gov"> Steve Penny
!</CONTACT>
!
!<OVERVIEW>
! This module conducts an assimilation cycle.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module applies corrections determined via an external
! data assimilation procedure.
!</DESCRIPTION>
!
!<NAMELIST NAME="godas_nml">
!
!  <DATA NAME="debug_godas" TYPE="logical">
!  For debugging the godas module.
!  </DATA>
!</NAMELIST>
!
use fms_mod,           only: open_namelist_file, check_nml_error, close_file, file_exist
use fms_mod,           only: read_data, write_data
use fms_mod,           only: FATAL, WARNING, NOTE, stdout, stdlog
use mpp_mod,           only: mpp_error, mpp_pe, mpp_npes, mpp_root_pe
use mpp_mod,           only: mpp_sync, ALL_PES
use mpp_mod,           only: mpp_broadcast, mpp_transmit
use mpp_domains_mod,   only: mpp_update_domains
use mpp_domains_mod,   only: mpp_global_sum, BITWISE_EXACT_SUM
use mpp_domains_mod,   only: mpp_get_compute_domains
use mpp_io_mod,        only: mpp_open, mpp_close
use mpp_io_mod,        only: MPP_WRONLY, MPP_RDONLY, MPP_IEEE32, MPP_DIRECT, MPP_SEQUENTIAL
use mpp_io_mod,        only: MPP_SINGLE, MPP_MULTI
use time_manager_mod,  only: time_type, set_date, get_date, set_time, get_time, print_date
use time_manager_mod,  only: increment_time, decrement_time, repeat_alarm
use time_manager_mod,  only: operator(-), operator(>), operator(<), operator(<=)
use diag_manager_mod,  only: register_diag_field, send_data !!
use constants_mod,     only: pi

use ocean_domains_mod,          only: get_local_indices, get_global_indices
use ocean_types_mod,            only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,            only: ocean_time_type
use ocean_parameters_mod,       only: missing_value
!use ocean_util_mod,             only: write_timestamp
!use godas_types_mod,            only: ocean_prog_tracer_type
use ocean_types_mod,            only: ocean_prog_tracer_type !STEVE: new
!use godas_types_mod,            only: ocean_external_mode_type
use ocean_types_mod,            only: ocean_external_mode_type
use godas_types_mod,            only: ocean_cor_tracer_type, ocean_rstr_tracer_type
use godas_types_mod,            only: ocean_obsz_type, ocean_obs0_type
use godas_data_mod,             only: id_cor, id_rstr
use godas_data_mod,             only: num_cor_tracers, num_rstr_tracers, rstr_time
use godas_data_mod,             only: sst_damp, sss_damp
use godas_data_mod,             only: kass, kass2, ksalt, nsgobs, nsgsobs, maxits, npits, jemx
use godas_data_mod,             only: no_asm_rep, asm_cnt, obs_trk_cnt
use godas_data_mod,             only: gds_step, scl_incr, tov0f, sov0f, tbv0f, sbv0f, tovrf, sovrf
use godas_data_mod,             only: tbvrf, sbvrf, hrscl, hrscl0, vcvn, vsclf, no_lat_mx, yscl, ys2
use godas_data_mod,             only: xcb, xce, xcsz, ycb, yce, ycsz
use godas_data_mod,             only: s2, s1, wgns, elipt
use godas_data_mod,             only: wcn, wea, wwe, wso, wno, wgta
use godas_data_mod,             only: wcn_s, wea_s, wwe_s, wso_s, wno_s, wgta_s
use godas_data_mod,             only: cvn, cvnsalt, vtmp, vsal, ev, wrkk, vtmp_s, vsal_s
 use godas_data_mod,             only: eta_clm, cdnz, cdnzs
 use godas_data_mod,             only: num_obsz, num_obs0, num_obsa
 use godas_data_mod,             only: temp_code, salt_code, sst_code, sss_code, altm_code, ts_code
 use godas_data_mod,             only: dtemp_max, dtemp_elm, dsalt_max, dsalt_elm
 use godas_data_mod,             only: dsst_max, dsst_elm, dsss_max, dsss_elm, daltm_max, daltm_elm
 use godas_data_mod,             only: tz_wndw_fwd, tz_wndw_bwd, rtzw
 use godas_data_mod,             only: sz_wndw_fwd, sz_wndw_bwd, rszw
 use godas_data_mod,             only: t0_wndw_fwd, t0_wndw_bwd, rt0w
 use godas_data_mod,             only: s0_wndw_fwd, s0_wndw_bwd, rs0w
 use godas_data_mod,             only: al_wndw_fwd, al_wndw_bwd, ralw, wndw_secs
 use godas_data_mod,             only: g_cg, d_cg, f_cg, e_cg, t_cg, h_cg
 use godas_data_mod,             only: g_cg_s, d_cg_s, f_cg_s, e_cg_s, t_cg_s, h_cg_s
use godas_data_mod,             only: gds_freq, alrm_dur
use godas_data_mod,             only: assrestrt, rstrestrt, restore_sfc, save_all_inv, debug_godas, ovr_alrm
use godas_data_mod,             only: single_incr, apply_incr, asm_ts_seq, asm_sfc_split, godas_at_end
!use godas_data_mod,             only: asm_Tz, asm_Sz, asm_T0, asm_S0, asm_Al
!use godas_data_mod,             only: spd

!use godas_obs_mod,              only: godas_obs_track
!use godas_rstr_mod,             only: godas_rstr_comp
!
implicit none

private

logical :: godas_module_initialized = .false.

character(len=256) :: version = '$Id: godas.F90,v 1.0 2014/04/25 01:47:00 gtn Exp $'
character(len=256) :: tagname = 'Tag $Name: gds4p1 increment only $'
character(len=48), parameter          :: mod_name = 'godas_mod'

#include <ocean_memory.h>

type(ocean_grid_type), pointer   :: Grd =>NULL()
type(ocean_domain_type), pointer :: Dom =>NULL()

integer         :: index_temp
integer         :: index_salt
integer         :: asm_code

type data_type
   character(len=3) :: gridname
   character(len=128) :: fieldname_code ! used in user's code (e.g. mdl_tvv, mdl_svv, etc.)
   character(len=128) :: fieldname_file ! fieldname used in the data file (not used)
   character(len=128) :: file_name      ! name of data file
   logical :: ongrid                    ! false, not relevant, here for compatibility
   real :: factor                       ! For unit conversion, default=1
end type data_type

integer, parameter :: max_table=10

type(data_type), dimension(max_table) :: data_table

real          :: aeval, dbsq

! for diagnostics
logical :: used

! for ascii output
! integer :: unit=6

public  godas_init
public  godas_increment
!public  godas_end

!STEVE: for debugging:
LOGICAL :: dodebug = .true.

! logical :: dbg = .true.

namelist /godas_nml/ num_cor_tracers, kass, nsgobs, nsgsobs, maxits, npits, &
                     gds_step, no_asm_rep, single_incr, tov0f, sov0f, tbv0f, &
                     sbv0f, tovrf, sovrf, tbvrf, sbvrf, hrscl, hrscl0, vcvn, vsclf, &
                     no_lat_mx, yscl, tz_wndw_fwd, tz_wndw_bwd, sz_wndw_fwd, &
                     sz_wndw_bwd, t0_wndw_fwd, t0_wndw_bwd, s0_wndw_fwd, &
                     s0_wndw_bwd, al_wndw_fwd, al_wndw_bwd, assrestrt, &
                     save_all_inv, asm_ts_seq, asm_sfc_split, godas_at_end, &
                     restore_sfc, num_rstr_tracers, rstr_time, sst_damp, sss_damp, &
                     rstrestrt, debug_godas

contains


!#######################################################################
! <FUNCTION NAME="godas_init">
!
! <DESCRIPTION>
! Initialization code for godas, returning a pointer to
! the T_cor array.
! </DESCRIPTION>
!
function godas_init (Grid, Domain, Time, T_prog, num_cor, debug) &
                    result (T_cor)

  type(ocean_grid_type), intent(in), target   :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  type(ocean_time_type), intent(in)           :: Time
  type(ocean_prog_tracer_type), intent(in)    :: T_prog(:)
  integer, intent(out)                        :: num_cor
  logical, intent(in), optional               :: debug
  ! return value
  type(ocean_cor_tracer_type), dimension(:), pointer :: T_cor

  integer :: n
  integer               :: pe, ierr, mtss, mtsd
  integer               :: year, month, day, hour, minute, second
  integer               :: i, j, k, ig, jg, kg
  integer               :: num_prog_tracers
  integer               :: ioun, io_status
  real, dimension(2)    :: range_array

  if (godas_module_initialized) then
!   call mpp_error(FATAL, trim(error_header) // ' GODAS already initialized')
    call mpp_error(FATAL, 'ERROR: GODAS already initialized')
  endif

  nullify(T_cor)

  write( stdlog(),'(/a/)') trim(version)

  num_prog_tracers = size(T_prog)
  do n=1, num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo

  pe = mpp_pe()

! set namelist defaults (see godas_data for descriptions)

  num_cor_tracers     = 2
  kass                = 30         ! standard for 40-level MOM; use 35 for deep
  nsgobs              = 1          ! standard
  nsgsobs             = 1          ! standard
  maxits              = 3
  npits               = 200
  gds_step            = 864000  ! 10 days; ! 432000, 5 days, if in seconds; ! 43200, 12 hours
  no_asm_rep          = 1
  single_incr         = .false.
  tov0f               = 1.0
  sov0f               = 0.1
  tbv0f               = 1.0
  sbv0f               = 0.1
  tovrf               = 1.0
  sovrf               = 1.0
  tbvrf               = 0.01
  sbvrf               = 0.01
  hrscl               = 3.99
  hrscl0              = 3.99
  vsclf               = 0.5
  vcvn                = 1.0
  no_lat_mx           = 63.0
  yscl                = 4.0
  tz_wndw_fwd         = 14
  tz_wndw_bwd         = 14
  sz_wndw_fwd         = 14
  sz_wndw_bwd         = 14
  t0_wndw_fwd         = 7
  t0_wndw_bwd         = 7
  s0_wndw_fwd         = 7
  s0_wndw_bwd         = 7
  al_wndw_fwd         = 14
  al_wndw_bwd         = 14
  assrestrt           = .true.
  restore_sfc         = .false.
  num_rstr_tracers    = 2
  rstr_time(:)        = -1
  rstrestrt           = .true.
  sst_damp            = 0.1
  sss_damp            = 0.1
  save_all_inv        = .false.
  asm_ts_seq          = .false.
  asm_sfc_split       = .false.
  godas_at_end        = .false.
  debug_godas         = .false.

! provide for namelist over-ride

  ioun = open_namelist_file()
  read  (ioun, godas_nml,iostat=io_status)
  ierr = check_nml_error(io_status,'godas_nml')
  call close_file (ioun)

  num_cor = num_cor_tracers

! do some namelist based adjustments
  ksalt = kass

  if (.not. restore_sfc) num_rstr_tracers = 0

  write (stdout(),'(/)')
  write (stdout(), godas_nml)
  write (stdlog(), godas_nml)

  if (PRESENT(debug) .and. .not. debug_godas) then
    debug_godas = debug
  endif

  gds_freq = set_time(gds_step, 0)
  call get_time(Time%Time_step, mtss, mtsd)
  alrm_dur = set_time(mtss, 0)

  ! allocate T_cor
  allocate( T_cor  (num_cor_tracers) )
  allocate( id_cor (num_cor_tracers) )

  id_cor(:) = -1

  do n=1,num_cor_tracers-1
    T_cor(n)%complete=.false.
  enddo
  T_cor(num_cor_tracers)%complete=.true.

  ! set local array indices
  Grd => Grid
  Dom => Domain

  call get_local_indices(Dom, isd, ied, jsd, jed, isc, iec, jsc, jec)
  call get_global_indices(Dom, isg, ieg, jsg, jeg)
  nk=Grd%nk
  nj=Grd%nj

  do j=1,nj
    if (Grd%grid_y_t(j) < no_lat_mx) jemx = j
  enddo
  ys2 = yscl * yscl

  do n=1,num_cor_tracers
#ifndef STATIC_MEMORY
!#ifndef MOM4_STATIC_ARRAYS
    allocate( T_cor(n)%fcor(isd:ied,jsd:jed,nk) )
#endif
    T_cor(n)%fcor(:,:,:)        = 0.0
  enddo

! ----------------------------------------------
! register diagnostics
! ----------------------------------------------
!
! if (dodebug) print *, "Register diagnostics..."
  do n=1,num_cor_tracers
    if (n == index_temp) then
      T_cor(n)%name='tcor'
      T_cor(n)%units='Deg_C'
      T_cor(n)%longname='potential temperature correction'
      T_cor(n)%min_range=-10.0
      T_cor(n)%max_range=100.0
      T_cor(n)%init=.false.
      T_cor(n)%file_in='INPUT/ocean_cor.res.nc'
      T_cor(n)%file_out='RESTART/ocean_cor.res.nc'
      T_cor(n)%name_in='tcor'
    else if (n == index_salt) then
      T_cor(n)%name='scor'
      T_cor(n)%units='psu'
      T_cor(n)%longname='salinity correction'
      T_cor(n)%min_range=-10.0
      T_cor(n)%max_range=100.0
      T_cor(n)%init=.false.
      T_cor(n)%file_in='INPUT/ocean_cor.res.nc'
      T_cor(n)%file_out='RESTART/ocean_cor.res.nc'
      T_cor(n)%name_in='scor'
    endif
  enddo

! register diagnostics  (only if godas_at_end=.false.)
! do n=1,num_cor_tracers
!   range_array(1) = T_cor(n)%min_range
!   range_array(2) = T_cor(n)%max_range
!   id_cor(n) = register_diag_field ('ocean_model', trim(T_cor(n)%name), &
!          Grd%tracer_axes(1:3),                                             &
!          Time%model_time, trim(T_cor(n)%longname), trim(T_cor(n)%units), &
!          missing_value=missing_value, range=range_array)
! enddo

! gds_step :: length_of_run in seconds
!STEVE: assuming mtss is time step in seconds?
  scl_incr = float(mtss) / float(gds_step)

  !STEVE: Read the existing increment
  do n=1,num_cor_tracers
    if (.not. T_cor(n)%init) then
      write (stdout(),'(/a,a)') 'Expecting to read a GODAS restart file, ', T_cor(n)%file_in
    endif

    T_cor(n)%fcor(:,:,:) = 0.0
  
    if (file_exist(T_cor(n)%file_in)) then

      ! READ IN THE CORRECTION FILE HERE:
      call read_data(T_cor(n)%file_in, T_cor(n)%name_in, T_cor(n)%fcor(:,:,:), Dom%domain2d, timelevel=1)

      ! Rescale the correction based on the model runtime:
      T_cor(n)%fcor(:,:,:) = scl_incr * T_cor(n)%fcor(:,:,:)

      write (stdout(),'(/a,1pe12.3)') 'GODAS restart increment rescaled: ',scl_incr

    else
      write (stdout(),'(/a)')'GODAS restart not found, increments set to zero.'
    endif
  enddo

  godas_module_initialized = .true.
! if (dodebug) print *, "DONE godas_init."

end function godas_init
! </FUNCTION> NAME="godas_init">


!#######################################################################
! <SUBROUTINE NAME="godas_increment">
!
! <DESCRIPTION>
! Apply corrections from analysis.  Update analysis at specified interval.
! </DESCRIPTION>
!
!subroutine godas_increment (Time, T_prog, Ext_mode, T_cor, obs_Z, obs_0, obs_A, T_rstr)
subroutine godas_increment (Time, T_prog, Ext_mode, T_cor, T_rstr) !STEVE: simplifying...

  type(ocean_time_type), intent(in)                    :: Time
  type(ocean_prog_tracer_type), intent(inout)          :: T_prog(:)
  type(ocean_external_mode_type), intent(inout)        :: Ext_mode
  type(ocean_cor_tracer_type), intent(inout)           :: T_cor(num_cor_tracers)
! type(ocean_obsz_type), intent(inout)                 :: obs_Z(:)
! type(ocean_obs0_type), intent(inout)                 :: obs_0(:)
! type(ocean_obs0_type), intent(inout)                 :: obs_A(:)
  type(ocean_rstr_tracer_type), intent(inout)          :: T_rstr(:)

  integer         :: n, taup1, pe
  integer         :: year, month, day, hour, minute, second
  integer         :: i, j, k

! if godas_at_end=.true. no action is taken here. instead a single analysis is called
!  in subroutine godas_end

! if (dodebug) print *, "STEVE: In godas_increment ======================================="

  pe = mpp_pe()

! if (dodebug) print *, "STEVE: call get_date"
  call get_date(Time%model_time, year, month, day, hour, minute, second)
! if (dodebug) print *, "STEVE: year, month, day, hour, minute, second = ", year, month, day, hour, minute, second

!! send increments/corrections to diag_manager
  do n=1,num_cor_tracers
    if (id_cor(n) > 0) used = send_data (id_cor(n), T_cor(n)%fcor(:,:,:), &
                                  Time%model_time,rmask=Grd%tmask(:,:,:), &
                                  is_in=isc, js_in=jsc, ks_in=1, ie_in=iec, je_in=jec, ke_in=nk)
  enddo

!! send restorations to diag_manager
!   do n=1,num_rstr_tracers
!     if (id_rstr(n) > 0) used = send_data (id_rstr(n), T_rstr(n)%frstr(:,:), &
!                                 Time%model_time,rmask=Grd%tmask(:,:,1), &
!                                 is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
!   enddo

!  Apply (scaled) Increments 
!  increments are applied in equal parts at each time step of the model, run between analyses. 
  taup1   = Time%taup1
  do n=1,num_cor_tracers
    where (T_prog(index_temp)%field(:,:,:,taup1) > 0.0) &
      T_prog(n)%field(:,:,:,taup1) = T_prog(n)%field(:,:,:,taup1) + T_cor(n)%fcor(:,:,:)
  enddo

  if (restore_sfc) then
    do n=1,num_rstr_tracers
      T_prog(n)%field(:,:,1,taup1) = T_prog(n)%field(:,:,1,taup1) + T_rstr(n)%frstr(:,:)
    enddo
  endif

end subroutine godas_increment
! </SUBROUTINE> NAME="godas_increment">

end module godas_mod
