! caoalpha 1/31/2006 monthly alpha & dsdt
! caosbias 11/05
! caoggbias  8/2005  adding "bias"
! caonst  1/19/05 for adding "nst" data
! caovert 1/6/05 for adding vertical corr
 
module params
 
!====================== Assimilation Parameters ========================
!
!     kmb    = number of levels analysis will be done
!     kdh    = number of levels used to calculate dynamic height
!     ksh    = number of levels of temp derived from altimetry data
! caoggbias
!     kmc    = number of levels having model bias correction
!
!     mbx    = maximum number of observations for a batch
!     maxd   = maximum number of observations for one level
!     nbatch = total number of batches for a level
!     bzone  = Batch size ti collect obs for each grid point (km)
!     maxar  = maximum possible obs. number in a single patch
!
    integer, parameter :: imt = 360, jmt = 180, km = 50  ! lgchen: km originally 40
 
!   integer, parameter :: kmb = 22, kdh = 22, ksh = 19, kmc = 22,  ! lgchen: this is the original setting
!    &                    kobs = 14, kvar = 3,
!    &                    mbx  = 100,
!    &                    nbatch = imt*jmt

    integer, parameter ::  &
        kmb  = 34, kdh = 34, ksh = 19, kmc = 34,  &
        kobs = 14, kvar = 3,                      &
        mbx  = 100,                               &
        nbatch = imt*jmt
 
    logical ::  &
        lgcl_vert, lgcl_bias,                     &
        lgcl_ctd, lgcl_sst, lgcl_nst, lgcl_altm,  &
        ! cao_altm
        lgcl_geosat, lgcl_ers1_c, lgcl_ers1_g,    &
        lgcl_ers2, lgcl_jason1, lgcl_topex

    logical :: do_deepocn_relaxation
    integer :: deepocn_relaxation_time
 
    ! integer, parameter :: maxar = 3000 
    real, parameter :: bzone = 750  !1000  !300.  
end module params
 
 
 

module constants
    real, parameter :: pi = 3.14157,  radius = 6367000. 
    real, parameter :: deg_to_radians = pi/180.
    real, parameter :: km_per_deg = 111.0
end module constants

 
 
 
module io_units
  ! stdin  = unit number for standard input.
  ! stdout = unit number for standard output.
  ! stderr = unit number for standard error.
 
    use params
 
    integer, parameter ::  stdin = 5, stdout = 6, stderr = 6
 
    integer :: io_list, io_model, io_output
    integer :: io_grid_xy, io_grid_z
 
    ! ccaoalpha
    integer :: io_vert, io_dsdt
    integer :: io_levitus_temp, io_levitus_salt 
    integer, dimension(kobs) :: io_obs_unit, io_obs_work_unit
    integer, dimension(kvar) ::              io_var_work_unit

    integer :: io_temp_bias, io_salt_bias 
    integer :: io_g5_t, io_g1_t

    ! ccaoalpha
    integer ::  &
        ioalpha,                   &
        ! caonew
        io_geo_msfh, io_geo_date,  &
        io_es1_msfh, io_es1_date,  &
        io_es2_msfh, io_es2_date,  &
        io_js1_msfh, io_js1_date,  &
        io_tpx_msfh, io_tpx_date
end module io_units
 

 
  
module model_grid
    use params
 
    integer, dimension(imt, jmt) :: kmt
    integer                      :: year, month, day
 
    ! model grid specification
    real, dimension(imt) :: xt 
    real, dimension(jmt) :: yt 
    real, dimension(km ) :: zt, dzt
 
    ! model fields for temperature and salinity
    real, dimension(imt, jmt, km) :: ts, temp_first_guess
    real, dimension(imt, jmt, km) :: ss, salt_first_guess
    real, dimension(imt, jmt, km) :: ctemp_first_guess
    real, dimension(imt, jmt, km) :: csalt_first_guess
    real, dimension(imt, jmt    ) :: sst_obs, sss_obs
    real, dimension(imt, jmt    ) :: sfh_first_guess
    
    ! caoggbias
    real, dimension(imt, jmt, kmc) :: temp_bias, salt_bias
    real                           :: rmssfh 
 
    real, dimension(imt, jmt, km)  :: levitus_temp, levitus_salt
    ! caovert
    real, dimension(kmb, kmb) :: cvrt
    ! ccaoalpha 1/2006 & ccaosbias  11/05
    ! old version where dsdt was used for bias correction and altimenty assimilation
!    real, dimension(imt, jmt, kmc) :: dsdt 
    ! new version where dsdt is used for dT -> dS updating
    real, dimension(imt, jmt, km) :: dsdt
  
    integer, dimension(imt, jmt)  :: mldk 
 
  ! data zt /  &
  !     5.   , 15.  , 25.  , 35.  , 46.  , 57.  , 70.  , 82.  , 96. , 112.,  &
  !     129. , 148. , 171. , 197. , 229. , 268. , 317. , 381. , 465.,        &
  !     579. , 729. , 918. , 1139., 1378., 1625., 1875., 2125.,              &
  !     2375., 2624., 2874., 3124., 3374., 3624., 3874., 4124.,              &
  !     4374., 4624., 4874., 5124., 5374.                    /

    ! lgchen modified on 20140905
    data zt /  &
        5.03354978561401, 15.1006498336792, 25.2193508148193,  &
        35.3584518432617, 45.5763511657715, 55.8532485961914,  &
        66.2617492675781, 76.802848815918 , 87.5769500732422,  &
        98.6232528686523, 110.096153259277, 122.106651306152,  &
        134.908554077148, 148.746643066406, 164.053756713867,  &
        181.312454223633, 201.262954711914, 224.777252197266,  &
        253.068145751953, 287.550842285156, 330.007751464844,  &
        382.365051269531, 446.726348876953, 524.982421875   ,  &
        618.703125      , 728.692077636719, 854.993530273438,  &
        996.715270996094, 1152.37609863281, 1319.99670410156,  &
        1497.56188964844, 1683.05725097656, 1874.78820800781,  &
        2071.251953125  , 2271.3232421875 , 2474.04296875   ,  &
        2678.75708007812, 2884.89794921875, 3092.1171875    ,  &
        3300.0859375    , 3508.63305664062, 3717.56713867188,  &
        3926.81274414062, 4136.25146484375, 4345.86376953125,  &
        4555.56591796875, 4765.36865234375, 4975.20947265625,  &
        5185.11083984375, 5395.02294921875                    /
end module model_grid
 

 
 
module obs_types
    use params
 
    integer                  :: obs_read, obs_index
    integer, dimension(kobs) :: obs_intt
    integer, dimension(kobs) :: num_obs_read
 
    character(len=3), dimension(kobs) :: obs_type
end module obs_types


 
 
! Variable types i.e. temperature, salinity, velocity, ice ......
module var_types
    integer, dimension(5) :: num_var_used
 
    character (len=3), dimension(5)  :: var_type
end module var_types
 
 

 
module one_d
    use params
 
    real, dimension(kmb) :: td, d, tp, sp, sd
end module one_d
 
 

 
module batches
    use params
 
    real   , dimension(nbatch, 2)   :: x_center, y_center
    real   , dimension(nbatch, 2)   :: x_west, x_east
    real   , dimension(nbatch, 2)   :: y_south, y_north
    integer, dimension(nbatch, 2)   :: i_model_west, i_model_east
    integer, dimension(nbatch, 2)   :: j_model_south, j_model_north
    integer, dimension(nbatch, 2)   :: i_model_center, j_model_center
    integer, dimension(kmb)         :: iobs, nbmax
  ! integer, dimension(nbatch, kmb) :: ibch
    integer, save, dimension(nbatch, kmb) :: ibch
    integer                         :: batch
end module batches

 
 
 
module altimetry
    use params
 
    ! ccaoalpha
    real, dimension(imt, jmt, ksh) :: alpha  ! used to get temp & salt profiles
end module altimetry
 

 
 
module loadbh
!================================================================
!
! to load in temp & salt profiles for one batch 
!
!    lv   = use (1) or not use (0) more than 1 level ob data
!    lvt  = nomber of total levels of obs data for OI on a batch
!    maxb = maximum number of observasions for one batch
!    bti  = dates of observasions
!    btp  = type of observasion data
!    bxg  = longitudal degree of observasions
!    byg  = latitudal degree of observasions
!ccaocovindp 2/01
!    nii  = nearest model i corresponding to bxg
!    njj  = nearest model j corresponding to byg

!    bxm  = longitudal distance of observasions in meter
!    bym  = latitudal distance of observasions in meter
!    btt  = temperature of observasions 
!    nlv  = depth corresponding to btt
!
!================================================================
 
    use params
 
    save
 
    ! caovert
  ! parameter(lv = 0, lvt = 1 + 2*lv)
  ! integer :: maxb
    integer :: lv, lvt, maxb
    integer, allocatable, dimension(:) :: nlv, nii, njj
    real   , allocatable, dimension(:) :: bti, bxg, byg, bxm, bym, btt
    character(len=3), allocatable, dimension(:) :: btp
  ! integer, dimension(maxar) :: nlv, nii, njj
  ! real   , dimension(maxar) :: bti, bxg, byg, bxm, bym, btt
  ! character(len=3), dimension(maxar) :: btp
 
    integer :: ithrnum
    real    :: elapsed1, elapsed, rtc
      
!$omp threadprivate(bti, btp, bxg, byg, bxm, nii, njj, bym, btt, nlv, maxb)
end module loadbh

 
 
 

module loadts
!================================================================
!
! to load in temp & salt profiles for analysis
!
!    maxd   = maximum number of observasions for one level
!    jdatg  = array of observasion dates
!    xlatg  = array of observasion latitudes
!    xlngg  = array of observasion longitudes
!    jlatg  = array of the nearest model j corresponding to xlatg
!    ilngg  = array of the nearest model i corresponding to xlngg
!
!    typeg  = array of observasion type (xbt, slt, slh, ...)
!    tempg  = array of observasion profiles   
!
!================================================================
 
    save
 
    integer :: maxd
    character(len=3) :: obtyp
    real            , parameter :: tmiss = -99.9
    integer         , allocatable, dimension(:  ) :: jdatg
    real            , allocatable, dimension(:  ) :: xlatg, xlngg
    integer         , allocatable, dimension(:  ) :: jlatg, ilngg
    character(len=3), allocatable, dimension(:  ) :: typeg
    real            , allocatable, dimension(:,:) :: tempg
end module loadts
 


 
! caoggbias  added 8/24/2005
module bias 
    use params
    use io_units   

    save
 
    integer, parameter :: np = 3, nx = 73, ny = 37
    integer :: ntu

    ! caosbias
    real, dimension(nx , ny , kmc, np) :: G5_t
    real, dimension(imt, jmt, kmc, np) :: G1_t
  ! real, dimension(nx , ny , kmc, np) :: G5_s
  ! real, dimension(imt, jmt, kmc, np) :: G1_s
end module bias
