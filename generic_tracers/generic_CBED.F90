module generic_CBED

   use g_tracer_utils, only : g_tracer_type, g_tracer_get_common, g_tracer_get_domain
   use g_tracer_utils, only : g_tracer_set_values, g_tracer_get_values
   use g_tracer_utils, only : g_tracer_get_pointer
   use g_tracer_utils, only : register_diag_field=>g_register_diag_field, g_send_data
   use cobalt_types,   only : generic_COBALT_type, phytoplankton, missing_value1, sperd, spery, epsln
   use cobalt_types,   only : SMALL, MEDIUM, LARGE, DIAZO, NUM_PHYTO
   use time_manager_mod,  only: time_type
   use field_manager_mod, only: fm_string_len, fm_path_name_len
   use mpp_domains_mod,  only : domain2D,mpp_define_io_domain
   use data_override_mod, only: data_override
   use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, read_restart, write_restart
   use fms2_io_mod, only: register_restart_field, register_axis, register_field
   use fms_mod, only: error_mesg, NOTE, WARNING, FATAL

   implicit none; private

   character(len=fm_string_len), parameter :: mod_name       = 'generic_CBED'
   character(len=fm_string_len), parameter :: package_name   = 'generic_cbed'

   public generic_CBED_sediments_update_from_source
   public generic_CBED_init, generic_CBED_end
   public generic_CBED_reg_diagnostics, generic_CBED_send_diagnostics

   integer, parameter :: nk_cbed = 20    ! Number of benthic layers

   type generic_CBED_type
      ! TODO: change read_porosity_from_file into a namelist variable
      logical :: read_porosity_from_file = .true.   ! flag to read porosity from file

      !real, dimension(:,:,:), allocatable :: f_tr1  ! tracer 1 concentration field
      real, dimension(:,:,:), allocatable :: f_o2   ! tracer o2 concentration field
      real, dimension(:,:,:), allocatable :: f_om1   ! tracer organic matter 1 (fast reacting) concentration field
      real, dimension(:,:,:), allocatable :: f_om2   ! tracer organic matter 2 (medium reacting) concentration field
      real, dimension(:,:,:), allocatable :: f_om3   ! tracer organic matter 3 (slow reacting) concentration field
      real, dimension(:,:,:), allocatable :: f_nh4   ! tracer nh4 (ammonium) concentration field
      real, dimension(:,:,:), allocatable :: f_no3   ! tracer no3 (nitrate) concentration field
      real, dimension(:,:,:), allocatable :: f_dic   ! tracer dic (dissolved inorganic carbon) concentration field
      real, dimension(:,:,:), allocatable :: f_odu   ! tracer odu (oxygen deficit unit) concentration field
      real, dimension(:,:,:), allocatable :: f_talk   ! tracer talk (total alkalinity) concentration field
      ! Diagnostics
      ! 3D diags
      real, dimension(:,:,:), allocatable :: TOC              ! total organic carbon (wt %) in sediment
      real, dimension(:,:,:), allocatable :: R_om_o2          ! OM resep via aerobic process
      real, dimension(:,:,:), allocatable :: R_om_no3         ! OM resep via denitrification
      real, dimension(:,:,:), allocatable :: R_om_anaerobic   ! OM resep via other anaerobic process
      real, dimension(:,:,:), allocatable :: R_dic            ! total remineralization

      ! 2D diags
      real, dimension(:,:), allocatable :: o2_flux !benthic o2 flux
      real, dimension(:,:), allocatable :: nh4_flux !benthic nh4 flux
      real, dimension(:,:), allocatable :: no3_flux !benthic no3 flux
      real, dimension(:,:), allocatable :: dic_flux !benthic dic flux
      real, dimension(:,:), allocatable :: burial_om !organic matter burial at the bottom of sediment column
      real, dimension(:,:), allocatable :: denit
      real, dimension(:,:), allocatable :: cbed_k1
      real, dimension(:,:), allocatable :: cbed_k2
      real, dimension(:,:), allocatable :: cbed_k3
      real, dimension(:,:), allocatable :: cbed_w
      !real, dimension(:,:), allocatable :: cbed_anammox
      !real, dimension(:,:), allocatable :: cbed_o2resp
      !real, dimension(:,:), allocatable :: cbed_no3resp
      ! 3D diags, CBED grid
      real, dimension(:,:,:), allocatable :: dz_cbed             ! cbed grid thickness
      real, dimension(:,:,:), allocatable :: z_cbed_mid          ! cbed layer mid points

      ! 3D diags (interfaces)(nk_cbed+1)
      real, dimension(:,:,:), allocatable :: cbed_Db
      real, dimension(:,:,:), allocatable :: cbed_bioirri
      real, dimension(:,:,:), allocatable :: cbed_D_o2
      real, dimension(:,:,:), allocatable :: cbed_D_dic
      real, dimension(:,:,:), allocatable :: cbed_D_nh4
      real, dimension(:,:,:), allocatable :: cbed_D_no3
      real, dimension(:,:,:), allocatable :: cbed_D_odu
      real, dimension(:,:,:), allocatable :: cbed_por
      real, dimension(:,:,:), allocatable :: cbed_svf



      !integer :: id_tr1                              ! tracer 1 diagnostics id
      integer :: id_o2                               ! tracer o2 diagnostics id
      integer :: id_om1                              ! tracer om1 diagnostics id
      integer :: id_om2                              ! tracer om2 diagnostics id
      integer :: id_om3                              ! tracer om3 diagnostics id
      integer :: id_nh4                              ! tracer nh4 diagnostics id
      integer :: id_no3                              ! tracer no3 diagnostics id
      integer :: id_dic                              ! tracer dic diagnostics id
      integer :: id_odu                              ! tracer odu diagnostics id
      integer :: id_talk                             ! tracer talk diagnostics id
      ! 3D diags
      integer :: id_TOC
      integer :: id_R_om_o2
      integer :: id_R_om_no3
      integer :: id_R_om_anaerobic
      integer :: id_R_dic
      ! 2D diags
      integer :: id_o2_flux
      integer :: id_nh4_flux
      integer :: id_no3_flux
      integer :: id_dic_flux
      integer :: id_burial_om
      integer :: id_denit
      integer :: id_cbed_k1
      integer :: id_cbed_k2
      integer :: id_cbed_k3
      integer :: id_cbed_w
      !integer :: id_cbed_anammox
      !integer :: id_cbed_o2resp
      !integer :: id_cbed_no3resp
      ! 3D diag, CBED grid
      integer :: id_dz_cbed
      integer :: id_z_cbed_mid

      ! 3D diags (interfaces)(nk_cbed+1)
      integer :: id_cbed_Db
      integer :: id_cbed_bioirri
      integer :: id_cbed_D_o2
      integer :: id_cbed_D_dic
      integer :: id_cbed_D_nh4
      integer :: id_cbed_D_no3
      integer :: id_cbed_D_odu
      integer :: id_cbed_por
      integer :: id_cbed_svf


   end type generic_CBED_type

   type(generic_CBED_type) :: cbed

   real, parameter :: pi = acos(-1.0)

! porosity. check with Niki
   !real, dimension(isc:iec,jsc:jec) :: por = 0.8  !niki
   !real :: por = 0.8

   ! grid
   ! local parameters
   real, parameter :: l_cbed = 0.20           ! length of sediment domain | sediment depth (m, 20 cm)
   real, parameter :: dz1_cbed = 0.003        ! thickness of the first layer (m). For increasing thickness
   real, parameter :: rho_s = 2.5             ! solid density (g/cm³)
   real, parameter :: Db_l = 0.08             ! bioturbation length scale (m) 8 cm.
   real, parameter :: bioirri_l = 0.018       ! bioirrigation length scale (m) 1.8 cm.

   ! sediment grid and state variables (to be allocated)
   !real, allocatable :: dz_cbed(:)                 ! sediment layer thickness (m)
   !real, allocatable :: z_cbed(:)              ! sediment depth points (m)
   real :: dz_cbed(nk_cbed)              ! thickness of each cbed layers (m)
   real :: z_cbed_int(nk_cbed+1)         ! layer interfaces (m)
   real :: z_cbed_mid(nk_cbed)           ! layer mid points (m)

   ! grid param end.

   real, dimension(:,:,:), allocatable :: por       !porosity
   real, dimension(:,:,:), allocatable :: svf       !solid volume fraction (1-porosity)

   real, dimension(:,:,:), allocatable :: w       !sedimentation rate
   real, dimension(:,:), allocatable :: Db_0    !max bioturbation rate
   real, dimension(:,:,:), allocatable :: Db    !bioturbation
   real, dimension(:,:), allocatable :: bioirri_0    !max bioirrigation rate
   real, dimension(:,:,:), allocatable :: bioirri    !bioirrigation

   real, dimension(:,:,:), allocatable :: D_o2    !diffusion coefficient for o2
   real, dimension(:,:,:), allocatable :: D_dic    !diffusion coefficient for DIC
   real, dimension(:,:,:), allocatable :: D_nh4    !diffusion coefficient for NH4
   real, dimension(:,:,:), allocatable :: D_no3    !diffusion coefficient for NO3
   real, dimension(:,:,:), allocatable :: D_odu    !diffusion coefficient for ODU (H2S)

   real, dimension(:,:), allocatable :: k1       ! k1 is the rate constant for the first order reaction of OM1 decomposition
   real, dimension(:,:), allocatable :: k2       ! k2 is the rate constant for the first order reaction of OM2 decomposition
   real, dimension(:,:), allocatable :: k3       ! k3 is the rate constant for the first order reaction of OM3 decomposition




!     call grid_cbed(nk_cbed, dz_cbed, z_cbed_mid, z_cbed_int)

!    ! define uniform sediment grid
!    dz_cbed = l_cbed / real(nk_cbed)
!    z_cbed_int(1) = 0.0   !this is likely the interface. dimention of z_cbed is nk_cbed+1. z_int_cbed. might need z_mid_cbed
!    do k = 1, nk_cbed
!        z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
!    end do
!
!    z_cbed_mid(1) = dz_cbed(1)/2   ! first layer mid point
!    do k = 1, nk_cbed-1
!        z_cbed_mid(k+1) = z_cbed_mid(k) + dz_cbed(k)
!    end do


contains

!! To make increasing thickness CBED grid
! Function to find the common ratio r using bisection method
   function find_r(nk_cbed, dz_first, total_height) result(r)
      integer, intent(in) :: nk_cbed
      real, intent(in) :: dz_first, total_height
      real :: r
      real :: r_low, r_high, r_mid
      real :: power, sum_geom
      integer :: i, j

      r_low = 1.0001
      r_high = 10.0
      do i = 1, 100
         r_mid = (r_low + r_high) / 2.0
         ! Compute r_mid**nk_cbed using iterative multiplication to avoid overflow
         power = 1.0
         do j = 1, nk_cbed
            power = power * r_mid
         end do
         sum_geom = dz_first * (power - 1.0) / (r_mid - 1.0)
         if (sum_geom < total_height) then
            r_low = r_mid
         else
            r_high = r_mid
         end if
      end do
      r = r_mid
   end function find_r

   subroutine generic_CBED_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk)
      integer,     intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk
      !Locals
      type(domain2D), pointer :: domain
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file

      integer :: i,j,k !for grid.
      real    :: r ! for grid

      !real,dimension(isc:iec,jsc:jec,nk_cbed)    :: cbed_tmask
      ! Make a cbed mask. Note: it seems grid_tmask(:,:,k) does not depend on k
      ! Note that grid_tmask is already on isc:iec, jsc:jec
      !do j = jsc, jec; do i = isc, iec; do k=1,nk_cbed ;
      !   cbed_tmask(i,j,k) = grid_tmask(i,j,nk) ; enddo; enddo; enddo

      !Allocate and initialize CBED arrays for tracer concentrations and other workarrays
      !allocate(cbed%f_tr1(isd:ied,jsd:jed,nk_cbed));cbed%f_tr1=0.0
      allocate(cbed%f_o2(isd:ied,jsd:jed,nk_cbed));cbed%f_o2=0.0
      allocate(cbed%f_om1(isd:ied,jsd:jed,nk_cbed));cbed%f_om1=0.0
      allocate(cbed%f_om2(isd:ied,jsd:jed,nk_cbed));cbed%f_om2=0.0
      allocate(cbed%f_om3(isd:ied,jsd:jed,nk_cbed));cbed%f_om3=0.0
      allocate(cbed%f_nh4(isd:ied,jsd:jed,nk_cbed));cbed%f_nh4=0.0
      allocate(cbed%f_no3(isd:ied,jsd:jed,nk_cbed));cbed%f_no3=0.0
      allocate(cbed%f_dic(isd:ied,jsd:jed,nk_cbed));cbed%f_dic=0.0
      allocate(cbed%f_odu(isd:ied,jsd:jed,nk_cbed));cbed%f_odu=0.0
      allocate(cbed%f_talk(isd:ied,jsd:jed,nk_cbed));cbed%f_talk=0.0
      !Diagnostics
      ! 3D diags
      allocate(cbed%TOC(isd:ied,jsd:jed,nk_cbed));cbed%TOC=0.0
      allocate(cbed%R_om_o2(isd:ied,jsd:jed,nk_cbed));cbed%R_om_o2=0.0
      allocate(cbed%R_om_no3(isd:ied,jsd:jed,nk_cbed));cbed%R_om_no3=0.0
      allocate(cbed%R_om_anaerobic(isd:ied,jsd:jed,nk_cbed));cbed%R_om_anaerobic=0.0
      allocate(cbed%R_dic(isd:ied,jsd:jed,nk_cbed));cbed%R_dic=0.0
      ! 2D diags
      allocate(cbed%o2_flux(isd:ied,jsd:jed)); cbed%o2_flux=0.0
      allocate(cbed%nh4_flux(isd:ied,jsd:jed)); cbed%nh4_flux=0.0
      allocate(cbed%no3_flux(isd:ied,jsd:jed)); cbed%no3_flux=0.0
      allocate(cbed%dic_flux(isd:ied,jsd:jed)); cbed%dic_flux=0.0
      allocate(cbed%burial_om(isd:ied,jsd:jed));cbed%burial_om=0.0
      allocate(cbed%denit(isd:ied,jsd:jed));cbed%denit=0.0
      allocate(cbed%cbed_k1(isd:ied,jsd:jed));cbed%cbed_k1=0.0
      allocate(cbed%cbed_k2(isd:ied,jsd:jed));cbed%cbed_k2=0.0
      allocate(cbed%cbed_k3(isd:ied,jsd:jed));cbed%cbed_k3=0.0
      allocate(cbed%cbed_w(isd:ied,jsd:jed));cbed%cbed_w=0.0
      !allocate(cbed%cbed_anammox(isd:ied,jsd:jed));cbed%cbed_anammox=0.0
      !allocate(cbed%cbed_o2resp(isd:ied,jsd:jed));cbed%cbed_o2resp=0.0
      !allocate(cbed%cbed_no3resp(isd:ied,jsd:jed));cbed%cbed_no3resp=0.0
      ! 3D diag, CBED grid
      allocate(cbed%dz_cbed(isd:ied,jsd:jed,nk_cbed));cbed%dz_cbed=0.0
      allocate(cbed%z_cbed_mid(isd:ied,jsd:jed,nk_cbed));cbed%z_cbed_mid=0.0

      ! 3D diags (interfaces)(nk_cbed+1)
      allocate(cbed%cbed_Db(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_Db=0.0
      allocate(cbed%cbed_bioirri(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_bioirri=0.0
      allocate(cbed%cbed_D_o2(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_o2=0.0
      allocate(cbed%cbed_D_dic(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_dic =0.0
      allocate(cbed%cbed_D_nh4(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_nh4=0.0
      allocate(cbed%cbed_D_no3(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_no3=0.0
      allocate(cbed%cbed_D_odu(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_D_odu=0.0
      allocate(cbed%cbed_por(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_por=0.0
      allocate(cbed%cbed_svf(isd:ied,jsd:jed,nk_cbed+1));cbed%cbed_svf=0.0


      allocate(por(isd:ied,jsd:jed,nk_cbed+1));        por=0.0  !initialize to zero
      allocate(svf(isd:ied,jsd:jed,nk_cbed+1));        svf=0.0  !solid volume fraction

      allocate(w(isd:ied,jsd:jed,nk_cbed+1));        w=0.0      !adding sedimentation rate initalize
      allocate(Db_0(isd:ied,jsd:jed));               Db_0=0.0   !bioturbation_0 init.
      allocate(Db(isd:ied,jsd:jed,nk_cbed+1));       Db=0.0     !bioturbation init.
      allocate(bioirri_0(isd:ied,jsd:jed));          bioirri_0=0.0   !bioirrigation_0 init.
      allocate(bioirri(isd:ied,jsd:jed,nk_cbed));    bioirri=0.0     !bioturbation init.

      allocate(D_o2(isd:ied,jsd:jed,nk_cbed+1)); D_o2=0.0     ! D_o2 init.
      allocate(D_dic(isd:ied,jsd:jed,nk_cbed+1)); D_dic=0.0
      allocate(D_nh4(isd:ied,jsd:jed,nk_cbed+1)); D_nh4=0.0
      allocate(D_no3(isd:ied,jsd:jed,nk_cbed+1)); D_no3=0.0
      allocate(D_odu(isd:ied,jsd:jed,nk_cbed+1)); D_odu=0.0

      allocate(k1(isd:ied,jsd:jed)); k1=0.0
      allocate(k2(isd:ied,jsd:jed)); k2=0.0
      allocate(k3(isd:ied,jsd:jed)); k3=0.0

      !if (cbed%read_porosity_from_file) then
      !   call data_override('OCN', 'por', por(isc:iec,jsc:jec,nk_cbed+1), model_time, override=.true.)
      !   svf(isc:iec,jsc:jec,nk_cbed+1) = 1.0 - por(isc:iec,jsc:jec,nk_cbed+1)
      !else
      !   por(isc:iec,jsc:jec,nk_cbed+1) = 0.8
      !   svf(isc:iec,jsc:jec,nk_cbed+1) = 0.2
      !endif


      ! Grid does not change with time, so can be define only once.

      !!!! define uniform sediment grid
      !dz_cbed = l_cbed / real(nk_cbed)
      !z_cbed_int(1) = 0.0   !this is likely the interface. dimention of z_cbed is nk_cbed+1. z_int_cbed. might need z_mid_cbed
      !do k = 1, nk_cbed
      !   z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
      !end do
      !
      !z_cbed_mid(1) = dz_cbed(1)/2   ! first layer mid point
      !do k = 1, nk_cbed-1
      !   z_cbed_mid(k+1) = z_cbed_mid(k) + dz_cbed(k)
      !end do

      !!!! grid with increasing thickness. Values are copied from R model (20 layers)
      !dz_cbed = (/ 0.001999965, 0.002295888, 0.002635596, 0.003025570, 0.003473245, 0.003987160, &
      !   0.004577116, 0.005254364, 0.006031821, 0.006924313, 0.007948862, 0.009125008, &
      !   0.010475180, 0.012025129, 0.013804415, 0.015846971, 0.018191752, 0.020883476, &
      !   0.023973478, 0.027520689 /)
      !
      !z_cbed_int = (/ 0.0, 0.001999965, 0.004295853, 0.006931449, 0.009957019, 0.013430264,&
      !   0.017417424, 0.021994540, 0.027248904, 0.033280725, 0.040205039, 0.048153901, &
      !   0.057278909, 0.067754089, 0.079779218, 0.093583633, 0.109430604, 0.127622356, &
      !   0.148505832, 0.172479311, 0.2 /)
      !
      !z_cbed_mid = (/ 0.0009999825, 0.0031479088, 0.0056136508, 0.0084442338, 0.0116936411, &
      !   0.0154238436, 0.0197059818, 0.0246217221, 0.0302648149, 0.0367428821, &
      !   0.0441794700, 0.0527164049, 0.0625164986, 0.0737666532, 0.0866814254, &
      !   0.1015071186, 0.1185264802, 0.1380640944, 0.1604925715, 0.1862396553 /)


      !!!! grid with increasing thickness. Values are copied from R model (20 layers) (works)
      !dz_cbed = (/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 /)

      !z_cbed_int = (/ 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20 /)

      !z_cbed_mid = (/ 0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.125, 0.135, 0.145, 0.155, 0.165, 0.175, 0.185, 0.195 /)


      ! generate grid automatically according to supplied values (length of sediment column, number of layers and thicness of the first layer)
      ! Compute the common ratio r
      r = find_r(nk_cbed, dz1_cbed, l_cbed)

      ! Generate dz_cbed as geometric progression
      dz_cbed(1) = dz1_cbed
      do k = 2, nk_cbed
         dz_cbed(k) = dz_cbed(k-1) * r
      end do

      ! Optional: Scale to ensure exact total sum (due to floating-point precision)
      ! real :: actual_sum
      ! actual_sum = sum(dz_cbed)
      ! if (abs(actual_sum - l_cbed) > 1e-10) then
      !   dz_cbed = dz_cbed * l_cbed / actual_sum
      ! end if

      ! Compute layer interfaces
      z_cbed_int(1) = 0.0
      do k = 1, nk_cbed
         z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
      end do

      ! Compute layer midpoints
      do k = 1, nk_cbed
         z_cbed_mid(k) = (z_cbed_int(k) + z_cbed_int(k+1)) / 2.0
      end do


   end subroutine generic_CBED_init

   subroutine generic_CBED_reg_diagnostics(axes,init_time)
      USE diag_manager_mod, ONLY: register_diag_field, diag_axis_init
      integer,         intent(in) :: axes(3)
      type(time_type), intent(in) :: init_time
      !Locals
      integer :: k, id_layer, id_layer_i
      real :: cbed_layers(1:nk_cbed)
      real :: cbed_layers_i(1:nk_cbed+1)

      !!BEGIN read_restart code block
      !Niki: This code block does not seem to belong here and should be in the _init routine instead.
      !      But the problem with that is due to MOM6 code flow, when generic_CBED_init is called
      !      the MOM "domain" is not yet created/updated and we cannot access it, which is needed for reading the restart here.
      type(domain2D), pointer :: domain
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file
      !Resgister restarts
      restart_file = 'INPUT/generic_CBED.res.nc'
      call g_tracer_get_domain(domain)
      file_open_success=open_file(fileobj, trim(restart_file),"read", domain, is_restart=.true.)
      if (file_open_success) then
         call register_axis(fileobj,'x','x')
         call register_axis(fileobj,'y','y')
         !!< Register the domain decomposed dimensions as variables so that the combiner can work correctly
         !call register_field(fileobj, "x", "double", (/"x"/))
         !call register_field(fileobj, "y", "double", (/"y"/))
         call register_axis(fileobj,'lev',nk_cbed)
         ! register the restart variables
         !call register_restart_field(fileobj, "cbed_tr1", cbed%f_tr1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_o2", cbed%f_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om1", cbed%f_om1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om2", cbed%f_om2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om3", cbed%f_om3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_nh4", cbed%f_nh4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_no3", cbed%f_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_dic", cbed%f_dic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_odu", cbed%f_odu, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_talk", cbed%f_talk, (/"x","y","lev"/))

         call read_restart(fileobj)
      endif
      !!END read_restart code block

      !!Register diagnostics
      !Niki: I am unsure of the diag axis thingy, ask Yi-Cheng
      !Define cbed layer axis, the x,y axes are the same as MOM6 since the horizontal grids are the same
      do k=1,nk_cbed; cbed_layers(k) = k; enddo
      id_layer = diag_axis_init('cbedlayer', cbed_layers, 'None', 'z', long_name='Benthos Layer', direction=-1)
      do k=1,nk_cbed+1; cbed_layers_i(k) = k; enddo
      id_layer_i = diag_axis_init('cbedlayer_i', cbed_layers_i, 'None', 'z', long_name='Benthos Layer Interface', direction=-1)

      !cbed%id_tr1 = register_diag_field(package_name, 'cbed_tr1_conc', (/axes(1),axes(2),id_layer/), init_time,&
      !   'cbed tracer1 concentration', 'unknown units', missing_value = missing_value1)
      cbed%id_o2 = register_diag_field(package_name, 'cbed_o2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed oxygen concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_om1 = register_diag_field(package_name, 'cbed_om1_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM1 concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_om2 = register_diag_field(package_name, 'cbed_om2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM2 concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_om3 = register_diag_field(package_name, 'cbed_om3_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM3 concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_nh4 = register_diag_field(package_name, 'cbed_nh4_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed ammonium concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_no3 = register_diag_field(package_name, 'cbed_no3_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed nitrate concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_dic = register_diag_field(package_name, 'cbed_dic_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed DIC concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_odu = register_diag_field(package_name, 'cbed_odu_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed ODU concentration', 'mol m-3', missing_value = missing_value1)
      cbed%id_talk = register_diag_field(package_name, 'cbed_talk_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed total alkalinity concentration', 'mol m-3', missing_value = missing_value1)
      ! diags
      ! 3D diags
      cbed%id_TOC = register_diag_field(package_name, 'cbed_TOC', (/axes(1),axes(2),id_layer/), init_time,&
         'Total Organic Carbon in sediment', 'wt %', missing_value = missing_value1)
      cbed%id_R_om_o2 = register_diag_field(package_name, 'cbed_R_om_o2', (/axes(1),axes(2),id_layer/), init_time,&
         'aerobic respiration in sediment 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_om_no3 = register_diag_field(package_name, 'cbed_R_om_no3', (/axes(1),axes(2),id_layer/), init_time,&
         'OM respiration via denitrification in sediment 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_om_anaerobic = register_diag_field(package_name, 'cbed_R_om_anaerobic', (/axes(1),axes(2),id_layer/), init_time,&
         'OM respiration via other anaerobic processes in sediment 3D field', 'mol C m-3 s-1', missing_value = missing_value1)
      cbed%id_R_dic = register_diag_field(package_name, 'cbed_R_dic', (/axes(1),axes(2),id_layer/), init_time,&
         'DIC produced in sediment via OM remineralization 3D field', 'mol C m-3 s-1', missing_value = missing_value1)

      ! 2D diags
      cbed%id_o2_flux = register_diag_field(package_name, 'cbed_o2_flux', (/axes(1),axes(2)/), init_time,&
         'benthic O2 flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_nh4_flux = register_diag_field(package_name, 'cbed_nh4_flux', (/axes(1),axes(2)/), init_time,&
         'benthic nh4 flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_no3_flux = register_diag_field(package_name, 'cbed_no3_flux', (/axes(1),axes(2)/), init_time,&
         'benthic no3 flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_dic_flux = register_diag_field(package_name, 'cbed_dic_flux', (/axes(1),axes(2)/), init_time,&
         'benthic dic flux', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_burial_om = register_diag_field(package_name, 'cbed_burial_om', (/axes(1),axes(2)/), init_time,&
         'cbed organic carbon burial', 'mol m-2 s-1', missing_value = missing_value1)
      cbed%id_denit = register_diag_field(package_name, 'cbed_denit', (/axes(1),axes(2)/), init_time,&
         'cbed total denitrification (denit + anammox)', 'mol N m-2 s-1', missing_value = missing_value1)
      cbed%id_cbed_k1 = register_diag_field(package_name, 'cbed_k1', (/axes(1),axes(2)/), init_time,&
         'OM1 decay rate constant', 's-1', missing_value = missing_value1)
      cbed%id_cbed_k2 = register_diag_field(package_name, 'cbed_k2', (/axes(1),axes(2)/), init_time,&
         'OM2 decay rate constant', 's-1', missing_value = missing_value1)
      cbed%id_cbed_k3 = register_diag_field(package_name, 'cbed_k3', (/axes(1),axes(2)/), init_time,&
         'OM3 decay rate constant', 's-1', missing_value = missing_value1)
      cbed%id_cbed_w = register_diag_field(package_name, 'cbed_w', (/axes(1),axes(2)/), init_time,&
         'sedimentation rate', 'm/s', missing_value = missing_value1)

      !cbed%id_cbed_anammox = register_diag_field(package_name, 'cbed_anammox', (/axes(1),axes(2)/), init_time,&
      !   'cbed anammox', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_o2resp = register_diag_field(package_name, 'cbed_o2resp', (/axes(1),axes(2)/), init_time,&
      !   'cbed OM respiration by O2', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_no3resp = register_diag_field(package_name, 'cbed_no3resp', (/axes(1),axes(2)/), init_time,&
      !   'cbed OM respiration by NO3', 'mol/m2/s', missing_value = missing_value1)
      ! 3D diags, CBED grid
      cbed%id_dz_cbed = register_diag_field(package_name, 'cbed_dz_cbed', (/axes(1),axes(2),id_layer/), init_time,&
         'CBED grid layer thickess', 'm', missing_value = missing_value1)
      cbed%id_z_cbed_mid = register_diag_field(package_name, 'cbed_z_cbed_mid', (/axes(1),axes(2),id_layer/), init_time,&
         'CBED grid layer midpoints', 'm', missing_value = missing_value1)

      ! 3D diags, CBED grid interfaces
      cbed%id_cbed_Db = register_diag_field(package_name, 'cbed_Db', (/axes(1),axes(2),id_layer_i/), init_time,&
         'bioturbation coefficient', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_bioirri = register_diag_field(package_name, 'cbed_bioirri', (/axes(1),axes(2),id_layer_i/), init_time,&
         'bioirrigation coefficient', 's-1', missing_value = missing_value1)
      cbed%id_cbed_D_o2 = register_diag_field(package_name, 'cbed_D_o2', (/axes(1),axes(2),id_layer_i/), init_time,&
         'O2 molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_dic = register_diag_field(package_name, 'cbed_D_dic', (/axes(1),axes(2),id_layer_i/), init_time,&
         'DIC molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_nh4 = register_diag_field(package_name, 'cbed_D_nh4', (/axes(1),axes(2),id_layer_i/), init_time,&
         'NH4 molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_no3 = register_diag_field(package_name, 'cbed_D_no3', (/axes(1),axes(2),id_layer_i/), init_time,&
         'NO3 molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_D_odu = register_diag_field(package_name, 'cbed_D_odu', (/axes(1),axes(2),id_layer_i/), init_time,&
         'ODU molecular diffusion coefficient in sediment', 'm2/s', missing_value = missing_value1)
      cbed%id_cbed_por = register_diag_field(package_name, 'cbed_por', (/axes(1),axes(2),id_layer_i/), init_time,&
         'sediment porosity', 'dimensionless', missing_value = missing_value1)
      cbed%id_cbed_svf = register_diag_field(package_name, 'cbed_svf', (/axes(1),axes(2),id_layer_i/), init_time,&
         'sediment solid volume fraction', 'dimensionless', missing_value = missing_value1)

   end subroutine generic_CBED_reg_diagnostics

   subroutine generic_CBED_send_diagnostics(model_time, isc,iec,jsc,jec, isd,ied,jsd,jed,nk, grid_tmask)
      USE diag_manager_mod, ONLY: send_data
      type(time_type),          intent(in) :: model_time
      real, dimension(:,:,:),    pointer   :: grid_tmask
      integer,                  intent(in) :: isc,iec,jsc,jec, isd,ied,jsd,jed,nk
      ! local
      logical :: used
      integer :: i,j,k
      real,dimension(isd:ied,jsd:jed,nk_cbed)    :: cbed_tmask
      real,dimension(isd:ied,jsd:jed,nk_cbed+1)    :: cbed_tmask_i
      ! Make a cbed mask. Note: it seems grid_tmask(:,:,k) does not depend on k
      !do k=1,nk_cbed ; cbed_tmask(:,:,k) = grid_tmask(:,:,nk) ; enddo
      do j = jsc, jec; do i = isc, iec; do k=1,nk_cbed ;
               cbed_tmask(i,j,k) = grid_tmask(i,j,nk) ; enddo; enddo; enddo
      do j = jsc, jec; do i = isc, iec; do k=1,nk_cbed+1 ;
               cbed_tmask_i(i,j,k) = grid_tmask(i,j,nk) ; enddo; enddo; enddo


      !used = send_data(cbed%id_tr1, cbed%f_tr1, model_time, rmask = cbed_tmask,&
      !   is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_o2, cbed%f_o2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_om1, cbed%f_om1, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_om2, cbed%f_om2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_om3, cbed%f_om3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_nh4, cbed%f_nh4, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_no3, cbed%f_no3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_dic, cbed%f_dic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_odu, cbed%f_odu, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_talk, cbed%f_talk, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)

      ! diags
      ! 3D diags
      used = send_data(cbed%id_TOC, cbed%TOC, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_o2, cbed%R_om_o2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_no3, cbed%R_om_no3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_anaerobic, cbed%R_om_anaerobic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_dic, cbed%R_dic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      ! 2D diags
      used = send_data(cbed%id_o2_flux, cbed%o2_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_nh4_flux, cbed%nh4_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_no3_flux, cbed%no3_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_dic_flux, cbed%dic_flux, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_burial_om, cbed%burial_om, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_denit, cbed%denit, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_k1, cbed%cbed_k1, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_k2, cbed%cbed_k2, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_k3, cbed%cbed_k3, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_cbed_w, cbed%cbed_w, model_time, rmask = cbed_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      ! 1D diags
      used = send_data(cbed%id_dz_cbed, cbed%dz_cbed, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_z_cbed_mid, cbed%z_cbed_mid, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)

      ! 3D diags, CBED grid interfaces
      used = send_data(cbed%id_cbed_Db, cbed%cbed_Db, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_bioirri, cbed%cbed_bioirri, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_o2, cbed%cbed_D_o2, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in= jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_dic, cbed%cbed_D_dic, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_nh4, cbed%cbed_D_nh4, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_no3, cbed%cbed_D_no3, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_D_odu, cbed%cbed_D_odu, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_por, cbed%cbed_por, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)
      used = send_data(cbed%id_cbed_svf, cbed%cbed_svf, model_time, rmask = cbed_tmask_i,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed+1)



   end subroutine generic_CBED_send_diagnostics

   subroutine generic_CBED_end()
      !Locals
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file
      type(domain2D),pointer :: domain

      !Resgister restarts
      call g_tracer_get_domain(domain)
      restart_file = 'RESTART/generic_CBED.res.nc'
      file_open_success=open_file(fileobj, trim(restart_file),"overwrite", domain, is_restart=.true.)
      if (file_open_success) then
         call register_axis(fileobj,'x','x')
         call register_axis(fileobj,'y','y')
         !< Register the domain decomposed dimensions as variables so that the combiner can work correctly
         call register_field(fileobj, "x", "double", (/"x"/))
         call register_field(fileobj, "y", "double", (/"y"/))
         call register_axis(fileobj,'lev',nk_cbed)
         ! register the restart variables
         !call register_restart_field(fileobj, "cbed_tr1", cbed%f_tr1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_o2", cbed%f_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om1", cbed%f_om1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om2", cbed%f_om2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om3", cbed%f_om3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_nh4", cbed%f_nh4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_no3", cbed%f_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_dic", cbed%f_dic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_odu", cbed%f_odu, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_talk", cbed%f_talk, (/"x","y","lev"/))

         call write_restart(fileobj)
         call close_file(fileobj)
      else
         call error_mesg( 'generic_CBED_end', 'Cannot open restarts for write.', FATAL )
      endif

      !Deallocate arrays
      !deallocate(cbed%f_tr1)
      deallocate(cbed%f_o2)
      deallocate(cbed%f_om1)
      deallocate(cbed%f_om2)
      deallocate(cbed%f_om3)
      deallocate(cbed%f_nh4)
      deallocate(cbed%f_no3)
      deallocate(cbed%f_dic)
      deallocate(cbed%f_odu)
      deallocate(cbed%f_talk)

      ! diags
      ! 3D diags
      deallocate(cbed%TOC)
      deallocate(cbed%R_om_o2)
      deallocate(cbed%R_om_no3)
      deallocate(cbed%R_om_anaerobic)
      deallocate(cbed%R_dic)
      !2D diags
      deallocate(cbed%o2_flux)
      deallocate(cbed%nh4_flux)
      deallocate(cbed%no3_flux)
      deallocate(cbed%dic_flux)
      deallocate(cbed%burial_om)
      deallocate(cbed%denit)
      deallocate(cbed%cbed_k1)
      deallocate(cbed%cbed_k2)
      deallocate(cbed%cbed_k3)
      deallocate(cbed%cbed_w)
      ! 3D diags, CBED grid
      deallocate(cbed%dz_cbed)
      deallocate(cbed%z_cbed_mid)

      deallocate(por)
      deallocate(svf)

      deallocate(w)
      deallocate(Db_0)
      deallocate(Db)
      deallocate(bioirri_0)
      deallocate(bioirri)

      deallocate(D_o2)
      deallocate(D_dic)
      deallocate(D_nh4)
      deallocate(D_no3)
      deallocate(D_odu)

      deallocate(k1)
      deallocate(k2)
      deallocate(k3)

      ! 3D diags, CBED grid interfaces
      deallocate(cbed%cbed_Db)
      deallocate(cbed%cbed_bioirri)
      deallocate(cbed%cbed_D_o2)
      deallocate(cbed%cbed_D_dic)
      deallocate(cbed%cbed_D_nh4)
      deallocate(cbed%cbed_D_no3)
      deallocate(cbed%cbed_D_odu)
      deallocate(cbed%cbed_por)
      deallocate(cbed%cbed_svf)


   end subroutine generic_CBED_end



   ! subroutine vertdiff_CBED(cobalt_tracer_list,cobalt, cbed_field, field_name, D, w, VF, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
   !    type(g_tracer_type),          pointer       :: cobalt_tracer_list
   !    type(generic_COBALT_type),    intent(inout) :: cobalt
   !    real, dimension(:,:,:),       intent(inout) :: cbed_field  ! cbed tracer concentration field
   !    character(len=*),             intent(in)    :: field_name   !Name of the cbed field (e.g., "f_o2" or "f_nh4")
   !    real, dimension(:,:,:),       intent(in)    :: D   ! diffustion
   !    real, dimension(:,:,:),       intent(in)    :: w   !sinking velocity or sedimentation rate
   !    real, dimension(:,:,:),       intent(in)    :: VF    ! volumn fraction
   !    integer, dimension(:,:),      intent(in)    :: grid_kmt
   !    real,                         intent(in)    :: dt
   !    integer,                      intent(in)    :: tau
   !    integer,                      intent(in)    :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed

   !    !Locals
   !    integer :: i, j, k

   !    real, dimension(nk_cbed) :: a,b,c,f_old,h_old
   !    real, dimension(isc:iec,jsc:jec,nk_cbed) :: ea,eb
   !    real, dimension(isc:iec,jsc:jec,nk_cbed+1) :: sink

   !    real :: sfc_src, btm_src

   !    ! local parameters for bgc reactions
   !    real, parameter :: frac_OM1 = 0.70
   !    real, parameter :: frac_OM2 = 0.20
   !    real, parameter :: frac_OM3 = 0.10




   !    ! h_old
   !    do k = 1, nk_cbed
   !       h_old(k) = dz_cbed(k)
   !    enddo

   !    ! ea , eb
   !    do j = jsc, jec; do i = isc, iec
   !          if (grid_kmt(i,j) .gt. 0) then
   !             do k=1,nk_cbed
   !                ea(i,j,k) = VF(i,j,k)*D(i,j,k)*dt/h_old(k)
   !                eb(i,j,k) = VF(i,j,k)*D(i,j,k+1)*dt/h_old(k)
   !                !print *, field_name, " ea(", i, ",", j, ",", k, ") = ", ea(i,j,k)
   !                !print *, field_name, " eb(", i, ",", j, ",", k, ") = ", eb(i,j,k)
   !             enddo
   !          endif
   !       enddo; enddo
   !    ! print ea eb for o2 for debugging
   !    if (trim(field_name) == "f_o2") then
   !       do j = jsc, jec; do i = isc, iec
   !             if (grid_kmt(i,j) .gt. 0) then
   !                if (j == 375 .and. i == 475) then
   !                   do k=1,nk_cbed
   !                      print *, field_name, " ea(", i, ",", j, ",", k, ") = ", ea(i,j,k)
   !                      print *, field_name, " eb(", i, ",", j, ",", k, ") = ", eb(i,j,k)
   !                   enddo
   !                endif
   !             endif
   !          enddo; enddo
   !    endif


   !    ! sink(k+1)
   !    do j = jsc, jec; do i = isc, iec
   !          if (grid_kmt(i,j) .gt. 0) then
   !             do k=1,nk_cbed+1
   !                sink(i,j,k) = max(0.0, VF(i,j,k)*w(i,j,k)*dt )
   !             enddo
   !          endif
   !       enddo; enddo

   !    ! a, b, c, f_old
   !    do j = jsc, jec; do i = isc, iec
   !          if (grid_kmt(i,j) .gt. 0) then

   !             !! NEED TO TAKE CARE OF THE TOP AND BOTTOM FLUXES
   !             ! sfc_src = 0.0 ; btm_src = 0.0
   !             ! if (_ALLOCATED(g_tracer%stf)) sfc_src = (g_tracer%stf(i,j)*dt)*kg_m2_to_H
   !             ! if (_ALLOCATED(g_tracer%btf)) btm_src = (-g_tracer%btf(i,j)*dt)*kg_m2_to_H
   !             ! g_tracer%field(i,j,1,tau)  = g_tracer%field(i,j,1,tau)  + sfc_src/h_old(i,j,1)
   !             ! g_tracer%field(i,j,nz,tau) = g_tracer%field(i,j,nz,tau) + btm_src/h_old(i,j,nz)

   !             sfc_src = 0.0

   !             if ((trim(field_name) == "f_o2")) then

   !                if (j == 100 .and. i == 42) then
   !                   print *, "before update in vertdiff top layer o2 cbed%f_o2(", i, ",", j, ",1) = ", cbed_field(i,j,1)
   !                endif

   !                ! if (j == 5 .and. i == 11) then
   !                !    print *, "before update in vertdiff top layer o2 cbed%f_o2(", i, ",", j, ",1) = ", cbed_field(i,j,1)
   !                !    print *, "cobalt btm_o2 mol/kg = ", cobalt%btm_o2(i,j)
   !                !    print *, "cobalt Rho_0 = ", cobalt%Rho_0
   !                ! endif

   !                sfc_src = (VF(i,j,1)*D(i,j,1)*((cobalt%btm_o2(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0)) + VF(i,j,1)*w(i,j,1)*(cobalt%btm_o2(i,j)*cobalt%Rho_0))*dt ! top flux (diffuvive flux + advective flux)

   !                if (j == 100 .and. i == 42) then
   !                   print *, "sfc_src o2 (", i, ",", j, ",1) = ", sfc_src
   !                   print *, "diff flux part o2 = ", VF(i,j,1)*D(i,j,1)*((cobalt%btm_o2(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt
   !                   print *, "adv flux part o2 = ", VF(i,j,1)*w(i,j,1)*(cobalt%btm_o2(i,j)*cobalt%Rho_0)*dt
   !                   print *, "cobalt btm_o2 = ", (cobalt%btm_o2(i,j)*cobalt%Rho_0)
   !                   print *, "cobalt btm_o2 mol/kg = ", cobalt%btm_o2(i,j)
   !                   print *, "cobalt Rho_0 = ", cobalt%Rho_0
   !                endif

   !                if (j == 400 .and. i == 385) then
   !                   print *, "sfc_src o2 (", i, ",", j, ",1) = ", sfc_src
   !                   print *, "diff flux part o2 = ", VF(i,j,1)*D(i,j,1)*((cobalt%btm_o2(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt
   !                   print *, "adv flux part o2 = ", VF(i,j,1)*w(i,j,1)*(cobalt%btm_o2(i,j)*cobalt%Rho_0)*dt
   !                   print *, "cobalt btm_o2 = ", (cobalt%btm_o2(i,j)*cobalt%Rho_0)
   !                   print *, "cobalt btm_o2 mol/kg = ", cobalt%btm_o2(i,j)
   !                   print *, "cobalt Rho_0 = ", cobalt%Rho_0
   !                endif
   !                !print *, "before update in vertdiff top layer o2 cbed%f_o2(", i, ",", j, ",1) = ", cbed_field(i,j,1)

   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !                if (j == 100 .and. i == 42) then
   !                   print *, "after update in vertdiff top layer o2 cbed%f_o2(", i, ",", j, ",1) = ", cbed_field(i,j,1)
   !                   print *, "dt = ", dt
   !                endif
   !                !print *, "after update in vertdiff top layer o2 cbed%f_o2(", i, ",", j, ",1) = ", cbed_field(i,j,1)

   !                ! if (j == 5 .and. i == 11) then
   !                !    print *, "sfc_src o2 (", i, ",", j, ",1) = ", sfc_src
   !                !    print *, "diff flux part o2 = ", VF(i,j,1)*D(i,j,1)*((cobalt%btm_o2(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt
   !                !    print *, "adv flux part o2 = ", VF(i,j,1)*w(i,j,1)*(cobalt%btm_o2(i,j)*cobalt%Rho_0)*dt
   !                !    print *, "cobalt btm_o2 = ", (cobalt%btm_o2(i,j)*cobalt%Rho_0)
   !                !    print *, "after update in vertdiff top layer o2 cbed%f_o2(", i, ",", j, ",1) = ", cbed_field(i,j,1)
   !                ! endif

   !             else if ((trim(field_name) == "f_nh4")) then
   !                sfc_src =  (VF(i,j,1)*D(i,j,1)*((cobalt%f_nh4(i,j,nk)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0)) + VF(i,j,1)*w(i,j,1)*(cobalt%f_nh4(i,j,nk)*cobalt%Rho_0))*dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             else if ((trim(field_name) == "f_no3")) then
   !                sfc_src = (VF(i,j,1)*D(i,j,1)*((cobalt%btm_no3(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0)) + VF(i,j,1)*w(i,j,1)*(cobalt%btm_no3(i,j)*cobalt%Rho_0))*dt  ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             else if ((trim(field_name) == "f_dic") ) then
   !                sfc_src = (VF(i,j,1)*D(i,j,1)*((cobalt%btm_dic(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0)) + VF(i,j,1)*w(i,j,1)*(cobalt%btm_dic(i,j)*cobalt%Rho_0))*dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !                ! if (j == 10 .and. i == 14) then
   !                !    print *, "sfc_src dic (", i, ",", j, ",1) = ", sfc_src
   !                !    print *, "diff flux part dic = ", VF(i,j,1)*D(i,j,1)*((cobalt%btm_dic(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt
   !                !    print *, "adv flux part dic = ", VF(i,j,1)*w(i,j,1)*(cobalt%btm_dic(i,j)*cobalt%Rho_0)*dt
   !                !    print *, "cobalt btm_dic = ", (cobalt%btm_dic(i,j)*cobalt%Rho_0)
   !                !    print *, "cobalt btm_dic mol/kg = ", cobalt%btm_dic(i,j)
   !                !    print *, "cobalt Rho_0 = ", cobalt%Rho_0
   !                ! endif

   !             else if ((trim(field_name) == "f_odu") ) then
   !                sfc_src = VF(i,j,1)*D(i,j,1)*((0.0-cbed_field(i,j,1))/(dz_cbed(1)/2.0)) *dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             else if ((trim(field_name) == "f_talk") ) then
   !                sfc_src = (VF(i,j,1)*D(i,j,1)*((cobalt%btm_alk(i,j)*cobalt%Rho_0 - cbed_field(i,j,1))/(dz_cbed(1)/2.0)) + VF(i,j,1)*w(i,j,1)*(cobalt%btm_alk(i,j)*cobalt%Rho_0))*dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             else if ((trim(field_name) == "f_om1")) then
   !                sfc_src = frac_OM1*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             else if ((trim(field_name) == "f_om2")) then
   !                sfc_src = frac_OM2*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             else if ((trim(field_name) == "f_om3")) then
   !                sfc_src = frac_OM3*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt ! top flux
   !                cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

   !             endif

   !             ! bottom flux
   !             btm_src = 0.0
   !             btm_src = -(cbed_field(i,j,nk_cbed)*VF(i,j,nk_cbed)*w(i,j,nk_cbed+1)*dt) ! bottom flux
   !             cbed_field(i,j,nk_cbed) = cbed_field(i,j,nk_cbed) + btm_src/h_old(nk_cbed)



   !             do k=1,nk_cbed
   !                if (k == 1) then

   !                   a(k)= 0.0
   !                   b(k)= (h_old(k)+eb(i,j,k)+sink(i,j,k+1))/h_old(k)
   !                   c(k)= -eb(i,j,k)/h_old(k)

   !                elseif (k == nk_cbed) then

   !                   a(k)= (-ea(i,j,k)-sink(i,j,k))/h_old(k)
   !                   b(k)= (h_old(k)+ea(i,j,k))/h_old(k)
   !                   c(k)= 0.0

   !                else

   !                   a(k)= (-ea(i,j,k)-sink(i,j,k))/h_old(k)
   !                   b(k)= (h_old(k)+eb(i,j,k)+ea(i,j,k)+sink(i,j,k+1))/h_old(k)
   !                   c(k)= -eb(i,j,k)/h_old(k)

   !                   !f_old(k)= cbed_field(i,j,k)

   !                endif

   !                f_old(k)= cbed_field(i,j,k)

   !             enddo

   !             call CBED_tridag_solver_Press_et_al(a,b,c,f_old,cbed_field(i,j,:),nk_cbed)
   !          endif
   !       enddo; enddo

   ! end subroutine vertdiff_CBED

   ! new
   subroutine vertdiff_CBED(cobalt_tracer_list,cobalt, cbed_field, field_name, D, w, VF, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      type(g_tracer_type),          pointer       :: cobalt_tracer_list
      type(generic_COBALT_type),    intent(inout) :: cobalt
      real, dimension(:,:,:),       intent(inout) :: cbed_field  ! cbed tracer concentration field
      character(len=*),             intent(in)    :: field_name   !Name of the cbed field (e.g., "f_o2" or "f_nh4")
      real, dimension(:,:,:),       intent(in)    :: D   ! diffustion
      real, dimension(:,:,:),       intent(in)    :: w   !sinking velocity or sedimentation rate
      real, dimension(:,:,:),       intent(in)    :: VF    ! volumn fraction
      integer, dimension(:,:),      intent(in)    :: grid_kmt
      real,                         intent(in)    :: dt
      integer,                      intent(in)    :: tau
      integer,                      intent(in)    :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed

      !Locals
      integer :: i, j, k

      real, dimension(nk_cbed) :: a,b,c,f_old,h_old
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: ea,eb
      real, dimension(isc:iec,jsc:jec,nk_cbed+1) :: sink

      ! CHANGE 1: Added variables for Robin boundary condition
      ! These implement the implicit boundary flux in the matrix system
      real :: delta_z, alpha, alpha_dt_over_h
      real, dimension(isc:iec,jsc:jec) :: btm_tracer_conc  ! Bottom water concentration for the tracer

      real :: sfc_src, btm_src

      ! local parameters for bgc reactions
      real, parameter :: frac_OM1 = 0.70
      real, parameter :: frac_OM2 = 0.20
      real, parameter :: frac_OM3 = 0.10




      ! h_old
      do k = 1, nk_cbed
         h_old(k) = dz_cbed(k)
      enddo

      ! ea , eb
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k=1,nk_cbed
                  ea(i,j,k) = VF(i,j,k)*D(i,j,k)*dt/h_old(k)
                  eb(i,j,k) = VF(i,j,k)*D(i,j,k+1)*dt/h_old(k)
               enddo
            endif
         enddo; enddo

      ! sink(k+1)
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k=1,nk_cbed+1
                  sink(i,j,k) = max(0.0, VF(i,j,k)*w(i,j,k)*dt )
               enddo
            endif
         enddo; enddo

      ! CHANGE 3: Get bottom water concentration based on tracer type
      ! This will be used in the Robin boundary condition
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               btm_tracer_conc(i,j) = 0.0
               if ((trim(field_name) == "f_o2")) then
                  btm_tracer_conc(i,j) = cobalt%btm_o2(i,j)*cobalt%Rho_0
               else if ((trim(field_name) == "f_nh4")) then
                  btm_tracer_conc(i,j) = cobalt%f_nh4(i,j,nk)*cobalt%Rho_0
               else if ((trim(field_name) == "f_no3")) then
                  btm_tracer_conc(i,j) = cobalt%btm_no3(i,j)*cobalt%Rho_0
               else if ((trim(field_name) == "f_dic")) then
                  btm_tracer_conc(i,j) = cobalt%btm_dic(i,j)*cobalt%Rho_0
               else if ((trim(field_name) == "f_odu")) then
                  btm_tracer_conc(i,j) = 0.0
               else if ((trim(field_name) == "f_talk")) then
                  btm_tracer_conc(i,j) = cobalt%btm_alk(i,j)*cobalt%Rho_0
               else if ((trim(field_name) == "f_om1")) then
                  ! CHANGE 4: For particulate organic matter, add source to RHS
                  ! This is handled separately below in f_old(1)
                  btm_tracer_conc(i,j) = 0.0
               else if ((trim(field_name) == "f_om2")) then
                  btm_tracer_conc(i,j) = 0.0
               else if ((trim(field_name) == "f_om3")) then
                  btm_tracer_conc(i,j) = 0.0
               endif
            endif
         enddo; enddo

      ! a, b, c, f_old
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then

               ! CHANGE 2: REMOVED explicit flux application before tridiagonal solve
               ! OLD CODE (INCORRECT):
               ! btm_src = -(cbed_field(i,j,nk_cbed)*VF(i,j,nk_cbed)*w(i,j,nk_cbed+1)*dt)
               ! cbed_field(i,j,nk_cbed) = cbed_field(i,j,nk_cbed) + btm_src/h_old(nk_cbed)
               !
               ! WHY REMOVED: Bottom advective loss is already handled by the sink term
               ! in the tridiagonal matrix (sink(nk_cbed+1) in the bottom layer equation).
               ! Applying it explicitly AND implicitly double-counts the loss.

               ! CHANGE 6: Build tridiagonal matrix with proper boundary conditions
               do k=1,nk_cbed

                  if (k == 1) then
                     ! CHANGE 7: TOP LAYER - Different boundary conditions for solutes vs solids

                     ! Matrix coefficients for top layer:
                     a(1) = 0.0  ! No layer above

                     ! Start with base formulation (internal diffusion + advection to layer 2)
                     b(1) = (h_old(1) + eb(i,j,1) + sink(i,j,2))/h_old(1)
                     c(1) = -eb(i,j,1)/h_old(1)

                     ! Determine if this is a SOLUTE (needs Robin BC) or SOLID (no top BC)
                     if ((trim(field_name) /= "f_om1") .and. &
                        (trim(field_name) /= "f_om2") .and. &
                        (trim(field_name) /= "f_om3")) then

                        ! SOLUTES: Add Robin boundary condition
                        ! Implements: -D*VF*(dC/dz)|_interface = D*VF*(C_bottom_water - C_sed(1))/(dz/2)
                        ! Plus advection: w*VF*C_bottom_water

                        delta_z = dz_cbed(1) / 2.0
                        alpha = VF(i,j,1) * D(i,j,1) / delta_z
                        alpha_dt_over_h = alpha * dt / h_old(1)

                        ! CHANGE 8: Add Robin boundary term to diagonal
                        ! This couples sediment surface to bottom water concentration
                        b(1) = b(1) + alpha_dt_over_h

                        ! CHANGE 9: Advection INTO top layer from bottom water
                        ! This represents w*VF*C_bottom_water entering from above
                        if (w(i,j,1) > 0.0) then
                           b(1) = b(1) + sink(i,j,1)/h_old(1)
                        endif
                     else
                        ! SOLIDS: No Robin BC at top (no exchange with bottom water)
                        ! Particles don't couple to bottom water concentration
                        ! Only internal bioturbation (eb term already included above)
                        ! No advection IN because particles are delivered as external source
                     endif

                     ! RHS: old concentration
                     f_old(1) = cbed_field(i,j,1)

                     ! CHANGE 10: Add boundary contributions based on tracer type
                     if ((trim(field_name) /= "f_om1") .and. &
                        (trim(field_name) /= "f_om2") .and. &
                        (trim(field_name) /= "f_om3")) then

                        ! SOLUTES: Add Robin boundary contribution
                        ! This brings in the influence of bottom water concentration
                        f_old(1) = f_old(1) + alpha_dt_over_h * btm_tracer_conc(i,j)
                        !f_old(1) = f_old(1) + alpha_dt_over_h * (btm_tracer_conc(i,j) - cbed_field(i,j,1))   ! This would be the full Robin BC contribution, but since cbed_field(i,j,1) is on the LHS, we only add the bottom water part to the RHS. The flux term is split in two parts in the above. 

                        ! Add advective flux from bottom water
                        if (w(i,j,1) > 0.0) then
                           f_old(1) = f_old(1) + (sink(i,j,1)/h_old(1)) * btm_tracer_conc(i,j)
                        endif

                     else
                        ! SOLIDS: Add particulate organic matter source terms
                        ! These are external rain fluxes from the water column
                        if ((trim(field_name) == "f_om1")) then
                           f_old(1) = f_old(1) + frac_OM1*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt/h_old(1)
                        else if ((trim(field_name) == "f_om2")) then
                           f_old(1) = f_old(1) + frac_OM2*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt/h_old(1)
                        else if ((trim(field_name) == "f_om3")) then
                           f_old(1) = f_old(1) + frac_OM3*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt/h_old(1)
                        endif
                     endif

                  else if (k == nk_cbed) then
                     ! CHANGE 13: BOTTOM LAYER - Zero flux at bottom boundary
                     ! Advective loss through bottom is handled by sink(nk_cbed+1)

                     a(k) = (-ea(i,j,k) - sink(i,j,k))/h_old(k)

                     ! CHANGE 14: Modified b(k) for bottom layer
                     ! OLD: b(k) = (h_old(k) + eb(k) + ea(k) + sink(k+1))/h_old(k)
                     ! NEW: No eb term (no diffusion out of bottom), but keep sink(k+1) for advection out
                     b(k) = (h_old(k) + ea(i,j,k) + sink(i,j,k+1))/h_old(k)

                     c(k) = 0.0  ! No layer below

                     f_old(k) = cbed_field(i,j,k)

                  else
                     ! CHANGE 15: INTERIOR LAYERS - Standard advection-diffusion
                     ! sink(k) represents advection FROM layer k-1 (appears in super-diagonal of k-1)
                     ! sink(k+1) represents advection TO layer k+1 (appears in diagonal of k)

                     a(k) = (-ea(i,j,k) - sink(i,j,k))/h_old(k)

                     b(k) = (h_old(k) + ea(i,j,k) + eb(i,j,k) + sink(i,j,k+1))/h_old(k)

                     c(k) = -eb(i,j,k)/h_old(k)

                     f_old(k) = cbed_field(i,j,k)
                  endif

               enddo

               call CBED_tridag_solver_Press_et_al(a,b,c,f_old,cbed_field(i,j,:),nk_cbed)
            endif
         enddo; enddo

   end subroutine vertdiff_CBED

!!!! Copied from Niki's code.
   subroutine CBED_tridag_solver_Press_et_al(a,b,c,r,u,n)
      integer, intent(in) :: n
      real,    intent(in) :: a(n),b(n),c(n),r(n)
      real,    intent(inout) :: u(n)
      real    :: bet,gam(n)
      integer :: k
      bet=b(1)
      u(1)=r(1)/bet
      do k=2,n
         gam(k)=c(k-1)/bet
         bet=b(k)-a(k)*gam(k)
         u(k)=(r(k)-a(k)*u(k-1))/bet
      enddo
      do k=n-1,1,-1
         u(k)=u(k)-gam(k+1)*u(k+1)
      enddo
   end subroutine CBED_tridag_solver_Press_et_al



   subroutine generic_CBED_sediments_update_from_source(cobalt_tracer_list, cobalt, phyto, ilb, jlb, mask_coast, &
      grid_tmask, grid_dat, grid_kmt, isc,iec, jsc,jec, isd,ied, jsd,jed, nk, r_dt, dt, tau, model_time, frunoff, rho_dzt, dzt, internal_heat)

      type(g_tracer_type),          pointer       :: cobalt_tracer_list
      type(generic_COBALT_type),    intent(inout) :: cobalt
      type(phytoplankton), dimension(NUM_PHYTO), intent(inout) :: phyto
      integer,                      intent(in)    :: ilb, jlb
      real, dimension(ilb:,jlb:),   intent(in)    :: grid_dat
      real, dimension(:,:,:),       intent(in)    :: grid_tmask
      integer, dimension(:,:),      intent(in)    :: mask_coast, grid_kmt
      integer,                      intent(in)    :: isc,iec, jsc,jec, isd,ied, jsd,jed, nk
      real,                         intent(in)    :: r_dt, dt
      integer,                      intent(in)    :: tau
      type(time_type),              intent(in)    :: model_time
      real, dimension(ilb:,jlb:),   intent(in)    :: frunoff
      real, dimension(ilb:,jlb:,:), intent(in)    :: rho_dzt, dzt
      real, dimension(ilb:,jlb:),   intent(in), optional :: internal_heat

      integer :: i, j, k
      real :: fpoc_btm, drho_dzt, log10_fpoc_btm
      integer, dimension(isc:iec,jsc:jec) :: k_bot
      real,    dimension(isc:iec,jsc:jec) :: rho_dzt_bot

      ! local parameters for bgc reactions
      real, parameter :: frac_OM1 = 0.70
      real, parameter :: frac_OM2 = 0.20
      real, parameter :: frac_OM3 = 0.10

      real, parameter :: k_adj_denit = 0.1
      real, parameter :: k_adj_anoxia = 0.005

      real, parameter :: ks_o2 = 0.008   ! O2 half saturation constant (mol/m3)
      real, parameter :: ks_no3 = 0.001  ! NO3 half saturation constant (mol/m3)

      real, parameter :: k_nox = (2.0*10.0**5.0)/spery   ! 2e5 ! mol-1 m3 s-1 (from the original: mmol-1 L yr-1) !nitrification rate constant
      real, parameter :: k_ana =  0.0 * (10.0**3.0) /spery     !   ! 1e5               !anammox rate constant
      real, parameter :: k_oduox = (10.0**5.0) /spery   !              1e6      !ODU oxidation rate constant

      real, parameter :: Q10 = 1.88
      real, dimension(isc:iec,jsc:jec) :: Q10_factor

      ! Reaction rates
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_o2, R_om2_o2, R_om3_o2
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_no3, R_om2_no3, R_om3_no3
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_anoxic, R_om2_anoxic, R_om3_anoxic
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_dic_om1, R_dic_om2, R_dic_om3
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_nox, R_ana, R_oduox, odu_depo
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_talk

      ! b terms
      real, dimension(isc:iec,jsc:jec) :: b_o2, b_dic, b_nh4, b_no3


      ! write the grid layers to register as diags
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed
                  cbed%dz_cbed(i,j,k) = dz_cbed(k)
                  cbed%z_cbed_mid(i,j,k) = z_cbed_mid(k)
               enddo
            endif
         enddo; enddo

      ! calculate Q10 factor
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               Q10_factor(i,j) = Q10**( (cobalt%btm_temp(i,j)-4.0)/10.0 )
            endif
         enddo; enddo

      ! Porosity and solid volume fraction read from file.
      if (cbed%read_porosity_from_file) then
         call data_override('OCN', 'por', por(isc:iec,jsc:jec,1), model_time)

         do j = jsc, jec; do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  do k = 2, nk_cbed+1
                     por(i,j,k) = por(i,j,1)
                  enddo
               endif
            enddo; enddo
         do j = jsc, jec; do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  do k = 1, nk_cbed+1
                     svf(i,j,k) = 1.0 - por(i,j,k)
                  enddo
               endif
            enddo; enddo
      else
         do j = jsc, jec; do i = isc, iec
               if (grid_kmt(i,j) .gt. 0) then
                  do k = 1, nk_cbed+1
                     por(i,j,k) = 0.8
                     svf(i,j,k) = 0.2
                  enddo
               endif
            enddo; enddo
         !por = 0.8
         !svf = 0.2
      endif

      ! Sedimentation rate calculation
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed+1

                  ! w (sedimentation rate, cm/year) m/s
                  w(i,j,k) = ( (cobalt%fcadet_arag_btm(i,j)*100.0/2.71 + &
                     cobalt%fcadet_calc_btm(i,j)*100.0/2.94 + &
                     cobalt%fsitot_btm(i,j)*60.0/2.65 + &
                     cobalt%flithdet_btm(i,j)/2.65 + &
                     cobalt%ffetot_btm(i,j)*160.0/5.24 + &
                     cobalt%fptot_btm(i,j)*120.0/2.3 + &
                     cobalt%fntot_btm(i,j)*cobalt%c_2_n*22.4/0.9)/10000.0*spery/svf(i,j,k) )/100.0/spery
               enddo
            endif
         enddo;enddo

      !Bioturbation
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               ! relation from Archer. POC flux unit in umol cm-2 y-1.
               Db_0(i,j) = ( 0.0232*((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**0.85) ) /1e4/spery ! in cobalt unit m2/s
            endif
         enddo;enddo

      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed+1
                  ! relation from Archer. POC flux unit in umol cm-2 y-1.
                  Db(i,j,k) = max(0.0, Db_0(i,j)*exp(-(z_cbed_int(k)/Db_l)**2)*(cobalt%btm_o2(i,j)*cobalt%Rho_0/(cobalt%btm_o2(i,j)*cobalt%Rho_0+(20/1e3))) )
               enddo
            endif
         enddo;enddo

      !Bioirrigation
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               ! relation from Archer. POC flux unit in umol cm-2 y-1.
               bioirri_0(i,j) = ( 11*(((atan((5*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery) -400)/400))/pi)+0.5) &
                  - 0.9 + 20*((cobalt%btm_o2(i,j)*cobalt%Rho_0)/(cobalt%btm_o2(i,j)*cobalt%Rho_0+0.01)) * exp(-cobalt%btm_o2(i,j)*cobalt%Rho_0/0.01) * &
                  ((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)/((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)+30)) )/spery   ! in cobalt unit s^-1

            endif
         enddo;enddo

      do j = jsc, jec; do i = isc, iec
            do k = 1, nk_cbed
               if (grid_kmt(i,j) .gt. 0) then
                  ! relation from Archer. POC flux unit in umol cm-2 y-1.
                  bioirri(i,j,k) = 0.0 ! max(0.0, bioirri_0(i,j)*exp(-(z_cbed_mid(k)/bioirri_l)**2) )
               endif
            enddo
         enddo;enddo


      !Calculate diffusion coefficients
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed+1
                  D_o2(i,j,k)  = ( (0.031558+0.001428*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)   ! m2/s
                  D_dic(i,j,k) = ( (0.015179+0.000795*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_nh4(i,j,k) = ( (0.030926+0.001225*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_no3(i,j,k) = ( (0.030863+0.001153*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_odu(i,j,k) = ( (0.028938+0.001314*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
               enddo
            endif
         enddo;enddo


      !Calculate k1,k2,k3
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               ! POC flux unit in umol cm-2 y-1. Unit of k is y-1
               k1(i,j) = ( 0.15*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**(0.85) )/spery
               k2(i,j) = ( 0.0023*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**(0.85) )/spery
               k3(i,j) = ( 0.00013*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**(0.85) )/spery
            endif
         enddo;enddo



      ! ! calculate reaction rates
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed

                  ! O₂ reaction rates
                  R_om1_o2(i,j,k) = k1(i,j)*cbed%f_om1(i,j,k)*(cbed%f_o2(i,j,k)/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  R_om2_o2(i,j,k) = k2(i,j)*cbed%f_om2(i,j,k)*(cbed%f_o2(i,j,k)/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  R_om3_o2(i,j,k) = k3(i,j)*cbed%f_om3(i,j,k)*(cbed%f_o2(i,j,k)/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  ! NO₃ reaction rates
                  R_om1_no3(i,j,k) = k_adj_denit*k1(i,j)*cbed%f_om1(i,j,k)*(cbed%f_no3(i,j,k)/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  R_om2_no3(i,j,k) = k_adj_denit*k2(i,j)*cbed%f_om2(i,j,k)*(cbed%f_no3(i,j,k)/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  R_om3_no3(i,j,k) = k_adj_denit*k3(i,j)*cbed%f_om3(i,j,k)*(cbed%f_no3(i,j,k)/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  ! ODU reaction rates
                  R_om1_anoxic(i,j,k) = k_adj_anoxia*k1(i,j)*cbed%f_om1(i,j,k)*(ks_no3/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  R_om2_anoxic(i,j,k) = k_adj_anoxia*k2(i,j)*cbed%f_om2(i,j,k)*(ks_no3/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)
                  R_om3_anoxic(i,j,k) = k_adj_anoxia*k3(i,j)*cbed%f_om3(i,j,k)*(ks_no3/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k))) * Q10_factor(i,j)

                  ! dic
                  R_dic_om1(i,j,k) = (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k))
                  R_dic_om2(i,j,k) = (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k))
                  R_dic_om3(i,j,k) = (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k))
                  ! nitrification
                  R_nox(i,j,k) = k_nox*cbed%f_nh4(i,j,k)*cbed%f_o2(i,j,k) * Q10_factor(i,j)
                  ! anammox
                  R_ana(i,j,k) = k_ana*cbed%f_nh4(i,j,k)*cbed%f_no3(i,j,k) * Q10_factor(i,j)
                  ! ODU oxidation
                  R_oduox(i,j,k) = k_oduox*cbed%f_odu(i,j,k)*cbed%f_o2(i,j,k) * Q10_factor(i,j)
                  odu_depo(i,j,k) = (R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))*min(1.0, 0.233*(w(i,j,k)*100.0*spery)**0.336)

                  ! TA calculation
                  R_talk(i,j,k) = svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) + &
                     svf(i,j,k)/por(i,j,k)*(0.8+1.0/cobalt%c_2_n)*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
                     svf(i,j,k)/por(i,j,k)*(1.0+1.0/cobalt%c_2_n)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
                     2.0*R_nox(i,j,k) - 1.0*R_oduox(i,j,k)

                  ! calculations for diagnostics
                  cbed%R_om_o2(i,j,k) = R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)
                  cbed%R_om_no3(i,j,k) = R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)
                  cbed%R_om_anaerobic(i,j,k) = R_om1_anoxic(i,j,k) + R_om2_anoxic(i,j,k) + R_om3_anoxic(i,j,k)
                  cbed%R_dic(i,j,k) = R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)

               enddo
            endif
         enddo;enddo

      ! Calculate the "b terms" to feed into cobalt.
      ! The b terms are calulated based on the t-1 time step. This is to maintain mass balance
      ! with cobalt-cbed. The surface flux calculated in the vertdiff_CBED subroutine
      ! is the flux at t-1, so we need to use the t-1 concentration field to calculate the b term,
      ! which will update the bottom flux as the verdiff_G is not yet called for cobalt.
      ! The CBED reaction rates however are calculated based on the t time step, so they will be
      ! calculated after source sink calculation and call to vertdiff_CBED.
      !
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               b_o2(i,j) = por(i,j,1)*D_o2(i,j,1)*((cobalt%btm_o2(i,j)*cobalt%Rho_0 - cbed%f_o2(i,j,1))/(dz_cbed(1)/2.0)) + por(i,j,1)*w(i,j,1)*(cobalt%btm_o2(i,j)*cobalt%Rho_0)
               b_nh4(i,j) = por(i,j,1)*D_nh4(i,j,1)*((cobalt%f_nh4(i,j,nk)*cobalt%Rho_0 - cbed%f_nh4(i,j,1))/(dz_cbed(1)/2.0)) + por(i,j,1)*w(i,j,1)*(cobalt%f_nh4(i,j,nk)*cobalt%Rho_0)
               b_no3(i,j) = por(i,j,1)*D_no3(i,j,1)*((cobalt%btm_no3(i,j)*cobalt%Rho_0 - cbed%f_no3(i,j,1))/(dz_cbed(1)/2.0)) + por(i,j,1)*w(i,j,1)*(cobalt%btm_no3(i,j)*cobalt%Rho_0)
               b_dic(i,j) = por(i,j,1)*D_dic(i,j,1)*((cobalt%btm_dic(i,j)*cobalt%Rho_0 - cbed%f_dic(i,j,1))/(dz_cbed(1)/2.0)) + por(i,j,1)*w(i,j,1)*(cobalt%btm_dic(i,j)*cobalt%Rho_0)

               if (j == 10 .and. i == 14) then
                  print *, "b_o2 (", i, ",", j, ",1) = ", b_o2(i,j)
                  print *, "b_nh4 (", i, ",", j, ",1) = ", b_nh4(i,j)
                  print *, "b_no3 (", i, ",", j, ",1) = ", b_no3(i,j)
                  print *, "b_dic (", i, ",", j, ",1) = ", b_dic(i,j)
                  print *, "cobalt btm_dic mol/kg = ", cobalt%btm_dic(i,j)
                  print *, "cobalt btm_o2 mol/kg = ", cobalt%btm_o2(i,j)
                  print *, "cobalt btm_nh4 mol/kg = ", (cobalt%f_nh4(i,j,nk))
                  print *, "cobalt btm_no3 mol/kg = ", (cobalt%btm_no3(i,j))
                  print *, "cobalt Rho_0 = ", cobalt%Rho_0
               endif

               if (j == 100 .and. i == 42) then
                  print *, "b_o2 (", i, ",", j, ",1) = ", b_o2(i,j)
                  print *, "b_nh4 (", i, ",", j, ",1) = ", b_nh4(i,j)
                  print *, "b_no3 (", i, ",", j, ",1) = ", b_no3(i,j)
                  print *, "b_dic (", i, ",", j, ",1) = ", b_dic(i,j)
                  print *, "cobalt btm_dic mol/kg = ", cobalt%btm_dic(i,j)
                  print *, "cobalt btm_o2 mol/kg = ", cobalt%btm_o2(i,j)
                  print *, "cobalt btm_nh4 mol/kg = ", (cobalt%f_nh4(i,j,nk))
                  print *, "cobalt btm_no3 mol/kg = ", (cobalt%btm_no3(i,j))
                  print *, "cobalt Rho_0 = ", cobalt%Rho_0
               endif


            endif
         enddo;enddo

      ! save the benthic fluxes as diagnostics. 2D diag
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               cbed%o2_flux(i,j)  = b_o2(i,j)
               cbed%nh4_flux(i,j) = b_nh4(i,j)
               cbed%no3_flux(i,j) = b_no3(i,j)
               cbed%dic_flux(i,j) = b_dic(i,j)
            endif
         enddo;enddo


      !TOC (total organic carbon) diagnostics
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed
                  cbed%TOC(i,j,k)  = (cbed%f_om1(i,j,k)+cbed%f_om2(i,j,k)+cbed%f_om3(i,j,k))*12.0/1e6/rho_s*100.0
               enddo
            endif
         enddo;enddo

      ! some other diags
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then

               cbed%burial_om(i,j)  = (cbed%f_om1(i,j,nk_cbed)+cbed%f_om2(i,j,nk_cbed)+cbed%f_om3(i,j,nk_cbed))*w(i,j,nk_cbed+1) * svf(i,j,nk_cbed+1) ! mol/m2/s

               cbed%denit(i,j) = sum(dz_cbed(:)*(svf(i,j,1:nk_cbed)*0.8*cbed%R_om_no3(i,j,:) + por(i,j,1:nk_cbed)*2.0*R_ana(i,j,:)))

               cbed%cbed_k1(i,j) = cobalt%f_o2(i,j,nk) !k1(i,j)
               cbed%cbed_k2(i,j) = cobalt%btm_o2(i,j) !k2(i,j)
               cbed%cbed_k3(i,j) = cobalt%f_no3(i,j,nk)*cobalt%Rho_0 !k3(i,j)
               cbed%cbed_w(i,j) = cobalt%btm_no3(i,j)   !w(i,j,1)

            endif
         enddo;enddo

      ! some more diags
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed+1
                  cbed%cbed_Db(i,j,k) = Db(i,j,k)
                  cbed%cbed_bioirri(i,j,k) = bioirri(i,j,k)
                  cbed%cbed_D_o2(i,j,k) = D_o2(i,j,k)
                  cbed%cbed_D_nh4(i,j,k) = D_nh4(i,j,k)
                  cbed%cbed_D_no3(i,j,k) = D_no3(i,j,k)
                  cbed%cbed_D_dic(i,j,k) = D_dic(i,j,k)
                  cbed%cbed_D_odu(i,j,k) = D_odu(i,j,k)
                  cbed%cbed_por(i,j,k) = por(i,j,k)
                  cbed%cbed_svf(i,j,k) = svf(i,j,k)
               enddo
            endif
         enddo;enddo



      ! Source-sink calculations
      !Test that we can change the value of concentration field of a CBED tracer
      do j = jsc, jec; do i = isc, iec  !{
            if (grid_kmt(i,j) .gt. 0) then
               do k = 1, nk_cbed

                  !cbed%f_tr1(i,j,k) = cbed%f_tr1(i,j,k) + 0.01 * k !fictitious dubious dynamics for testing purposes

                  if (j == 100 .and. i == 42) then
                     print *, "before update o2 cbed%f_o2(", i, ",", j, ",", k, ") = ", cbed%f_o2(i,j,k)
                  endif

                  cbed%f_o2(i,j,k)  = max(0.0, cbed%f_o2(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)) - &
                     (2.0*R_nox(i,j,k)+R_oduox(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_o2(i,j) - cbed%f_o2(i,j,k)) )*dt )

                  if (j == 100 .and. i == 42) then
                     print *, "after update o2 cbed%f_o2(", i, ",", j, ",", k, ") = ", cbed%f_o2(i,j,k)
                  endif
                  !print *, "after update o2 cbed%f_o2(", i, ",", j, ",", k, ") = ", cbed%f_o2(i,j,k)

                  cbed%f_om1(i,j,k) = max(0.0, cbed%f_om1(i,j,k) + ( - (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k)) )*dt )

                  cbed%f_om2(i,j,k) = max(0.0, cbed%f_om2(i,j,k) + ( - (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k)) )*dt )

                  cbed%f_om3(i,j,k) = max(0.0, cbed%f_om3(i,j,k) + ( - (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k)) )*dt )

                  cbed%f_nh4(i,j,k) = max(0.0, cbed%f_nh4(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(1.0/cobalt%c_2_n)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                     ( - R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%f_nh4(i,j,nk) - cbed%f_nh4(i,j,k)) )*dt )

                  !print *, "before update no3 cbed%f_no3(", i, ",", j, ",", k, ") = ", cbed%f_no3(i,j,k)

                  cbed%f_no3(i,j,k) = max(0.0, cbed%f_no3(i,j,k) + ( - svf(i,j,k)/por(i,j,k)*0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)) + &
                     (R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_no3(i,j) - cbed%f_no3(i,j,k)) )*dt )

                  !print *, "after update no3 cbed%f_no3(", i, ",", j, ",", k, ") = ", cbed%f_no3(i,j,k)

                  !print *, "before update dic cbed%f_dic(", i, ",", j, ",", k, ") = ", cbed%f_dic(i,j,k)

                  cbed%f_dic(i,j,k) = max(0.0, cbed%f_dic(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)) + &
                     bioirri(i,j,k)*(cobalt%btm_dic(i,j) - cbed%f_dic(i,j,k)) )*dt )

                  !print *, "after update dic cbed%f_dic(", i, ",", j, ",", k, ") = ", cbed%f_dic(i,j,k)

                  cbed%f_odu(i,j,k) = max(0.0, cbed%f_odu(i,j,k) + ( + svf(i,j,k)/por(i,j,k)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k)) - &
                     R_oduox(i,j,k) - odu_depo(i,j,k)  + bioirri(i,j,k)*(0.0 - cbed%f_odu(i,j,k)) )*dt )

                  cbed%f_talk(i,j,k) = max(0.0, cbed%f_talk(i,j,k) + ( + R_talk(i,j,k) + bioirri(i,j,k)*(cobalt%btm_alk(i,j) - cbed%f_talk(i,j,k)) )*dt )


               enddo
            endif
         enddo;enddo


      ! call vertdiff_CBED. This updates the fields.
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om1, "f_om1", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om2, "f_om2", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_om3, "f_om3", Db,    w, svf, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_o2,  "f_o2", D_o2,   w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_nh4, "f_nh4", D_nh4, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_no3, "f_no3", D_no3, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_dic, "f_dic", D_dic, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_odu, "f_odu", D_odu, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)
      call vertdiff_CBED(cobalt_tracer_list,cobalt, cbed%f_talk, "f_talk", D_dic, w, por, grid_kmt, dt, tau, isc,iec,jsc,jec,isd,ied,jsd,jed,nk, nk_cbed)




      !!==================================================================================================================
      !!The rest of this subrouine that follows is a copy of the COBALT code.
      !!It must be replaced by CBED calculations for 'btf' fluxes which are "set" for COBALT at the end of this subroutine.
      !!==================================================================================================================

      ! Calculate the bottom conditions and the fluxes to the bottom for diagnostics and benthic flux calculations.
      ! MOM4/5 used the bottom grid cell, but MOM6 often has a number of vanishingly thin layers overlying the bottom.
      ! Grid scale noise in these layers can occur, particularly for quantities with large bottom fluxes.  COBALT thus
      ! uses conditions over a specified bottom layer thickness (cobalt%bottom_thickness, default = 1m) for bottom calcs.

      do j = jsc, jec; do i = isc, iec  !{
            if (grid_kmt(i,j) .gt. 0) then !{

               ! Add the phytoplankton fluxes to the detritus fluxes to get total flux to benthos
               cobalt%fntot_btm(i,j) = cobalt%f_ndet_btf(i,j,1) + cobalt%f_ndi_btf(i,j,1) + &
                  cobalt%f_nsm_btf(i,j,1) + cobalt%f_nmd_btf(i,j,1) + cobalt%f_nlg_btf(i,j,1)
               cobalt%fptot_btm(i,j) = cobalt%f_pdet_btf(i,j,1) + cobalt%f_pdi_btf(i,j,1) + &
                  cobalt%f_psm_btf(i,j,1) + cobalt%f_pmd_btf(i,j,1) + cobalt%f_plg_btf(i,j,1)
               cobalt%ffetot_btm(i,j) = cobalt%f_fedet_btf(i,j,1) + cobalt%f_fedi_btf(i,j,1) + &
                  cobalt%f_fesm_btf(i,j,1) + cobalt%f_femd_btf(i,j,1) + cobalt%f_felg_btf(i,j,1)
               cobalt%fsitot_btm(i,j) = cobalt%f_sidet_btf(i,j,1) + cobalt%f_silg_btf(i,j,1) + &
                  cobalt%f_simd_btf(i,j,1)

               ! Calculate the values of tracers influencing the sedimentary transformations
               ! and fluxes over a layer defined by "bottom_thickess".
               rho_dzt_bot(i,j) = 0.0
               cobalt%btm_o2(i,j) = 0.0
               cobalt%btm_no3(i,j) = 0.0
               cobalt%btm_co3_sol_calc(i,j) = 0.0
               cobalt%btm_co3_ion(i,j) = 0.0
               cobalt%btm_omega_calc(i,j) = 0.0
               k_bot(i,j) = 0
               ! Note that grid_kmt is always the total number of layers in MOM6
               do k = grid_kmt(i,j),1,-1   !{
                  ! Check if the top of layer k is within the bottom thickness.  If so, include its properties in the bottom
                  ! layer averages.  Overshoots will be subtracted off later.
                  if (rho_dzt_bot(i,j).lt.(cobalt%Rho_0*cobalt%bottom_thickness)) then
                     k_bot(i,j) = k
                     rho_dzt_bot(i,j) = rho_dzt_bot(i,j) + rho_dzt(i,j,k)
                     cobalt%btm_o2(i,j) = cobalt%btm_o2(i,j) + cobalt%f_o2(i,j,k)*rho_dzt(i,j,k)
                     cobalt%btm_no3(i,j) = cobalt%btm_no3(i,j) + cobalt%f_no3(i,j,k)*rho_dzt(i,j,k)
                     cobalt%btm_co3_sol_calc(i,j) = cobalt%btm_co3_sol_calc(i,j) + cobalt%co3_sol_calc(i,j,k)*rho_dzt(i,j,k)
                     cobalt%btm_co3_ion(i,j) = cobalt%btm_co3_ion(i,j) + cobalt%f_co3_ion(i,j,k)*rho_dzt(i,j,k)
                  endif
               enddo
               ! Subtract off overshoot
               drho_dzt = rho_dzt_bot(i,j) - cobalt%Rho_0*cobalt%bottom_thickness
               cobalt%btm_o2(i,j)=cobalt%btm_o2(i,j)-cobalt%f_o2(i,j,k_bot(i,j))*drho_dzt
               cobalt%btm_no3(i,j)=cobalt%btm_no3(i,j)-cobalt%f_no3(i,j,k_bot(i,j))*drho_dzt
               cobalt%btm_co3_sol_calc(i,j)=cobalt%btm_co3_sol_calc(i,j)-cobalt%co3_sol_calc(i,j,k_bot(i,j))*drho_dzt
               cobalt%btm_co3_ion(i,j)=cobalt%btm_co3_ion(i,j)-cobalt%f_co3_ion(i,j,k_bot(i,j))*drho_dzt
               ! convert back to moles kg-1
               cobalt%btm_o2(i,j)=cobalt%btm_o2(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
               cobalt%btm_no3(i,j)=cobalt%btm_no3(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
               cobalt%btm_co3_sol_calc(i,j)=cobalt%btm_co3_sol_calc(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
               cobalt%btm_co3_ion(i,j)=cobalt%btm_co3_ion(i,j)/(cobalt%bottom_thickness*cobalt%Rho_0)
               ! calculate the saturation state with respect to calcite for subsequent calculations
               cobalt%btm_omega_calc(i,j)=cobalt%btm_co3_ion(i,j)/cobalt%btm_co3_sol_calc(i,j)

               ! Calculate the processing of organic matter in the sediment.  The fate of organic matter is partitioned
               ! between burial (i.e., removal from the system), aerobic remineralization, remineralization via
               ! denitrification, and remineralization via sulfate reduction.  Note that the latter pathway is effectively
               ! a "catch all" for any other anaerobic pathway and the sulfate cycle is not explicitly modeled.
               k = grid_kmt(i,j)
               if (cobalt%fntot_btm(i,j) .gt. 0.0) then !{

                  ! The Burial flux estimates are based on Dunne et al., 2007. A synthesis of global particle export from
                  ! the surface ocean and cycling through the ocean interior and on the seafloor.  Global Biogeochemical
                  ! Cycles. Vol. 21, GB4006, doi:10.1029/2006GB002907.  See Figure 2, eq. (3).  The default units of this
                  ! relationship are mmoles C m-2 day-1, and the local variable "fpoc_btm" is used to create a bottom flux
                  ! in these units.
                  !
                  ! As described in Dunne et al., (2007) this relationship was generally developed for deeper ocean areas
                  ! and its validity in shallow areas is unclear.  Past experiments suggest that it may overestimate burial
                  ! in shallow areas, resulting in large nutrient losses that are inconsistent with observations.  The
                  ! parameter "z_burial" thus provides a depth scale (an effective "half-saturation") for ramping up burial
                  ! from 0 to its full value.
                  !
                  ! Since burial is highly uncertain and often used in global earth system simulations to balance inputs and
                  ! outputs, a dimensionless scaling factor (cobalt%scale_burial) has also been included.
                  fpoc_btm = cobalt%fntot_btm(i,j)*cobalt%c_2_n*sperd*1000.0
                  cobalt%frac_burial(i,j) = 0.013 + 0.53*fpoc_btm**2.0/((7.0+fpoc_btm)**2.0) * &
                     cobalt%zt(i,j,k) / (cobalt%z_burial + cobalt%zt(i,j,k))
                  cobalt%frac_burial(i,j) = cobalt%scale_burial*cobalt%frac_burial(i,j)
                  cobalt%fn_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%fntot_btm(i,j)
                  cobalt%fp_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%fptot_btm(i,j)

                  ! Denitrification follows Middelburg et al., 1996. Denitrification in marine sediments: a modeling study
                  ! Global Biogeochemical Cycles 10(4).  pp. 661-673.  https://doi.org/10.1029/96GB02562. COBALT uses the
                  ! carbon flux-based relationship based on Middelburg's first extraction of his metamodel (the first
                  ! equation in Section 3.4 of the paper).  This relationship requires a flux to the benthos in micromoles C
                  ! cm-2 day-1.  This means that fpoc_btm defined for the burial calculation above must be multiplied by:
                  !
                  ! 1e3 micromoles/millimole*1e-4 cm2/m2 = 0.1
                  !
                  ! to get the proper units.  The Middelburg relationship yields a rate at which arriving particulate organic
                  ! carbon is denitrified in micromoles C cm-2 day-1.  This is converted to a rate at which arriving
                  ! particulate organic nitrogen denitrified in moles N m-2 sec-1 by dividing by:
                  !
                  ! c_2_n*sperd*1e6 micromoles/mole*1e-4 cm2/m2 = c_2_n*sperd*100
                  !
                  ! The nitrate demand associated with this denitrification (fno3denit_sed) is obtained by multiplying the
                  ! resulting value by the moles of NO3 required to denitrify each mole of organic N (n_2_n_denit).
                  !
                  ! A number of limiters are applied to support global application.  First, the C flux used in the
                  ! Middelburg relationship is capped at 43.0 micromoles C cm-2 day-1 to avoid anomalous extrapolation.
                  ! Second, denitrification is slowed when bottom nitrate is low by a) scaling rates with a nitrate
                  ! half-saturation constant with (k_no3_denit), b) preventing the exhaustion of bottom nitrate over
                  ! single time step, and c) limiting the total amount of organic carbon denitrified to that arriving at
                  ! the sediment minus that which was buried. Finally, to prevent excessive denitrification in very shallow
                  ! areas, a depth scale (z_denit) was included to ramp up rates to full Middelburg values only in deeper
                  ! waters.
                  log10_fpoc_btm = log10(min(43.0,0.1*fpoc_btm))
                  cobalt%fno3denit_sed(i,j) = min(cobalt%btm_no3(i,j)*cobalt%bottom_thickness*cobalt%Rho_0*r_dt,  &
                     min((cobalt%fntot_btm(i,j)-cobalt%fn_burial(i,j))*cobalt%n_2_n_denit, &
                     10.0**(-0.9543+0.7662*log10_fpoc_btm - 0.235*log10_fpoc_btm**2.0)/(cobalt%c_2_n*sperd*100.0)* &
                     cobalt%n_2_n_denit*cobalt%btm_no3(i,j)/(cobalt%k_no3_denit + cobalt%btm_no3(i,j)))) * &
                     cobalt%zt(i,j,k) / (cobalt%z_denit + cobalt%zt(i,j,k))

                  ! Calculate the rate of organic matter degradation in the sediment after accounting for burial
                  ! and denitrification.  Two pathways are tracked:
                  !
                  ! fnoxic_sed (moles N m-2 sec-1) accounts for organic material remineralized by processes using oxygen
                  ! *within the sediment*, resulting in an oxygen demand at the sediment-water interface.  These
                  ! include direct aerobic remineralization and sulfate reduction/HS- oxidation (see stoichiometry
                  ! for details).  Note that the partitioning between these two pathways is not calculated, just the
                  ! combined effect.
                  !
                  ! fnso4_sed (moles N m-2 sec-1) accounts for organic material that only undergoes sulfate reduction in
                  ! the sediment, but not HS- oxidation.  This results in HS- released from the sediment.  The latent O2
                  ! demand from the HS- is tracked if "do_fnso4red_sed = true".  The amount of material falling into this
                  ! category is equal to the remainder after all other pathways are accounted for and it generally only
                  ! occurs under anoxic conditions in COBALT.  The added O2 demand from HS- can result in negative O2
                  ! concentrations that should be interpreted as 0 moles O2 kg-1 + additional O2 demand from HS-
                  !
                  ! NOTE: The O2 demand from the sediment is calculated by multiplying the organic matter remineralized
                  !       by the O2 per N (i.e., fnoxic_sed*o2_2_nh4 or (fnoxic_sed+fnso4red_sed)*o2_2_nh4 if
                  !       do_fnso4red_sed is true).
                  ! NOTE: fnso4red_sed is not the total sulfate reduction in the sediment, only that which is not paired
                  !       with subsequent HS- oxidation.  COBALT does not calculate the total sulfate reduction in seds.
                  ! NOTE: The maximum organic remin supported by local O2 is:
                  ! btm_o2(moles O2 kg-1)*bottom_thickness(m)*density(kg m-3)* 1/dt(s-1)*molN/molO2 = moles N m-2 s-1
                  !
                  ! The thickness of the bottom boundary layer (cobalt%bottom_thickness) impacts this upper bound.
                  ! Efforts are underway to implement a more dynamic bottom boundary layer scheme.
                  !
                  if (cobalt%btm_o2(i,j) .gt. cobalt%o2_min) then  !{
                     cobalt%fnoxic_sed(i,j) = max(0.0, min(cobalt%btm_o2(i,j)*cobalt%bottom_thickness* &
                        cobalt%Rho_0*r_dt*(1.0/cobalt%o2_2_nh4), &
                        cobalt%fntot_btm(i,j) - cobalt%fn_burial(i,j) - &
                        cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit))
                  else
                     cobalt%fnoxic_sed(i,j) = 0.0
                  endif !}
                  cobalt%fnso4red_sed(i,j) = max(0.0, cobalt%fntot_btm(i,j)-cobalt%fnoxic_sed(i,j)- &
                     cobalt%fn_burial(i,j)-cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)
               else
                  cobalt%fnso4red_sed(i,j) = 0.0
                  cobalt%fno3denit_sed(i,j) = 0.0
                  cobalt%fnoxic_sed(i,j) = 0.0
               endif !}

               !
               ! Iron flux from the sediment
               !

               ! Iron from sediment (Dale, 2015).  The maximum release from the sediment is set by ffe_sed_max.  The
               ! hyperbolic tangent requires the flux of carbon to the sediments (as mmoles m-2 day-1) in the numerator
               ! and the bottom water oxygen concentration (in microMolar units) in the denominator. Note that ffe_sed_max
               ! was converted to moles Fe m-2 sec-1 during parameter input, so ffe_sed is in moles Fe m-2 sec-1
               cobalt%ffe_sed(i,j) = cobalt%ffe_sed_max * tanh( (cobalt%fntot_btm(i,j)*cobalt%c_2_n*sperd*1.0e3)/ &
                  max(cobalt%btm_o2(i,j)*1.0e6,epsln) )

               ! Additional coastal iron (Optional, default fe_coast = 0)
               !
               ! Coarse resolution models and/or intermediate resolution models in areas with exceptionally steep bathymetry
               ! can under-represent coastal iron because they don't resolve shallow regions. An option to add iron through
               ! the vertical face of the land mass has thus been included.  The flux is posed as a fraction (fe_coast) of
               ! the sediment Fe flux (moles Fe m-2 sec-1) that would have resulted from the sinking organic matter flux and
               ! O2 level of the adjacent waters.  This is then spread across the layer mass (rho_dzt(i,j,k)) to give an input
               ! in moles Fe kg-1 sec-1. Conceptually, this can be thought of as a net iron flux resulting from the fraction
               ! of the sinking flux that would have been intercepted at shallower depths were the model resolution finer.
               ! The default value of fe_coast is 0 (i.e., only the explicitly resolved benthic flux is included).
               !
               ! Old Expression:
               ! cobalt%jfe_coast(i,j,1) = cobalt%fe_coast * mask_coast(i,j) * grid_tmask(i,j,1) / &
               !     sqrt(grid_dat(i,j))
               !
               do k = 1, nk !{
                  if (cobalt%fe_coast == 0.0) then
                     cobalt%jfe_coast(i,j,k) = 0.0
                  else
                     cobalt%jfe_coast(i,j,k) = cobalt%fe_coast*dzt(i,j,k)*mask_coast(i,j)*grid_tmask(i,j,k)* &
                        cobalt%ffe_sed_max*tanh( ( (cobalt%f_ndet(i,j,k)*cobalt%wsink+ &
                        phyto(SMALL)%f_n(i,j,k)*phyto(SMALL)%vmove(i,j,k)+ &
                        phyto(MEDIUM)%f_n(i,j,k)*phyto(MEDIUM)%vmove(i,j,k)+ &
                        phyto(LARGE)%f_n(i,j,k)*phyto(LARGE)%vmove(i,j,k)+ &
                        phyto(DIAZO)%f_n(i,j,k)*phyto(DIAZO)%vmove(i,j,k))*cobalt%c_2_n*sperd*1.0e3 )/ &
                        max(cobalt%f_o2(i,j,k)*1.0e6,epsln) )/rho_dzt(i,j,k)
                  endif
               enddo  !} k

               ! Have ffe_geotherm default to zero if the internal_heat variable
               ! needed to calculate it is not available (if geothermal heating is disabled).
               if(present(internal_heat)) then
                  cobalt%ffe_geotherm(i,j) = cobalt%ffe_geotherm_ratio*internal_heat(i,j)*4184.0/dt
               else
                  cobalt%ffe_geotherm(i,j) = 0.0
               endif

               !
               ! Calcium carbonate flux and burial, based on Dunne et al., 2012
               !
               ! phi_surfresp_cased = 0.14307   ! const for enhanced diss., surf sed respiration (dimensionless)
               ! phi_deepresp_cased = 4.1228    ! const for enhanced diss., deep sed respiration (dimensionless)
               ! alpha_cased = 2.7488 ! exponent controlling non-linearity of deep dissolution
               ! beta_cased = -2.2185 ! exponent controlling non-linearity of effective thickness
               ! gamma_cased = 0.03607/spery   ! dissolution rate constant
               ! Co_cased = 8.1e3        ! moles CaCo3 m-3 for pure calcite sediment with porosity = 0.7
               !
               ! if cased_steady is true, burial is calculated from Dunne's eq. (2) assuming dcased/dt = 0.
               ! This ensures that all the calcite bottom flux is partitioned between burial and redissolution.
               ! The steady state cased value of cased is calculated to reflect the changing bottom conditions.
               ! This influences the the partitioning of burial and redissolution over time, but there are
               ! no alkalinity changes/drifts associated with the long-term evolution of cased
               !
               ! If cased_steady is false, calcite is partitioned between dissolution, burial and evolving
               ! cased as described in Dunne et al. (2012).  The multi-century scale evolution of cased
               ! impacts alkalinity, but care must to ensure that cased starts in equilibrium with the
               ! mean ocean state to avoid unrealistic drifts.

               k = grid_kmt(i,j)

               ! Enhanced dissolution by fast respiration near the sediment surface, proportional
               ! to organic flux, moles Ca m-2 s-1, limited to a max 1/2 the instantaneous calcite flux
               cobalt%fcased_redis_surfresp(i,j)=min(0.5*cobalt%f_cadet_calc_btf(i,j,1), &
                  cobalt%phi_surfresp_cased*cobalt%fntot_btm(i,j)*cobalt%c_2_n)

               ! Ca-specific dissolution coeficient, depends on calcite saturation state and is enhanced by
               ! respiration deep in the sediment (s-1), non-linearity controlled by alpha_cased
               cobalt%cased_redis_coef(i,j) = cobalt%gamma_cased*max(0.0,1.0-cobalt%btm_omega_calc(i,j)+ &
                  cobalt%phi_deepresp_cased*cobalt%fntot_btm(i,j)*cobalt%c_2_n*spery)**cobalt%alpha_cased

               ! Effective thickness term that enhances burial of calcite when total sediment accumulation is high
               ! dimensionless value between 0 and 1
               cobalt%cased_redis_delz(i,j) = max(1.0, &
                  cobalt%f_lithdet_btf(i,j,1)*spery+cobalt%f_cadet_calc_btf(i,j,1)*100.0*spery)**cobalt%beta_cased

               ! calculate the sediment redissolution rate (moles Ca m-2 sec-1). This calculation is subject to
               ! three limiters: a) a maximum of 1/2 of the total cased over one time step; b) a maximum of 0.01
               ! moles Ca per day; and c) a minimum of 0.0
               cobalt%fcased_redis(i,j) = max(0.0, min(0.01/sperd, min(0.5*cobalt%f_cased(i,j,1)*r_dt,  &
                  cobalt%fcased_redis_surfresp(i,j)+cobalt%cased_redis_coef(i,j)*cobalt%cased_redis_delz(i,j)*cobalt%f_cased(i,j,1))) )

               !
               ! Old expression
               !
               !cobalt%fcased_redis(i,j) = max(0.0, min(0.01/sperd,min(0.5 * cobalt%f_cased(i,j,1) * r_dt, min(0.5 *       &
               !   cobalt%f_cadet_calc_btf(i,j,1), 0.14307 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n) +        &
               !   0.03607 / spery * max(0.0, 1.0 - cobalt%omega_calc(i,j,k) +   &
               !   4.1228 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n * spery)**(2.7488) *                        &
               !   max(1.0, cobalt%f_lithdet_btf(i,j,1) * spery + cobalt%f_cadet_calc_btf(i,j,1) * 100.0 *  &
               !   spery)**(-2.2185) * cobalt%f_cased(i,j,1))))*grid_tmask(i,j,k)

               if (cobalt%cased_steady) then
                  cobalt%fcased_burial(i,j) = cobalt%f_cadet_calc_btf(i,j,1) - cobalt%fcased_redis(i,j)
                  cobalt%f_cased(i,j,1) = cobalt%fcased_burial(i,j)*cobalt%Co_cased/cobalt%f_cadet_calc_btf(i,j,1)
               else
                  cobalt%fcased_burial(i,j) = max(0.0, cobalt%f_cadet_calc_btf(i,j,1) * cobalt%f_cased(i,j,1) / &
                     cobalt%Co_cased)
                  cobalt%f_cased(i,j,1) = cobalt%f_cased(i,j,1) + (cobalt%f_cadet_calc_btf(i,j,1) -            &
                     cobalt%fcased_redis(i,j) - cobalt%fcased_burial(i,j)) / cobalt%z_sed * dt *                &
                     grid_tmask(i,j,k)
               endif

               !
               ! Bottom flux boundaries passed to the vertical mixing routine
               ! (negative values are fluxes into the ocean)
               !
               cobalt%b_dic(i,j) =  - cobalt%fcased_redis(i,j) - cobalt%f_cadet_arag_btf(i,j,1) -       &
                  (cobalt%fntot_btm(i,j) - cobalt%fn_burial(i,j)) * cobalt%c_2_n
               cobalt%b_fed(i,j) = - cobalt%ffe_sed(i,j) - cobalt%ffe_geotherm(i,j)
               cobalt%b_nh4(i,j) = - cobalt%fntot_btm(i,j) + cobalt%fn_burial(i,j)
               cobalt%b_no3(i,j) = cobalt%fno3denit_sed(i,j)
               ! Include latent O2 demand and alkalinity effects of HS- (see stoichiometry)
               if (cobalt%do_fnso4red_sed) then
                  cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * (cobalt%fnoxic_sed(i,j) + cobalt%fnso4red_sed(i,j))
                  cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
                     cobalt%fnoxic_sed(i,j) - cobalt%fno3denit_sed(i,j)*cobalt%alk_2_n_denit - cobalt%fnso4red_sed(i,j)
               else
                  cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * cobalt%fnoxic_sed(i,j)
                  cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+cobalt%f_cadet_arag_btf(i,j,1)) -    &
                     cobalt%fnoxic_sed(i,j) - cobalt%fno3denit_sed(i,j)*cobalt%alk_2_n_denit
               endif
               cobalt%b_po4(i,j) = - cobalt%fptot_btm(i,j) + cobalt%fp_burial(i,j)
               cobalt%b_sio4(i,j)= - cobalt%fsitot_btm(i,j)

            endif !}
         enddo; enddo  !} i, j

      do k = 2, nk ; do j = jsc, jec ; do i = isc, iec   !{
               cobalt%f_cased(i,j,k) = 0.0
            enddo; enddo ; enddo  !} i,j,k



      ! set the cobalt%b_* terms. These are used in some other places in COBALT as well. So just "b_o2" might not work.
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               !cobalt%b_dic(i,j) = b_dic(i,j)
               !cobalt%b_o2(i,j)  = b_o2(i,j)
               !cobalt%b_nh4(i,j) = b_nh4(i,j)
               !cobalt%b_no3(i,j) = b_no3(i,j)
            endif
         enddo; enddo


      call g_tracer_set_values(cobalt_tracer_list,'alk',  'btf', cobalt%b_alk ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'dic',  'btf', cobalt%b_dic ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'fed',  'btf', cobalt%b_fed ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'nh4',  'btf', cobalt%b_nh4 ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'no3',  'btf', cobalt%b_no3 ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'o2',   'btf', cobalt%b_o2  ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'po4',  'btf', cobalt%b_po4 ,isd,jsd)
      call g_tracer_set_values(cobalt_tracer_list,'sio4', 'btf', cobalt%b_sio4,isd,jsd)

   end subroutine generic_CBED_sediments_update_from_source

end module generic_CBED
