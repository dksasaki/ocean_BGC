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
   use fms2_io_mod, only: FmsNetcdfDomainFile_t, open_file, close_file, read_restart, write_restart
   use fms2_io_mod, only: register_restart_field, register_axis
   use fms_mod, only: error_mesg, NOTE, WARNING, FATAL

   implicit none; private

   character(len=fm_string_len), parameter :: mod_name       = 'generic_CBED'
   character(len=fm_string_len), parameter :: package_name   = 'generic_cbed'

   public generic_CBED_sediments_update_from_source
   public generic_CBED_init, generic_CBED_end
   public generic_CBED_reg_diagnostics, generic_CBED_send_diagnostics

   integer, parameter :: nk_cbed = 20    ! Number of benthic layers

   type generic_CBED_type
      real, dimension(:,:,:), allocatable :: f_tr1  ! tracer 1 concentration field
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
      real, dimension(:,:,:), allocatable :: R_om_o2          ! OM resep via aerobic process
      real, dimension(:,:,:), allocatable :: R_om_no3         ! OM resep via denitrification
      real, dimension(:,:,:), allocatable :: R_om_anaerobic   ! OM resep via other anaerobic process
      real, dimension(:,:,:), allocatable :: R_dic            ! total remineralization

      ! 2D diags
      real, dimension(:,:), allocatable :: o2_flux !benthic o2 flux
      real, dimension(:,:), allocatable :: nh4_flux !benthic nh4 flux
      real, dimension(:,:), allocatable :: no3_flux !benthic no3 flux
      real, dimension(:,:), allocatable :: dic_flux !benthic dic flux
      !real, dimension(:,:), allocatable :: burial_om !organic matter burial at the bottom of sediment column
      !real, dimension(:,:), allocatable :: cbed_denit
      !real, dimension(:,:), allocatable :: cbed_anammox
      !real, dimension(:,:), allocatable :: cbed_o2resp
      !real, dimension(:,:), allocatable :: cbed_no3resp
      integer :: id_tr1                              ! tracer 1 diagnostics id
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
      integer :: id_R_om_o2
      integer :: id_R_om_no3
      integer :: id_R_om_anaerobic
      integer :: id_R_dic
      ! 2D diags
      integer :: id_o2_flux
      integer :: id_nh4_flux
      integer :: id_no3_flux
      integer :: id_dic_flux
      !integer :: id_burial_om
      !integer :: id_cbed_denit
      !integer :: id_cbed_anammox
      !integer :: id_cbed_o2resp
      !integer :: id_cbed_no3resp



   end type generic_CBED_type

   type(generic_CBED_type) :: cbed

   real, parameter :: pi = acos(-1.0)

! porosity. check with Niki
   !real, dimension(isc:iec,jsc:jec) :: por = 0.8  !niki
   !real :: por = 0.8

   ! grid
   ! local parameters
   real, parameter :: l_cbed = 0.20           ! length of sediment domain | sediment depth (m, 20 cm)
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

   subroutine generic_CBED_init(isc,iec,jsc,jec,isd,ied,jsd,jed,nk)
      integer,     intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk
      !Locals
      type(domain2D), pointer :: domain
      type(FmsNetcdfDomainFile_t) :: fileobj ! netCDF file object returned by call to fms2_open_file
      character(len=64)           :: restart_file
      logical                     :: file_open_success ! result returned by call to fms2_open_file

      integer :: k !for grid.

      !Allocate and initialize CBED arrays for tracer concentrations and other workarrays
      allocate(cbed%f_tr1(isd:ied,jsd:jed,nk_cbed));cbed%f_tr1=0.0
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
      allocate(cbed%R_om_o2(isd:ied,jsd:jed,nk_cbed));cbed%R_om_o2=0.0
      allocate(cbed%R_om_no3(isd:ied,jsd:jed,nk_cbed));cbed%R_om_no3=0.0
      allocate(cbed%R_om_anaerobic(isd:ied,jsd:jed,nk_cbed));cbed%R_om_anaerobic=0.0
      allocate(cbed%R_dic(isd:ied,jsd:jed,nk_cbed));cbed%R_dic=0.0
      ! 2D diags
      allocate(cbed%o2_flux(isd:ied,jsd:jed)); cbed%o2_flux=0.0
      allocate(cbed%nh4_flux(isd:ied,jsd:jed)); cbed%nh4_flux=0.0
      allocate(cbed%no3_flux(isd:ied,jsd:jed)); cbed%no3_flux=0.0
      allocate(cbed%dic_flux(isd:ied,jsd:jed)); cbed%dic_flux=0.0
      !allocate(cbed%burial_om(isd:ied,jsd:jed));cbed%burial_om=0.0
      !allocate(cbed%cbed_denit(isd:ied,jsd:jed));cbed%cbed_denit=0.0
      !allocate(cbed%cbed_anammox(isd:ied,jsd:jed));cbed%cbed_anammox=0.0
      !allocate(cbed%cbed_o2resp(isd:ied,jsd:jed));cbed%cbed_o2resp=0.0
      !allocate(cbed%cbed_no3resp(isd:ied,jsd:jed));cbed%cbed_no3resp=0.0

      allocate(por(isc:iec,jsc:jec,nk_cbed+1));        por=0.8  !porosity=0.8 assumed constant for whole seafloor.
      allocate(svf(isc:iec,jsc:jec,nk_cbed+1));        svf=0.2 !solid volume fraction

      allocate(w(isc:iec,jsc:jec,nk_cbed+1));        w=0.0      !adding sedimentation rate initalize
      allocate(Db_0(isc:iec,jsc:jec));               Db_0=0.0   !bioturbation_0 init.
      allocate(Db(isc:iec,jsc:jec,nk_cbed+1));       Db=0.0     !bioturbation init.
      allocate(bioirri_0(isc:iec,jsc:jec));          bioirri_0=0.0   !bioirrigation_0 init.
      allocate(bioirri(isc:iec,jsc:jec,nk_cbed));    bioirri=0.0     !bioturbation init.

      allocate(D_o2(isc:iec,jsc:jec,nk_cbed+1)); D_o2=0.0     ! D_o2 init.
      allocate(D_dic(isc:iec,jsc:jec,nk_cbed+1)); D_dic=0.0
      allocate(D_nh4(isc:iec,jsc:jec,nk_cbed+1)); D_nh4=0.0
      allocate(D_no3(isc:iec,jsc:jec,nk_cbed+1)); D_no3=0.0
      allocate(D_odu(isc:iec,jsc:jec,nk_cbed+1)); D_odu=0.0

      allocate(k1(isc:iec,jsc:jec)); k1=0.0
      allocate(k2(isc:iec,jsc:jec)); k2=0.0
      allocate(k3(isc:iec,jsc:jec)); k3=0.0


      ! Grid does not change with time, so can be define only once.
      ! define uniform sediment grid
      dz_cbed = l_cbed / real(nk_cbed)
      z_cbed_int(1) = 0.0   !this is likely the interface. dimention of z_cbed is nk_cbed+1. z_int_cbed. might need z_mid_cbed
      do k = 1, nk_cbed
         z_cbed_int(k+1) = z_cbed_int(k) + dz_cbed(k)
      end do

      z_cbed_mid(1) = dz_cbed(1)/2   ! first layer mid point
      do k = 1, nk_cbed-1
         z_cbed_mid(k+1) = z_cbed_mid(k) + dz_cbed(k)
      end do


   end subroutine generic_CBED_init

   subroutine generic_CBED_reg_diagnostics(axes,init_time)
      USE diag_manager_mod, ONLY: register_diag_field, diag_axis_init
      integer,         intent(in) :: axes(3)
      type(time_type), intent(in) :: init_time
      !Locals
      integer :: k, id_layer
      real :: cbed_layers(1:nk_cbed)

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
         call register_axis(fileobj,'lev',nk_cbed)
         ! register the restart variables
         call register_restart_field(fileobj, "cbed_tr1", cbed%f_tr1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_o2", cbed%f_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om1", cbed%f_om1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om2", cbed%f_om2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om3", cbed%f_om3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_nh4", cbed%f_nh4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_no3", cbed%f_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_dic", cbed%f_dic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_odu", cbed%f_odu, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_talk", cbed%f_talk, (/"x","y","lev"/))
         !diags
         ! 3D diags
         call register_restart_field(fileobj, "cbed_R_om_o2", cbed%R_om_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_R_om_no3", cbed%R_om_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_R_om_anaerobic", cbed%R_om_anaerobic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_R_dic", cbed%R_dic, (/"x","y","lev"/))
         ! 2D diags
         call register_restart_field(fileobj, "cbed_o2_flux", cbed%o2_flux, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_nh4_flux", cbed%nh4_flux, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_no3_flux", cbed%no3_flux, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_dic_flux", cbed%dic_flux, (/"x","y"/))
         !call register_restart_field(fileobj, "cbed_burial_om", cbed%burial_om, (/"x","y"/))
         !call register_restart_field(fileobj, "cbed_denit", cbed%cbed_denit, (/"x","y"/))
         !call register_restart_field(fileobj, "cbed_anammox", cbed%cbed_anammox, (/"x","y"/))
         !call register_restart_field(fileobj, "cbed_o2resp", cbed%cbed_o2resp, (/"x","y"/))
         !call register_restart_field(fileobj, "cbed_no3resp", cbed%cbed_no3resp, (/"x","y"/))

         call read_restart(fileobj)
      endif
      !!END read_restart code block

      !!Register diagnostics
      !Niki: I am unsure of the diag axis thingy, ask Yi-Cheng
      !Define cbed layer axis, the x,y axes are the same as MOM6 since the horizontal grids are the same
      do k=1,nk_cbed; cbed_layers(k) = k; enddo
      id_layer = diag_axis_init('cbedlayer', cbed_layers, 'None', 'z', long_name='Benthos Layer', direction=-1)

      cbed%id_tr1 = register_diag_field(package_name, 'cbed_tr1_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed tracer1 concentration', 'unknown units', missing_value = missing_value1)
      cbed%id_o2 = register_diag_field(package_name, 'cbed_o2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed oxygen concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_om1 = register_diag_field(package_name, 'cbed_om1_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM1 concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_om2 = register_diag_field(package_name, 'cbed_om2_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM2 concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_om3 = register_diag_field(package_name, 'cbed_om3_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed OM3 concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_nh4 = register_diag_field(package_name, 'cbed_nh4_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed ammonium concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_no3 = register_diag_field(package_name, 'cbed_no3_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed nitrate concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_dic = register_diag_field(package_name, 'cbed_dic_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed DIC concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_odu = register_diag_field(package_name, 'cbed_odu_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed ODU concentration', 'mol/kg', missing_value = missing_value1)
      cbed%id_talk = register_diag_field(package_name, 'cbed_talk_conc', (/axes(1),axes(2),id_layer/), init_time,&
         'cbed total alkalinity concentration', 'mol/kg', missing_value = missing_value1)
      ! diags
      ! 3D diags
      cbed%id_R_om_o2 = register_diag_field(package_name, 'cbed_R_om_o2', (/axes(1),axes(2),id_layer/), init_time,&
         'aerobic respiration in sediment 3D field', 'mol C/kg', missing_value = missing_value1)
      cbed%id_R_om_no3 = register_diag_field(package_name, 'cbed_R_om_no3', (/axes(1),axes(2),id_layer/), init_time,&
         'OM respiration via denitrification in sediment 3D field', 'mol C/kg', missing_value = missing_value1)
      cbed%id_R_om_anaerobic = register_diag_field(package_name, 'cbed_R_om_anaerobic', (/axes(1),axes(2),id_layer/), init_time,&
         'OM respiration via other anaerobic processes in sediment 3D field', 'mol C/kg', missing_value = missing_value1)
      cbed%id_R_dic = register_diag_field(package_name, 'cbed_R_dic', (/axes(1),axes(2),id_layer/), init_time,&
         'DIC produced in sediment via OM remineralization 3D field', 'mol C/kg', missing_value = missing_value1)

      ! 2D diags
      cbed%id_o2_flux = register_diag_field(package_name, 'cbed_o2_flux', (/axes(1),axes(2)/), init_time,&
         'benthic O2 flux', 'mol/m2/s', missing_value = missing_value1)
      cbed%id_nh4_flux = register_diag_field(package_name, 'cbed_nh4_flux', (/axes(1),axes(2)/), init_time,&
         'benthic nh4 flux', 'mol/m2/s', missing_value = missing_value1)
      cbed%id_no3_flux = register_diag_field(package_name, 'cbed_no3_flux', (/axes(1),axes(2)/), init_time,&
         'benthic no3 flux', 'mol/m2/s', missing_value = missing_value1)
      cbed%id_dic_flux = register_diag_field(package_name, 'cbed_dic_flux', (/axes(1),axes(2)/), init_time,&
         'benthic dic flux', 'mol/m2/s', missing_value = missing_value1)

      !cbed%id_burial_om = register_diag_field(package_name, 'cbed_burial_om', (/axes(1),axes(2)/), init_time,&
      !   'cbed organic carbon burial', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_denit = register_diag_field(package_name, 'cbed_denit', (/axes(1),axes(2)/), init_time,&
      !   'cbed denitrification', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_anammox = register_diag_field(package_name, 'cbed_anammox', (/axes(1),axes(2)/), init_time,&
      !   'cbed anammox', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_o2resp = register_diag_field(package_name, 'cbed_o2resp', (/axes(1),axes(2)/), init_time,&
      !   'cbed OM respiration by O2', 'mol/m2/s', missing_value = missing_value1)
      !cbed%id_cbed_no3resp = register_diag_field(package_name, 'cbed_no3resp', (/axes(1),axes(2)/), init_time,&
      !   'cbed OM respiration by NO3', 'mol/m2/s', missing_value = missing_value1)

   end subroutine generic_CBED_reg_diagnostics

   subroutine generic_CBED_send_diagnostics(model_time,grid_tmask, isc,iec,jsc,jec, isd,ied,jsd,jed,nk)
      USE diag_manager_mod, ONLY: send_data
      type(time_type),          intent(in) :: model_time
      real, dimension(:,:,:),    pointer   :: grid_tmask
      integer,                  intent(in) :: isc,iec,jsc,jec, isd,ied,jsd,jed,nk
      ! local
      logical :: used
      integer :: k
      real,dimension(isd:ied,jsd:jed,nk_cbed)    :: cbed_tmask
      !Make a cbed mask. Note: it seems grid_tmask(:,:,k) does not depend on k
      do k=1,nk_cbed ; cbed_tmask(:,:,k) = grid_tmask(:,:,nk) ; enddo

      used = send_data(cbed%id_tr1, cbed%f_tr1, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
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
      used = send_data(cbed%id_R_om_o2, cbed%R_om_o2, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_no3, cbed%R_om_no3, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_om_anaerobic, cbed%R_om_anaerobic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      used = send_data(cbed%id_R_dic, cbed%R_dic, model_time, rmask = cbed_tmask,&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec, ks_in=1, ke_in=nk_cbed)
      ! 2D diags
      used = send_data(cbed%id_o2_flux, cbed%o2_flux, model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_nh4_flux, cbed%nh4_flux, model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_no3_flux, cbed%no3_flux, model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
      used = send_data(cbed%id_dic_flux, cbed%dic_flux, model_time, rmask = grid_tmask(:,:,1),&
         is_in=isc, js_in=jsc,ie_in=iec, je_in=jec)
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
         call register_axis(fileobj,'lev',nk_cbed)
         ! register the restart variables
         call register_restart_field(fileobj, "cbed_tr1", cbed%f_tr1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_o2", cbed%f_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om1", cbed%f_om1, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om2", cbed%f_om2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_om3", cbed%f_om3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_nh4", cbed%f_nh4, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_no3", cbed%f_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_dic", cbed%f_dic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_odu", cbed%f_odu, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_talk", cbed%f_talk, (/"x","y","lev"/))
         ! diags
         ! 3D diags
         call register_restart_field(fileobj, "cbed_R_om_o2", cbed%R_om_o2, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_R_om_no3", cbed%R_om_no3, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_R_om_anaerobic", cbed%R_om_anaerobic, (/"x","y","lev"/))
         call register_restart_field(fileobj, "cbed_R_dic", cbed%R_dic, (/"x","y","lev"/))
         ! 2D diags
         call register_restart_field(fileobj, "cbed_o2_flux", cbed%o2_flux, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_nh4_flux", cbed%nh4_flux, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_no3_flux", cbed%no3_flux, (/"x","y"/))
         call register_restart_field(fileobj, "cbed_dic_flux", cbed%dic_flux, (/"x","y"/))

         call write_restart(fileobj)
         call close_file(fileobj)
      else
         call error_mesg( 'generic_CBED_end', 'Cannot open restarts for write.', FATAL )
      endif

      !Deallocate arrays
      deallocate(cbed%f_tr1)
      deallocate(cbed%f_o2)
      deallocate(cbed%f_om1)
      deallocate(cbed%f_om2)
      deallocate(cbed%f_om3)
      deallocate(cbed%f_nh4)
      deallocate(cbed%f_no3)
      deallocate(cbed%f_dic)
      deallocate(cbed%f_odu)
      deallocate(cbed%f_talk)

      deallocate(cbed%R_om_o2)
      deallocate(cbed%R_om_no3)
      deallocate(cbed%R_om_anaerobic)
      deallocate(cbed%R_dic)
      deallocate(cbed%o2_flux)
      deallocate(cbed%nh4_flux)
      deallocate(cbed%no3_flux)
      deallocate(cbed%dic_flux)

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


   end subroutine generic_CBED_end



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
      real, dimension(isd:ied,jsd:jed,nk_cbed) :: ea,eb
      real, dimension(isc:iec,jsc:jec,nk_cbed+1) :: sink

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
                  ea(i,j,k) = D(i,j,k)*dt/h_old(k)
                  eb(i,j,k) = D(i,j,k+1)*dt/h_old(k)
               enddo
            endif
         enddo; enddo

      ! sink(k+1)
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               do k=1,nk_cbed+1
                  sink(i,j,k) = max(0.0, w(i,j,k)*dt )
               enddo
            endif
         enddo; enddo

      ! a, b, c, f_old
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then

               !! NEED TO TAKE CARE OF THE TOP AND BOTTOM FLUXES
               ! sfc_src = 0.0 ; btm_src = 0.0
               ! if (_ALLOCATED(g_tracer%stf)) sfc_src = (g_tracer%stf(i,j)*dt)*kg_m2_to_H
               ! if (_ALLOCATED(g_tracer%btf)) btm_src = (-g_tracer%btf(i,j)*dt)*kg_m2_to_H
               ! g_tracer%field(i,j,1,tau)  = g_tracer%field(i,j,1,tau)  + sfc_src/h_old(i,j,1)
               ! g_tracer%field(i,j,nz,tau) = g_tracer%field(i,j,nz,tau) + btm_src/h_old(i,j,nz)

               sfc_src = 0.0

               if ((trim(field_name) == "f_o2")) then
                  sfc_src = VF(i,j,1)*D(i,j,1)*((cobalt%btm_o2(i,j)-cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_nh4")) then
                  sfc_src =  VF(i,j,1)*D(i,j,1)*((cobalt%f_nh4(i,j,nk)-cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_no3")) then
                  sfc_src = VF(i,j,1)*D(i,j,1)*((cobalt%btm_no3(i,j)-cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_dic") ) then
                  sfc_src = VF(i,j,1)*D(i,j,1)*((cobalt%btm_dic(i,j)-cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_odu") ) then
                  sfc_src = VF(i,j,1)*D(i,j,1)*((0.0-cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_talk") ) then
                  sfc_src = VF(i,j,1)*D(i,j,1)*((cobalt%btm_alk(i,j)-cbed_field(i,j,1))/(dz_cbed(1)/2.0))*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_om1")) then
                  sfc_src = frac_OM1*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_om2")) then
                  sfc_src = frac_OM2*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               else if ((trim(field_name) == "f_om3")) then
                  sfc_src = frac_OM3*cobalt%fntot_btm(i,j)*cobalt%c_2_n*dt ! top flux
                  cbed_field(i,j,1)  = cbed_field(i,j,1)  + sfc_src/h_old(1)

               endif

               ! bottom flux
               btm_src = 0.0
               btm_src = -(cbed_field(i,j,nk_cbed)*VF(i,j,nk_cbed)*w(i,j,nk_cbed+1)*dt) ! bottom flux
               cbed_field(i,j,nk_cbed) = cbed_field(i,j,nk_cbed) + btm_src/h_old(nk_cbed)



               do k=1,nk_cbed

                  !a(k)= -(ea(i,j,k)+sink(i,j,k))/(VF(i,j,k)*h_old(k))

                  !b(k)=  (VF(i,j,k)*h_old(k)+eb(i,j,k)+ea(i,j,k)+sink(i,j,k+1))/(VF(i,j,k)*h_old(k))

                  !c(k)= -eb(i,j,k)/(VF(i,j,k)*h_old(k))

                  a(k)= (-ea(i,j,k)-sink(i,j,k))/(h_old(k))

                  b(k)= (h_old(k)+eb(i,j,k)+ea(i,j,k)+sink(i,j,k+1))/(h_old(k))

                  c(k)= -eb(i,j,k)/(h_old(k))

                  f_old(k)= cbed_field(i,j,k)

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
      grid_tmask, grid_dat, grid_kmt, isc,iec, jsc,jec, isd,ied, jsd,jed, nk, r_dt, dt, tau, frunoff, rho_dzt, dzt, internal_heat)

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

      real, parameter :: ks_o2 = 0.008 /1e3   ! O2 half saturation constant (mol/kg)
      real, parameter :: ks_no3 = 0.001 /1e3  ! NO3 half saturation constant (mol/kg)

      real, parameter :: k_nox = 1e6 /1e3/spery    ! mol/kg/s (mol/L/s) !nitrification rate constant
      real, parameter :: k_ana = 1e5 /1e3/spery    !                    !anammox rate constant
      real, parameter :: k_oduox = 1e6 /1e3/spery  !                    !ODU oxidation rate constant

      real, parameter :: Q10 = 1.88

      ! Reaction rates
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_o2, R_om2_o2, R_om3_o2
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_no3, R_om2_no3, R_om3_no3
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_om1_anoxic, R_om2_anoxic, R_om3_anoxic
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_dic_om1, R_dic_om2, R_dic_om3
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_nox, R_ana, R_oduox, odu_depo
      real, dimension(isc:iec,jsc:jec,nk_cbed) :: R_talk

      ! b terms
      real, dimension(isc:iec,jsc:jec) :: b_o2, b_dic, b_nh4, b_no3


      ! Sedimentation rate calculation
      do j = jsc, jec; do i = isc, iec
            do k = 1, nk_cbed+1
               if (grid_kmt(i,j) .gt. 0) then
                  ! w (sedimentation rate, cm/year) m/s
                  w(i,j,k) = ( (cobalt%fcadet_arag_btm(i,j)*100.0/2.71 + &
                     cobalt%fcadet_calc_btm(i,j)*100.0/2.94 + &
                     cobalt%fsitot_btm(i,j)*60.0/2.65 + &
                     cobalt%flithdet_btm(i,j)/2.65 + &
                     cobalt%ffetot_btm(i,j)*160.0/5.24 + &
                     cobalt%fptot_btm(i,j)*120.0/2.3 + &
                     cobalt%fntot_btm(i,j)*cobalt%c_2_n*22.4/0.9)/10000.0*spery/(1-por(i,j,k)) )/100.0/spery
               endif
            enddo
         enddo;enddo

      !Bioturbation
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               ! relation from Archer. POC flux unit in umol cm-2 y-1.
               Db_0(i,j) = ( 0.0232*((cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery)**0.85) ) /1e4/spery ! in cobalt unit m2/s
            endif
         enddo;enddo

      do j = jsc, jec; do i = isc, iec
            do k = 1, nk_cbed+1
               if (grid_kmt(i,j) .gt. 0) then
                  ! relation from Archer. POC flux unit in umol cm-2 y-1.
                  Db(i,j,k) = max(0.0, Db_0(i,j)*exp(-(z_cbed_int(k)/Db_l)**2)*(cobalt%btm_o2(i,j)/(cobalt%btm_o2(i,j)+(20/1e6))) )
               endif
            enddo
         enddo;enddo

      !Bioirrigation
      do j = jsc, jec; do i = isc, iec
            if (grid_kmt(i,j) .gt. 0) then
               ! relation from Archer. POC flux unit in umol cm-2 y-1.
               bioirri_0(i,j) = ( 11*(((atan((5*(cobalt%fntot_btm(i,j)*cobalt%c_2_n *1e6/1e4*spery) -400)/400))/pi)+0.5) &
                  - 0.9 + 20*((cobalt%btm_o2(i,j)*1e6)/(cobalt%btm_o2(i,j)*1e6+10)) * exp(-cobalt%btm_o2(i,j)*1e6/10) * &
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
            do k = 1, nk_cbed+1
               if (grid_kmt(i,j) .gt. 0) then
                  D_o2(i,j,k)  = ( (0.031558+0.001428*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)   ! m2/s
                  D_dic(i,j,k) = ( (0.015179+0.000795*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_nh4(i,j,k) = ( (0.030926+0.001225*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_no3(i,j,k) = ( (0.030863+0.001153*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
                  D_odu(i,j,k) = ( (0.028938+0.001314*cobalt%btm_temp(i,j))/(1-2*log(por(i,j,k))) )/spery + Db(i,j,k)
               endif
            enddo
         enddo;enddo


      ! calculate k1,k2,k3

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
            do k = 1, nk_cbed
               if (grid_kmt(i,j) .gt. 0) then

                  ! O₂ reaction rates
                  R_om1_o2(i,j,k) = k1(i,j)*cbed%f_om1(i,j,k)*(cbed%f_o2(i,j,k)/(ks_o2 + cbed%f_o2(i,j,k)))
                  R_om2_o2(i,j,k) = k2(i,j)*cbed%f_om2(i,j,k)*(cbed%f_o2(i,j,k)/(ks_o2 + cbed%f_o2(i,j,k)))
                  R_om3_o2(i,j,k) = k3(i,j)*cbed%f_om3(i,j,k)*(cbed%f_o2(i,j,k)/(ks_o2 + cbed%f_o2(i,j,k)))
                  ! NO₃ reaction rates
                  R_om1_no3(i,j,k) = k_adj_denit*k1(i,j)*cbed%f_om1(i,j,k)*(cbed%f_no3(i,j,k)/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                  R_om2_no3(i,j,k) = k_adj_denit*k2(i,j)*cbed%f_om2(i,j,k)*(cbed%f_no3(i,j,k)/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                  R_om3_no3(i,j,k) = k_adj_denit*k3(i,j)*cbed%f_om3(i,j,k)*(cbed%f_no3(i,j,k)/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                  ! ODU reaction rates
                  R_om1_anoxic(i,j,k) = k_adj_anoxia*k1(i,j)*cbed%f_om1(i,j,k)*(ks_no3/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                  R_om2_anoxic(i,j,k) = k_adj_anoxia*k2(i,j)*cbed%f_om2(i,j,k)*(ks_no3/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))
                  R_om3_anoxic(i,j,k) = k_adj_anoxia*k3(i,j)*cbed%f_om3(i,j,k)*(ks_no3/(ks_no3 + cbed%f_no3(i,j,k)))*(ks_o2/(ks_o2 + cbed%f_o2(i,j,k)))

                  ! dic
                  R_dic_om1(i,j,k) = (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k))
                  R_dic_om2(i,j,k) = (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k))
                  R_dic_om3(i,j,k) = (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k))
                  ! nitrification
                  R_nox(i,j,k) = k_nox*cbed%f_nh4(i,j,k)*cbed%f_o2(i,j,k)
                  ! anammox
                  R_ana(i,j,k) = k_ana*cbed%f_nh4(i,j,k)*cbed%f_no3(i,j,k)
                  ! ODU oxidation
                  R_oduox(i,j,k) = k_oduox*cbed%f_odu(i,j,k)*cbed%f_o2(i,j,k)
                  odu_depo(i,j,k) = (R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))*min(1, 0.233*(w(i,j,k)*100.0*spery)**0.336)

                  ! TA calculation
                  R_talk(i,j,k) = (1.0/cobalt%c_2_n)*(R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k))/por(i,j,k) + &
                     (0.8+1.0/cobalt%c_2_n)*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k))/por(i,j,k) + &
                     (1.0+1.0/cobalt%c_2_n)*(R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))/por(i,j,k) - &
                     2.0*R_nox(i,j,k) - 2.0*R_oduox(i,j,k)

                  ! calculations for diagnostics
                  cbed%R_om_o2(i,j,k) = R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k)
                  cbed%R_om_no3(i,j,k) = R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k)
                  cbed%R_om_anaerobic(i,j,k) = R_om1_anoxic(i,j,k) + R_om2_anoxic(i,j,k) + R_om3_anoxic(i,j,k)
                  cbed%R_dic(i,j,k) = R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k)


               endif
            enddo
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
               b_o2(i,j) = por(i,j,1)*D_o2(i,j,1)*((cobalt%btm_o2(i,j)-cbed%f_o2(i,j,1))/(dz_cbed(1)/2.0))
               b_nh4(i,j) = por(i,j,1)*D_nh4(i,j,1)*((cobalt%f_nh4(i,j,nk)-cbed%f_nh4(i,j,1))/(dz_cbed(1)/2.0))
               b_no3(i,j) = por(i,j,1)*D_no3(i,j,1)*((cobalt%btm_no3(i,j)-cbed%f_no3(i,j,1))/(dz_cbed(1)/2.0))
               b_dic(i,j) = por(i,j,1)*D_dic(i,j,1)*((cobalt%btm_dic(i,j)-cbed%f_dic(i,j,1))/(dz_cbed(1)/2.0))
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


      !Test that we can change the value of concentration field of a CBED tracer
      do j = jsc, jec; do i = isc, iec  !{
            do k=1,nk_cbed
               if (grid_kmt(i,j) .gt. 0) then
                  cbed%f_tr1(i,j,k) = cbed%f_tr1(i,j,k) + 0.01 * k !fictitious dubious dynamics for testing purposes

                  cbed%f_o2(i,j,k)  = cbed%f_o2(i,j,k) - (R_om1_o2(i,j,k) + R_om2_o2(i,j,k) + R_om3_o2(i,j,k))/por(i,j,k) - &
                     (2.0*R_nox(i,j,k)+R_oduox(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_o2(i,j) - cbed%f_o2(i,j,k))

                  cbed%f_om1(i,j,k) = cbed%f_om1(i,j,k) - (R_om1_o2(i,j,k) + R_om1_no3(i,j,k) + R_om1_anoxic(i,j,k))

                  cbed%f_om2(i,j,k) = cbed%f_om2(i,j,k) - (R_om2_o2(i,j,k) + R_om2_no3(i,j,k) + R_om2_anoxic(i,j,k))

                  cbed%f_om3(i,j,k) = cbed%f_om3(i,j,k) - (R_om3_o2(i,j,k) + R_om3_no3(i,j,k) + R_om3_anoxic(i,j,k))

                  cbed%f_nh4(i,j,k) = cbed%f_nh4(i,j,k) + (1.0/cobalt%c_2_n)*(R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k))/por(i,j,k) + &
                     ( - R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%f_nh4(i,j,nk) - cbed%f_nh4(i,j,k))

                  cbed%f_no3(i,j,k) = cbed%f_no3(i,j,k) - 0.8*(R_om1_no3(i,j,k) + R_om2_no3(i,j,k) + R_om3_no3(i,j,k))/por(i,j,k) + &
                     (R_nox(i,j,k) - R_ana(i,j,k)) + bioirri(i,j,k)*(cobalt%btm_no3(i,j) - cbed%f_no3(i,j,k))

                  cbed%f_dic(i,j,k) = cbed%f_dic(i,j,k) + (R_dic_om1(i,j,k) + R_dic_om2(i,j,k) + R_dic_om3(i,j,k))/por(i,j,k) + &
                     bioirri(i,j,k)*(cobalt%btm_dic(i,j) - cbed%f_dic(i,j,k))

                  cbed%f_odu(i,j,k) = cbed%f_odu(i,j,k) + (R_om1_anoxic(i,j,k)+R_om2_anoxic(i,j,k)+R_om3_anoxic(i,j,k))/por(i,j,k) - &
                     R_oduox(i,j,k) - odu_depo(i,j,k)/por(i,j,k)  + bioirri(i,j,k)*(0.0 - cbed%f_odu(i,j,k))

                  cbed%f_talk(i,j,k) = cbed%f_talk(i,j,k) + R_talk(i,j,k) + bioirri(i,j,k)*(cobalt%btm_alk(i,j) - cbed%f_talk(i,j,k))


               endif

            enddo
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
               cobalt%b_dic(i,j) = b_dic(i,j)
               cobalt%b_o2(i,j)  = b_o2(i,j)
               cobalt%b_nh4(i,j) = b_nh4(i,j)
               cobalt%b_no3(i,j) = b_no3(i,j)
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
