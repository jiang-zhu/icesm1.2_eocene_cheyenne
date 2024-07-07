!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module estuary_mod

! !DESCRIPTION:
!  This module computes estuary parameterization
!
!   estuary_type 
!     = 'top-layer'   :  all sw absorbed in top ocean layer
!     = 'equ-depth'   :  equivalent depth
!     = 'sca-weigh'   :  weighted function
! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use kinds_mod
   use domain_size
   use domain
   use constants
   use io
   use io_types
   use grid
   use forcing_fields,only:ROFF_F
   use time_management
   use prognostic,only:TRACER,curtime   
!TS
   use global_reductions, only: global_sum_prod
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_estuary, &
             add_estuary_param, &
!TS
             set_estuary_vsf_forcing,&
             set_ep

! !PUBLIC DATA MEMBERS:
   real (r8), dimension(1:km) :: ep_frac
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),public ::& 
      ep_ratio                     

   character (char_len), public ::       &
      estuary_type          ! type of estuary_parameterization
   integer(int_kind), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      k_h_c                         ! height     mg/m^3
   integer (int_kind), dimension(nx_block,ny_block,max_blocks_clinic,12) :: &
      k_h_c_DATA    ! k index of global h_c data

   character (char_len) ::       &
      h_c_option,                &! h_c option ['file', 'model']
      h_c_filename,              &! h_c data file name
      h_c_file_fmt,              &! h_c data file format
      h_c_data_name               ! h_c short name (eg, 'H_C')

   integer (int_kind) ::         &
      h_c_bndy_loc,              &! location and field types for ghost
      h_c_bndy_type               !    cell update routines
!TS
   real(r8), public :: vsf_river_correction
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::&
        FLUX_ROFF_VSF_SRF

!-----------------------------------------------------------------------
!
!  tavg ids for tavg diagnostics related to EBM
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
         tavg_FLUX_ROFF_VSF_SRF !Surf virt salt flux from river

!EOC
!***********************************************************************

   contains

!***********************************************************************

   subroutine init_estuary
   
!-----------------------------------------------------------------------
!
!     initialize estuary parameterization
!
!-----------------------------------------------------------------------0

   implicit none

   integer (int_kind) :: &
      i,j,k,n,           &! level index
      nu,                &! unit for input dataset
      nml_error,         &! namelist error flag
      bid,               &! local block address for this block
      iblock
   real (r8), allocatable, dimension(:,:,:) :: &
      TEMP_DATA   ! temporary data array
   real (r8), allocatable, dimension(:,:,:,:) :: &
      h_c_DATA    ! global h_c data 
   type (datafile) ::        &
      forcing_file            !data file structure for input forcing file

   type (io_field_desc) ::   &
      io_h_c                  ! io field descriptor for input sss field

   type (io_dim) ::          &
      i_dim, j_dim,          &! dimension descriptors for horiz dimensions
      month_dim               ! dimension descriptor  for monthly data

   namelist /estuary_nml/ &
        estuary_type,     &! 'top-layer', 'equ-depth', 'sca_weigh'
        h_c_option,       &! h_c option ['file', 'model']
        h_c_filename,     &! include local filepath in name
        h_c_file_fmt       ! include local filepath in name

   estuary_type   = 'top-layer'
   h_c_option     = 'model'
!!TS   h_c_option     = 'file'
   h_c_file_fmt = 'nc'
   h_c_filename = 'runoff_hc.nc'

   IF(.FALSE.)THEN
   if (my_task == master_task) then
   open (nml_in, file=nml_filename, status='old',iostat=nml_error)
   if (nml_error /= 0) then
      nml_error = -1
   else
      nml_error =  1
   endif
   do while (nml_error > 0)
      read(nml_in, nml=estuary_nml,iostat=nml_error)
   end do
   if (nml_error == 0) close(nml_in)
   endif
 
   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading estuary_nml')
   endif

   if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' Estuary parameterization:'
       write(stdout,blank_fmt)
       write(stdout,*) ' Estuary parametierzation namelist  '
       write(stdout,blank_fmt)
       write(stdout, estuary_nml)
       write(stdout,blank_fmt)
   endif

   call broadcast_scalar(estuary_type,    master_task)
   call broadcast_scalar(h_c_option,            master_task)
   call broadcast_scalar(h_c_filename,          master_task)
   ENDIF!!!SET FALSE

   IF(my_task.EQ.0)WRITE(stdout,*)'Hahaha-estuary_type',estuary_type
   if (estuary_type .ne. 'top-layer'.AND.&
       estuary_type .ne. 'equ-depth'.AND.&
       estuary_type .ne. 'sca-weigh' ) then
     call exit_POP(sigAbort,'ERROR estuary_type unknown')
   endif
   allocate( h_c_DATA (nx_block,ny_block,max_blocks_clinic,12))
   allocate( TEMP_DATA(nx_block,ny_block,max_blocks_clinic))

   if (h_c_option == 'file') then

    h_c_bndy_loc  = field_loc_center
    h_c_bndy_type = field_type_scalar
    h_c_data_name = 'H_C'

    h_c_DATA = c0;TEMP_DATA = c0

    forcing_file = construct_file(h_c_file_fmt,                 &
                                   full_name=trim(h_c_filename), &
                                   record_length = rec_type_dbl,  &
                                   recl_words=nx_global*ny_global)

    call data_set(forcing_file,'open_read')

    i_dim     = construct_io_dim('nx',nx_global)
    j_dim     = construct_io_dim('ny',ny_global)
    io_h_c = construct_io_field( &
                    trim(h_c_data_name),                        &
                    dim1=i_dim, dim2=j_dim,                     &
                    field_loc  = h_c_bndy_loc,                  &
                    field_type = h_c_bndy_type,                 &
                    d2d_array=TEMP_DATA(:,:,:))
    call data_set(forcing_file,'define',io_h_c)
    call data_set(forcing_file,'read'  ,io_h_c)
    IF(mod(MY_TASK,4).EQ.0)WRITE(*,*)'Haha-ini ',MY_TASK,TEMP_DATA(10,10,1)
    call destroy_io_field(io_h_c)

      !*** re-order data

      !$OMP PARALLEL DO PRIVATE(iblock, n)
!TEMP_DATA unit="m", need to convert it into "cm"
    do iblock=1,nblocks_clinic
         do n=1,12
            h_c_DATA(:,:,iblock,n) = TEMP_DATA(:,:,iblock)*100.d0
         do j=1,ny_block
         do i=1,nx_block
            do k=1,KMT(i,j,iblock) 
             k_h_c_data(:,:,iblock,n)=k
             IF(h_c_DATA(i,j,iblock,n).GE. zw(k))exit
            end do
         enddo
         enddo
         end do
      end do
      !$OMP END PARALLEL DO
      call data_set(forcing_file,'close')
      call destroy_file(forcing_file)
      if (my_task.eq.master_task) then
       write(stdout,blank_fmt)
       write(stdout,*) ' Estuary Parameterization h_c monthly file read: ',h_c_filename
      endif
    elseif(h_c_option == 'model') then
      !$OMP PARALLEL DO PRIVATE(iblock, n)
    do iblock=1,nblocks_clinic
         do n=1,12
            h_c_DATA(:,:,iblock,n) = 10.d0*100.d0
         end do
      end do
      !$OMP END PARALLEL DO
    else
      call exit_POP(sigAbort,'ERROR h_c_option unknown')
    endif
    k_h_c_data=0
!      !$OMP PARALLEL DO PRIVATE(iblock, n)
    do n=1,12
      do iblock=1,nblocks_clinic
         do j=1,ny_block
         do i=1,nx_block         
           do k=1,KMT(i,j,iblock)
             k_h_c_data(i,j,iblock,n)=k
             if(h_c_DATA(i,j,iblock,n).LE. zw(k))exit
           end do
         end do
         end do
      end do
    end do
!      !$OMP END PARALLEL DO
    deallocate(TEMP_DATA)
    deallocate(h_c_DATA)
    k_h_c(:,:,:) = k_h_c_DATA(:,:,:,imonth) ! constant across month  
    IF(mod(MY_TASK,4).EQ.0)WRITE(*,*)'Haha-h_h_c',my_task,KMT(10,10,1),k_h_c_DATA(10,10,1,imonth)
!-----------------------------------------------------------------------
!
!  define tavg fields related to EBM
!
!-----------------------------------------------------------------------
   call define_tavg_field(tavg_FLUX_ROFF_VSF_SRF,'FLUX_ROFF_VSF_SRF',2,                     &
                    long_name='Surface Virtual Salt Flux Associated with Rivers',&
                          units='g/kg*cm/s', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')

   end subroutine init_estuary


!***********************************************************************
!BOP
! !IROUTINE: add_estuary_param
! !INTERFACE:

 subroutine add_estuary_param(T_SOURCE, k, this_block)

! !DESCRIPTION:
!  If surface short wave heat flux is available, this routine caculates
!  the flux which passes through the top layer and enters lower vertical
!  depths in the ocean.  This flux is added as a source term in the
!  baroclinic equations.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      T_SOURCE     ! source terms for all tracers (to avoid copies)
                   ! sw absorption added only to other potential
                   ! temperature tracers

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k            ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block info for this block
   integer(int_kind)::i,j
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid                ! local block index
 
   real (r8), dimension(nx_block,ny_block) :: &
      WORK,WORK1    ! temporary work space

   real (r8)::tlatd,tlond
      bid = this_block%local_id
!!TS
      ep_ratio=0.d0
      where(K.LE.k_h_c(:,:,bid))
      ep_ratio(:,:,bid)=1.d0/DBLE(k_h_c(:,:,bid))
      endwhere
!-----------------------------------------------------------------------
!EOC
!TS      WORK = RCALCT(:,:,bid)*ROFF_F(:,:,bid)*salinity_factor
!TS use the local surface salinity because the bottom salinity may be unrealistically large 
        WORK=RCALCT(:,:,bid)*ROFF_F(:,:,bid)*&
        (-MAX(TRACER(:,:,1,2,curtime,bid)*c1000*fwflux_factor,0._r8))

      WORK1= WORK*dzr(k)*ep_ratio(:,:,bid)
      do j=1,ny_block
      do i=1,nx_block
         tlatd = TLAT(i,j,bid)*radian
         tlond = TLON(i,j,bid)*radian
      if ( tlatd.gt.45.5d0.and.tlatd.lt.47.5d0.and.&
           tlond.gt.233.d0.and.tlond.lt.235.7d0.and.&
           k.le.1) then
         WRITE(*,*)'Haha-A',real(tlond),real(tlatd),k_h_c(i,j,bid),real(ep_ratio(i,j,bid))
         WRITE(*,*)'Haha-B',real(tlond),real(tlatd),i,j,real(WORK(i,j)),real(WORK1(i,j))

!TS       WRITE(*,*)'Haha-A',real(tlond),dzr(k),k_h_c(i,j,bid),ep_ratio(i,j,bid),WORK(i,j),WORK1(i,j)
      endif
      enddo
      enddo
      T_SOURCE(:,:,2) = T_SOURCE(:,:,2)+WORK1
 end subroutine add_estuary_param

!***********************************************************************

!***********************************************************************

   subroutine set_ep

   logical (kind=log_kind) :: first_call = .true.

   integer (int_kind) ::     &
      iblock
 
!-----------------------------------------------------------------------
!
!   Update equavalent depth h_c; update every new month
!
!-----------------------------------------------------------------------
      if (imonth .ne. imonth_last) then      ! update every timestep for 'model'
         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock=1,nblocks_clinic
               k_h_c(:,:,iblock) = k_h_c_DATA(:,:,iblock,imonth)            ! constant across month
         end do
         !$OMP END PARALLEL DO
         IF(MY_TASK.EQ.0)WRITE(stdout,*)'Update EP month',imonth
      endif

   first_call = .false.

   end subroutine set_ep

!***********************************************************************
!BOP
! !IROUTINE: set_estuary_vsf_forcing
! !INTERFACE:

 subroutine set_estuary_vsf_forcing

! !DESCRIPTION:
!  This routine calucates the salinity flux through the sea surface
!  for the virtual salt flux forcing option. It uses the local SSS.
!
! Need to add code to compute correction for global conservation
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: iblock

   real(r8) :: gsum_vsf_sref, gsum_vsf_sloc

!-----------------------------------------------------------------------

   ! Compute surface virtual salt flux using local salinity
   do iblock=1,nblocks_clinic
      FLUX_ROFF_VSF_SRF(:,:,iblock) = -fwmass_to_fwflux*ROFF_F(:,:,iblock) &
          * MAX(TRACER(:,:,1,2,curtime,iblock), c0) * RCALCT(:,:,iblock)
!-----------------------------------------------------------------------
!
!  compute diagnostics if necessary.
!
!-----------------------------------------------------------------------
   call accumulate_tavg_field(FLUX_ROFF_VSF_SRF(:,:,iblock), &
                                    tavg_FLUX_ROFF_VSF_SRF,iblock,1)
      !QS: diag output
      IF (my_task.EQ.47) THEN
        WRITE(stdout,71)iblock,ROFF_F(34,30,iblock), &
            &  TRACER(34,30,1,2,curtime,iblock),     &
            &  FLUX_ROFF_VSF_SRF(34,30,iblock)
71      FORMAT('QS: tot_ROFF_VSF (iblock,ROFF_F,SSS,ROFF_VSF_SRF): ', &
            &  1(i3,1x),3(1pe11.4,1x))
      END IF

   enddo

   ! Compute correction for global conservation
   gsum_vsf_sref = global_sum_prod(ROFF_F,TAREA,distrb_clinic,&
        & field_loc_center, RCALCT)*salinity_factor

   gsum_vsf_sloc = global_sum_prod(FLUX_ROFF_VSF_SRF,TAREA,distrb_clinic,&
        & field_loc_center, RCALCT)

   vsf_river_correction = (gsum_vsf_sref - gsum_vsf_sloc)/area_t

!!$   if ( my_task == master_task) then
!!$      write(stdout,*) ' ### ESTUARY FORCING ###'
!!$      write(stdout,*) '   vsf_sref = ', gsum_vsf_sref
!!$      write(stdout,*) '   vsf_sloc = ', gsum_vsf_sloc
!!$      write(stdout,*) '   vsf correction = ', vsf_river_correction
!!$   endif

   return
 end subroutine set_estuary_vsf_forcing
   
   end module estuary_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
