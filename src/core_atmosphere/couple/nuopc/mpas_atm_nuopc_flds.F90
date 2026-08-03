module mpas_atm_nuopc_flds

  !-----------------------------------------------------------------------------
  ! Import and export fields related routines 
  !-----------------------------------------------------------------------------

  use ESMF, only: operator(==)
  use ESMF, only: ESMF_GridComp, ESMF_LOGMSG_ERROR, ESMF_FAILURE
  use ESMF, only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
  use ESMF, only: ESMF_State, ESMF_StateGet, ESMF_STATEITEM_FIELD
  use ESMF, only: ESMF_Finalize, ESMF_END_ABORT, ESMF_StateItem_Flag
  use ESMF, only: ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_FieldWriteVTK
  use ESMF, only: ESMF_TYPEKIND_R8, ESMF_KIND_R8, ESMF_MAXSTR
  use ESMF, only: ESMF_Field, ESMF_FieldGet, ESMF_Mesh, ESMF_StateRemove
  use ESMF, only: ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
  use ESMF, only: ESMF_FieldWrite, ESMF_UtilString2Double

  use NUOPC, only: NUOPC_Advertise, NUOPC_Realize, NUOPC_IsConnected
  use NUOPC_Model, only: NUOPC_ModelGet

  use mpas_atm_nuopc_shr, only: ChkErr
  use mpas_atm_nuopc_types, only: mpas_cpl_type, mpas_cpl

  use mpas_kind_types, only: rkind
  use mpas_derived_types, only: block_type, mpas_pool_type
  use mpas_derived_types, only: domain_type
  use mpas_pool_routines, only: mpas_pool_get_array
  use mpas_pool_routines, only: mpas_pool_get_subpool
  use mpas_pool_routines, only: mpas_pool_get_dimension

  implicit none
  private

  !-----------------------------------------------------------------------------
  ! Public module routines
  !-----------------------------------------------------------------------------

  public :: ChkErr
  public :: advertise_fields
  public :: realize_fields
  public :: state_diagnose
  public :: export_fields
  public :: import_fields

  !-----------------------------------------------------------------------------
  ! Public module data 
  !-----------------------------------------------------------------------------

  type fldListType
     character(len=128) :: stdname
     character(len=128) :: internalgroup
     character(len=128) :: internalname
     integer :: level = 0
     real(ESMF_KIND_R8) :: scale_factor = 1.0d0
     real(ESMF_KIND_R8) :: add_offset = 0.0d0
     real(ESMF_KIND_R8) :: valid_min = -1.0d20
     real(ESMF_KIND_R8) :: valid_max = 1.0d20
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
     logical :: connected = .false.
  end type fldListType

  integer, parameter :: fldsMax = 30
  integer :: fldsToMPAS_num = 0
  integer :: fldsFrMPAS_num = 0
  type(fldListType) :: fldsToMPAS(fldsMax)
  type(fldListType) :: fldsFrMPAS(fldsMax)

  !-----------------------------------------------------------------------------
  ! Private module data
  !-----------------------------------------------------------------------------

  character(*), parameter :: modName = "(mpas_atm_fields)"

  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    integer,             intent(out) :: rc

    ! local variables
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    integer           :: n, itemCount
    character(len=*), parameter :: subname=trim(modName)//':(advertise_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    ! export scalar
    ! set as constant in here but actually set in export_fields() routine later
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_z', 'diag', 'hgt_lowest', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_u', 'diag', 'uReconstructZonal', level=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_v', 'diag', 'uReconstructMeridional', level=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_tbot', 'diag', 'tmp_lowest', level=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_pbot', 'diag', 'pressure', level=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_shum', 'diag', 'qv_lowest', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_dens', 'diag', 'rho', level=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_ptem', 'diag', 'theta', level=1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Sa_pslv', 'diag', 'surface_pressure', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! export flux
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swnet', 'diag_physics', 'gsw', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_lwdn' , 'diag_physics', 'lwdnb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_lwup' , 'diag_physics', 'lwupb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_lwnet', 'diag', 'lwnetb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swdn' , 'diag_physics', 'swdnb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swup' , 'diag_physics', 'swupb', rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_sen', 'diag_physics', 'hfx', scale_factor=-1.0d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_lat', 'diag_physics', 'lh', scale_factor=-1.0d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_evap', 'diag_physics', 'qfx', scale_factor=-1.0d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_rain', 'diag', 'rain_total', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_snow', 'diag', 'snow_total', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! The ratios used to split net shortwave radiation is taken from CMEPS mediator
    ! Ref: https://github.com/NOAA-EMC/CMEPS/blob/fc8b9140e08465dcb5eab48056d4d5636c0e1716/mediator/med_phases_prep_ocn_mod.F90#L504
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swndr', 'diag_physics', 'swdnb', scale_factor=0.285d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swvdr', 'diag_physics', 'swdnb', scale_factor=0.285d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swndf', 'diag_physics', 'swdnb', scale_factor=0.215d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_swvdf', 'diag_physics', 'swdnb', scale_factor=0.215d0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_taux', 'diag', 'taux', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldlist_add(fldsFrMPAS_num, fldsFrMPAS, 'Faxa_tauy', 'diag', 'tauy', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now advertise above export fields
    do n = 1, fldsFrMPAS_num
       call NUOPC_Advertise(exportState, standardName=fldsFrMPAS(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    ! import from ocn 
    call fldlist_add(fldsToMPAS_num, fldsToMPAS, 'So_t', 'coupling', 'sst_c', valid_min=100.0d0, valid_max = 1.0d20, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Now advertise import fields
    do n = 1, fldsToMPAS_num
       call ESMF_StateGet(importState, itemSearch=trim(fldsToMPAS(n)%stdname), itemCount=itemCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (itemCount == 0) then
          call NUOPC_Advertise(importState, standardName=fldsToMPAS(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine advertise_fields

  !=============================================================================

  subroutine realize_fields(importState, exportState, mesh, rc)

    ! input/output variables
    type(ESMF_State), intent(inout) :: importState
    type(ESMF_State), intent(inout) :: exportState
    type(ESMF_mesh) , intent(in)    :: mesh 
    integer         , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call fldlist_realize( &
         state=exportState, &
         fldList=fldsFrMPAS, &
         numflds=fldsFrMPAS_num, &
         tag=subname//':MPAS_Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToMPAS, &
         numflds=fldsToMPAS_num, &
         tag=subname//':MPAS_Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine realize_fields

  !=============================================================================

  subroutine fldlist_add(num, fldlist, stdname, intgrp, intname, level, &
                scale_factor, add_offset, &
                valid_min, valid_max, &
                ungridded_lbound, ungridded_ubound, rc)

    ! input/output variables
    integer,           intent(inout) :: num
    type(fldListType), intent(inout) :: fldlist(:)
    character(len=*),  intent(in)    :: stdname
    character(len=*),  intent(in)    :: intgrp
    character(len=*),  intent(in)    :: intname
    integer, optional, intent(in)    :: level
    real(ESMF_KIND_R8), optional, intent(in) :: scale_factor 
    real(ESMF_KIND_R8), optional, intent(in) :: add_offset
    real(ESMF_KIND_R8), optional, intent(in) :: valid_min
    real(ESMF_KIND_R8), optional, intent(in) :: valid_max
    integer, optional, intent(in)    :: ungridded_lbound
    integer, optional, intent(in)    :: ungridded_ubound
    integer, optional, intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(fldlist_add)'
    !---------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Set up a list of field information
    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
         ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    endif

    fldlist(num)%stdname = trim(stdname)
    fldlist(num)%internalgroup = trim(intgrp)
    fldlist(num)%internalname = trim(intname)

    if (present(level)) then
       fldlist(num)%level = level
    end if

    if (present(scale_factor)) then
       fldlist(num)%scale_factor = scale_factor
    end if

    if (present(add_offset)) then
       fldlist(num)%add_offset = add_offset
    end if

    if (present(valid_min)) then
       fldlist(num)%valid_min = valid_min
    end if

    if (present(valid_max)) then
       fldlist(num)%valid_max = valid_max
    end if

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine fldlist_add

  !=============================================================================

  subroutine fldlist_realize(state, fldList, numflds, mesh, tag, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: state
    type(fldListType), intent(inout) :: fldList(:)
    integer          , intent(in)    :: numflds
    character(len=*) , intent(in)    :: tag
    type(ESMF_Mesh)  , intent(in)    :: mesh
    integer          , intent(inout) :: rc

    ! local variables
    integer :: n
    type(ESMF_Field) :: field
    character(len=80) :: stdname
    character(len=*), parameter :: subname=trim(modName)//':fldlist_realize)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    do n = 1, numflds
       stdname = trim(fldList(n)%stdname)
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          ! Create the field
          if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                  ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                  ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                  gridToFieldMap=(/2/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
               ESMF_LOGMSG_INFO)

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Set flag for connected fields
          fldList(n)%connected = .true.
       else
          call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
               ESMF_LOGMSG_INFO)
          call ESMF_StateRemove(state, (/stdname/), relaxedflag=.true., rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine fldlist_realize

  !===============================================================================

  subroutine export_fields(exportState, domain, rc)

    ! input/output variables
    type(ESMF_State), intent(inout) :: exportState
    type(domain_type), intent(in), pointer :: domain
    integer, intent(out) :: rc

    ! local variables
    integer :: n, iCell, gCell, nCells, cell_offset
    logical :: apply_conversion
    type(ESMF_Field) :: lfield
    type(ESMF_StateItem_Flag) :: itemType
    type(block_type), pointer :: block => null()
    type(mpas_pool_type), pointer :: meshPool
    type(mpas_pool_type), pointer :: mpasPtrPool
    real(kind=rkind), dimension(:), pointer :: fldPtr
    real(kind=rkind), dimension(:,:), pointer :: fldPtr2d
    real(ESMF_KIND_R8), dimension(:), pointer :: fldPtrExport
    integer, dimension(:), pointer :: nCellsArray
    character(len=*), parameter :: subname=trim(modName)//':(export_fields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! -----------------------
    ! Loop over export fields and update them 
    ! -----------------------

    do n = 1, fldsFrMPAS_num
       ! Check field
       call ESMF_StateGet(exportState, itemName=trim(fldsFrMPAS(n)%stdname), itemType=itemType, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (itemType == ESMF_STATEITEM_FIELD) then
          ! Get field
          call ESMF_StateGet(exportState, itemName=trim(fldsFrMPAS(n)%stdname), field=lfield, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Query field pointer and initialize
          call ESMF_FieldGet(lfield, farrayPtr=fldPtrExport, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          fldPtrExport(:) = 1.0d20

          ! In case of filling with constant value
          if (trim(fldsFrMPAS(n)%internalgroup) == 'const') then
             fldPtrExport(:) = fldsFrMPAS(n)%add_offset
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             cycle
          end if

          ! Check if we need to apply conversion
          apply_conversion = .false.
          if (fldsFrMPAS(n)%scale_factor /= 1.0d0 .or. fldsFrMPAS(n)%add_offset /= 0.0d0) then
             apply_conversion = .true.
          end if

          ! Query internal pointer and fill export field 
          cell_offset = 0
          block => domain % blocklist
          do while (associated(block))
             ! Get number of cells in decomposition block
             call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
             call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
             nCells = nCellsArray(1)

             ! Access internal pointer pool
             call mpas_pool_get_subpool(block % structs, trim(fldsFrMPAS(n)%internalgroup), mpasPtrPool)

             ! Pass data 
             if (fldsFrMPAS(n)%level > 0) then ! 2d field
                ! Access internal field pointer
                call mpas_pool_get_array(mpasPtrPool, trim(fldsFrMPAS(n)%internalname), fldptr2d)

                ! Put data to export field 
                if (.not. associated(fldptr2d)) then
                   ! TODO: Throw error and exit
                   call ESMF_LogWrite(subname//' '//trim(fldsFrMPAS(n)%internalname)//&
                      ' is not found in '//trim(fldsFrMPAS(n)%internalgroup), ESMF_LOGMSG_INFO)
                else
                   if (apply_conversion) then
                      do iCell = 1, nCells
                         gCell = iCell + cell_offset
                         fldPtrExport(gCell) = dble(fldptr2d(fldsFrMPAS(n)%level,iCell))*fldsFrMPAS(n)%scale_factor+fldsFrMPAS(n)%add_offset
                      end do
                   else
                      do iCell = 1, nCells
                         gCell = iCell + cell_offset
                         fldPtrExport(gCell) = dble(fldptr2d(fldsFrMPAS(n)%level,iCell))
                      end do
                   end if
                end if

                ! Nullify pointer
                nullify(fldptr2d)

             else ! 1d field
                ! Access internal field pointer
                call mpas_pool_get_array(mpasPtrPool, trim(fldsFrMPAS(n)%internalname), fldptr)

                ! Put data to export field 
                if (.not. associated(fldptr)) then
                   ! TODO: Throw error and exit
                   call ESMF_LogWrite(subname//' '//trim(fldsFrMPAS(n)%internalname)//&
                      ' is not found in '//trim(fldsFrMPAS(n)%internalgroup), ESMF_LOGMSG_INFO)
                else
                   if (apply_conversion) then
                      do iCell = 1, nCells
                         gCell = iCell + cell_offset
                         fldPtrExport(gCell) = dble(fldptr(iCell))*fldsFrMPAS(n)%scale_factor+fldsFrMPAS(n)%add_offset
                      end do
                   else
                      do iCell = 1, nCells
                         gCell = iCell + cell_offset
                         fldPtrExport(gCell) = dble(fldptr(iCell))
                      end do
                   end if
                end if

                ! Nullify pointer
                nullify(fldptr)

             end if

             ! Increment cell offset
             cell_offset = cell_offset + nCells

             ! Go to next block
             block => block % next

          end do

          ! Init pointers
          nullify(fldPtrExport)
       else
          call ESMF_LogWrite(subname//' '//trim(fldsFrMPAS(n)%stdname)//' is not in the state!', ESMF_LOGMSG_INFO)
       end if
    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine export_fields

  !===============================================================================

  subroutine import_fields(importState, domain, rc)

    ! input/output variables
    type(ESMF_State), intent(inout) :: importState
    type(domain_type), intent(in), pointer :: domain
    integer, intent(out) :: rc

    ! local variables
    integer :: n, iCell, gCell, nCells, cell_offset
    logical :: apply_conversion
    type(ESMF_Field) :: lfield
    type(ESMF_StateItem_Flag) :: itemType
    type(block_type), pointer :: block => null()
    type(mpas_pool_type), pointer :: meshPool
    type(mpas_pool_type), pointer :: mpasPtrPool
    type(mpas_pool_type), pointer :: sfcInputPool
    type(mpas_pool_type), pointer :: couplingPool
    integer, dimension(:), pointer :: mask_c
    real(kind=rkind), dimension(:), pointer :: xland
    real(kind=rkind), dimension(:), pointer :: fldPtr
    real(ESMF_KIND_R8), dimension(:), pointer :: fldPtrImport
    integer, dimension(:), pointer :: nCellsArray
    character(len=*), parameter :: subname=trim(modName)//':(import_fields)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! -----------------------
    ! Loop over export fields and update them 
    ! -----------------------

    do n = 1, fldsToMPAS_num
       ! Check field
       call ESMF_StateGet(importState, itemName=trim(fldsToMPAS(n)%stdname), itemType=itemType, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (itemType == ESMF_STATEITEM_FIELD) then
          ! Get field
          call ESMF_StateGet(importState, itemName=trim(fldsToMPAS(n)%stdname), field=lfield, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Query field pointer and initialize
          call ESMF_FieldGet(lfield, farrayPtr=fldPtrImport, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Check if we need to apply conversion
          apply_conversion = .false.
          if (fldsToMPAS(n)%scale_factor /= 1.0d0 .or. fldsToMPAS(n)%add_offset /= 0.0d0) then
             apply_conversion = .true.
          end if

          ! Query internal pointer and fill export field 
          cell_offset = 0
          block => domain % blocklist
          do while (associated(block))
             ! Access internal pointer
             call mpas_pool_get_subpool(block % structs, trim(fldsToMPAS(n)%internalgroup), mpasPtrPool)
             call mpas_pool_get_array(mpasPtrPool, trim(fldsToMPAS(n)%internalname), fldptr)

             ! Get land sea mask
             call mpas_pool_get_subpool(block % structs, 'sfc_input', sfcInputPool)
             call mpas_pool_get_subpool(block % structs, 'coupling', couplingPool)
             call mpas_pool_get_array(sfcInputPool, 'xland', xland)
             call mpas_pool_get_array(couplingPool, 'mask_c', mask_c)
             
             ! Get number of cells in decomposition block
             call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
             call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
             nCells = nCellsArray(1)

             ! Loop over cells and fill pointer of export field
             mask_c(:) = 1
             if (.not. associated(fldptr)) then
                ! TODO: Throw error and exit
                call ESMF_LogWrite(subname//' '//trim(fldsToMPAS(n)%internalname)//&
                   ' is not found in '//trim(fldsToMPAS(n)%internalgroup), ESMF_LOGMSG_INFO)
             else
                if (apply_conversion) then
                   do iCell = 1, nCells
                      gCell = iCell + cell_offset
                      if(xland(iCell) .gt. 1.5 .and. &
                         fldPtrImport(gCell) >= fldsToMPAS(n)%valid_min .and. &
                         fldPtrImport(gCell) <= fldsToMPAS(n)%valid_max) then
                         fldptr(iCell) = fldPtrImport(gCell)*fldsToMPAS(n)%scale_factor+fldsToMPAS(n)%add_offset
                         mask_c(iCell) = 0
                      end if
                   end do
                else
                   do iCell = 1, nCells
                      gCell = iCell + cell_offset
                      if(xland(iCell) .gt. 1.5 .and. &
                         fldPtrImport(gCell) >= fldsToMPAS(n)%valid_min .and. &
                         fldPtrImport(gCell) <= fldsToMPAS(n)%valid_max) then
                         fldptr(iCell) = fldPtrImport(gCell)
                         mask_c(iCell) = 0
                      end if
                   end do
                end if
             end if 

             ! Increment cell offset
             cell_offset = cell_offset + nCells

             ! Go to next block
             block => block % next

             ! Nullify pointer
             nullify(fldptr)
          end do

          ! Init pointers
          nullify(fldPtrImport)
       else
          call ESMF_LogWrite(subname//' '//trim(fldsToMPAS(n)%stdname)//' is not in the state!', ESMF_LOGMSG_INFO)
       end if
    end do

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine import_fields

  !===============================================================================

  subroutine state_diagnose(state, string, rc)

    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    integer                         :: n
    type(ESMF_Field)                :: lfield
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    real(ESMF_KIND_R8), pointer     :: dataPtr1d(:)
    character(len=1024)             :: msgString
    character(len=*), parameter     :: subname='(state_diagnose)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
       call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldGet(lfield, farrayPtr=dataPtr1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
             minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    deallocate(lfieldnamelist)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine state_diagnose

end module mpas_atm_nuopc_flds
