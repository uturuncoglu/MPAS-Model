module mpas_atm_nuopc

  !-----------------------------------------------------------------------------
  ! This is the NUOPC cap for MPAS-Atmosphere 
  !-----------------------------------------------------------------------------

  use ESMF, only: operator(+), operator(/=)
  use ESMF, only: ESMF_GridComp, ESMF_GridCompSetEntryPoint, ESMF_GridCompGet
  use ESMF, only: ESMF_VM, ESMF_VMGet
  use ESMF, only: ESMF_State, ESMF_StateGet, ESMF_Field
  use ESMF, only: ESMF_Clock, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockPrint
  use ESMF, only: ESMF_Time, ESMF_TimeGet, ESMF_TimePrint
  use ESMF, only: ESMF_TimeInterval, ESMF_TimeIntervalGet
  use ESMF, only: ESMF_LogWrite, ESMF_LogSetError, ESMF_MeshWriteVTK
  use ESMF, only: ESMF_Mesh, ESMF_MeshCreate, ESMF_FILEFORMAT_ESMFMESH
  use ESMF, only: ESMF_DistGrid, ESMF_DistGridCreate
  use ESMF, only: ESMF_SUCCESS, ESMF_FAILURE, ESMF_RC_VAL_WRONG
  use ESMF, only: ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
  use ESMF, only: ESMF_METHOD_INITIALIZE
  use ESMF, only: ESMF_KIND_R8

  use NUOPC, only: NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC, only: NUOPC_SetAttribute, NUOPC_IsUpdated
  use NUOPC, only: NUOPC_CompDerive, NUOPC_CompFilterPhaseMap
  use NUOPC, only: NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize

  use NUOPC_Model, only: SetVM
  use NUOPC_Model, only: NUOPC_ModelGet
  use NUOPC_Model, only: model_routine_SS => SetServices
  use NUOPC_Model, only: model_label_SetClock => label_SetClock
  use NUOPC_Model, only: model_label_Advance => label_Advance

  use mpas_atm_nuopc_shr, only: ChkErr
  use mpas_atm_nuopc_types, only: mpas_cpl_type, mpas_cpl
  use mpas_atm_nuopc_flds, only: advertise_fields
  use mpas_atm_nuopc_flds, only: realize_fields
  use mpas_atm_nuopc_flds, only: export_fields
  use mpas_atm_nuopc_flds, only: import_fields
  use mpas_atm_nuopc_flds, only: state_diagnose

  use mpas_derived_types, only: block_type, mpas_pool_type
  use mpas_pool_routines, only: mpas_pool_get_subpool
  use mpas_pool_routines, only: mpas_pool_get_dimension
  use mpas_pool_routines, only: mpas_pool_get_array
  use mpas_pool_routines, only: mpas_pool_get_config
  use atm_core, only: atm_core_run_start, atm_core_run_advance
  use mpas_subdriver, only: mpas_init, mpas_finalize

  implicit none
  private

  !-----------------------------------------------------------------------------
  ! Public module routines
  !-----------------------------------------------------------------------------

  public SetVM
  public SetServices

  !-----------------------------------------------------------------------------
  ! Module private routines
  !-----------------------------------------------------------------------------

  private :: InitializeP0        ! Phase zero of initialization
  private :: InitializeAdvertise ! Advertise the fields that can be passed

  !-----------------------------------------------------------------------------
  ! Private module data
  !-----------------------------------------------------------------------------

  character(len=*), parameter :: modName = "(mpas_atm_nuopc)"

  character(len=*), parameter :: u_FILE_u = &
     __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname = trim(modName)//':(SetServices) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !------------------
    ! register the generic methods
    !------------------

    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! switching to IPD versions
    !------------------
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !------------------
    ! set entry point for methods that require specific implementation
    !------------------
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !call NUOPC_CompSpecialize(gcomp, specLabel=label_SetClock, &
    !     specRoutine=SetClock, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Advertise coupling fields
    ! ---------------------

    call advertise_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    type(ESMF_Mesh) :: mesh
    type(ESMF_DistGrid) :: distGrid
    integer :: n, nCells, iCell
    integer :: ierr, petCount, localPet, comm
    integer , allocatable   :: gindex(:)
    character(len=255) :: cvalue, mesh_atm
    logical :: isSet, isPresent
    type(block_type), pointer :: block => null()
    type(mpas_pool_type), pointer :: meshPool
    integer, dimension(:), pointer :: indexToCellID
    integer, dimension(:), pointer :: nCellsArray
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! ---------------------
    ! Query VM
    ! ---------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, mpiCommunicator=comm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Initialize MPAS
    ! ---------------------

    call mpas_init(mpas_cpl%corelist, mpas_cpl%domain, external_comm=comm)

    ! ---------------------
    ! Prepare MPAS to run
    ! ---------------------

    ierr = atm_core_run_start( &
       mpas_cpl%domain, &
       mpas_cpl%block_ptr, & 
       mpas_cpl%dt, &
       mpas_cpl%config_do_restart, &
       mpas_cpl%config_restart_timestamp_name, &
       mpas_cpl%diag_start_time, &
       mpas_cpl%diag_stop_time, &
       mpas_cpl%state, &
       mpas_cpl%diag, &
       mpas_cpl%diag_physics, &
       mpas_cpl%mesh, &
       mpas_cpl%output_start_time, &
       mpas_cpl%output_stop_time, &
       mpas_cpl%config_apply_lbcs, &
       mpas_cpl%currTime, &
       mpas_cpl%timestamp, &
       mpas_cpl%itimestep)
    if (ierr /= 0) then
       call ESMF_LogWrite(trim(subname)//": "//' MPAS atm_core_run_start() is failed!. Exiting ...', ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    end if

    ! ---------------------
    ! Query MPAS configuration
    ! ---------------------

    call mpas_pool_get_config(mpas_cpl % domain % blocklist % configs, 'config_enable_import', mpas_cpl % enable_import)

    ! ---------------------
    ! Get coupling specific options
    ! ---------------------

    call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (isPresent .and. isSet) then
      mesh_atm = trim(cvalue)
      call ESMF_LogWrite(trim(subname)//': MPAS ESMF mesh file = '//trim(mesh_atm), ESMF_LOGMSG_INFO)
    else
      call ESMF_LogWrite(trim(subname)//': MPAS ESMF mesh is required! Please set <mesh_atm>. Exiting ....', ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    end if

    ! ---------------------
    ! Determine the global index space needed for the distgrid
    ! ---------------------

    n = 0
    block => mpas_cpl % domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray(1)
       n = n + nCells
       block => block % next
    end do
    allocate(gindex(n))

    n = 0
    block => mpas_cpl % domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       nCells = nCellsArray(1)
       do iCell = 1, nCells
          gindex(n+iCell) = indexToCellID(iCell)
       enddo
       n = n + nCells
       block => block % next
    end do

    ! ---------------------
    ! Create distGrid from global index array
    ! ---------------------

    distGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Create MPAS mesh
    ! ---------------------

    mesh = ESMF_MeshCreate(filename=trim(mesh_atm), fileformat=ESMF_FILEFORMAT_ESMFMESH, elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !NOTE: ESMF framework has bug to write high-order meshes in VTK format but will work >= 9.0.0b13
    !call ESMF_MeshWriteVTK(mesh, filename="mpas_mesh", rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Realize coupling fields
    ! ---------------------

    call realize_fields(importState, exportState, mesh, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine SetClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(SetClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetClock

  !===============================================================================

  subroutine DataInitialize(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
  
    ! local variables
    integer :: n, fieldCount
    character(len=64), allocatable :: fieldNameList(:)
    type(ESMF_Field) :: field
    type(ESMF_Clock) :: clock
    type(ESMF_State) :: importState, exportState
    character(len=*), parameter :: subname=trim(modName)//':(DataInitialize) '
    !-------------------------------------------------------------------------------
  
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------
    ! Query the Component for its clock, importState and exportState
    !-----------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------
    ! Update export state
    !-----------------------

    call export_fields(exportState, mpas_cpl%domain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_diagnose(exportState, 'export', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------
    ! Update attribute of the fields in export state
    !-----------------------

    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount
       call ESMF_StateGet(exportState, itemName=fieldNameList(n), field=field, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    deallocate(fieldNameList)

    !-----------------------
    ! Check whether all Fields in the exportState are "Updated"
    !-----------------------

    if (NUOPC_IsUpdated(exportState)) then
       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite("MPAS - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine DataInitialize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: n, ierr
    integer, save :: nSteps = 0
    logical, save :: first_time = .true.
    real(ESMF_KIND_R8) :: dt_cpl
    type(ESMF_Clock) :: clock
    type(ESMF_Time) :: startTime, currTime
    type(ESMF_TimeInterval) :: timestep
    type(ESMF_State) :: importState, exportState
    character(len=255) :: msgString 
    character(len=*), parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------
    ! Query the Component for its clock, importState and exportState
    !-----------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ---------------------
    ! Calculate number of MPAS advance steps
    ! ---------------------

    if (first_time) then
       ! Query model clock
       call ESMF_ClockGet(clock, timestep=timestep, rc=rc) 
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeIntervalGet(timestep, s_r8=dt_cpl, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Check if coupling time step is evenly divisible by MPAS time step or not 
       if (mpas_cpl%dt /= 0.0d0 .and. abs(dt_cpl/mpas_cpl%dt - nint(dt_cpl/mpas_cpl%dt)) < 1.0d-12) then
          nSteps = nint(dt_cpl/mpas_cpl%dt)
       else
          call ESMF_LogWrite(trim(subname)//": "//&
             "Coupling time step must be evenly divisible by MPAS time step!", ESMF_LOGMSG_ERROR)
          write(msgString, fmt="(A,F8.1,A,F8.1,A)") &
             "Coupling time step is ", dt_cpl, "s but MPAS time step is ", mpas_cpl%dt, "s" 
          call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
       write(msgString, fmt="(A,F8.1,A)") "Coupling time step =", dt_cpl, "s"
       call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_INFO)
       write(msgString, fmt="(A,F8.1,A)") "MPAS time step = ", mpas_cpl%dt, "s"
       call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_INFO)
       write(msgString, fmt="(A,I5)") "MPAS will advance #step = ", nSteps
       call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_INFO)

       first_time = .false.
    end if

    !----------------------
    ! Ingest data from import state
    !----------------------

    if (mpas_cpl % enable_import) then
       call import_fields(importState, mpas_cpl%domain, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_diagnose(importState, 'import', rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)// &
         ": config_enable_import set to False in MPAS configuration. Skip importing", ESMF_LOGMSG_INFO)
    end if

    ! ---------------------
    ! Run MPAS
    ! ---------------------

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    call ESMF_ClockPrint(clock, options="currTime", &
       preString="------>Advancing ATM from: ", unit=msgString, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
       timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimePrint(currTime + timeStep, &
       preString="--------------------------------> to: ", unit=msgString, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_INFO)

    do n = 1, nSteps
       write(msgString, fmt="(A,I8)") "Advancing MPAS, itimestep = ", mpas_cpl%itimestep
       call ESMF_LogWrite(trim(subname)//": "//trim(msgString), ESMF_LOGMSG_INFO) 
       ierr = atm_core_run_advance( &
          mpas_cpl%domain, &
          mpas_cpl%timestamp, &
          mpas_cpl%block_ptr, &
          mpas_cpl%config_apply_lbcs, &
          mpas_cpl%input_start_time, &
          mpas_cpl%input_stop_time, &
          mpas_cpl%output_start_time, &
          mpas_cpl%output_stop_time, &
          mpas_cpl%input_stream, &
          mpas_cpl%read_time, &
          mpas_cpl%stream_dir, &
          mpas_cpl%integ_start_time, &
          mpas_cpl%integ_stop_time, &
          mpas_cpl%diag_start_time, &
          mpas_cpl%diag_stop_time, &
          mpas_cpl%dt, &
          mpas_cpl%itimestep, &
          mpas_cpl%state, &
          mpas_cpl%mesh, &
          mpas_cpl%diag, &
          mpas_cpl%diag_physics, &
          mpas_cpl%tend, &
          mpas_cpl%tend_physics, &
          mpas_cpl%config_restart_timestamp_name)
       if (ierr /= 0) then
          call ESMF_LogWrite(trim(subname)//": "//' MPAS atm_core_run_advance() is failed!. Exiting ...', ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       end if
    end do

    !----------------------
    ! Put updated data to export state
    !----------------------

    call export_fields(exportState, mpas_cpl%domain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_diagnose(exportState, 'export', rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock) :: dclock, mclock
    type(ESMF_Time) :: dcurrtime, dstoptime
    type(ESMF_Time) :: mcurrtime, mstoptime
    type(ESMF_TimeInterval) :: dtimestep, mtimestep
    character(len=128) :: dtimestring, mtimestring
    character(len=*), parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------
  
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------
    ! Query the component for its clock
    !----------------------

    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=dstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Check that the current time in the model and driver are the same
    !--------------------------------

    if (mcurrtime /= dcurrtime) then
      call ESMF_TimeGet(dcurrtime, timeString=dtimestring, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_TimeGet(mcurrtime, timeString=mtimestring, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
           msg=subname//": ERROR in time consistency: "//trim(dtimestring)//" != "//trim(mtimestring),  &
           line=__LINE__, file=__FILE__, rcToReturn=rc)
      return
    endif

    !--------------------------------
    ! Force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelCheckImport(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
  
    ! local variables
    character(len=*), parameter :: subname=trim(modName)//':(ModelCheckImport) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

  
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelCheckImport

  subroutine ModelFinalize(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: ierr
    character(len=*), parameter :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Finalize model
    !--------------------------------

    call mpas_finalize(mpas_cpl%corelist, mpas_cpl%domain)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

end module mpas_atm_nuopc
