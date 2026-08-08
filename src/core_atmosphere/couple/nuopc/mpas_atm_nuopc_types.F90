module mpas_atm_nuopc_types

  use mpas_kind_types, only: rkind, r8kind, strkind
  use mpas_derived_types, only: core_type, domain_type
  use mpas_derived_types, only: block_type, mpas_pool_type, mpas_time_type

  public

  type mpas_cpl_type
     type(core_type), pointer :: corelist => null()
     type(domain_type), pointer :: domain => null()
     type(block_type), pointer :: block_ptr => null()
     real(kind=rkind) :: dt
     real(kind=rkind) :: dt_cpl
     logical, pointer :: config_do_restart => null()
     character(len=strkind), pointer :: config_restart_timestamp_name => null()
     real(kind=r8kind) :: diag_start_time
     real(kind=r8kind) :: diag_stop_time
     type(mpas_pool_type), pointer :: state => null()
     type(mpas_pool_type), pointer :: diag => null()
     type(mpas_pool_type), pointer :: diag_physics => null()
     type(mpas_pool_type), pointer :: mesh => null()
     real(kind=r8kind) :: input_start_time
     real(kind=r8kind) :: input_stop_time
     real(kind=r8kind) :: output_start_time
     real(kind=r8kind) :: output_stop_time
     real(kind=r8kind) :: integ_start_time
     real(kind=r8kind) :: integ_stop_time
     logical, pointer :: config_apply_lbcs => null()
     type(mpas_time_type), pointer :: currTime => null()
     type (mpas_pool_type), pointer :: tend => null()
     type (mpas_pool_type), pointer :: tend_physics => null()
     character(len=strkind) :: timestamp
     integer :: itimestep
     character(len=strkind) :: input_stream
     character(len=strkind) :: read_time
     integer :: stream_dir
     type(mpas_pool_type), pointer :: sfc_input
     logical, pointer :: enable_import => null()
  end type mpas_cpl_type

  type(mpas_cpl_type), target :: mpas_cpl

end module mpas_atm_nuopc_types
