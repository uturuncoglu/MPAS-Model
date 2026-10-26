module mpas_coupling_diagnostics

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
    use mpas_kind_types, only : RKIND, STRKIND
    use mpas_constants
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: diag_physics

    type (MPAS_clock_type), pointer :: clock

    real (kind=RKIND), pointer :: config_dt
    real (kind=RKIND), pointer :: config_bucket_rainc
    real (kind=RKIND), pointer :: config_bucket_rainnc
    character (len=STRKIND), pointer :: config_bucket_update

    ! Alarm used to advance the coupling-period precipitation accumulation/rate bookkeeping
    ! exactly once per coupling period, regardless of how many internal MPAS timesteps
    ! (atm_core_run_advance / coupling_diagnostics_compute calls) occur within that period.
    character (len=*), parameter :: cplDiagAlarmID = 'coupling_diagnostics'
    logical :: cplDiagAlarmCreated = .false.

    public :: coupling_diagnostics_setup, &
              coupling_diagnostics_compute

    contains

    subroutine coupling_diagnostics_setup(configs, all_pools, simulation_clock)

        use mpas_pool_routines, only : mpas_pool_get_subpool
        use mpas_pool_routines, only : mpas_pool_get_config

        implicit none

        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
        call mpas_pool_get_config(configs, 'config_dt', config_dt)
        call mpas_pool_get_config(configs, 'config_bucket_update', config_bucket_update)
        call mpas_pool_get_config(configs, 'config_bucket_rainc', config_bucket_rainc)
        call mpas_pool_get_config(configs, 'config_bucket_rainnc', config_bucket_rainnc)

        clock => simulation_clock

    end subroutine coupling_diagnostics_setup

    subroutine coupling_diagnostics_compute(configs)

        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
        use mpas_pool_routines, only : mpas_pool_get_config
        use mpas_timekeeping, only : MPAS_Time_Type, MPAS_TimeInterval_Type, MPAS_START_TIME, &
                                      mpas_get_clock_time, mpas_set_timeInterval, &
                                      mpas_add_clock_alarm, mpas_is_alarm_ringing, mpas_reset_clock_alarm, &
                                      operator(+)

        implicit none

        type (MPAS_pool_type), pointer :: configs

        integer :: iCell
        integer, pointer :: nCells
        integer, pointer :: varIndx
        real (kind=RKIND) :: wspd_sfc
        type (MPAS_Time_Type) :: cplAlarmStartTime
        type (MPAS_TimeInterval_Type) :: cplAlarmTimeStep
        logical :: cplPeriodComplete
        integer :: alarm_ierr

        integer, dimension(:), pointer:: i_rainc
        integer, dimension(:), pointer:: i_rainnc

        real (kind=RKIND), pointer :: config_cpl_dt
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: zGrid
        real (kind=RKIND), dimension(:,:), pointer :: rho
        real (kind=RKIND), dimension(:,:), pointer :: uReconstructZonal, uReconstructMeridional
        real (kind=RKIND), dimension(:,:), pointer :: theta
        real (kind=RKIND), dimension(:,:), pointer :: exner
        real (kind=RKIND), dimension(:), pointer :: hgt_lowest
        real (kind=RKIND), dimension(:), pointer :: qv_lowest
        real (kind=RKIND), dimension(:), pointer :: tmp_lowest
        real (kind=RKIND), dimension(:), pointer :: ust, wspd
        real (kind=RKIND), dimension(:), pointer :: taux, tauy
        real (kind=RKIND), dimension(:), pointer :: lwnetb, lwdnb, lwupb
        real (kind=RKIND), dimension(:), pointer :: rainc, rainnc
        real (kind=RKIND), dimension(:), pointer :: snownc, graupelnc
        real (kind=RKIND), dimension(:), pointer :: rainnc_total, rainc_total
        real (kind=RKIND), dimension(:), pointer :: snownc_total, snowc_total
        real (kind=RKIND), dimension(:), pointer :: rain_total, snow_total
        real (kind=RKIND), dimension(:), pointer :: liq_acc_conv
        real (kind=RKIND), dimension(:), pointer :: liq_acc_nonconv, frz_acc_nonconv
        real (kind=RKIND), dimension(:), pointer :: liq_acc_conv_prev
        real (kind=RKIND), dimension(:), pointer :: liq_acc_nonconv_prev, frz_acc_nonconv_prev

        call mpas_log_write('')
        call mpas_log_write('Calculating coupling specific diagnostics:')

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)

        ! querying config_cpl_dt in here to get correct value.
        ! config_cpl_dt is set in the cap just after calling
        ! mpas_init(), which already calls coupling_diagnostics_setup()
        ! so, we could get correct value in coupling_diagnostics_setup() 
        call mpas_pool_get_config(configs, 'config_cpl_dt', config_cpl_dt)

        ! coupling_diagnostics_compute is called once per internal MPAS timestep (config_dt), but a
        ! coupling period (config_cpl_dt) can span several internal timesteps. The precipitation
        ! accumulation/rate bookkeeping below (the liq_acc_*_prev update and the corresponding rate
        ! calculation) must only advance once per coupling period, not once per internal timestep,
        ! otherwise it ends up differencing against a baseline only one internal timestep old while
        ! still dividing by the full (larger) config_cpl_dt, which silently discards most of the
        ! period's precipitation and badly misscales what little it does capture. Use a clock alarm,
        ! created here (once config_cpl_dt has its real coupling-period value; it is not yet set to
        ! that value when coupling_diagnostics_setup() runs) to gate that part of the computation.
        if (.not. cplDiagAlarmCreated) then
           call mpas_set_timeInterval(cplAlarmTimeStep, dt=config_cpl_dt, ierr=alarm_ierr)
           cplAlarmStartTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr=alarm_ierr) + cplAlarmTimeStep
           call mpas_add_clock_alarm(clock, cplDiagAlarmID, cplAlarmStartTime, cplAlarmTimeStep, ierr=alarm_ierr)
           cplDiagAlarmCreated = .true.
        end if

        cplPeriodComplete = mpas_is_alarm_ringing(clock, cplDiagAlarmID, ierr=alarm_ierr)
        if (cplPeriodComplete) then
           call mpas_reset_clock_alarm(clock, cplDiagAlarmID, ierr=alarm_ierr)
        end if

        call mpas_pool_get_array(mesh, 'zgrid', zGrid)
        call mpas_pool_get_dimension(state, 'index_qv', varIndx)
        call mpas_pool_get_array(state, 'scalars', scalars, 1)
        call mpas_pool_get_array(diag, 'theta', theta)
        call mpas_pool_get_array(diag, 'exner', exner)
        call mpas_pool_get_array(diag, 'rho', rho)
        call mpas_pool_get_array(diag_physics, 'ust', ust)
        call mpas_pool_get_array(diag_physics, 'wspd', wspd)
        call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
        call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
        call mpas_pool_get_array(diag_physics, 'lwdnb', lwdnb)
        call mpas_pool_get_array(diag_physics, 'lwupb', lwupb)
        call mpas_pool_get_array(diag_physics, 'rainc', rainc)
        call mpas_pool_get_array(diag_physics, 'rainnc', rainnc)
        call mpas_pool_get_array(diag_physics, 'snownc', snownc)
        call mpas_pool_get_array(diag_physics, 'graupelnc', graupelnc)
        call mpas_pool_get_array(diag_physics, 'i_rainc', i_rainc)
        call mpas_pool_get_array(diag_physics, 'i_rainnc', i_rainnc)

        call mpas_pool_get_array(diag, 'hgt_lowest', hgt_lowest)
        call mpas_pool_get_array(diag, 'qv_lowest', qv_lowest)
        call mpas_pool_get_array(diag, 'tmp_lowest', tmp_lowest)
        call mpas_pool_get_array(diag, 'taux', taux)
        call mpas_pool_get_array(diag, 'tauy', tauy)
        call mpas_pool_get_array(diag, 'lwnetb', lwnetb)
        call mpas_pool_get_array(diag, 'snownc_total', snownc_total)
        call mpas_pool_get_array(diag, 'snowc_total', snowc_total)
        call mpas_pool_get_array(diag, 'rainnc_total', rainnc_total)
        call mpas_pool_get_array(diag, 'rainc_total', rainc_total)
        call mpas_pool_get_array(diag, 'snow_total', snow_total)
        call mpas_pool_get_array(diag, 'rain_total', rain_total)
        call mpas_pool_get_array(diag, 'liq_acc_conv', liq_acc_conv)
        call mpas_pool_get_array(diag, 'liq_acc_nonconv', liq_acc_nonconv)
        call mpas_pool_get_array(diag, 'frz_acc_nonconv', frz_acc_nonconv)
        call mpas_pool_get_array(diag, 'liq_acc_conv_prev', liq_acc_conv_prev)
        call mpas_pool_get_array(diag, 'liq_acc_nonconv_prev', liq_acc_nonconv_prev)
        call mpas_pool_get_array(diag, 'frz_acc_nonconv_prev', frz_acc_nonconv_prev)

        do iCell = 1,nCells
            ! height at lowest model layer, m
            hgt_lowest(iCell) = 0.5_RKIND*(zGrid(2,iCell) - zGrid(1,iCell))

            ! specific humidity at lowest model layer, kg kg-1
            ! scalars(index_qv,...) is the water vapor MIXING RATIO w (kg vapor / kg dry air);
            ! convert to specific humidity q (kg vapor / kg moist air) via q = w / (1 + w).
            qv_lowest(iCell) = scalars(varIndx,1,iCell) / (1.0_RKIND + scalars(varIndx,1,iCell))

            ! temperature at lowest model layer, K
            tmp_lowest(iCell) = theta(1,iCell)*exner(1,iCell)

            ! surface momentum flux components, N m-2
            ! use a local, floored copy of wspd (min 0.1 m/s, same floor the surface-layer
            ! scheme itself uses, e.g. sf_sfclayrev.F90) so we don't divide by zero/near-zero
            wspd_sfc = max(min(wspd(iCell), 1.e+10_RKIND), 0.1_RKIND)
            taux(iCell) = rho(1,iCell) * ust(iCell)**2 * ( uReconstructZonal(1,iCell)/wspd_sfc )
            tauy(iCell) = rho(1,iCell) * ust(iCell)**2 * ( uReconstructMeridional(1,iCell)/wspd_sfc )

            ! net longwave radiation, W m-2
            lwnetb(iCell) = lwdnb(iCell)-lwupb(iCell)

            ! calculate total accumulations for rainnc and rainc, accounting for bucket accumulation if present
            ! NOTE: this is always recomputed from the raw physics bucket accumulators (rainc, rainnc,
            ! i_rainc, i_rainnc), which are themselves restart fields, so this total is reproduced exactly
            ! across a restart.
            if (config_bucket_update /= "none") then
               if (config_bucket_rainnc .gt. 0._RKIND) then
                  rainnc_total(iCell) = i_rainnc(iCell) * config_bucket_rainnc + rainnc(iCell)
               else
                  rainnc_total(iCell) = rainnc(iCell)
               end if
               if (config_bucket_rainc .gt. 0._RKIND) then
                  rainc_total(iCell) = i_rainc(iCell) * config_bucket_rainc + rainc(iCell)
               else
                  rainc_total(iCell) = rainc(iCell)
               end if
            else
               rainnc_total(iCell) = rainnc(iCell)
               rainc_total(iCell) = rainc(iCell)
            end if

            ! calculate accumulations of liquid and frozen precipitation since start of simulation
            ! convective
            liq_acc_conv(iCell) = rainc_total(iCell)

            ! non-convective
            frz_acc_nonconv(iCell) = snownc(iCell) + graupelnc(iCell)
            liq_acc_nonconv(iCell) = rainnc_total(iCell) - frz_acc_nonconv(iCell)

            ! calculate average liquid and frozen rates over the previous coupling period, and
            ! advance the liq_acc_*_prev baseline -- but only once per coupling period (see note
            ! above), not on every internal-timestep call to this routine.
            if (cplPeriodComplete) then
               ! convective
               snowc_total(iCell) = 0._RKIND
               rainc_total(iCell) = (liq_acc_conv(iCell) - liq_acc_conv_prev(iCell)) / config_cpl_dt

               ! non-convective
               snownc_total(iCell) = (frz_acc_nonconv(iCell) - frz_acc_nonconv_prev(iCell)) / config_cpl_dt
               rainnc_total(iCell) = (liq_acc_nonconv(iCell) - liq_acc_nonconv_prev(iCell)) / config_cpl_dt

               ! total
               snow_total(iCell) = snowc_total(iCell) + snownc_total(iCell)
               rain_total(iCell) = rainc_total(iCell) + rainnc_total(iCell)

               ! save current accumulation values as previous values for next coupling time
               ! convective
               liq_acc_conv_prev(iCell) = liq_acc_conv(iCell)

               ! non-convective
               frz_acc_nonconv_prev(iCell) = frz_acc_nonconv(iCell)
               liq_acc_nonconv_prev(iCell) = liq_acc_nonconv(iCell)
            end if
        end do

    end subroutine coupling_diagnostics_compute

    !subroutine coupling_diagnostics_reset()
    !end coupling_diagnostics_reset

    !subroutine coupling_diagnostics_cleanup()
    !end coupling_diagnostics_cleanup

end module mpas_coupling_diagnostics
