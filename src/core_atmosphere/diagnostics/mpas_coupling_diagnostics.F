module mpas_coupling_diagnostics

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
    use mpas_kind_types, only : RKIND
    use mpas_constants
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: diag_physics

    type (MPAS_clock_type), pointer :: clock

    real (kind=RKIND), pointer :: config_dt

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

        clock => simulation_clock

    end subroutine coupling_diagnostics_setup

    !subroutine coupling_diagnostics_update()
    !end subroutine coupling_diagnostics_update

    subroutine coupling_diagnostics_compute()

        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none

        integer :: iCell
        integer, pointer :: nCells

        real (kind=RKIND), dimension(:,:), pointer :: rho
        real (kind=RKIND), dimension(:,:), pointer :: uReconstructZonal, uReconstructMeridional
        real (kind=RKIND), dimension(:), pointer :: ust, wspd
        real (kind=RKIND), dimension(:), pointer :: taux, tauy
        real (kind=RKIND), dimension(:), pointer :: lwnetb, lwdnb, lwupb
        real (kind=RKIND), dimension(:), pointer :: raincv, rainncv
        real (kind=RKIND), dimension(:), pointer :: snowncv, graupelncv 
        real (kind=RKIND), dimension(:), pointer :: rain_total, snow_total

        call mpas_log_write('')
        call mpas_log_write('Calculating coupling specific diagnostics:')

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)

        call mpas_pool_get_array(diag, 'taux', taux)
        call mpas_pool_get_array(diag, 'tauy', tauy)
        call mpas_pool_get_array(diag, 'rho', rho)
        call mpas_pool_get_array(diag_physics, 'ust', ust)
        call mpas_pool_get_array(diag_physics, 'wspd', wspd)
        call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
        call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)

        call mpas_log_write('   surface momentum flux')
        do iCell = 1,nCells
            wspd = min(wspd, 1.e+10)
            taux(iCell) = rho(1,iCell) * ust(iCell)**2 * ( uReconstructZonal(1,iCell)/wspd(iCell) )
            tauy(iCell) = rho(1,iCell) * ust(iCell)**2 * ( uReconstructMeridional(1,iCell)/wspd(iCell) )
        enddo

        call mpas_pool_get_array(diag, 'lwnetb', lwnetb)
        call mpas_pool_get_array(diag_physics, 'lwdnb', lwdnb)
        call mpas_pool_get_array(diag_physics, 'lwupb', lwupb)

        call mpas_log_write('   net longwave radiation')
        do iCell = 1,nCells
            lwnetb(iCell) = lwdnb(iCell)-lwupb(iCell)
        enddo

        call mpas_pool_get_array(diag, 'snow_total', snow_total)
        call mpas_pool_get_array(diag, 'rain_total', rain_total)
        call mpas_pool_get_array(diag_physics, 'raincv', raincv)
        call mpas_pool_get_array(diag_physics, 'rainncv', rainncv)
        call mpas_pool_get_array(diag_physics, 'snowncv', snowncv)
        call mpas_pool_get_array(diag_physics, 'graupelncv', graupelncv)

        call mpas_log_write('   instantaneous precipitation')
        do iCell = 1,nCells
            snow_total(iCell) = (snowncv(iCell) + graupelncv(iCell))/config_dt
            rain_total(iCell) = (rainncv(iCell) + raincv(iCell))/config_dt - snow_total(iCell)
        enddo

    end subroutine coupling_diagnostics_compute

    !subroutine coupling_diagnostics_reset()
    !end coupling_diagnostics_reset

    !subroutine coupling_diagnostics_cleanup()
    !end coupling_diagnostics_cleanup

end module mpas_coupling_diagnostics
