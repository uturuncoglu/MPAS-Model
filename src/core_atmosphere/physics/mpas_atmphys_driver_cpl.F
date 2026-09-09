! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
module mpas_atmphys_driver_cpl
 
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_atmphys_constants

 implicit none
 private
 public :: driver_cpl

!MPAS driver for coupling 
!Ufuk Turuncoglu (send comments to turuncu@ucar.edu).
!2026-05-18.

contains

!=================================================================================================================
 subroutine driver_cpl(configs,mesh,diag,diag_physics,sfc_input,coupling)
!=================================================================================================================

 !input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag

 !inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: coupling

 integer :: iCell
 integer, pointer:: nCellsSolve
 real(kind=RKIND), dimension(:), pointer:: sst
 real(kind=RKIND), dimension(:), pointer:: skintemp
 real(kind=RKIND), dimension(:), pointer:: xland
 real(kind=RKIND), dimension(:), pointer:: hfx
 real(kind=RKIND), dimension(:), pointer:: lh
 real(kind=RKIND), dimension(:), pointer:: ust
 integer, dimension(:), pointer:: mask_cpl
 real(kind=RKIND), dimension(:), pointer:: sst_cpl
 real(kind=RKIND), dimension(:), pointer:: uo_cpl
 real(kind=RKIND), dimension(:), pointer:: vo_cpl
 real(kind=RKIND), dimension(:), pointer:: hfx_cpl
 real(kind=RKIND), dimension(:), pointer:: lh_cpl
 real(kind=RKIND), dimension(:), pointer:: ust_cpl
!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine driver_cpl:')

 call mpas_pool_get_dimension(mesh    , 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_array(sfc_input   , 'sst'        , sst)
 call mpas_pool_get_array(sfc_input   , 'skintemp'   , skintemp)
 call mpas_pool_get_array(sfc_input   , 'xland'      , xland)
 call mpas_pool_get_array(diag_physics, 'hfx'        , hfx)
 call mpas_pool_get_array(diag_physics, 'lh'         , lh)
 call mpas_pool_get_array(coupling    , 'mask_cpl'   , mask_cpl)
 call mpas_pool_get_array(coupling    , 'sst_cpl'    , sst_cpl)
 call mpas_pool_get_array(coupling    , 'uo_cpl'     , uo_cpl)
 call mpas_pool_get_array(coupling    , 'vo_cpl'     , vo_cpl)
 call mpas_pool_get_array(coupling    , 'hfx_cpl'    , hfx_cpl)
 call mpas_pool_get_array(coupling    , 'lh_cpl'     , lh_cpl)

 do iCell = 1, nCellsSolve
    ! only update over ocean and active coupling region
    if(xland(iCell) .gt. 1.5 .and. mask_cpl(iCell) .eq. 0) then
       sst(iCell) = sst_cpl(iCell)
       skintemp(iCell) = sst(iCell)
       hfx(iCell) = hfx_cpl(iCell)
       lh(iCell) = lh_cpl(iCell)
    end if
 end do

 call mpas_log_write('--- end subroutine driver_cpl.')

 end subroutine driver_cpl

!=================================================================================================================
end module mpas_atmphys_driver_cpl
!=================================================================================================================
