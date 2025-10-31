#include "../include/dprec.fh"
#include "../include/assert.fh"

module pol_gauss_mdin
  implicit none
  private

  integer :: ipgm                                   ! control if pGM is used
  integer :: pol_gauss_verbose = 0                  ! pGM verbose mode for detailed printing
  integer :: use_average = 0                        ! Use average error for convergence check
  integer :: dipole_print = 0                       ! specify the number of steps as interval to print total dipole
  integer :: saved_nprint = 0                       ! For detailed dipole printing
  integer :: dipole_scf_init = 3                    ! Initial dipole options: 1=alpha E, 2=last step, 3=extrapolation
  integer :: dipole_scf_init_order = 3              ! Extrapolation order(when dipole_scf_init=3)
  integer :: dipole_scf_init_step = 2               ! Number of steps for extrapolation(when dipole_scf_init=3)
  integer :: scf_solv_opt = 3                       ! SCF solver option: 1=CG, 2=PCG, 3=PKCG, 4,5=SOR
  integer :: scf_sor_niter = 100                    ! Maximum number of iterations in SCF SOR
  integer :: scf_cg_niter = 50                      ! Maximum number of iterations in SCF CG
  integer :: scf_local_niter = 3                    ! Maximum number of iterations in local preconditioner
  integer :: pol_gauss_ips = 0                      ! use IPS instead of PME
  _REAL_ :: dipole_scf_tol = 0.01d0                 ! SCF convergence criterion
  _REAL_ :: scf_sor_coefficient = 0.65d0            ! SOR coefficient (Peek step is one additional SOR)
  _REAL_ :: scf_local_cut = 4.0d0                   ! local preconditioner cutoff (when scf_solv_opt = 3)
  _REAL_ :: ee_dsum_cut = 9.0d0                     ! Cutoff value for elec interactions
  _REAL_ :: ee_damped_cut=4.5d0                     ! Cutoff for damped elec interactions (not used)
  _REAL_ :: erfc_tol = 1.0e-5                       ! Erfc tolerance (not used)
  _REAL_ :: ee_gauss_cut = 3.12342d0                ! Max error of 10^-5 in erfc (not used)

  ! pack these into common blocks for easy bcast
  integer, parameter :: pgm_ctrl_int_cnt = 12
  integer, parameter :: pgm_ctrl_dbl_cnt = 7
  common /pgm_ctrl_int/ pol_gauss_verbose, use_average, dipole_print, saved_nprint, &
                        dipole_scf_init, dipole_scf_init_order, dipole_scf_init_step, &
                        scf_solv_opt, scf_sor_niter, scf_cg_niter, scf_local_niter, &
                        pol_gauss_ips
  save :: /pgm_ctrl_int/
  common /pgm_ctrl_dbl/ dipole_scf_tol, scf_sor_coefficient, scf_local_cut, &
                        ee_dsum_cut, ee_damped_cut, erfc_tol, ee_gauss_cut
  save :: /pgm_ctrl_dbl/

  public :: POL_GAUSS_read_mdin, ipgm, pol_gauss_verbose, dipole_print, saved_nprint, &
            dipole_scf_tol, use_average, dipole_scf_init, dipole_scf_init_order, dipole_scf_init_step, &
            scf_solv_opt, scf_cg_niter, scf_local_niter, &
            scf_sor_coefficient, scf_sor_niter, &
            ee_dsum_cut, ee_damped_cut, scf_local_cut, pgm_ctrl_int_cnt, pgm_ctrl_dbl_cnt, &
            pol_gauss_ips

#ifdef MPI
  public :: pGM_mdin_bcast
#endif

  contains
!-------------------------------------------------------------------------------
#ifdef MPI
subroutine pGM_mdin_bcast()
  include 'mpif.h'
# include "extra.h"
# include "parallel.h"

  integer :: ibuf(pgm_ctrl_int_cnt)
  common /pgm_ctrl_int/ ibuf
  save :: /pgm_ctrl_int/
  _REAL_ :: buf(pgm_ctrl_dbl_cnt)
  common /pgm_ctrl_dbl/ buf
  save :: /pgm_ctrl_dbl/

  integer :: ier

  call mpi_bcast(ibuf, pgm_ctrl_int_cnt, MPI_INTEGER, 0, commsander, ier)
  call mpi_bcast(buf, pgm_ctrl_dbl_cnt, MPI_DOUBLE_PRECISION, 0, commsander, ier)

end subroutine pGM_mdin_bcast
#endif
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_read_mdin(nf)
  use file_io_dat

  integer,intent(in) :: nf

  namelist/pol_gauss/pol_gauss_verbose,use_average,dipole_print, &
                     dipole_scf_tol,dipole_scf_init,dipole_scf_init_order,dipole_scf_init_step, &
                     scf_solv_opt,scf_sor_coefficient,scf_sor_niter, &
                     scf_cg_niter,scf_local_niter,scf_local_cut, &
                     ee_dsum_cut,ee_damped_cut, &
                     erfc_tol,ee_gauss_cut, pol_gauss_ips

  read(nf,nml=pol_gauss)

  scf_cg_niter = max(2, scf_cg_niter)
end subroutine POL_GAUSS_read_mdin
!-------------------------------------------------------------------------------
end module pol_gauss_mdin
