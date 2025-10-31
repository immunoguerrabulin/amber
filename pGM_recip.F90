#include "../include/dprec.fh"
#include "../include/assert.fh"

!---------------------------------------------------------
module pol_gauss_recip

  implicit none
  private

  integer, parameter :: Max_Bspline_order=25
  logical, save :: perm_field_done = .FALSE.
  integer, save :: num_ks
  integer, save :: nmine ! number of atoms in each thread, = numatoms in serial
  _REAL_, allocatable, save :: G_func(:,:,:),Qperm(:,:,:), &
                             Q1(:,:,:), &
                             theta1(:,:,:),theta2(:,:,:),theta3(:,:,:), &
                             fractional_multipole(:,:),perm_F_field(:,:)
  integer,allocatable, save :: init_grid_ind(:,:)
#ifdef MPI
  integer, allocatable, save :: my_cg(:) ! atom index for each thread in MPI
#endif
# include "pol_gauss_mpole_index.h"
  public pGM_RECIP_perm_field,pGM_RECIP_dipole_field,pGM_RECIP_ene_frc, &
         pGM_RECIP_allocate,pGM_RECIP_deallocate

  contains
!---------------------------------------------------------
subroutine pGM_RECIP_allocate(numatoms)
  integer,intent(in) :: numatoms

# include "ew_pme_recip.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
  include 'mpif.h'
#endif

  integer :: nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw, &
             dr_order,Bspline_order


  call get_fftdims(nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)

  if (.not. allocated(Qperm)) allocate(Qperm(2*nfftdim1,nfftdim2,nfftdim3))
  if (.not. allocated(Q1)) allocate(Q1(2*nfftdim1,nfftdim2,nfftdim3))
  if (.not. allocated(G_func)) allocate(G_func(nfft3,nfftdim1,nfft2))

#ifdef MPI
  if(num_ks == 0) &
         num_ks = min((((nxyslab(0)+order-1)*numatoms*4)/ &
         (3*nfft3)),numatoms)
#else
  num_ks = numatoms
#endif

  Bspline_order = order
  dr_order = 3
  if ( Bspline_order < dr_order + 2 )then
    write(6,*)'Spline order too small. Must be at least ',dr_order + 2
    call mexit(6,1) 
  endif
  if ( Bspline_order > Max_Bspline_order )then
    write(6,*)'Bspline_order too big! Max = ',Max_Bspline_order
    call mexit(6,1)
  endif

  if (.not. allocated(theta1)) allocate(theta1(0:dr_order,1:Bspline_order,num_ks))
  if (.not. allocated(theta2)) allocate(theta2(0:dr_order,1:Bspline_order,num_ks))
  if (.not. allocated(theta3)) allocate(theta3(0:dr_order,1:Bspline_order,num_ks))
  if (.not. allocated(init_grid_ind)) allocate(init_grid_ind(3,num_ks))
#ifdef MPI
  if (.not. allocated(my_cg)) allocate(my_cg(num_ks))
#endif
  if (.not. allocated(fractional_multipole)) allocate(fractional_multipole(10,numatoms))
  if (.not. allocated(perm_F_field)) allocate(perm_F_field(10,numatoms))
  
end subroutine pGM_RECIP_allocate
!---------------------------------------------------------
subroutine pGM_RECIP_deallocate

   implicit none

  if ( allocated(Qperm)) deallocate(Qperm)
  if ( allocated(Q1)) deallocate(Q1)
  if ( allocated(G_func)) deallocate(G_func)
  if ( allocated(theta1)) deallocate(theta1)
  if ( allocated(theta2)) deallocate(theta2)
  if ( allocated(theta3)) deallocate(theta3)
  if ( allocated(init_grid_ind)) deallocate(init_grid_ind)
#ifdef MPI
  if ( allocated(my_cg)) deallocate(my_cg)
#endif
  if ( allocated(fractional_multipole)) deallocate(fractional_multipole)
  if ( allocated(perm_F_field)) deallocate(perm_F_field)

end subroutine pGM_RECIP_deallocate
!---------------------------------------------------------
subroutine pGM_RECIP_perm_field(numatoms,crd,cart_dipole_field,x)
  use pol_gauss_multipoles, only : global_multipole
  use nblist, only: recip,volume
  use stack
  implicit none

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,*)
  _REAL_,intent(inout) :: cart_dipole_field(3,*)
  _REAL_,intent(in) :: x(*)

# include "ew_pme_recip.h"
# include "def_time.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
#endif

  character(kind=1,len=20) :: routine="pGM_RECIP_perm_field"
  integer :: nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
  integer :: Bspline_order,dr_order,p_tmpy,p_alpha,p_beta,p_fftw,p_FdipF

  Bspline_order = order
  dr_order = 3
  call timer_start(TIME_BSPL)
  call get_fftdims(nfft1,nfft2,nfft3, &
                   nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)

  ! fill theta1-3 and their derivatives
  ! these are saved for use in induction scf & final energy, frc
  call pGM_RECIP_Bspline_fill(numatoms,crd,recip, &
                 nfft1,nfft2,nfft3,Bspline_order,dr_order, &
#ifdef MPI
                 num_ks,my_cg, &
#endif
                 nmine)

  call pGM_RECIP_perm_to_fractional(numatoms,recip,nfft1,nfft2,nfft3,  &
                 global_multipole,fractional_multipole, &
#ifdef MPI
                 my_cg, &
#endif
                 nmine)

  call timer_stop_start(TIME_BSPL,TIME_FILLG)
  call pGM_RECIP_perm_fillgrid(numatoms,Bspline_order,dr_order, &
                 fractional_multipole,nfft1,nfft2,nfft3, &
                 nfftdim1,nfftdim2,nfftdim3,Q1, &
#ifdef MPI
                 my_cg, &
#endif
                 nmine)

  ! we are ready for backward fft of Q1
  call timer_stop_start(TIME_FILLG,TIME_FFT)
  call get_stack(p_tmpy,2*nfftdim1,routine)
  call get_stack(p_alpha,nfft1,routine)
  call get_stack(p_beta,nfft1,routine)
  call get_stack(p_fftw,sffw,routine)
  if (.not. rstack_ok) then
    deallocate(r_stack)
    allocate(r_stack(1:lastrst),stat=alloc_ier)
    call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  call fft_backrc(Q1,x(lfftable),r_stack(p_fftw), &
                  nfft1,nfft2,nfft3,nfftdim1,nfftdim2, &
                  r_stack(p_tmpy),r_stack(p_alpha),r_stack(p_beta))

  ! make a copy of Q1 for final virial, force calculation
  ! compute G function and save for later, in the same time do scalar sum
  call array_copy(Q1,Qperm,siz_q)
  call timer_stop_start(TIME_FFT,TIME_SCSUM)
  call pGM_RECIP_get_G_Q(x(lprefac1),x(lprefac2),x(lprefac3), &
                 recip,volume,ew_coeff,  &
                 nfft1,nfftdim1,nfft2,nfft3,G_func,Q1)

  ! now forward fft of Q1.
  call timer_stop_start(TIME_SCSUM,TIME_FFT)
  call fft_forwardrc(Q1,x(lfftable),r_stack(p_fftw), &
                     nfft1,nfft2,nfft3,nfftdim1,nfftdim2, &
                     r_stack(p_tmpy),r_stack(p_alpha),r_stack(p_beta))
  call free_stack(p_fftw,routine)
  call free_stack(p_beta,routine)
  call free_stack(p_alpha,routine)
  call free_stack(p_tmpy,routine)

  ! finally we are ready for grad sum and field.
  call timer_stop_start(TIME_FFT,TIME_GRADS)
  call get_stack(p_FdipF,3*nmine,routine)
  if (.not. rstack_ok) then
    deallocate(r_stack)
    allocate(r_stack(1:lastrst),stat=alloc_ier)
    call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)

  call pGM_RECIP_get_perm_F_field(numatoms, &
                 dr_order,Bspline_order, &
                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                 Q1,perm_F_field,r_stack(p_FdipF), &
#ifdef MPI
                 my_cg, &
#endif
                 nmine)

  call pGM_RECIP_Fdip_to_Cdip_field(numatoms,nfft1,nfft2,nfft3, &
                 recip,r_stack(p_FdipF),cart_dipole_field,&
#ifdef MPI
                 my_cg, &
#endif
                 nmine)
  call free_stack(p_FdipF,routine)
  
  perm_field_done = .true.
  call timer_stop(TIME_GRADS)

end subroutine pGM_RECIP_perm_field
!---------------------------------------------------------
subroutine pGM_RECIP_dipole_field(numatoms,x,ind_dip,dip_field)
  use nblist, only: recip,volume
  use stack
  implicit none

  character(kind=1,len=22) :: routine="pGM_RECIP_dipole_field"
  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(out) :: dip_field(3,*)

# include "ew_pme_recip.h"
# include "def_time.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
#endif

  integer :: nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
  integer :: Bspline_order,dr_order,p_tmpy,p_alpha,p_beta,p_fftw, &
             p_fdip,p_frac_field

  Bspline_order = order
  dr_order = 3
  call timer_start(TIME_BSPL)
  call get_fftdims(nfft1,nfft2,nfft3, &
                   nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)

  call get_stack(p_fdip,3*numatoms,routine)
  call get_stack(p_frac_field,3*numatoms,routine)
  if(.not. rstack_ok)then
      deallocate(r_stack)
      allocate(r_stack(1:lastrst),stat=alloc_ier)
      call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  call pGM_RECIP_dipole_to_fractional(numatoms,nfft1,nfft2,nfft3, &
                 recip,ind_dip,r_stack(p_fdip),&
#ifdef MPI
                 my_cg,&
#endif
                 nmine)

  call timer_stop_start(TIME_BSPL,TIME_FILLG)
  call pGM_RECIP_dipole_fillgrid(numatoms,Bspline_order,dr_order, &
                 r_stack(p_fdip),  &
                 nfft1,nfft2,nfft3, &
                 nfftdim1,nfftdim2,nfftdim3,Q1,&
#ifdef MPI
                 my_cg,&
#endif
                 nmine)

  call timer_stop_start(TIME_FILLG,TIME_FFT)
  call get_stack(p_tmpy,2*nfftdim1,routine)
  call get_stack(p_alpha,nfft1,routine)
  call get_stack(p_beta,nfft1,routine)
  call get_stack(p_fftw,sffw,routine)
  if (.not. rstack_ok) then
      deallocate(r_stack)
      allocate(r_stack(1:lastrst),stat=alloc_ier)
      call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  call fft_backrc(Q1,x(lfftable),r_stack(p_fftw), &
                 nfft1,nfft2,nfft3, nfftdim1,nfftdim2, &
                 r_stack(p_tmpy),r_stack(p_alpha),r_stack(p_beta))

  call timer_stop_start(TIME_FFT,TIME_SCSUM)
  call pGM_RECIP_G_times_Q(nfft1,nfftdim1,nfft2,nfft3,G_func,Q1)

  call timer_stop_start(TIME_SCSUM,TIME_FFT)
  call fft_forwardrc(Q1,x(lfftable),r_stack(p_fftw), &
                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2, &
                 r_stack(p_tmpy),r_stack(p_alpha),r_stack(p_beta))
  call free_stack(p_fftw,routine)
  call free_stack(p_beta,routine)
  call free_stack(p_alpha,routine)
  call free_stack(p_tmpy,routine)

  call timer_stop_start(TIME_FFT,TIME_GRADS)
  call pGM_RECIP_get_dip_F_field(numatoms, &
                 dr_order,Bspline_order, &
                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                 Q1,r_stack(p_frac_field),&
#ifdef MPI
                 my_cg, &
#endif
                 nmine)

  call pGM_RECIP_Fdip_to_Cdip_field(numatoms,nfft1,nfft2,nfft3, &
                 recip,r_stack(p_frac_field),dip_field,&
#ifdef MPI
                 my_cg, &
#endif
                 nmine)

  call free_stack(p_frac_field,routine)
  call free_stack(p_fdip,routine)
  call timer_stop(TIME_GRADS)

end subroutine pGM_RECIP_dipole_field
!---------------------------------------------------------
subroutine pGM_RECIP_ene_frc(numatoms,crd,x,ind_dip,virial,phi)
  use pol_gauss_multipoles, only : global_multipole
  use pol_gauss_multipoles, only : coulomb_const_kcal_per_mole
  use nblist, only: recip,volume
  use stack

  character(kind=1,len=17) :: routine="pGM_RECIP_ene_frc"
  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,*),x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(inout) :: virial(3,3)
  _REAL_ :: phi(10,numatoms)
  integer :: i, kq

#include "ew_pme_recip.h"
#include "def_time.h"
#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
  integer kbot0,kbot,ktop
#endif
  integer :: nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw,j,k
  integer :: Bspline_order,dr_order,p_tmpy,p_alpha,p_beta,p_fftw, &
             p_fdip,p_frac_field,p_cdf

  if ( .not. perm_field_done )then ! this occurs if pGM_RECIP_perm_field
                                   ! was not called since last call to here
                                   ! i.e. no induced dipoles
                                   ! we need permanent field for ene_perm
    call get_stack(p_cdf,3*numatoms,routine)
    if (.not. rstack_ok) then
      deallocate(r_stack)
      allocate(r_stack(1:lastrst),stat=alloc_ier)
      call reassign_rstack(routine)
    endif
    REQUIRE(rstack_ok)
    call pGM_RECIP_perm_field(numatoms,crd,r_stack(p_cdf),x)
    call free_stack(p_cdf,routine)
  endif
  perm_field_done = .false.  ! reset for next force call

  Bspline_order = order
  dr_order = 3
  call timer_start(TIME_BSPL)
  call get_fftdims(nfft1,nfft2,nfft3, &
                 nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)

  call get_stack(p_fdip,3*numatoms,routine)
  call get_stack(p_frac_field,10*numatoms,routine)
  if (.not. rstack_ok) then
    deallocate(r_stack)
    allocate(r_stack(1:lastrst),stat=alloc_ier)
    call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  call pGM_RECIP_dipole_to_fractional(numatoms,nfft1,nfft2,nfft3, &
                 recip,ind_dip,r_stack(p_fdip),&
#ifdef MPI
                 my_cg,&
#endif
                 nmine)

  call timer_stop_start(TIME_BSPL,TIME_FILLG)
  call pGM_RECIP_dipole_fillgrid(numatoms,Bspline_order,dr_order, &
                 r_stack(p_fdip),  &
                 nfft1,nfft2,nfft3, &
                 nfftdim1,nfftdim2,nfftdim3,Q1,&
#ifdef MPI
                 my_cg,&
#endif
                 nmine)

  call timer_stop_start(TIME_FILLG,TIME_FFT)
  call get_stack(p_tmpy,2*nfftdim1,routine)
  call get_stack(p_alpha,nfft1,routine)
  call get_stack(p_beta,nfft1,routine)
  call get_stack(p_fftw,sffw,routine)
  if (.not. rstack_ok) then
    deallocate(r_stack)
    allocate(r_stack(1:lastrst),stat=alloc_ier)
    call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  call fft_backrc(Q1,x(lfftable),r_stack(p_fftw), &
                 nfft1,nfft2,nfft3, nfftdim1,nfftdim2, &
                 r_stack(p_tmpy),r_stack(p_alpha),r_stack(p_beta))
  call timer_stop_start(TIME_FFT,TIME_SCSUM)

  call pGM_RECIP_scalar_sum(recip,ew_coeff,nfft1,nfftdim1,nfft2,nfft3, &
                 G_func,Qperm,Q1,virial)

  call timer_stop_start(TIME_SCSUM,TIME_FFT)
  call fft_forwardrc(Q1,x(lfftable),r_stack(p_fftw), &
                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2, &
                 r_stack(p_tmpy),r_stack(p_alpha),r_stack(p_beta))

  call free_stack(p_fftw,routine)
  call free_stack(p_beta,routine)
  call free_stack(p_alpha,routine)
  call free_stack(p_tmpy,routine)

  call timer_stop_start(TIME_FFT,TIME_GRADS)

  ! get the field due to the induced dipoles
  ! (up to grad of quadrupole terms i.e. octupolar order)
  call pGM_RECIP_get_ind_F_field(numatoms, &
                 dr_order,Bspline_order, &
                 nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                 Q1,r_stack(p_frac_field),&
#ifdef MPI
                 my_cg,&
#endif
                 nmine)

  ! Haixin: convert to cart field from fractional grid field
  ! And obtain total phi and its derivatives
  call pGM_RECIP_tot_phi_grad(numatoms,nfft1,nfft2,nfft3,recip,perm_F_field,&
                 r_stack(p_frac_field),phi, &
#ifdef MPI
                 my_cg,&
#endif
                 nmine)

  call free_stack(p_frac_field,routine)
  call free_stack(p_fdip,routine)
  call timer_stop(TIME_GRADS)

end subroutine pGM_RECIP_ene_frc
!-------------------------------------------------------------------------------
subroutine pGM_RECIP_perm_fillgrid(numatoms,Bspline_order,dr_order, &
                     Fmpole,nfft1,nfft2,nfft3, &
                     nfftdim1,nfftdim2,nfftdim3,Q, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,Bspline_order,dr_order,nmine
  _REAL_,intent(in) :: Fmpole(10,*)
  integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  _REAL_,intent(out) :: Q(2*nfftdim1,nfftdim2,nfftdim3)
#ifdef MPI
  integer my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"

  integer kbot0,kbot,ktop
#endif
  integer :: im,n,ith1,ith2,ith3,i,i0,i00,j,j0,j00,k,k0,kq
  _REAL_ :: t0,t1,t2,u0,u1,u2,v0,v1,v2,term0,term1,term2

  call zero_array(Q,2*nfftdim1*nfftdim2*nfftdim3)

#ifdef MPI
  kbot0 = mxystart(mytaskid)
  kbot = kbot0 + 1
  ktop = kbot0 + mxyslabs
#endif

  do im = 1, nmine
#ifdef MPI
    n = my_cg(im)
#else
    n = im
#endif
    ! begin index for each direction
    i00 = int(init_grid_ind(1,im))
    j00 = int(init_grid_ind(2,im))
    k0  = int(init_grid_ind(3,im))

    do ith3 = 1, Bspline_order
      k0 = k0 + 1
      k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
         
#ifdef MPI
      if ( k >= kbot .and. k <= ktop ) then
        kq = k - kbot0
#else
        kq = k
#endif

        ! everytime jump iqk numbers in array
        ! (as there are imag part in the array),
        ! to  fill in the real parts of the charge only
            
        v0 = theta3(0,ith3,im) !theta3
        v1 = theta3(1,ith3,im) !1st deriv of theta3

        j0 = j00
        do ith2 = 1,Bspline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2

          u0 = theta2(0,ith2,im) !theta2
          u1 = theta2(1,ith2,im) !1st deriv of theta2
 
          ! hardwire our knowledge of layout of theta1,2,3 to pre-assemble
          ! factors
          term0 = Fmpole(Ind_000,im)*u0*v0 + &
                  Fmpole(Ind_010,im)*u1*v0 + &
                  Fmpole(Ind_001,im)*u0*v1
          term1 = Fmpole(Ind_100,im)*u0*v0
          
          i0 = i00
          do ith1 = 1, Bspline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2

            t0 = theta1(0,ith1,im) !theta1
            t1 = theta1(1,ith1,im) !1st deriv of theta1

            Q(i,j,kq) = Q(i,j,kq) + term0*t0 + term1*t1
          end do ! ith1 = 1, Bspline_order
        end do ! ith2 = 1, Bspline_order
#ifdef MPI
      end if
#endif
    end do ! ith3 = 1, Bspline_order
  end do !  im = 1, nmine

end subroutine pGM_RECIP_perm_fillgrid
!---------------------------------------------------------
subroutine pGM_RECIP_dipole_fillgrid(numatoms,Bspline_order,dr_order, &
                     fdip,  &
                     nfft1,nfft2,nfft3, &
                     nfftdim1,nfftdim2,nfftdim3,Q1, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,Bspline_order,dr_order, nmine
  _REAL_,intent(in) :: fdip(3,*)
  integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  _REAL_,intent(out) :: Q1(2*nfftdim1,nfftdim2,nfftdim3)
#ifdef MPI
  integer my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"

  integer kbot0,kbot,ktop
#endif
  integer :: n,igrd0,jgrd0,kgrd0,ith1,ith2,ith3,i,i0,j,j0,k,k0,ntot,im,kq
  integer :: i00, j00
  _REAL_ :: t0,t1,u0,u1,v0,v1,term1_0,term0,term1,term2

  ntot = 2*nfftdim1*nfftdim2*nfftdim3
  call zero_array(Q1,ntot)

#ifdef MPI
  kbot0 = mxystart(mytaskid)
  kbot = kbot0 + 1
  ktop = kbot0 + mxyslabs
#endif

  do im = 1, nmine
#ifdef MPI
    n = my_cg(im)
#else
    n = im
#endif
    ! begin index for each direction
    i00 = int(init_grid_ind(1,im))
    j00 = int(init_grid_ind(2,im))
    k0  = int(init_grid_ind(3,im))

    do ith3 = 1, Bspline_order
      k0 = k0 + 1

      k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
         
#ifdef MPI
      if ( k >= kbot .and. k <= ktop ) then
        kq = k - kbot0
#else
        kq = k
#endif

        ! everytime jump iqk numbers in array
        ! (as there are imag part in the array),
        ! to  fill in the real parts of the charge only
            
        v0 = theta3(0,ith3,im) !theta3
        v1 = theta3(1,ith3,im) !1st deriv of theta3

        j0 = j00
        do ith2 = 1,Bspline_order
          j0 = j0 + 1

          j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2

          u0 = theta2(0,ith2,im) !theta2
          u1 = theta2(1,ith2,im) !1st deriv of theta2
 
          ! hardwire our knowledge of layout of theta1,2,3 to pre-assemble
          ! factors
          term0 = fdip(3,im)*u0*v1 + &
                  fdip(2,im)*u1*v0 
          term1 = fdip(1,im)*u0*v0
          
          i0 = i00
          do ith1 = 1, Bspline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2

            t0 = theta1(0,ith1,im) !theta1
            t1 = theta1(1,ith1,im) !1st deriv of theta1

            Q1(i,j,kq) = Q1(i,j,kq) + term0*t0 + term1*t1
          end do ! ith1 = 1, Bspline_order
        end do ! ith2 = 1, Bspline_order
#ifdef MPI
      end if
#endif
    end do ! ith3 = 1, Bspline_order
  end do !  im = 1, nmine

end subroutine pGM_RECIP_dipole_fillgrid
!---------------------------------------------------------
subroutine pGM_RECIP_Bspline_fill(numatoms,crd,recip, &
                     nfft1,nfft2,nfft3,Bspline_order,dr_order, &
#ifdef MPI
                     num_ks,my_cg, &
#endif
                     nmine)
  implicit none               

  integer,intent(in)  :: numatoms,nfft1,nfft2,nfft3,Bspline_order,dr_order
  _REAL_,intent(in)   :: crd(3,*),recip(3,3)
  integer,intent(out) :: nmine
#ifdef MPI
  integer num_ks,my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"
#endif

  integer n,imine
  _REAL_ fr3n,fr2n,fr1n,w,w1,w2,w3
  _REAL_ anint
#ifdef MPI
  integer kbot0,kbot,ktop1,ktop,ido,k00
#endif

  ! this is for generating a list, my_cg, to save the atoms needed to generate the grid at current thread
#ifdef MPI
  kbot0 = mxystart(mytaskid)
  kbot = kbot0 + 1
  ktop = kbot0 + mxyslabs
  ktop1 = ktop + Bspline_order - 2
  imine = 0
   
  do n = 1,numatoms
    w = crd(1,n)*recip(1,3) &
          +crd(2,n)*recip(2,3)+crd(3,n)*recip(3,3)
    fr3n = nfft3*(w - (anint(w) - 0.5d0))
    k00 = int(fr3n)
    ! code for filtering atoms. In single proc mode, do all atoms
    ido = 0
    if ( ktop1 >= nfft3)then
       if ( k00 >= kbot0 .or. k00 <= ktop1 - nfft3 )ido = 1
    else
       if ( k00 >= kbot0 .and. k00 <= ktop1)ido = 1
    end if
    if ( ido == 1)then
       imine = imine + 1
       !           ----- do not fill my_cg past num_ks ------
       if( imine <= num_ks) my_cg(imine)=n
    end if
  end do
  nmine=imine
  !   ----- ERROR condition met --------------------
  if(imine > num_ks) return
  !   ----------------------------------------------
#else
  nmine=numatoms
#endif

  imine = 0
  do imine = 1,nmine
#ifdef MPI
    n=my_cg(imine)
#else
    n=imine
#endif
    w = crd(1,n)*recip(1,3) &
          +crd(2,n)*recip(2,3)+crd(3,n)*recip(3,3)
    fr3n = nfft3*(w - (anint(w) - 0.5d0))
    w = crd(1,n)*recip(1,1) &
          +crd(2,n)*recip(2,1)+crd(3,n)*recip(3,1)
    fr1n = nfft1*(w - (anint(w) - 0.5d0))
    w = crd(1,n)*recip(1,2) &
          +crd(2,n)*recip(2,2)+crd(3,n)*recip(3,2)
    fr2n = nfft2*(w - (anint(w) - 0.5d0))
    w1 = fr1n-int(fr1n)
    w2 = fr2n-int(fr2n)
    w3 = fr3n-int(fr3n)

    init_grid_ind(1,imine) = int(fr1n) - Bspline_order
    init_grid_ind(2,imine) = int(fr2n) - Bspline_order
    init_grid_ind(3,imine) = int(fr3n) - Bspline_order
    call pGM_RECIP_bspline_fill_gen1(w1,Bspline_order,dr_order, &
                   theta1(0,1,imine))

    call pGM_RECIP_bspline_fill_gen1(w2,Bspline_order,dr_order, &
                   theta2(0,1,imine))

    call pGM_RECIP_bspline_fill_gen1(w3,Bspline_order,dr_order, &
                   theta3(0,1,imine))
  end do

end subroutine pGM_RECIP_Bspline_fill
!---------------------------------------------------------
! ZH/UCI: Replace previous bspline_fill_gen by adding nderiv=3 for pGM
!---------------------------------------------------------
subroutine pGM_RECIP_bspline_fill_gen1(w,spline_order,dr_order,array)
  use ew_bspline, only : fill_bspline_0,fill_bspline_1,fill_bspline_2,fill_bspline_3
  implicit none

  integer,intent(in) :: spline_order, dr_order
  _REAL_,intent(in) :: w
  _REAL_,intent(out) :: array(0:dr_order,spline_order)

  integer :: k
  _REAL_ :: theta(spline_order)
  _REAL_ :: dtheta(spline_order)
  _REAL_ :: d2theta(spline_order)
  _REAL_ :: d3theta(spline_order)

  if      ( dr_order == 0 )then
    call fill_bspline_0(w,spline_order,theta)
  else if ( dr_order == 1 )then
    call fill_bspline_1(w,spline_order,theta,dtheta)
  else if ( dr_order == 2 )then
    call fill_bspline_2(w,spline_order,theta,dtheta,d2theta)
  else if ( dr_order == 3 )then ! ZH: added nderiv==3 for pGM
    call fill_bspline_3(w,spline_order,theta,dtheta,d2theta,d3theta)
  end if
  do k = 1, spline_order
    array(0,k) = theta(k)
    array(1,k) = dtheta(k)
    array(2,k) = d2theta(k)
    array(3,k) = d3theta(k)
  end do

end subroutine pGM_RECIP_bspline_fill_gen1
!---------------------------------------------------------
subroutine pGM_RECIP_get_G_Q( prefac1,prefac2,prefac3,  &
                     recip,volume,ewald_coeff,  &
                     nfft1,nfftdim1,nfft2,nfft3,G,Q)
  use constants, only : pi
  implicit none

  integer,intent(in) :: nfft1,nfftdim1,nfft2,nfft3
  _REAL_,intent(in) :: prefac1(nfft1),prefac2(nfft2),prefac3(nfft3)
  _REAL_,intent(in) :: recip(3,3),volume,ewald_coeff
  _REAL_,intent(out) :: G(nfft3,nfftdim1,nfft2)
  _REAL_,intent(out) :: Q(2,nfft3,nfftdim1,nfft2)

# include "extra.h"
#ifdef MPI
# include "ew_parallel.h"
# include "parallel.h"
#endif

  integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3, k1q, k2q, k3q
  _REAL_  :: fac,mhat1,mhat2,mhat3,msq,denom

  fac = pi**2/ewald_coeff**2
  nf1 = nfft1/2
  if ( 2*nf1 < nfft1 )nf1 = nf1+1
  nf2 = nfft2/2
  if ( 2*nf2 < nfft2 )nf2 = nf2+1
  nf3 = nfft3/2
  if ( 2*nf3 < nfft3 )nf3 = nf3+1

  ! Insist that Q(1,1,1,1) is set to 0 (true already for neutral)
  if (master) then
    Q(1,1,1,1) = 0.d0
    Q(2,1,1,1) = 0.d0
  end if

#ifdef MPI
  do k2q = 1, mxzslabs
    if (master) then
      k2 = k2q
    else
      k2 = k2q + mxzstart(mytaskid)
    end if
#else
  do k2q = 1, nfft2
    k2 = k2q
#endif
    m2 = k2 - 1
    if ( k2 > nf2 )m2 = k2 - 1 - nfft2

    do k3 = 1,nfft3
      m3 = k3 - 1
      if ( k3 > nf3 )m3 = k3 - 1 - nfft3
      k10 = 1
      if(master)then
        if (k3+k2 == 2) k10 = 2
      end if
      do k1 = k10, nf1+1
        m1 = k1 - 1
        if ( k1 > nf1 )m1 = k1 - 1 - nfft1
        mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
        mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
        mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
        msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
        denom = 1.0d0/(pi*volume*msq)
        G(k3,k1,k2q) = exp(-fac*msq)*prefac1(k1)* &
                prefac2(k2)*prefac3(k3)*denom
        Q(1,k3,k1,k2q) = G(k3,k1,k2q)*Q(1,k3,k1,k2q)
        Q(2,k3,k1,k2q) = G(k3,k1,k2q)*Q(2,k3,k1,k2q)
      end do !k1 = k10, nf1+1
    end do !k3 = 1, nfft3
  end do !k2 = 1, nfft2

end subroutine pGM_RECIP_get_G_Q
!--------------------------------------------------------------
subroutine pGM_RECIP_G_times_Q(nfft1,nfftdim1,nfft2,nfft3,G,Q)
  implicit none

  integer,intent(in) :: nfft1,nfftdim1,nfft2,nfft3
  _REAL_,intent(in) ::  G(nfft3,nfftdim1,nfft2)
  _REAL_,intent(out) :: Q(2,nfft3,nfftdim1,nfft2)

# include "extra.h"
#ifdef MPI
# include "ew_parallel.h"
# include "parallel.h"
#endif

  integer :: k1,k2,k3,k10,nf1,k2q

  nf1 = nfft1/2
  if ( 2*nf1 < nfft1 )nf1 = nf1+1
  ! why here so simple?

  ! Insist that Q(1,1,1,1) is set to 0 (true already for neutral)
  if (master) then
    Q(1,1,1,1) = 0.d0
    Q(2,1,1,1) = 0.d0
  end if

#ifdef MPI
  do k2q = 1, mxzslabs
    if (master) then
      k2 = k2q
    else
      k2 = k2q + mxzstart(mytaskid)
    end if
#else
  do k2q = 1, nfft2
    k2 = k2q
#endif
    do k3 = 1,nfft3
      k10 = 1
      if(master)then
        if (k3+k2 == 2) k10 = 2
      end if
      do k1 = k10, nf1+1
        Q(1,k3,k1,k2q) = G(k3,k1,k2q)*Q(1,k3,k1,k2q)
        Q(2,k3,k1,k2q) = G(k3,k1,k2q)*Q(2,k3,k1,k2q)
      end do !k1 = k10, nf1+1
    end do !k3 = 1, nfft3
  end do !k2 = 1, nfft2

end subroutine pGM_RECIP_G_times_Q
!--------------------------------------------------------------
subroutine pGM_RECIP_scalar_sum(recip,ewald_coeff, &
                     nfft1,nfftdim1,nfft2,nfft3, &
                     G,Q,Q1,virial)
  use constants, only : pi
  implicit none

  _REAL_,intent(in) :: recip(3,3),ewald_coeff
  integer,intent(in) :: nfft1,nfftdim1,nfft2,nfft3
  _REAL_,intent(in) :: G(nfft3,nfftdim1,nfft2)
  _REAL_,intent(in) :: Q(2,nfft3,nfftdim1,nfft2)
  _REAL_,intent(inout) :: Q1(2,nfft3,nfftdim1,nfft2)
  _REAL_,intent(out) :: virial(3,3)

# include "extra.h"
#ifdef MPI
# include "ew_parallel.h"
# include "parallel.h"
#endif

  integer :: k1,k2,k3,m1,m2,m3,k10,nf1,nf2,nf3, k1q, k2q, k3q
  _REAL_  :: fac,mhat1,mhat2,mhat3,msq,struc2,eterm,vterm, &
             q11,q12,q21,q22,tmp1,tmp2,mult,vxx,vxy,vxz,vyy,vyz,vzz

  fac = pi**2/ewald_coeff**2
  nf1 = nfft1/2
  if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
  nf2 = nfft2/2
  if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
  nf3 = nfft3/2
  if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

#ifdef MPI
  do k2q = 1, mxzslabs
    if ( master ) then
      k2 = k2q
    else
      k2 = k2q + mxzstart(mytaskid)
    end if
#else
  do k2q = 1, nfft2
    k2 = k2q
#endif
    m2 = k2 - 1
    if ( k2 > nf2 )m2 = k2 - 1 - nfft2

    do k3 = 1,nfft3
      m3 = k3 - 1
      if ( k3 > nf3 )m3 = k3 - 1 - nfft3
      k10 = 1
      if ( master ) then
        if (k3+k2 == 2) k10 = 2
      end if
      do k1 = k10, nf1+1
        if ( k1 > 1 ) then
          mult = 2.d0
        else
          mult = 1.d0
        end if
        m1 = k1 - 1
        if ( k1 > nf1 )m1 = k1 - 1 - nfft1
        mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
        mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
        mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
        msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
        eterm = mult*G(k3,k1,k2q)
        vterm = 2.d0*(fac*msq + 1.d0)/msq
        q11 = Q(1,k3,k1,k2q) + Q1(1,k3,k1,k2q)
        q12 = Q(2,k3,k1,k2q) + Q1(2,k3,k1,k2q)
        struc2 = q11**2+q12**2
        tmp1 = eterm*struc2
        tmp2 = tmp1*vterm
        vxx = vxx + tmp2*mhat1*mhat1 - tmp1
        vxy = vxy + tmp2*mhat1*mhat2
        vxz = vxz + tmp2*mhat1*mhat3
        vyy = vyy + tmp2*mhat2*mhat2 - tmp1
        vyz = vyz + tmp2*mhat2*mhat3
        vzz = vzz + tmp2*mhat3*mhat3 - tmp1

        Q1(1,k3,k1,k2q) = G(k3,k1,k2q)*Q1(1,k3,k1,k2q)
        Q1(2,k3,k1,k2q) = G(k3,k1,k2q)*Q1(2,k3,k1,k2q)
      end do !k1 = k10, nf1+1
    end do !k3 = 1, nfft3
  end do !k2 = 1, nfft2

  virial(1,1) = virial(1,1) + 0.5d0*vxx
  virial(1,2) = virial(1,2) + 0.5d0*vxy
  virial(2,1) = virial(2,1) + 0.5d0*vxy
  virial(1,3) = virial(1,3) + 0.5d0*vxz
  virial(3,1) = virial(3,1) + 0.5d0*vxz
  virial(2,2) = virial(2,2) + 0.5d0*vyy
  virial(2,3) = virial(2,3) + 0.5d0*vyz
  virial(3,2) = virial(3,2) + 0.5d0*vyz
  virial(3,3) = virial(3,3) + 0.5d0*vzz

end subroutine pGM_RECIP_scalar_sum
!--------------------------------------------------------------
subroutine pGM_RECIP_get_perm_F_field(numatoms, &
                     dr_order,Bspline_order, &
                     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                     Q_p,F_perm_field,Fperm_field, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,dr_order,Bspline_order
  integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  _REAL_,intent(in) :: Q_p(2*nfftdim1,nfftdim2,nfftdim3)
  _REAL_,intent(out) :: F_perm_field(10,*),Fperm_field(3,*)
  integer :: nmine
#ifdef MPI
  integer my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"

  integer kbot0,kbot,ktop
#endif
! Note order of field elements
! 1 Potential F corresponds to 000 in spline indices
! 2-4 -> F_x, F_y, F_z respectively or to 100,010,001 in indices
! 5-10 -> F_xx,F_yy,F_zz,F_xy,F_xz, F_yz resp. or to 
!            200, 020, 002, 110, 101, 011 in indices
! RL: The following terms are deleted
! 11-20 -> F_xxx,F_yyy,F_zzz,F_xxy,F_xxz,F_xyy,F_yyz,F_xzz,F_yzz,F_xyz or to
!          300,030,003,210,201,120,021,102,012,111  in indices
  integer :: m,n,i,i0,j,j0,k,k0,igrd0,jgrd0,kgrd0,ith1,ith2,ith3,im,kq
  _REAL_ :: tq_p,t_p(0:3),u(0:3),v(0:3)
  _REAL_ :: tu_p(6),tuv_p(10)

#ifdef MPI
  kbot0 = mxystart(mytaskid)
  kbot = kbot0 + 1
  ktop = kbot0 + mxyslabs
#endif

  do im = 1, nmine
#ifdef MPI
    n = my_cg(im)
#else
    n = im
#endif

    igrd0 = init_grid_ind(1,im) !begin index in 1st direction
    jgrd0 = init_grid_ind(2,im) !begin index in 2nd direction
    kgrd0 = init_grid_ind(3,im) !begin index in 3rd direction
    do m = 1, 10
      tuv_p(m) = 0.d0
    end do
    k0 = kgrd0
    do ith3 = 1, Bspline_order
      do m = 1, 6
        tu_p(m) = 0.d0
      end do
      do m = 0, 2
        v(m) = theta3(m,ith3,im)
      end do
      k0 = k0 + 1
      k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
#ifdef MPI
      if ( k >= kbot .and. k <= ktop ) then
        kq = k - kbot0
#else
        kq = k
#endif
       
      j0 = jgrd0
      do ith2 = 1, Bspline_order
        j0 = j0 + 1
        j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
        i0 = igrd0
        do m = 0, 2
          u(m) = theta2(m,ith2,im)
          t_p(m) = 0.d0
        end do
        do ith1 = 1, Bspline_order
          i0 = i0 + 1
          i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
          tq_p = Q_p(i,j,kq)
          t_p(0) = t_p(0) + tq_p*theta1(0,ith1,im)
          t_p(1) = t_p(1) + tq_p*theta1(1,ith1,im)
          t_p(2) = t_p(2) + tq_p*theta1(2,ith1,im)
        end do !ith1 = 1, Bspline_order
        tu_p(Ind_00) = tu_p(Ind_00) + t_p(0)*u(0)
        tu_p(Ind_10) = tu_p(Ind_10) + t_p(1)*u(0)
        tu_p(Ind_01) = tu_p(Ind_01) + t_p(0)*u(1)
        tu_p(Ind_20) = tu_p(Ind_20) + t_p(2)*u(0)
        tu_p(Ind_02) = tu_p(Ind_02) + t_p(0)*u(2)
        tu_p(Ind_11) = tu_p(Ind_11) + t_p(1)*u(1)
      end do !ith2 = 1, Spline_order
      tuv_p(Ind_000) = tuv_p(Ind_000) + tu_p(Ind_00)*v(0)
      tuv_p(Ind_100) = tuv_p(Ind_100) + tu_p(Ind_10)*v(0)
      tuv_p(Ind_010) = tuv_p(Ind_010) + tu_p(Ind_01)*v(0)
      tuv_p(Ind_001) = tuv_p(Ind_001) + tu_p(Ind_00)*v(1)
      tuv_p(Ind_200) = tuv_p(Ind_200) + tu_p(Ind_20)*v(0)
      tuv_p(Ind_020) = tuv_p(Ind_020) + tu_p(Ind_02)*v(0)
      tuv_p(Ind_002) = tuv_p(Ind_002) + tu_p(Ind_00)*v(2)
      tuv_p(Ind_110) = tuv_p(Ind_110) + tu_p(Ind_11)*v(0)
      tuv_p(Ind_101) = tuv_p(Ind_101) + tu_p(Ind_10)*v(1)
      tuv_p(Ind_011) = tuv_p(Ind_011) + tu_p(Ind_01)*v(1)
#ifdef MPI
      end if
#endif
    end do !ith3 = 1, Spline_order
    do m = 1, 10
      F_perm_field(m,im) = tuv_p(m)
    end do
    Fperm_field(1,im) = tuv_p(Ind_100)
    Fperm_field(2,im) = tuv_p(Ind_010)
    Fperm_field(3,im) = tuv_p(Ind_001)
  end do !im = 1, nmine

end subroutine pGM_RECIP_get_perm_F_field
!--------------------------------------------------------------
! RL: same as pGM_RECIP_get_ind_F_field except it is for induction cycle
subroutine pGM_RECIP_get_dip_F_field(numatoms, &
                     dr_order,Bspline_order, &
                     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q1, &
                     F_dip_field, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)                     
  implicit none

  integer,intent(in) :: numatoms,dr_order,Bspline_order
  integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  _REAL_,intent(in) :: Q1(2*nfftdim1,nfftdim2,nfftdim3)
  _REAL_,intent(out) :: F_dip_field(3,*)
#ifdef MPI
  integer,intent(in) :: my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"

  integer :: kbot0,kbot,ktop
#endif
  integer :: n,i,i0,j,j0,k,k0,igrd0,jgrd0,kgrd0,ith1,ith2,ith3, nmine,kq,im
  _REAL_ :: tq1,t0_1,t1_1,u0,u1,v0,v1,tu00_1,tu10_1,tu01_1,&
            tuv001_1,tuv010_1,tuv100_1

#ifdef MPI
  kbot0 = mxystart(mytaskid)
  kbot = kbot0 + 1
  ktop = kbot0 + mxyslabs
#endif

  do im = 1, nmine
#ifdef MPI
    n = my_cg(im)
#else
    n = im
#endif

    igrd0 = init_grid_ind(1,im) !begin index in 1st direction
    jgrd0 = init_grid_ind(2,im) !begin index in 2nd direction
    kgrd0 = init_grid_ind(3,im) !begin index in 3rd direction

    tuv001_1 = 0.d0
    tuv010_1 = 0.d0
    tuv100_1 = 0.d0
    k0 = kgrd0
    do ith3 = 1,Bspline_order
      v0 = theta3(0,ith3,im) !theta3
      v1 = theta3(1,ith3,im) !1st deriv of theta3
      tu00_1 = 0.d0
      tu10_1 = 0.d0
      tu01_1 = 0.d0
      k0 = k0 + 1
      k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
#ifdef MPI
      if ( k >= kbot .and. k <= ktop ) then
        kq = k - kbot0
#else
        kq = k
#endif
      j0 = jgrd0
      do ith2 = 1,Bspline_order
        j0 = j0 + 1
        j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
        i0 = igrd0
        u0 = theta2(0,ith2,im) !theta2
        u1 = theta2(1,ith2,im) !1st deriv of theta2
        t0_1 = 0.d0
        t1_1 = 0.d0
        do ith1 = 1,Bspline_order
          i0 = i0 + 1
          i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
          tq1 = Q1(i,j,kq)
          t0_1 = t0_1 + tq1*theta1(0,ith1,im)
          t1_1 = t1_1 + tq1*theta1(1,ith1,im)
        end do !ith1 = 1,Bspline_order
        tu00_1 = tu00_1 + t0_1*u0
        tu10_1 = tu10_1 + t1_1*u0
        tu01_1 = tu01_1 + t0_1*u1
      end do !ith2 = 1,Spline_order
      tuv100_1 = tuv100_1 + tu10_1*v0
      tuv010_1 = tuv010_1 + tu01_1*v0
      tuv001_1 = tuv001_1 + tu00_1*v1
#ifdef MPI
      end if
#endif
    end do !ith3 = 1,Spline_order
    F_dip_field(1,im) = tuv100_1
    F_dip_field(2,im) = tuv010_1
    F_dip_field(3,im) = tuv001_1
  end do !n = 1,numatoms

end subroutine pGM_RECIP_get_dip_F_field
!--------------------------------------------------------------
! RL: same as pGM_RECIP_get_dip_F_field except it is for final ene/frc
subroutine pGM_RECIP_get_ind_F_field(numatoms, &
                     dr_order,Bspline_order, &
                     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
                     Q1,F_dip_field, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,dr_order,Bspline_order,nmine
  integer,intent(in) :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
  _REAL_,intent(in) :: Q1(2*nfftdim1,nfftdim2,nfftdim3)
  _REAL_,intent(out) :: F_dip_field(10,*)
#ifdef MPI
  integer my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"

  integer kbot0,kbot,ktop
#endif
! Note order of field elements
! 1 Potential F corresponds to 000 in spline indices
! 2-4 -> F_x, F_y, F_z respectively or to 100,010,001 in indices
! 5-10 -> F_xx,F_yy,F_zz,F_xy,F_xz, F_yz resp. or to 
!            200, 020, 002, 110, 101, 011 in indices
! RL: The following terms are deleted
! 11-20 -> F_xxx,F_yyy,F_zzz,F_xxy,F_xxz,F_xyy,F_yyz,F_xzz,F_yzz,F_xyz or to
!          300,030,003,210,201,120,021,102,012,111  in indices
  integer :: im, kq
  integer :: m,n,i,i0,j,j0,k,k0,igrd0,jgrd0,kgrd0,ith1,ith2,ith3
  _REAL_ :: tq_1,t_1(0:3),u(0:3),v(0:3)
  _REAL_ :: tu_1(6),tuv_1(10)

#ifdef MPI
  kbot0 = mxystart(mytaskid)
  kbot = kbot0 + 1
  ktop = kbot0 + mxyslabs
#endif

  do im = 1, nmine
#ifdef MPI
    n = my_cg(im)
#else
    n = im
#endif
     igrd0 = init_grid_ind(1,im) !begin index in 1st direction
     jgrd0 = init_grid_ind(2,im) !begin index in 2nd direction
     kgrd0 = init_grid_ind(3,im) !begin index in 3rd direction
     do m = 1, 10
       tuv_1(m) = 0.d0
     end do
     k0 = kgrd0
     do ith3 = 1, Bspline_order
       do m = 1, 6
         tu_1(m) = 0.d0
       end do
       do m = 0, 2
         v(m) = theta3(m,ith3,im)
       end do
       k0 = k0 + 1
       k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
#ifdef MPI
      if ( k >= kbot .and. k <= ktop ) then
        kq = k - kbot0
#else
        kq = k
#endif
       j0 = jgrd0
       do ith2 = 1, Bspline_order
         j0 = j0 + 1
         j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
         i0 = igrd0
         do m = 0, 2
           u(m) = theta2(m,ith2,im)
           t_1(m) = 0.d0
         end do
         do ith1 = 1, Bspline_order
           i0 = i0 + 1
           i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
           tq_1 = Q1(i,j,kq)
           t_1(0) = t_1(0) + tq_1*theta1(0,ith1,im)
           t_1(1) = t_1(1) + tq_1*theta1(1,ith1,im)
           t_1(2) = t_1(2) + tq_1*theta1(2,ith1,im)
         end do !ith1 = 1, Bspline_order
         tu_1(Ind_00) = tu_1(Ind_00) + t_1(0)*u(0)
         tu_1(Ind_10) = tu_1(Ind_10) + t_1(1)*u(0)
         tu_1(Ind_01) = tu_1(Ind_01) + t_1(0)*u(1)
         tu_1(Ind_20) = tu_1(Ind_20) + t_1(2)*u(0)
         tu_1(Ind_02) = tu_1(Ind_02) + t_1(0)*u(2)
         tu_1(Ind_11) = tu_1(Ind_11) + t_1(1)*u(1)
       end do !ith2 = 1, Spline_order
       tuv_1(Ind_000) = tuv_1(Ind_000) + tu_1(Ind_00)*v(0)
       tuv_1(Ind_100) = tuv_1(Ind_100) + tu_1(Ind_10)*v(0)
       tuv_1(Ind_010) = tuv_1(Ind_010) + tu_1(Ind_01)*v(0)
       tuv_1(Ind_001) = tuv_1(Ind_001) + tu_1(Ind_00)*v(1)
       tuv_1(Ind_200) = tuv_1(Ind_200) + tu_1(Ind_20)*v(0)
       tuv_1(Ind_020) = tuv_1(Ind_020) + tu_1(Ind_02)*v(0)
       tuv_1(Ind_002) = tuv_1(Ind_002) + tu_1(Ind_00)*v(2)
       tuv_1(Ind_110) = tuv_1(Ind_110) + tu_1(Ind_11)*v(0)
       tuv_1(Ind_101) = tuv_1(Ind_101) + tu_1(Ind_10)*v(1)
       tuv_1(Ind_011) = tuv_1(Ind_011) + tu_1(Ind_01)*v(1)
#ifdef MPI
      end if
#endif
     end do !ith3 = 1, Spline_order
     do m = 1, 10
       F_dip_field(m,im) = tuv_1(m)
     end do
  end do !n = 1,numatoms

end subroutine pGM_RECIP_get_ind_F_field
!--------------------------------------------------------------
subroutine pGM_RECIP_tot_phi_grad(numatoms,nfft1,nfft2,nfft3,recip, &
                     F_perm_field,F_dip_field,phi, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,nfft1,nfft2,nfft3
  _REAL_,intent(in) :: recip(3,3)
  _REAL_,intent(in) :: F_perm_field(10,*)
  _REAL_,intent(inout) :: F_dip_field(10,*)
  _REAL_,intent(out) :: phi(10,numatoms)
#ifdef MPI
  integer my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"
#endif

  _REAL_ :: f1,f2,f3,dfx,dfy,dfz
  _REAL_ :: mpole_xform_3x3(3,3),field_xform_3x3(3,3)
  integer :: k,n,j1,j2,j3, im, nmine

  ! obtain total potentials and derivatives
  F_dip_field(1:10,1:nmine) = F_dip_field(1:10,1:nmine) + F_perm_field(1:10,1:nmine)

  call pGM_RECIP_xform_matrices(nfft1,nfft2,nfft3,recip, &
                 mpole_xform_3x3,field_xform_3x3)

  do im = 1, nmine
#ifdef MPI
    n = my_cg(im)
#else
    n = im
#endif
    phi(1,n) = F_dip_field(1,im)

    phi(2,n) = field_xform_3x3(1,1)*F_dip_field(2,im) +&
               field_xform_3x3(1,2)*F_dip_field(3,im) +&
               field_xform_3x3(1,3)*F_dip_field(4,im)
    phi(3,n) = field_xform_3x3(2,1)*F_dip_field(2,im) +&
               field_xform_3x3(2,2)*F_dip_field(3,im) +&
               field_xform_3x3(2,3)*F_dip_field(4,im)
    phi(4,n) = field_xform_3x3(3,1)*F_dip_field(2,im) +&
               field_xform_3x3(3,2)*F_dip_field(3,im) +&
               field_xform_3x3(3,3)*F_dip_field(4,im)

    phi(5,n) = field_xform_3x3(1,1)*field_xform_3x3(1,1)*F_dip_field(Ind_200,im) +&
               field_xform_3x3(1,1)*field_xform_3x3(1,2)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(1,1)*field_xform_3x3(1,3)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(1,1)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(1,2)*F_dip_field(Ind_020,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(1,3)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(1,1)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(1,2)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(1,3)*F_dip_field(Ind_002,im)

    phi(6,n) = field_xform_3x3(2,1)*field_xform_3x3(2,1)*F_dip_field(Ind_200,im) +&
               field_xform_3x3(2,1)*field_xform_3x3(2,2)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(2,1)*field_xform_3x3(2,3)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(2,2)*field_xform_3x3(2,1)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(2,2)*field_xform_3x3(2,2)*F_dip_field(Ind_020,im) +&
               field_xform_3x3(2,2)*field_xform_3x3(2,3)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(2,3)*field_xform_3x3(2,1)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(2,3)*field_xform_3x3(2,2)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(2,3)*field_xform_3x3(2,3)*F_dip_field(Ind_002,im)

    phi(7,n) = field_xform_3x3(3,1)*field_xform_3x3(3,1)*F_dip_field(Ind_200,im) +&
               field_xform_3x3(3,1)*field_xform_3x3(3,2)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(3,1)*field_xform_3x3(3,3)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(3,2)*field_xform_3x3(3,1)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(3,2)*field_xform_3x3(3,2)*F_dip_field(Ind_020,im) +&
               field_xform_3x3(3,2)*field_xform_3x3(3,3)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(3,3)*field_xform_3x3(3,1)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(3,3)*field_xform_3x3(3,2)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(3,3)*field_xform_3x3(3,3)*F_dip_field(Ind_002,im)

    phi(8,n) = field_xform_3x3(1,1)*field_xform_3x3(2,1)*F_dip_field(Ind_200,im) +&
               field_xform_3x3(1,1)*field_xform_3x3(2,2)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(1,1)*field_xform_3x3(2,3)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(2,1)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(2,2)*F_dip_field(Ind_020,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(2,3)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(2,1)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(2,2)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(2,3)*F_dip_field(Ind_002,im)

    phi(9,n) = field_xform_3x3(1,1)*field_xform_3x3(3,1)*F_dip_field(Ind_200,im) +&
               field_xform_3x3(1,1)*field_xform_3x3(3,2)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(1,1)*field_xform_3x3(3,3)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(3,1)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(3,2)*F_dip_field(Ind_020,im) +&
               field_xform_3x3(1,2)*field_xform_3x3(3,3)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(3,1)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(3,2)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(1,3)*field_xform_3x3(3,3)*F_dip_field(Ind_002,im)

    phi(10,n)= field_xform_3x3(2,1)*field_xform_3x3(3,1)*F_dip_field(Ind_200,im) +&
               field_xform_3x3(2,1)*field_xform_3x3(3,2)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(2,1)*field_xform_3x3(3,3)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(2,2)*field_xform_3x3(3,1)*F_dip_field(Ind_110,im) +&
               field_xform_3x3(2,2)*field_xform_3x3(3,2)*F_dip_field(Ind_020,im) +&
               field_xform_3x3(2,2)*field_xform_3x3(3,3)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(2,3)*field_xform_3x3(3,1)*F_dip_field(Ind_101,im) +&
               field_xform_3x3(2,3)*field_xform_3x3(3,2)*F_dip_field(Ind_011,im) +&
               field_xform_3x3(2,3)*field_xform_3x3(3,3)*F_dip_field(Ind_002,im)
  end do

end subroutine pGM_RECIP_tot_phi_grad
!--------------------------------------------------------------
subroutine pGM_RECIP_perm_to_fractional(numatoms,recip, &
                     nfft1,nfft2,nfft3,glob_mpole,frac_mpole, &
#ifdef MPI
                     imy_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,nfft1,nfft2,nfft3
  _REAL_,intent(in) :: recip(3,3)
  _REAL_,intent(in) :: glob_mpole(10,*)
  _REAL_,intent(out) :: frac_mpole(10,*)
  integer,intent(in) :: nmine
#ifdef MPI
  integer imy_cg(*)

# include "parallel.h"
# include "ew_parallel.h"
#endif

  integer order,dimxy,n
  _REAL_ :: mpole_xform_3x3(3,3),field_xform_3x3(3,3)
  _REAL_ Mpole_xy(35*35) ! the maximum is up to hexadecapole

  ! ZH/UCI: only include charge and dipole
  order = 4
  dimxy = 4

  ! first get mpole_xform_3x3
  call pGM_RECIP_xform_matrices(nfft1,nfft2,nfft3,recip, &
                      mpole_xform_3x3,field_xform_3x3)
  call XFORM_MPOLE_matrix(mpole_xform_3x3,Mpole_xy,order)

  ! current index in each thread should be imy_cg + n - 1,
  ! while we save in the same place in the frac mpole
  !
  ! note nmine = numatoms for serial job
  do n = 1, nmine
#ifdef MPI
    call XFORM_MPOLE(Mpole_xy,dimxy,glob_mpole(:,imy_cg(n)),frac_mpole(:,n),order)
#else
    call XFORM_MPOLE(Mpole_xy,dimxy,glob_mpole(:,n),frac_mpole(:,n),order)
#endif
  end do

  return
end subroutine pGM_RECIP_perm_to_fractional
!---------------------------------------------------------
subroutine pGM_RECIP_Fdip_to_Cdip_field( &
                     numatoms,nfft1,nfft2,nfft3, &
                     recip,frac_dipole_field,cart_dipole_field, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,nfft1,nfft2,nfft3
  _REAL_,intent(in) :: recip(3,3),frac_dipole_field(3,*)
  _REAL_,intent(out) :: cart_dipole_field(3,*)
  integer :: nmine,imine
#ifdef MPI
  integer my_cg(*)

# include "parallel.h"
# include "ew_parallel.h"
#endif

  integer :: n,i,j
  _REAL_ :: mpole_xform_3x3(3,3),field_xform_3x3(3,3)

  call pGM_RECIP_xform_matrices(nfft1,nfft2,nfft3,recip, &
                 mpole_xform_3x3,field_xform_3x3)

  ! add recip contribution to cart_dipole_field
  do imine = 1,nmine
#ifdef MPI
    n = my_cg(imine)
#else
    n = imine
#endif
    do i = 1,3
      do j = 1,3
        cart_dipole_field(i,n) = cart_dipole_field(i,n) + &
             field_xform_3x3(i,j)*frac_dipole_field(j,imine)
      end do
    end do
  end do

end subroutine pGM_RECIP_Fdip_to_Cdip_field
!--------------------------------------------------------------
subroutine pGM_RECIP_dipole_to_fractional( &
                     numatoms,nfft1,nfft2,nfft3, &
                     recip,ind_dip,f_dip, &
#ifdef MPI
                     my_cg, &
#endif
                     nmine)
  implicit none

  integer,intent(in) :: numatoms,nfft1,nfft2,nfft3, nmine
  _REAL_,intent(in)  :: recip(3,3),ind_dip(3,*)
  _REAL_,intent(out) :: f_dip(3,*)
#ifdef MPI
   integer my_cg(*)

#  include "parallel.h"
#  include "ew_parallel.h"
#endif

  integer :: n,i,j, imine
  _REAL_ :: mpole_xform_3x3(3,3),field_xform_3x3(3,3)

  call pGM_RECIP_xform_matrices(nfft1,nfft2,nfft3,recip, &
                      mpole_xform_3x3,field_xform_3x3)

  do imine = 1, nmine
#ifdef MPI
    n = my_cg(imine)
#else
    n = imine
#endif
    do i = 1,3
      f_dip(i,imine) = 0.d0
      do j = 1,3
        f_dip(i,imine) = f_dip(i,imine) + mpole_xform_3x3(i,j)*ind_dip(j,n)
      end do
    end do
  end do

end subroutine pGM_RECIP_dipole_to_fractional
!-----------------------------------------------------------------------
! RL: This routine is called repeatedly. Should save it for each step
subroutine pGM_RECIP_xform_matrices(nfft1,nfft2,nfft3,recip, &
                     mpole_xform_3x3,field_xform_3x3)
  implicit none

  integer,intent(in) :: nfft1,nfft2,nfft3
  _REAL_,intent(in) :: recip(3,3)
  _REAL_,intent(out) :: mpole_xform_3x3(3,3),field_xform_3x3(3,3)

  _REAL_ du1_dx,du1_dy,du1_dz,du2_dx,du2_dy,du2_dz,du3_dx,du3_dy,du3_dz

  du1_dx = nfft1*recip(1,1)
  du1_dy = nfft1*recip(2,1)
  du1_dz = nfft1*recip(3,1)
  du2_dx = nfft2*recip(1,2)
  du2_dy = nfft2*recip(2,2)
  du2_dz = nfft2*recip(3,2)
  du3_dx = nfft3*recip(1,3)
  du3_dy = nfft3*recip(2,3)
  du3_dz = nfft3*recip(3,3)

  field_xform_3x3(1,1) = du1_dx
  field_xform_3x3(1,2) = du2_dx
  field_xform_3x3(1,3) = du3_dx
  field_xform_3x3(2,1) = du1_dy
  field_xform_3x3(2,2) = du2_dy
  field_xform_3x3(2,3) = du3_dy
  field_xform_3x3(3,1) = du1_dz
  field_xform_3x3(3,2) = du2_dz
  field_xform_3x3(3,3) = du3_dz

  mpole_xform_3x3(1,1) = du1_dx
  mpole_xform_3x3(1,2) = du1_dy
  mpole_xform_3x3(1,3) = du1_dz
  mpole_xform_3x3(2,1) = du2_dx
  mpole_xform_3x3(2,2) = du2_dy
  mpole_xform_3x3(2,3) = du2_dz
  mpole_xform_3x3(3,1) = du3_dx
  mpole_xform_3x3(3,2) = du3_dy
  mpole_xform_3x3(3,3) = du3_dz

end subroutine pGM_RECIP_xform_matrices
!--------------------------------------------------------------
end module pol_gauss_recip
