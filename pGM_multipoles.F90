#include "../include/dprec.fh"
#include "../include/assert.fh"

!---------------------------------------------------------
module pol_gauss_multipoles

  implicit none
  private

  integer, parameter :: MAXMP=10
  _REAL_, parameter :: coulomb_const_kcal_per_mole = 332.05382d0 !Tinker

  integer, save :: num_multipoles
  integer, save :: start_multipoles,end_multipoles
  integer, save, allocatable :: rcvcntMP(:), displMP(:)
  _REAL_, save, allocatable :: global_multipole(:,:)

  !Gaussian radii
  _REAL_, parameter :: betamax = sqrt(2.0d0)/0.7147597d0 !Duan different
  _REAL_, dimension(:), allocatable  :: r_gauss, b_polar

  !Covalent local frame lists
  integer, save :: ncovalent ! add to global var
  integer, dimension(:), allocatable :: covalent_ptr, covalent_atm
  _REAL_, dimension(:), allocatable :: covalent_monop, covalent_dip

  public global_multipole,r_gauss,betamax,covalent_ptr,covalent_atm,covalent_dip, &
         MAXMP,coulomb_const_kcal_per_mole,start_multipoles,end_multipoles, &
         pGM_MPOLE_readparm,pGM_MPOLE_deallocate, &
         pGM_MPOLE_global_to_local,pGM_MPOLE_local_to_global
#ifdef MPI
  public :: pGM_MPOLE_bcast
#endif

  contains
!---------------------------------------------------------
#ifdef MPI
subroutine pGM_MPOLE_bcast()
  include 'mpif.h'
# include "extra.h"
# include "parallel.h"

  integer :: ier

  call mpi_bcast(num_multipoles,1,MPI_INTEGER,0,commsander,ier)
  call mpi_bcast(ncovalent,1,MPI_INTEGER,0,commsander,ier)
  if (.not.master) then
    allocate(covalent_ptr(0:num_multipoles),    &
             covalent_monop(num_multipoles),    &  ! charges
             b_polar(num_multipoles),           &  ! gaussian multipole radius
             r_gauss(num_multipoles),           &  ! stores 2*b_polar**2
             covalent_atm(1:ncovalent),         &
             covalent_dip(1:ncovalent),         &
             global_multipole(MAXMP,num_multipoles), &
             stat=ier)
    REQUIRE(ier==0)
  end if
     
  call mpi_bcast(covalent_ptr, num_multipoles + 1, mpi_integer, 0, commsander, ier)
  call mpi_bcast(covalent_monop, num_multipoles, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  call mpi_bcast(b_polar, num_multipoles, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  call mpi_bcast(r_gauss, num_multipoles, MPI_DOUBLE_PRECISION, 0, commsander, ier)
  call mpi_bcast(covalent_atm, ncovalent, MPI_INTEGER, 0, commsander, ier)
  call mpi_bcast(covalent_dip, ncovalent, MPI_DOUBLE_PRECISION, 0, commsander, ier)
 
end subroutine pGM_MPOLE_bcast
#endif 
!---------------------------------------------------------
function pGM_MPOLE_readparm(nf,num_atoms)
  use pol_gauss_mdin, only : pol_gauss_verbose

  integer :: pGM_MPOLE_readparm
  integer, intent(in) :: nf,num_atoms

  integer :: n,ier,dim1

  pGM_MPOLE_readparm = 0

  num_multipoles = num_atoms

  ! first. read in covalent local frame pointer
  allocate(covalent_ptr(0:num_atoms), stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call pGM_read_list_data('POL_GAUSS_COVALENT_POINTERS_',nf, &
                 dim1,num_atoms,covalent_ptr(1))
  ncovalent = sum (covalent_ptr(1:num_atoms))
  if ( pol_gauss_verbose == 3 ) then
    write(116,'(///)')
    write(116,*) 'Finding covalent pointers', num_atoms
    write(116,'(10i8)') covalent_ptr(1:num_atoms)
    write(116,'(///)')
  end if
  covalent_ptr(0) = 0
  do n = 1, num_atoms
    covalent_ptr(n) = covalent_ptr(n) + covalent_ptr(n-1)
  end do

  ! only read dipoles if present
  if ( ncovalent > 0 ) then

  ! second. read in covalent local frame atoms
  allocate(covalent_atm(1:ncovalent),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call pGM_read_list_data('POL_GAUSS_COVALENT_ATOMS_',nf, &
                 dim1,ncovalent,covalent_atm)
  if ( pol_gauss_verbose == 3 ) then
    write(116,'(///)')
    write(116,*) 'Finding covalent atoms', ncovalent
    write(116,'(10i8)') covalent_atm(1:ncovalent)
    write(116,'(///)')
  end if

  ! third. read in throughbond moments (in electron*Angstrom)
  allocate(covalent_dip(1:ncovalent),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call pGM_read_real_list_data('POL_GAUSS_COVALENT_DIPOLES_',nf, &
                 dim1,ncovalent,covalent_dip)
  if ( pol_gauss_verbose == 3 ) then
    write(116, '(///)')
    write(116,*) 'Finding throughbond moments',ncovalent 
    write(116,'(5E16.8)') covalent_dip(1:ncovalent)
    write(116, '(///)')
  end if

  end if

  ! fourth. read in gaussian monopole
  allocate(covalent_monop(num_atoms),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call pGM_read_real_list_data('POL_GAUSS_MONOPOLES_',nf, &
                 dim1,num_atoms,covalent_monop)
  if ( pol_gauss_verbose == 3 ) then
    write(116,'(///)')
    write(116,*) 'Finding gaussian monopoles', num_atoms
    write(116,'(5E16.8)') covalent_monop(1:num_atoms)
    write(116,'(///)')
  end if

  ! allocate global multipoles
  allocate(global_multipole(MAXMP,num_atoms),stat=ier)
  REQUIRE(ier==0)

  ! fifth. read in gaussian multipole radius
  allocate(b_polar(num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(r_gauss(num_atoms),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call pGM_read_real_list_data('POL_GAUSS_RADII_',nf, &
                 dim1,num_atoms,b_polar)
  if ( pol_gauss_verbose == 3 ) then
    write(116,'(///)')
    write(116,*) 'Finding gaussian radii', num_atoms
    write(116,'(5E16.8)') b_polar(1:num_atoms)
    write(116,'(///)')
  end if
  ! b_polar is read in as radius
  ! preprocess radius and beta for efficient combination
  ! this is to be consistent with pyresp
  r_gauss(1:num_atoms) = 2.0d0*b_polar(1:num_atoms)**2
  ! this is the original condition used for MD testing.
  !r_gauss(1:num_atoms) = 1.0d0*b_polar(1:num_atoms)**2

  pGM_MPOLE_readparm = 1
end function pGM_MPOLE_readparm
!---------------------------------------------------------
subroutine pGM_MPOLE_deallocate()
  if ( allocated(global_multipole) ) deallocate(global_multipole)
  if ( allocated(covalent_ptr)     ) deallocate(covalent_ptr)
  if ( allocated(covalent_monop)   ) deallocate(covalent_monop)
  if ( allocated(b_polar)          ) deallocate(b_polar)
  if ( allocated(r_gauss)          ) deallocate(r_gauss)
  if ( allocated(covalent_atm)     ) deallocate(covalent_atm)
  if ( allocated(covalent_dip)     ) deallocate(covalent_atm)
end subroutine pGM_MPOLE_deallocate
!---------------------------------------------------------
! RL: There is allocation/deallocation called every step
subroutine pGM_MPOLE_local_to_global(crd)
  use pol_gauss_mdin, only : pol_gauss_verbose
  use constants, only : ZERO, ONE

#ifdef MPI
  include 'mpif.h'
# include "parallel.h"
#endif

  _REAL_, intent(in) :: crd(3,*)

  integer ier, i
  integer :: iatm, jbond, jatm, start_bonds, end_bonds
  _REAL_, allocatable :: vector(:,:), norm(:)
#ifdef MPI
  integer :: id, startlist, endlist
  _REAL_, allocatable :: buffer(:,:)

  ! initialize the MPI send/receive data
  if (.not. allocated(rcvcntMP)) then
    allocate(rcvcntMP(0:numtasks-1), stat=ier)
    REQUIRE(ier==0)
  end if
  if (.not. allocated(displMP)) then
    allocate(displMP(0:numtasks-1), stat=ier)
    REQUIRE(ier==0)
  end if

  do id = 0, numtasks-1
    startlist = iparpt(id) + 1
    endlist = iparpt(id+1)
    rcvcntMP(id) = MAXMP*(endlist-startlist+1)
    displMP(id) = MAXMP*(startlist-1)+1
  end do

  ! initialize the MPI send/receive partition
  allocate(buffer(0:MAXMP-1,num_multipoles), stat=ier)
  REQUIRE(ier==0)
  buffer(0:MAXMP-1,1:num_multipoles) = ZERO

  ! set up multipole partitions to be consistent with sander routines 
  start_multipoles = iparpt(mytaskid) + 1
  end_multipoles = iparpt(mytaskid+1)
#else
  start_multipoles = 1
  end_multipoles = num_multipoles
#endif

  global_multipole(1,start_multipoles:end_multipoles) = &
  covalent_monop(start_multipoles:end_multipoles)
  global_multipole(2:4,start_multipoles:end_multipoles) = ZERO
  do iatm = start_multipoles, end_multipoles
    start_bonds = covalent_ptr(iatm-1) + 1
    end_bonds = covalent_ptr(iatm)

    ! if there is no dipole on this atom, skip it
    if ( end_bonds < start_bonds ) exit

    ! compute through-bond unit vectors
    allocate(vector(1:3,start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    allocate(norm(start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    do jbond = start_bonds, end_bonds
      jatm = covalent_atm(jbond)
      vector(1:3,jbond) = crd(1:3,jatm) - crd(1:3,iatm)
      norm(jbond) = ONE/sqrt(sum(vector(1:3,jbond)**2))
    enddo
    do jbond = start_bonds, end_bonds
      vector(1:3,jbond) = vector(1:3,jbond)*norm(jbond)
    enddo

    ! recover global dipoles
    do jbond = start_bonds, end_bonds
      global_multipole(2:4,iatm) = global_multipole(2:4,iatm) + covalent_dip(jbond)*vector(1:3,jbond)
    enddo
    deallocate(vector, stat=ier)
    REQUIRE(ier==0)
    deallocate(norm, stat=ier)
    REQUIRE(ier==0)
  enddo
 
#ifdef MPI
  call mpi_barrier(commsander, ier)
  call mpi_allgatherv(global_multipole(1,(displMP(mytaskid)-1)/MAXMP+1), rcvcntMP(mytaskid),&
          MPI_DOUBLE_PRECISION, buffer, rcvcntMP, displMP, MPI_DOUBLE_PRECISION, commsander, ier)
  global_multipole(1:4,1:num_multipoles) = buffer(1:4,1:num_multipoles)
  deallocate(buffer, stat=ier)
  REQUIRE(ier==0)
#endif

  if ( pol_gauss_verbose == 3 ) then
    write(6,*) 'converted global multipoles from through-bond dipoles (D)'
    write(6,'(f15.6)') global_multipole(1:1,1:num_multipoles)
    write(6,'(3f15.6)') global_multipole(2:4,1:num_multipoles)*4.803239377d0
  end if

end subroutine pGM_MPOLE_local_to_global
!---------------------------------------------------------
! if global multipoles are read from prmtop file in different units
! do the conversions here.
!---------------------------------------------------------
subroutine pGM_MPOLE_rescale_multipoles()

  integer j,n
  _REAL_ bohr

  ! change dipole unit to electron*Angstrom
  !bohr = 0.5291772083d0 ! R. Luo: if input uses Bohr as the length unit
  !bohr = 1.0d0 ! R. Luo: if input dipoles are read in as electron*Angstrom
  bohr = 1.0d0/4.803239377d0 ! R. Luo: if input dipoles are read in as debyes

  !the following is changed to scale input global moments only
  do n = 1,num_multipoles
    do j = 2,4
      global_multipole(j,n) = bohr*global_multipole(j,n)
    end do
  end do
end subroutine pGM_MPOLE_rescale_multipoles
!---------------------------------------------------------
! if multipoles are read from prmtop file in the global frame
! the input multipoles can be converted into covalent dipoles here
!---------------------------------------------------------
subroutine pGM_MPOLE_global_to_local(crd)
  use pol_gauss_mdin, only : pol_gauss_verbose
  use constants, only : ZERO, ONE
  use stack

  _REAL_, intent(in) :: crd(3,*)

  integer :: iatm, jbond, jatm, start_bonds, end_bonds
  integer, parameter :: mp = 3, np = 4
  integer i, j, info, ier
  integer m, n
  _REAL_, parameter :: TOL = 1d-6
  _REAL_ a(mp,np), u(mp,np), w(np), v(np,np), vt(np,np), b(mp), t(np)
  _REAL_ wmax, thresh
  _REAL_ work(512)
  _REAL_, allocatable :: vector(:,:), norm(:)

  start_bonds = 1
  do iatm = start_multipoles, end_multipoles
    start_bonds = start_bonds + covalent_ptr(iatm-1)
    end_bonds = start_bonds + covalent_ptr(iatm) - 1

    ! compute through-bond basis vectors
    allocate(vector(1:3,start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    allocate(norm(start_bonds:end_bonds), stat=ier)
    REQUIRE(ier==0)
    do jbond = start_bonds, end_bonds
      jatm = covalent_atm(jbond)
      vector(1:3,jbond) = crd(1:3,jatm) - crd(1:3,iatm)
      norm(jbond) = ONE/sqrt(sum(vector(1:3,jbond)**2))
    end do
    do jbond = start_bonds, end_bonds
      vector(1:3,jbond) = vector(1:3,jbond)*norm(jbond)
    end do

    ! set up linear system for least square fitting.
    ! Here the goal is to obtain the coefficients T for X, the basis vectors
    ! so that X * T = b. Thus linear system matrix A is X as written when
    ! solving for T.
    a = ZERO; u = ZERO; w = ZERO; vt = ZERO; b = ZERO; t = ZERO
    m = 3
    n = covalent_ptr(iatm)
    if (n > np) then
      write(6,*) ' pGM_MPOLE_global_to_local(): Too many bonds for SVD'
      call mexit(6,1)
    end if
    do jbond = start_bonds, end_bonds
      a(1:3, jbond-start_bonds+1) = vector(1:3, jbond)
    end do
    b(1:3) = global_multipole(2:4, iatm)
    deallocate(vector, stat=ier)
    REQUIRE(ier==0)
    deallocate(norm, stat=ier)
    REQUIRE(ier==0)

    ! svd call for U, W, V as in A = U * 1/W * V^T 
    ! DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call dgesvd('A','A',m,n,a,mp,w,u,mp,vt,np,work,512,info)
    do i = 1, np
      do j = 1, np
        v(i,j) = vt(j,i)
      end do
    end do
    wmax = ZERO
    do j = 1, n
      if ( w(j) > wmax ) wmax = w(j)
    end do
    thresh = TOL * wmax
    do j = 1, n
      if ( w(j) < thresh ) w(j) = ZERO
    end do

    ! back substitution call for A^-1 * b
    call svbksb(u,w,v,m,n,mp,np,b,t)

    ! store through-bond dipole moment in T for later
    if ( pol_gauss_verbose == 3 ) write(6,*) ' Atom info', iatm, info
    do jbond = start_bonds, end_bonds
      covalent_dip(jbond) = t(jbond-start_bonds+1)
      if ( pol_gauss_verbose == 3 ) write(6,*) &
        ' through-bond moments', covalent_dip(jbond)
    end do
  end do ! end of iatm = start_multipoles, end_multipoles

  if ( pol_gauss_verbose == 3 ) then
    write(6,*) 'global multipoles usd to compute through-bond moments (D)'
    write(6,'(f15.6)') global_multipole(1:1,1:num_multipoles)
    write(6,'(3f15.6)') global_multipole(2:4,1:num_multipoles)*4.803239377d0
  end if
end subroutine pGM_MPOLE_global_to_local
!---------------------------------------------------------
end module pol_gauss_multipoles
