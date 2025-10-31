#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"


module qm_ewald_comm
  implicit none

  !
  ! Some external semiempirical QM packages are interfaced to sander;
  ! however, they are typically parallelized with OpenMP rather than
  ! MPI. Only the rank 1 sander process is used within the external
  ! program; therefore the distributed memory of the QM-Ewald is
  ! not readily available in the external program's SCF procedure.
  ! The QMQMEwCommT will use MPI to gather the required k-space
  ! quantities and store them within a QMQMEwCommT object.
  !
  ! The class has 2 public routines
  !
  !  qmqmewcomm_setup -
  !          Each mpi rank calls this at every time step before
  !          the scf procedure is started. It will precompute the
  !          QM/MM potential and the necessary quantities for
  !          the QM/QM Ewald. The QM/QM Ewald quantities are
  !          then copied to rank 1.
  !
  !  qmqmewcomm_calcene -
  !          This should be called only by the rank 1 process
  !          within the external program SCF procedure to
  !          get the potential experienced by each QM atom.
  !          It also returns a correction to the energy 
  !          corresponding to the QM-QM Ewald interaction.
  !

  
  private

  public :: QMQMEwCommT
  public :: QMQMEwComm_Setup
  public :: QMQMEwComm_CalcEne

  
  type QMQMEwCommT
     
     integer, pointer :: kcounts(:) => null()
     integer, pointer :: kdisps(:) => null()
     integer, pointer :: ktsizes(:) => null()
     integer, pointer :: ktdisps(:) => null()
     _REAL_,pointer :: qmktable(:,:,:) => null()
     _REAL_,pointer :: kvec(:) => null()

  end type QMQMEwCommT

  
  
contains

  subroutine QMQMEwComm_Setup(self,first_call,natom,scaled_mm_charges)

    use qmmm_module, only: qmmm_struct
    use qmmm_module, only: qmmm_nml
    use qmmm_module, only: qm2_struct
    use qmmm_module, only: qm2_rij_eqns
    use qmmm_module, only: qmewald

#ifdef MPI
    use qmmm_module, only : qmmm_scratch
    use qmmm_module, only : qmmm_mpi
#endif

    implicit none

    type(QMQMEwCommT),intent(inout) :: self
    logical,intent(in) :: first_call
    integer,intent(in) :: natom
    _REAL_,intent(in) :: scaled_mm_charges(natom)

#ifdef MPI
    integer :: ier
    include 'mpif.h'
#endif

    if ( qmmm_nml%qm_ewald>0 ) then

       if ( first_call ) then
          call qm2_load_params_and_allocate(.false.)
       end if

       call qm2_calc_rij_and_eqns(qmmm_struct%qm_coords, &
            qmmm_struct%nquant_nlink,qmmm_struct%qm_xcrd, &
            natom, qmmm_struct%qm_mm_pairs)

       qm2_struct%hmatrix = 0.d0  !Note this zeros the entire array
       qmmm_struct%enuclr_qmqm = 0.d0
       qmmm_struct%enuclr_qmmm = 0.d0

       call qm_ewald_mm_pot(qmmm_struct%qm_xcrd, qmmm_struct%qm_mm_pairs, &
            & qmmm_struct%qm_coords, natom, &
            & qmmm_nml%qmmmrij_incore, qm2_rij_eqns%qmmmrijdata, &
            & scaled_mm_charges, qmewald%kvec )

#ifdef MPI
#ifdef USE_MPI_IN_PLACE
       if (qmmm_mpi%commqmmm_master) then
          call mpi_reduce(MPI_IN_PLACE,qmewald%mmpot,qmmm_struct%nquant_nlink, &
               MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
       else
          call mpi_reduce(qmewald%mmpot,0,qmmm_struct%nquant_nlink, &
               MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
       end if
#else
       call mpi_reduce(qmewald%mmpot,qmmm_scratch%matsize_red_scratch, &
            & qmmm_struct%nquant_nlink, &
            & MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier )
       if (qmmm_mpi%commqmmm_master) &
            qmewald%mmpot(1:qmmm_struct%nquant_nlink) = &
            & qmmm_scratch%matsize_red_scratch(1:qmmm_struct%nquant_nlink)
#endif
#endif

#ifdef MPI
       if ( qmmm_mpi%numthreads > 1 ) then
          if ( first_call ) then
             call QMQMEwComm_init(self)
          end if
          ! This collections kvec and qmktable on the lead rank
          call QMQMEwComm_gather(self)
       end if
#endif

    end if ! qm_ewald > 0
  end subroutine QMQMEwComm_Setup


  subroutine QMQMEwComm_init(self)
    use qmmm_module, only : qmmm_struct, qmmm_nml, qmewald
#ifdef MPI
    use qmmm_module, only : qmmm_mpi
#endif
    implicit none
    type(QMQMEwCommT), intent(inout) :: self
    integer :: ier,i

#ifdef MPI
    !integer,allocatable :: sbuf(:)
    include 'mpif.h'

    if ( associated(self%kcounts) ) then
       deallocate(self%kcounts)
    end if
    allocate( self%kcounts( qmmm_mpi%numthreads ) )
    self%kcounts = 0

    if ( associated(self%kdisps) ) then
       deallocate(self%kdisps)
    end if
    allocate( self%kdisps( qmmm_mpi%numthreads ) )
    self%kdisps = 0

    if ( associated(self%ktsizes) ) then
       deallocate(self%ktsizes)
    end if
    allocate( self%ktsizes( qmmm_mpi%numthreads ) )
    self%ktsizes = 0

    if ( associated(self%ktdisps) ) then
       deallocate(self%ktdisps)
    end if
    allocate( self%ktdisps( qmmm_mpi%numthreads ) )
    self%ktdisps = 0

    if ( associated( self%qmktable ) ) then
       deallocate( self%qmktable )
    end if
    allocate( self%qmktable(6,qmmm_struct%nquant_nlink,qmewald%totkq) )
    self%qmktable = 0.d0

    if ( associated( self%kvec ) ) then
       deallocate( self%kvec )
    end if
    allocate( self%kvec(qmewald%totkq) )
    self%kvec = 0.d0

    ier = 0
    call mpi_gather( qmmm_mpi%totkq_count, 1, mpi_integer, &
         & self%kcounts(1), 1, mpi_integer, 0, qmmm_mpi%commqmmm, ier )

    self%ktsizes = self%kcounts * 6*qmmm_struct%nquant_nlink
    do i=2,qmmm_mpi%numthreads
       self%ktdisps(i) = self%ktdisps(i-1) + self%ktsizes(i-1)
       self%kdisps(i) = self%kdisps(i-1) + self%kcounts(i-1)
    end do

    call mpi_bcast(self%kcounts(1),qmmm_mpi%numthreads,mpi_integer,&
         & 0, qmmm_mpi%commqmmm, ier )

    call mpi_bcast(self%kdisps(1),qmmm_mpi%numthreads,mpi_integer,&
         & 0, qmmm_mpi%commqmmm, ier )

    call mpi_bcast(self%ktsizes(1),qmmm_mpi%numthreads,mpi_integer,&
         & 0, qmmm_mpi%commqmmm, ier )

    call mpi_bcast(self%ktdisps(1),qmmm_mpi%numthreads,mpi_integer,&
         & 0, qmmm_mpi%commqmmm, ier )


#endif
  end subroutine QMQMEwComm_init


  
  subroutine QMQMEwComm_gather(self)
    !
    ! Copy the distributed memory from the
    ! qmewald data structure to the QMQMEwCommT
    ! object on rank 1
    !
    use qmmm_module, only : qmewald
#ifdef MPI
    use qmmm_module, only : qmmm_mpi
#endif
    implicit none
    type(QMQMEwCommT), intent(inout) :: self
#ifdef MPI
    integer :: ier
    include 'mpif.h'

    ier = 0

    self%qmktable = 0.d0

    call mpi_gatherv( qmewald%qmktable(1,1,1), &
         & self%ktsizes(qmmm_mpi%mytaskid+1), &
         & mpi_double_precision, &
         & self%qmktable(1,1,1), &
         & self%ktsizes, &
         & self%ktdisps, &
         & mpi_double_precision, &
         & 0, qmmm_mpi%commqmmm, ier )


    self%kvec = 0.d0

    call mpi_gatherv( qmewald%kvec(1), &
         & self%kcounts(qmmm_mpi%mytaskid+1), &
         & mpi_double_precision, &
         & self%kvec(1), &
         & self%kcounts, &
         & self%kdisps, &
         & mpi_double_precision, &
         & 0, qmmm_mpi%commqmmm, ier )

#endif

  end subroutine QMQMEwComm_gather





  subroutine QMQMEwComm_qm_pot(self,scf_mchg,qm_coords)
    ! This is basically a copy-paste of
    ! qm_ewald_qm_pot in qm_ewald.F90
    ! written by Ross Walker (TSRI, 2005)
    ! HOWEVER, the calculation is forced to be
    ! serial using the MPI-gathered data
    ! stored in the QMQMEwCommT object
    ! rather than the distributed memory in
    ! qmmm_struct.

    !Calculates QM potential in qmewald%qmpot
    !Requires the Mulliken Charges (or equivalent).

    use qmmm_module, only : qmmm_struct, qmewald
    use constants, only : PI, INVSQRTPI, zero, one, two
    use nblist, only : volume  
    implicit none

    !Passed in

    !integer, intent(in) :: nquant, nlink
    type(QMQMEwCommT),intent(in) :: self
    _REAL_, intent(in) :: scf_mchg(qmmm_struct%nquant_nlink)
    _REAL_, intent(in) :: qm_coords(3,qmmm_struct%nquant_nlink) 
    !_REAL_, intent(in) :: kvec(qmmm_mpi%totkq_count)

    !Local
    integer :: i, j
    integer :: loop_count
    integer :: iminus, iplus
    _REAL_ :: esfact, wfact, qmi(1:3), vec(1:3), kvectmp
    _REAL_ :: r2, oneRIJ, RIJ, erfcx, chg_qm
    _REAL_ :: ktgs(8), x_cos, y_cos, z_cos, x_sin, y_sin, z_sin, ksum

#include "qm2_array_locations.h"

    !For ew_coeff
    !#include "ew_pme_recip.h"  -replaced by qmewald%kappa: TL

    !We only need to do most of this work if we are changing the mulliken charges on the QM image atoms
    !at every step (qm_ewald=1). If qm_ewald=2 then the charges on the QM image atoms are fixed at the SCF
    !charges from the previous MD step and were included in mmpot. But if we are not updating the image charges
    !we still need to if this is our first MD step.

    !===================================
    !    Calculate K Space Potential 
    !===================================
    !Calculate the K space potential between only QM atoms
    loop_count = 0
    qmewald%qmpot(1:qmmm_struct%nquant_nlink) = zero !Zero's entire array
    do i = 1, qmewald%totkq
       !do i = qmmm_mpi%kvec_start, qmmm_mpi%kvec_end
       loop_count = loop_count+1
       ktgs(1:8) = zero
       kvectmp = two*self%kvec(loop_count)
       do j = 1, qmmm_struct%nquant_nlink
          x_cos =self%qmktable(1,j,loop_count)
          x_sin =self%qmktable(2,j,loop_count)
          y_cos =self%qmktable(3,j,loop_count)
          y_sin =self%qmktable(4,j,loop_count)
          z_cos =self%qmktable(5,j,loop_count)
          z_sin =self%qmktable(6,j,loop_count)
          chg_qm = scf_mchg(j)
          ktgs(1) = ktgs(1) + chg_qm*x_cos*y_cos*z_cos
          ktgs(2) = ktgs(2) + chg_qm*x_sin*y_sin*z_sin
          ktgs(3) = ktgs(3) + chg_qm*x_cos*y_cos*z_sin
          ktgs(4) = ktgs(4) + chg_qm*x_cos*y_sin*z_cos
          ktgs(5) = ktgs(5) + chg_qm*x_cos*y_sin*z_sin
          ktgs(6) = ktgs(6) + chg_qm*x_sin*y_sin*z_cos
          ktgs(7) = ktgs(7) + chg_qm*x_sin*y_cos*z_sin
          ktgs(8) = ktgs(8) + chg_qm*x_sin*y_cos*z_cos
       end do
       do j = 1, qmmm_struct%nquant_nlink
          x_cos =self%qmktable(1,j,loop_count)
          x_sin =self%qmktable(2,j,loop_count)
          y_cos =self%qmktable(3,j,loop_count)
          y_sin =self%qmktable(4,j,loop_count)
          z_cos =self%qmktable(5,j,loop_count)
          z_sin =self%qmktable(6,j,loop_count)
          ksum = zero
          ksum = ksum + x_cos*y_cos*z_cos*ktgs(1)
          ksum = ksum + x_sin*y_sin*z_sin*ktgs(2)
          ksum = ksum + x_cos*y_cos*z_sin*ktgs(3)
          ksum = ksum + x_cos*y_sin*z_cos*ktgs(4)
          ksum = ksum + x_cos*y_sin*z_sin*ktgs(5)
          ksum = ksum + x_sin*y_sin*z_cos*ktgs(6)
          ksum = ksum + x_sin*y_cos*z_sin*ktgs(7)
          ksum = ksum + x_sin*y_cos*z_cos*ktgs(8)
          qmewald%qmpot(j) = qmewald%qmpot(j) + ksum*kvectmp
       end do
    end do
    !===================================
    !  End Calculate K Space Potential 
    !===================================

    !======================================
    !    Calculate Real Space Potential
    !======================================

    !We have to calculate QM-QM distances on the fly for qmewald
    ! do i = 1, nquant
    !Note: Ross Walker - Should really split this loop at some point
    !                    but be careful to do balanced in parallel

    esfact = -two*qmewald%kappa*INVSQRTPI
    wfact = -PI/(qmewald%kappa*qmewald%kappa)/volume

    do i =1,qmmm_struct%nquant_nlink
       !do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
       chg_qm = scf_mchg(i)
       qmi(1:3) = qm_coords(1:3,i)
       iminus = i-1
       do j = 1,iminus
          vec(1:3) = qmi(1:3) - qm_coords(1:3,j)
          r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
          onerij = one / sqrt(r2)
          rij = r2*oneRIJ !one/oneRIJ
          call erfcfun(rij*qmewald%kappa, erfcx)
          qmewald%qmpot(j) = qmewald%qmpot(j) - chg_qm*(one - erfcx)*onerij + wfact*chg_qm
       end do
       !i=j
       qmewald%qmpot(i) = qmewald%qmpot(i) + esfact*chg_qm + wfact*chg_qm
       iplus = i+1
       do j = iplus, qmmm_struct%nquant_nlink
          vec(1:3) = qmi(1:3) - qm_coords(1:3,j)
          r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
          onerij = one / sqrt(r2)
          rij = r2*oneRIJ !one/oneRIJ
          call erfcfun(rij*qmewald%kappa, erfcx)
          qmewald%qmpot(j) = qmewald%qmpot(j) - chg_qm*(one - erfcx)*onerij + wfact*chg_qm
       end do
    end do

    !======================================
    !  End Calculate Real Space Potential
    !======================================

    return

  end subroutine QMQMEwComm_qm_pot


  subroutine QMQMEwComm_CalcEne(self,scf_mchg,emmew,eqmcor,shift)

    !
    ! scf_mchg : the current mulliken charges (au)
    ! shift : the net electrostatic potential at the qm atoms (au)
    ! emmew : the QM charge interaction with the MM atoms, excluding
    !         the nearby real-space QM/MM interactions
    ! eqmew : the QM charge interaction with the QM images
    !
    ! The net QM/MM electrostatic energy is
    !
    !  E = eqmcor + emmew + E_{qm,coulomb} + E_{qm/mm,within cutoff}
    !
    !  E_{qm/mm,within cutoff} is the gas-phase-like real-space
    !   interactions within a cutoff
    !
    !  E_{qm,coulomb} is the Coulomb self-interaction of the
    !  primary cell's QM region.
    !
    !  eqmcor is the interaction of the primary cell's QM region
    !  with the distant QM images (the reciprocal space potential
    !  corrected with the short range Gaussian potential)
    !
    ! emmew similarly the interaction of the primary cell's QM region
    ! with all of the MM region (including its images) EXCEPT for
    ! the gas-phase-like real-space interactions with the nearby
    ! point charges within a cutoff.
    !
    ! Some QM implementations will use the QM potentials (the shift)
    ! array to construct a Fock matrix and compute an energy term
    ! Eext = \sum_{ij} F_{ij} P_{ij}
    ! where P_{ij} is the density matrix
    !
    ! In this case, one is double counting the interaction of the
    ! QM region with the QM images, so one would need to correct
    ! this by subtracting Eext = Eext - eqmcor
    !
    ! In other QM implementations, the shift values polarize the
    ! Fock matrix, but no energy term is directly computed from it.
    ! In this case, one must manually calculate
    ! Eext = eqmcor + qmmew
    !  
    !
    
#ifdef MPI
    use qmmm_module, only : qmmm_scratch
#endif
    use qmmm_module, only : qmmm_struct, qmmm_nml, qmewald, qmmm_mpi
    use constants, only: BOHRS_TO_A

    implicit none

    type(QMQMEwCommT),intent(in) :: self
    _REAL_,intent(in) :: scf_mchg(:)
    _REAL_,intent(out) :: emmew,eqmcor
    _REAL_,intent(inout) :: shift(:)

    _REAL_,allocatable :: netshift(:)

    integer :: i
#ifdef MPI
    include 'mpif.h'
    integer :: ier
#endif

    emmew = 0.d0
    eqmcor = 0.d0

    if ( qmmm_nml%qm_ewald /= 1 ) then
       return
    end if

    allocate( netshift( qmmm_struct%nquant_nlink ) )
    netshift = 0.d0

    ! Calculates the Ewald potential at the position of the atoms
    if (qmmm_nml%qm_ewald==1 .OR. qmewald%ewald_startup) then
       qmewald%qmpot = 0.d0

#ifdef MPI
       if ( qmmm_mpi%numthreads == 1 ) then
          call qm_ewald_qm_pot(qmmm_struct%nquant, qmmm_struct%nlink, scf_mchg,&
               qmmm_struct%qm_coords,qmewald%kvec)
       else
          call QMQMEwComm_qm_pot(self,scf_mchg,qmmm_struct%qm_coords)
       end if
#else
       call qm_ewald_qm_pot(qmmm_struct%nquant, qmmm_struct%nlink, scf_mchg,&
            qmmm_struct%qm_coords,qmewald%kvec)
#endif


       ! The potential comes in au/A, and needs to be converted to au/Bohrs.
       do i = 1, qmmm_struct%nquant_nlink
          netshift(i) = ( qmewald%mmpot(i) + qmewald%qmpot(i) ) * BOHRS_TO_A
          !netshift(i) = netshift(i) + qmewald%mmpot(i) * BOHRS_TO_A
          !netshift(i) = netshift(i) + qmewald%qmpot(i) * BOHRS_TO_A
          emmew  = emmew + &
               & scf_mchg(i) * qmewald%mmpot(i) * BOHRS_TO_A
          eqmcor = eqmcor + 0.5d0 * &
               & scf_mchg(i) * qmewald%qmpot(i) * BOHRS_TO_A
       end do
       
       ! This removes 1/2 of the long-range qm/qm interactions
       ! because the self term should have a leading factor
       ! of 1/2, but energy contribution calculated from
       ! the potential is missing that factor, e.g.,
       ! Eext = dot(Qqm,Pqm) + dot(Qqm,Pmm) (INCORRECT)
       !   versus
       ! Eext = 0.5 * dot(Qqm,Pqm) + dot(Qqm,Pmm) (CORRECT)
       ! What's actually computed is
       !   Eext = dot(Qqm,Pqm) + dot(Qqm,Pmm) - 0.5*dot(Qqm,Pqm)
       ! The energy being returned is is -0.5*dot(Qqm,Pqm)
       
       !ene = - 0.5d0 * ene

       ! ... actually, it depends on the implementation
       ! Some codes will compute Eext from the shift array
       ! whereas others will only polarize the Fock matrix
       
    end if

    do i=1,qmmm_struct%nquant_nlink
       shift(i) = shift(i) + netshift(i)
    end do

  end subroutine QMQMEwComm_CalcEne

end module qm_ewald_comm
