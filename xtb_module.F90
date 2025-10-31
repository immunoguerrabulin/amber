#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

! QM/MM interface to xtb library
!
! Written by Timothy J. Giese
! Date: December 2023
!
! Example Installation (07 Dec 2023)
! ---------------------------------
!
! 1. Download xtb from https://github.com/grimme-lab/xtb.git
!
! mkdir ~/software
! cd ~/software
! git clone https://github.com/grimme-lab/xtb.git
! cd xtb
!
! 2. Configure xtb with API shared library and XTB support.
!
! export XTBHOME=${HOME}/software/xtb/local
! export FC=gfortran
! export CC=gcc
! cmake -B_build \
!        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
!        -DBUILD_SHARED_LIBS=TRUE \
!        -DINSTALL_MODULES=TRUE \
!        -DCMAKE_EXE_LINKER_FLAGS="-lm" \
!        -DBLAS_LIBRARIES="-lopenblas" \
!        -DLAPACK_LIBRARIES="-lopenblas" \
!        -DCMAKE_INSTALL_PREFIX=${XTBHOME}
!
! 3. Build and install dftbplus
!
! make -C _build
! make -C _build install
!
! 4. You should now see the following files:
!
! ${XTBHOME}/lib64/libxtb.so
! ${XTBHOME}/include/xtb.h
! ${XTBHOME}/include/xtb_xtb_calculator.mod
!
! If you are missing the shared library, then
! you may need to DOWNGRADE your cmake to version 3.26.0
! because the xtb developers have disabled support for
! more recent versions due to a bug in cmake.
!
! 5. Change to the amber build directory
!
! cd /path/to/amber/build
!
! 6. Edit run_cmake to define USE_XTB and XTB_DIR; e.g.,
!
!   -DUSE_XTB=ON
!   -DXTB_LIBRARY=${XTBHOME}/lib64/libxtb.so
!   -DXTB_INCLUDE_DIR="${XTBHOME}/include"
!
! 7. Compile and install amber
!
! bash ./run_cmake
! make install
!
!

module xtb_module

  !use constants, only: A_TO_BOHRS, AU_TO_KCAL
  use constants, only: CODATA08_A_TO_BOHRS, CODATA08_AU_TO_KCAL
  use qm_ewald_comm
  
#ifdef XTB
  use xtb_type_solvation
  use xtb_type_environment
  use xtb_type_restart
  use xtb_type_molecule
  use xtb_xtb_calculator
  use xtb_type_data, only : scc_results
#endif
  
  implicit none
  private

  public :: get_xtb_qmmm_forces
  
#ifdef XTB

  
  ! XTB namelist contains user definable settings for XTB
  type xtb_nml_type
     
     ! GFN1-XTB,GFN2-XTB
     character(len=80) :: qm_level = "GFN2-XTB"
     
     ! electronic temperature (K)
     _REAL_  :: tfermi = 300.d0

     ! The scvconv value is 1.e-6 * "accuracy",
     ! where accuracy is in the range 1.e-6 to 1.e+3
     _REAL_  :: accuracy = 1.0d-3

     ! max scf cycles
     integer :: maxiter = 250

     ! chemical hardness (au) of the MM atoms
     ! when computing electrostatics
     _REAL_  :: mmhardness = 0.d0
     
     logical :: debug = .false.
     
  end type xtb_nml_type

  !
  ! qm_level, default = "GFN2-XTB"
  !           valid options: GFN2-XTB, GFN1-XTB
  !           
  ! tfermi, default = 300.0
  !
  ! accuracy, default = 1.e-3
  !           valid: 1.e-5 <= accuracy <= 1.e+3
  !           Smaller values mean more accuracy.
  !
  !           XTB will use a scfconv value of
  !           1.e-6 * accuracy.
  !
  !   The accuracy value roughly reflects the
  !   desired relative force error. I've found
  !   that 1.e-5 gives very good agreement between
  !   analytic and numerical forces. A value of
  !   1.e-4 appears to give accurate forces relative
  !   to 1.e-5, but in comparison to the 1.e-4
  !   numerical values the force error is in the
  !   range 1.e-4 to 1.e-3.
  !
  !
  ! mmhardness, default = 0.0
  !
  !    If mmhardness > 0, then ALL MM atoms use the same
  !    hardness with the specified value.
  !
  !    If mmhardness = 0.0, then ALL MM atoms se the
  !    hardness of hydrogen. This is the same as setting
  !    mmhardness=0.405771
  !
  !    If mmhardness = -1.0, then each MM atom's hardness
  !    is chosen from its atomic number, and extra point
  !    particles use a hardness of hydrogen.
  !
  !    More generally, if mmhardness < 0, then
  !           mmgam(i) = abs(mmhardness) * gam(z(i))
  !    That is, the hardness is chosen from the MM atomic
  !    number, but it is scaled by abs(mmhardness).
  !

  
  integer, parameter :: wp=kind(1.0d0)


!  _REAL_,parameter :: AU_PER_AMBER_CHARGE = 1.d0 / AMBER_ELECTROSTATIC
!  _REAL_,parameter :: AU_PER_AMBER_ANGSTROM = 1.d0 / BOHRS_TO_A
!  _REAL_,parameter :: AU_PER_AMBER_KCAL_PER_MOL = &
!       & AU_PER_AMBER_CHARGE *  AU_PER_AMBER_CHARGE / AU_PER_AMBER_ANGSTROM
!  _REAL_,parameter :: AU_PER_AMBER_FORCE = &
!       & AU_PER_AMBER_KCAL_PER_MOL / AU_PER_AMBER_ANGSTROM


  !_REAL_,parameter :: AU_PER_ANG = A_TO_BOHRS
  !_REAL_,parameter :: AU_PER_KCAL = 1.d0 / AU_TO_KCAL
  
  _REAL_,parameter :: AU_PER_ANG = CODATA08_A_TO_BOHRS
  _REAL_,parameter :: AU_PER_KCAL = 1.d0 / CODATA08_AU_TO_KCAL

  
  
  type, extends(TSolvation) :: XTBCallbackT
     type(QMQMEwCommT) :: ewald
     _REAL_ :: ewald_qmqm
     _REAL_ :: ewald_qmmm
     _REAL_ :: ewald_ene
   contains
     procedure :: update => xtb_sander_update
     procedure :: addShift => xtb_sander_addShift
     procedure :: getEnergy => xtb_sander_getEnergy
     procedure :: addGradient => xtb_sander_addGradient
  end type XTBCallbackT


  type XTBIface
     type(xtb_nml_type) :: xtb_nml
     type(TEnvironment) :: env
     type(TxTBCalculator) :: calc
     type(TRestart) :: res
     type(TMolecule) :: mol
     type(scc_results) :: spRes
     logical :: exitRun
     _REAL_  :: sigma(3,3)
     _REAL_  :: egap
  end type xtbiface

  
  type(XTBIface),save :: XTBstate
  type(XTBCallbackT),save :: XTBcallback

#endif

  
contains

  
#ifdef XTB


  subroutine xtb_sander_update(self, env, num, xyz)

    implicit none
    
    !> Instance of the solvation model
    class(XTBCallbackT), intent(inout) :: self
    
    !> Computation environment
    type(TEnvironment), intent(inout) :: env
    
    !> Atomic numbers
    integer, intent(in) :: num(:)
    
    !> Cartesian coordinates
    real(wp), intent(in) :: xyz(:, :)

    
  end subroutine xtb_sander_update


  
  subroutine xtb_sander_addShift(self, env, qat, &
       & qsh, atomicshift, shellshift)

    implicit none
    
    !> Instance.
    class(XTBCallbackT), intent(inout) :: self
    
    !> Computation environment
    type(TEnvironment), intent(inout) :: env

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(wp), intent(in) :: qat(:)

    !> Shell-resolved net number of electrons
    real(wp), intent(in) :: qsh(:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(wp), intent(inout) :: atomicshift(:)

    !> Shell-resolved external potential shift
    real(wp), intent(inout) :: shellshift(:)

    call QMQMEwComm_CalcEne( XTBcallback%ewald, qat, &
         & XTBcallback%ewald_qmmm, &
         & XTBcallback%ewald_qmqm, &
         & atomicshift )
    
    XTBcallback%ewald_ene = XTBcallback%ewald_qmmm + XTBcallback%ewald_qmqm
    
    !write(6,*)"shift=",atomicshift
    
  end subroutine xtb_sander_addshift


  subroutine xtb_sander_getEnergy(self, env, qat, qsh, energy)
    implicit none

    !> Instance of the solvation model
    class(XTBCallbackT), intent(inout) :: self
    
    !> Computation environment
    type(TEnvironment), intent(inout) :: env
    
    !> Atomic partial charges
    real(wp), intent(in) :: qat(:)
    
    !> Shell-resolved partial charges
    real(wp), intent(in) :: qsh(:)
    
    !> Total solvation energy
    real(wp), intent(out) :: energy

    energy = XTBcallback%ewald_ene
    
  end subroutine xtb_sander_getEnergy

  
  subroutine xtb_sander_addGradient(self, env, num, xyz, qat, qsh, gradient)
    implicit none
    
    !> Instance of the solvation model
    class(XTBCallbackT), intent(inout) :: self
    
    !> Computation environment
    type(TEnvironment), intent(inout) :: env
    
    !> Atomic numbers
    integer, intent(in) :: num(:)
    
    !> Cartesian coordinates
    real(wp), intent(in) :: xyz(:, :)
    
    !> Atomic partial charges
    real(wp), intent(in) :: qat(:)
    
    !> Shell-resolved partial charges
    real(wp), intent(in) :: qsh(:)
    
    !> Molecular gradient
    real(wp), intent(inout) :: gradient(:, :)
    
  end subroutine xtb_sander_addGradient


  subroutine newXTBIface(self,nqm,atnums,qm_coords_bohr)
    
    use qmmm_module, only: qmmm_nml

    type(XTBIface),intent(inout) :: self
    integer,intent(in) :: nqm
    integer,intent(in) :: atnums(nqm)
    _REAL_,intent(in) :: qm_coords_bohr(3,nqm)

    logical :: exitRun
    
    call read_xtb_nml(self%xtb_nml)
         
    call mctc_init('xtb/sander',11,.true.)
    
    call init(self%env,.true.)
  
    call init(self%mol, atnums, qm_coords_bohr, chrg=real(qmmm_nml%qmcharge,wp) )
    
    call self%env%check(exitRun)
    if (exitRun) then
       call self%env%error("Could not initialize xtb molecule", "sander/xtb_module.F90")
       return
    end if

    if ( trim(self%xtb_nml%qm_level) == "GFN1-XTB" ) then
       call newXTBCalculator(self%env, self%mol, self%calc, method=1)
       call self%env%check(exitRun)
       if (exitRun) then
          call self%env%error("Could not construct GFN1-xTB calculator", "sander/xtb_module.F90")
          return
       end if
    else if ( trim(self%xtb_nml%qm_level) == "GFN2-XTB" ) then
       call newXTBCalculator(self%env, self%mol, self%calc, method=2)
       call self%env%check(exitRun)
       if (exitRun) then
          call self%env%error("Could not construct GFN2-xTB calculator", "sander/xtb_module.F90")
          return
       end if
    else
       write(6,*)"Invalid xtb_nml%qm_level: ",self%xtb_nml%qm_level
       stop
    end if

    self%calc%etemp = max(self%xtb_nml%tfermi,1.d-6)
    self%calc%maxiter = max(self%xtb_nml%maxiter,1)
    self%calc%accuracy = max(1.d-5,min(1.e+3,self%xtb_nml%accuracy))
    
    call newWavefunction(self%env, self%mol, self%calc, self%res)
    call self%env%check(exitRun)
    if (exitRun) then
       call self%env%error("Could not initialize xtb wavefunction", "sander/xtb_module.F90")
       return
    end if

    if ( qmmm_nml%qm_ewald>0 ) then
       self%calc%lSolv = .true.
    end if
         
  end subroutine newXTBIface
  
#endif



  
  
  ! Calculate the XTB qmmm force and energy
  subroutine get_xtb_qmmm_forces( &
       & natom, scaled_mm_charges, &
       & nqmatoms, atnums, qm_coords, &
       & nclatoms, clcoords, &
       & escf, dxyzqm, dxyzcl, scf_mchg )

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

    integer, intent(in) :: natom
    _REAL_, intent(in) :: scaled_mm_charges(natom)
    integer, intent(in) :: nqmatoms
    integer, intent(in) :: atnums(nqmatoms)
    _REAL_, intent(in)  :: qm_coords(3,nqmatoms) 
    integer, intent(in) :: nclatoms
    _REAL_, intent(in)  :: clcoords(4,nclatoms)     ! MM atom coordinates
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    _REAL_, intent(out) :: scf_mchg(nqmatoms)

#ifdef XTB
    ! local variables
    integer :: i,z
    integer :: ier = 0
    logical :: exitRun
    logical, save :: first_call = .true.
    _REAL_ :: qm_coords_bohr(3, nqmatoms)
    _REAL_ :: clcoords_bohr(3, nclatoms)
    _REAL_ :: clcharges(nclatoms)
    _REAL_ :: extPot(nqmatoms)
    _REAL_ :: extPotGrad(3, nqmatoms)
    _REAL_ :: atomCharges(nqmatoms)
    _REAL_ :: ecor
    integer :: mmatnums(nclatoms)
#endif
    
#ifdef MPI
    !_REAL_,allocatable :: buf(:)
    include 'mpif.h'
#endif


#ifndef XTB
    
       call sander_bomb('get_xtb_qmmm_forces','XTB is not enabled', &
            'Check your installation or reconfigure with the XTB support.')
    
#else

    dxyzqm = 0.d0
    dxyzcl = 0.d0
    escf = 0.d0
    scf_mchg = 0.d0
    ecor = 0.d0
    
    !sigma = 0.d0
    !egap = 0.d0

    mmatnums = 1
    qm_coords_bohr(:,:) = qm_coords * AU_PER_ANG
    if (nclatoms > 0) then
       clcoords_bohr(:,:) = clcoords(1:3,:) * AU_PER_ANG
       clcharges(:) = clcoords(4,:)
       do i=1,nclatoms
          mmatnums(i) = max(1,qmmm_struct%qm_mm_pair_atom_numbers(i))
       end do
    end if

    !write(6,*)"nqmatoms",nqmatoms
    !write(6,*)"nclatoms",nclatoms
    !write(6,*)"natom",natom
    !write(6,*)clcharges
    !write(6,*)mmatnums

#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif
       if (XTBstate%xtb_nml%debug) then 
          write(6,*)
          write(6,*) ' XTB input coordinates:'
          do i = 1, size(qm_coords_bohr, dim=2)
             write(6,'(3(f16.8,2x))') qm_coords_bohr(:,i)
          end do
          if ( nclatoms > 0 ) then
             write(6,*)
             write(6,*) ' XTB external point coordinates and charges:'
             do i = 1, size(clcoords_bohr, dim=2)
                write(6,'(4(f16.8,2x))') clcoords_bohr(:,i), clcharges(i)
             end do
          end if
       end if
#ifdef MPI
    end if
#endif

    !
    ! Create interface state if this is the first call
    !

#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif
       if (first_call) then
          call newXTBIFace(XTBstate,nqmatoms,atnums,qm_coords_bohr)
       end if
#ifdef MPI
    end if
#endif

    !
    ! Use all ranks to compute ewald quantities and then
    ! gather them on the first rank
    !

    call QMQMEwComm_Setup(XTBcallback%ewald,first_call,natom,scaled_mm_charges)
    XTBcallback%ewald_ene = 0.d0
    XTBcallback%ewald_qmqm = 0.d0
    XTBcallback%ewald_qmmm = 0.d0

    
#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif
       
       !
       ! Update the QM coordinates
       !

       XTBstate%mol%xyz = qm_coords_bohr


       !
       ! Set the external MM point charges
       !

       if (nclatoms > 0) then

          call XTBstate%calc%pcem%allocate(nclatoms)
          do i=1,nclatoms
             XTBstate%calc%pcem%xyz(:,i) = clcoords_bohr(:,i)
             XTBstate%calc%pcem%q(i) = clcharges(i)
          end do
          if ( XTBstate%xtb_nml%mmhardness > 1.d-6 ) then
             do i=1,nclatoms
                XTBstate%calc%pcem%gam(i) = XTBstate%xtb_nml%mmhardness
             end do
          else if ( XTBstate%xtb_nml%mmhardness > -1.d-6 ) then
             do i=1,nclatoms
                XTBstate%calc%pcem%gam(i) = &
                     & XTBstate%calc%xtbData%coulomb%chemicalHardness(1)
             end do
          else
             do i=1,nclatoms
                z = max(1,mmatnums(i))
                XTBstate%calc%pcem%gam(i) = abs(XTBstate%xtb_nml%mmhardness) * &
                     & XTBstate%calc%xtbData%coulomb%chemicalHardness(z)
             end do
          end if

          XTBstate%calc%pcem%grd = 0.d0
       end if

       !
       ! Perform a single point calculation
       !

       dxyzqm = 0.d0
       escf = 0.d0
       call xtbsinglept(XTBstate%calc, XTBstate%env, XTBstate%mol, &
            & XTBstate%res, XTBcallback, &
            & qmmm_nml%verbosity, .true., &
            & escf, dxyzqm, XTBstate%sigma, XTBstate%egap, XTBstate%spRes)

       ! check SCF convergence; if failed, throw error
       call XTBstate%env%check(exitRun)
       if (exitRun) then
          call XTBstate%env%show('QMMM: Unable to achieve self '// &
               & 'consistency to the tolerances specified')
           call sander_bomb('get_xtb_qmmm_forces','xtb SCF failed', &
               'Check the above message from xtb.')
       end if

       !
       ! Extract the atomic charges
       !

       scf_mchg = XTBstate%res%wfn%q


       !
       ! Extract the energy
       !

       ! Convert the Energy unit from Hartree to (kcal/mol)
       !write(6,*)"escf=",escf," out",QMQMCOR,QMMMCOR
       escf = escf / AU_PER_KCAL

       ! This gets added-in via xtb_sander_getEnergy, so it is unnecessary here
       !
       !ecor = XTBcallback%ewald_ene / AU_PER_KCAL
       !write(6,*)"escf=",escf, "kcal",ecor
       !escf = escf + ecor


       !
       ! Extract the QM gradients
       !

       ! Convert the forces unit from Hartree/Bohr to kcal/mol A  
       dxyzqm(:,:) = dxyzqm(:,:) * AU_PER_ANG / AU_PER_KCAL


       !
       ! Extract the MM gradients
       !

       if (nclatoms > 0 ) then         
          ! Convert the forces unit from Hartree/Bohr to kcal/mol A  
          dxyzcl(:,:) = XTBstate%spRes%pcem%grd * AU_PER_ANG / AU_PER_KCAL
       end if


       if (XTBstate%xtb_nml%debug) then 
          write(6,*)
          write(6,*) ' XTB energy:', escf
          write(6,*) ' XTB forces and mulliken charges'
          do i = 1, size(dxyzqm, dim=2)
             write(6,'(4es14.5)') dxyzqm(:,i), scf_mchg(i)
          end do
          write(6,*)
          write(6,*) ' XTB MM forces'
          do i=1, size(dxyzcl, dim=2)
             write(6,'(3es14.5)') dxyzcl(:,i)
          end do
       end if

       !
       ! bcast the final charges so we can compute the ewald part
       ! of the forces in qm_mm
       !

#ifdef MPI
       call mpi_barrier(qmmm_mpi%commqmmm,ier)
       if ( qmmm_nml%qm_ewald>0 ) then
          if ( qmmm_mpi%numthreads > 1 ) then
             call mpi_bcast(scf_mchg,qmmm_struct%nquant_nlink, &
                  & mpi_double_precision,0,qmmm_mpi%commqmmm,ier)
          end if
       end if
    else
       call mpi_barrier(qmmm_mpi%commqmmm,ier)
       if ( qmmm_nml%qm_ewald>0 ) then
          if ( qmmm_mpi%numthreads > 1 ) then
             call mpi_bcast(scf_mchg,qmmm_struct%nquant_nlink, &
                  & mpi_double_precision,0,qmmm_mpi%commqmmm,ier)
          end if
       end if
    endif
#endif

    first_call=.false.

    call XTBstate%calc%pcem%deallocate()

#endif
    
  end subroutine get_xtb_qmmm_forces



  
#ifdef XTB

  ! Read XTB namelist from mdin file
  subroutine read_xtb_nml(xtb_nml)
#ifdef MPI
    use qmmm_module, only : qmmm_mpi
#endif
    implicit none
    type(xtb_nml_type) :: xtb_nml
    ! local variables
    integer, parameter :: iu_mdin = 5  ! assume mdin file connected to unit 5
    integer :: ierr
    logical :: is_open

    ! namelist variables
    character(len=80) :: qm_level,tmpstr
    _REAL_  :: tfermi,accuracy,mmhardness
    integer :: maxiter
    logical debug

    namelist /xtb/ qm_level, tfermi, accuracy, maxiter, mmhardness, debug


    tmpstr=""

    ! Default namelist variable values
    qm_level = 'GFN2-xTB'
    tfermi = 300.d0
    accuracy = 1.d-3
    maxiter = 250
    mmhardness = 0.d0
    debug = .false.


#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif

       ! Read namelist
       inquire(unit=iu_mdin, opened=is_open)
       if ( .not. is_open) then
          call sander_bomb('read_xtb_settings', &
               'mdin file not connected to unit 5', &
               'Stopping now.')
       end if
       rewind(unit=iu_mdin)
       read(unit=iu_mdin, nml=xtb, iostat=ierr)
       if ( ierr /= 0 ) then
          call sander_bomb('read_xtb_settings', &
               '&xtb namelist read error', &
               'Please check your input.')
       end if

       call str2upper(qm_level,tmpstr)
       qm_level = tmpstr

       !write(6,'(3A)')"qm_level='",trim(qm_level),"'"
       ! Assign values to xtb_nml based on input namelist values
       if ( trim(tmpstr) == 'GFN1-XTB') then
          xtb_nml%qm_level = "GFN1-XTB"
       else if ( trim(tmpstr) == 'GFN2-XTB') then
          xtb_nml%qm_level = "GFN2-XTB"
       else
          call sander_bomb('read_xtb_settings', &
               '&xtb namelist read error', &
               'Please check your qm_level.')
       end if

       xtb_nml%tfermi = tfermi
       xtb_nml%accuracy = accuracy
       xtb_nml%maxiter = maxiter
       xtb_nml%mmhardness = mmhardness
       xtb_nml%debug = debug


#ifdef MPI
    end if
#endif

  end subroutine read_xtb_nml

  
  ! Print the XTB namelist (Debug)
  subroutine print_xtb_nml(self)

    implicit none
    type(xtb_nml_type), intent(in) :: self

    write(6,'(/,a)')     '     ======== XTB settings ======== '
    write(6,'(a,a)')     ' qm_level          : ', trim(self%qm_level)
    write(6,'(a,e10.2)') ' tfermi             : ', self%tfermi
    write(6,'(a,e10.2)') ' accuracy          : ', self%accuracy
    write(6,'(a,i4)')    ' maxiter           : ', self%maxiter
    write(6,'(a,l)')     ' debug             : ', self%debug

  end subroutine print_xtb_nml

  
  subroutine str2upper(str,upper)
    implicit none
    character(len=*), intent(in) :: str
    character(len=*), intent(out) :: upper
    integer :: i

    ! not sure if this will work on all kinds of machines
    ! a more general way to convert to uppercase may be needed
    integer, parameter :: offset = ichar('A') - ichar('a')
    do i=1,len_trim(str)
       if (lge(str(i:i),'a') .and. lle(str(i:i),'z')) then
          upper(i:i) =  char(ichar(str(i:i)) + offset)
       else
          upper(i:i) = str(i:i)
       end if
    end do
  end subroutine str2upper



  subroutine xtbsinglept(self, env, mol, chk, mysolvation, &
       & printlevel, restart, &
       & energy, gradient, sigma, hlgap, results)

    use xtb_scf, only : scf


    !> Source of the generated errors
    character(len=*), parameter :: source = 'xtb_module/xtbsinglept'

    !> Calculator instance
    class(TxTBCalculator), intent(inout) :: self

    !> Computational environment
    type(TEnvironment), intent(inout) :: env

    !> Molecular structure data
    type(TMolecule), intent(inout) :: mol

    !> Wavefunction data
    type(TRestart), intent(inout) :: chk

    class(TSolvation), intent(inout) :: mysolvation

    !> Print level for IO
    integer, intent(in) :: printlevel

    !> Restart from previous results
    logical, intent(in) :: restart

    !> Total energy
    real(wp), intent(out) :: energy

    !> Molecular gradient
    real(wp), intent(out) :: gradient(:, :)

    !> Strain derivatives
    real(wp), intent(out) :: sigma(:, :)

    !> HOMO-LUMO gap
    real(wp), intent(out) :: hlgap

    !> Detailed results
    type(scc_results), intent(out) :: results

    class(TSolvation), allocatable :: solvation

    integer :: i,ich
    integer :: mode_sp_run = 1
    real(wp) :: efix
    logical :: inmol
    logical, parameter :: ccm = .true.
    logical :: exitRun
    character(len=*),parameter :: outfmt = &
         '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

    if ( self%lSolv ) then
       solvation = mysolvation
    end if

    call mol%update

    energy = 0.0_wp
    gradient(:, :) = 0.0_wp
    sigma(:, :) = 0.0_wp
    hlgap = 0.0_wp
    efix = 0.0_wp

    call scf(env,mol,chk%wfn,self%basis,self%pcem,self%xtbData,solvation, &
         &   hlgap,self%etemp,self%maxiter,printlevel,restart,.true., &
         &   self%accuracy,energy,gradient,results)

    call env%check(exitRun)
    if (exitRun) then
       call env%error("Electronic structure method terminated", source)
       return
    end if

    ! ------------------------------------------------------------------------
    !  post processing of gradient and energy

    ! point charge embedding gradient file
    ! if (allocated(set%pcem_grad) .and. self%pcem%n > 0) then
    !    call open_file(ich,set%pcem_grad,'w')
    !    do i=1,self%pcem%n
    !       write(ich,'(3f12.8)')self%pcem%grd(1:3,i)
    !    enddo
    !    call close_file(ich)
    ! endif

    ! ------------------------------------------------------------------------
    !  various external potentials
    !call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
    !call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
    !call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
    !call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
    !call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)

    ! ------------------------------------------------------------------------
    !  fixing of certain atoms
    !energy = energy + efix
    results%e_total = energy
    results%gnorm = norm2(gradient)
    !if (fixset%n.gt.0) then
    !   do i=1, fixset%n
    !      gradient(1:3,fixset%atoms(i))=0
    !   enddo
    !endif

    ! save point charge gradients in results
    if (self%pcem%n > 0) then
       results%pcem = self%pcem
    endif

    if (printlevel.ge.2) then
       ! start with summary header
       !if (.not.set%silent) then
       write(env%unit,'(9x,53(":"))')
       write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
       !endif
       write(env%unit,'(9x,53(":"))')
       write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
       !if (.not.set%silent.and.allocated(self%solvation)) then
       !   write(env%unit,outfmt) "total w/o Gsasa/hb", &
       !      &  results%e_total-results%g_sasa-results%g_hb-results%g_shift, "Eh   "
       !endif
       write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
       write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
       !if (.not.set%silent) then
       !   if (set%verbose) then
       write(env%unit,'(9x,"::",49("."),"::")')
       write(env%unit,outfmt) "HOMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo),  "eV   "
       write(env%unit,outfmt) "LUMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo+1),"eV   "
       !   endif
       write(env%unit,'(9x,"::",49("."),"::")')
       if (self%xtbData%level.eq.2) then
          call print_gfn2_results(env%unit,results,.true.,allocated(self%solvation))
       else if (self%xtbData%level.eq.1) then
          call print_gfn1_results(env%unit,results,.true.,allocated(self%solvation))
       end if
       write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
       write(env%unit,outfmt) "total charge      ", sum(chk%wfn%q), "e    "
       ! if (set%verbose) then
       write(env%unit,'(9x,"::",49("."),"::")')
       write(env%unit,outfmt) "atomisation energy", results%e_atom, "Eh   "
       ! endif
       !endif
       write(env%unit,'(9x,53(":"))')
       write(env%unit,'(a)')
    endif

  end subroutine xtbsinglept


  subroutine print_gfn1_results(iunit,res,verbose,lsolv)
    use xtb_type_data

    character(len=*),parameter :: outfmt = &
         '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

    integer, intent(in) :: iunit ! file handle (usually output_unit=6)
    type(scc_results),    intent(in) :: res
    logical,intent(in) :: verbose,lsolv
    write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
    write(iunit,outfmt) "-> electrostatic  ", res%e_es,   "Eh   "
    !if (lsolv) then
    !write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
    !write(iunit,outfmt) "   -> Gelec       ", res%g_born, "Eh   "
    !write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
    !write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
    !write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
    !endif
    write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
    write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
    write(iunit,outfmt) "halogen bond corr.", res%e_xb,   "Eh   "
  end subroutine print_gfn1_results

  subroutine print_gfn2_results(iunit,res,verbose,lsolv)
    use xtb_type_data

    character(len=*),parameter :: outfmt = &
         '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

    integer, intent(in) :: iunit ! file handle (usually output_unit=6)
    type(scc_results),    intent(in) :: res
    logical,intent(in) :: verbose,lsolv
    write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
    write(iunit,outfmt) "-> isotropic ES   ", res%e_es,   "Eh   "
    write(iunit,outfmt) "-> anisotropic ES ", res%e_aes,  "Eh   "
    write(iunit,outfmt) "-> anisotropic XC ", res%e_axc,  "Eh   "
    write(iunit,outfmt) "-> dispersion     ", res%e_disp, "Eh   "
    !if (lsolv) then
    !write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
    !write(iunit,outfmt) "   -> Gelec       ", res%g_born, "Eh   "
    !write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
    !write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
    !write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
    !endif
    write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
  end subroutine print_gfn2_results


#endif

end module xtb_module

