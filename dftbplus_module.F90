#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

! QM/MM interface to DFTBPLUS library
!
! Written by Timothy J. Giese
! Date: December 2023
!
! Original draft written by Yu-Ching Hu, Andreas W. Goetz
! Date: August 2016
!       March 2018
!
! Example Installation (13 Dec 2023)
! ----------------------------------
!
! 1. Download DFTB+ from https://github.com/dftbplus/dftbplus
!
! mkdir ~/software
! cd ~/software
! git clone https://github.com/dftbplus/dftbplus.git
! cd dftbplus
!
! 2. Configure dftbplus with API and shared library support.
!
! export DFTBPLUSHOME=${HOME}/software/dftbplus/local
! export FC=gfortran
! export CC=gcc
! cmake -B_build \
!        -DCMAKE_INSTALL_PREFIX=${DFTBPLUSHOME} \
!        -DBLAS_LIBRARY=-lopenblas \
!        -DWITH_API=TRUE \
!        -DWITH_SDFTD3=TRUE \
!        -DBUILD_SHARED_LIBS=TRUE \
!        -B_build .
!
! 3. Build and install dftbplus
!
! cmake --build _build -- -j
! cmake --install _build
!
!
! 4. You should now see the following files:
!
!    ~/software/dftbplus/local/lib64/libdftbplus.so
!    ~/software/dftbplus/local/lib64/cmake/dftbplus/dftbplus-config.cmake
!    ~/software/dftbplus/local/lib64/cmake/dftbplus/dftbplus-targets.cmake
!    ~/software/dftbplus/local/include/dftbplus/modfiles/
!
! 5. Change to the amber build directory
!
! cd /path/to/amber/build
!
! 6. Edit run_cmake to define USE_DFTBPLUS and CMAKE_PREFIX_PATH; e.g.,
!
!   -DUSE_DFTBPLUS=ON
!   -DCMAKE_PREFIX_PATH="${DFTBPLUSHOME}"
!
!    Note: To define multiple paths in CMAKE_PREFIX_PATH, use ';'
!    as the delimiter rather than ':'.
!
! 7. Compile and install amber
!
! bash ./run_cmake
! make install
!
!

module dftbplus_module

#ifdef DFTBPLUS
  use, intrinsic :: iso_c_binding
  use dftbplus
#endif
  use, intrinsic ::  iso_fortran_env, only : output_unit

  !use constants, only: A_TO_BOHRS, AU_TO_KCAL
  use constants, only: CODATA08_A_TO_BOHRS, CODATA08_AU_TO_KCAL
  use qm_ewald_comm
  
  implicit none
  private

  public :: get_dftbplus_qmmm_forces
  
#ifdef DFTBPLUS

  
  ! DFTBPLUS namelist contains user definable settings for DFTBPLUS
  type dftbplus_nml_type
     
     character(len=80) :: qm_level = "DFTB3"
     _REAL_  :: scftol ! SCF convergence threshold
     _REAL_  :: tfermi ! Fermi smearing, temperature in Kelvin
     integer :: maxiter 
     logical :: hcorrection
     logical :: thirdorderfull
     logical :: orbitalresolvedscc
     logical :: debug ! print debug output
     logical :: silent ! suppress all DFTB+ output
     logical :: reqconv ! if true, then program exits when convergence not met

     character(len=80) :: mixer = "" ! BROYDEN, DIIS, or SIMPLE
     _REAL_ :: broyden_mixpar = 0.2d0
     _REAL_ :: broyden_invjacwt = 0.01d0
     _REAL_ :: broyden_minwt = 1.0d0
     _REAL_ :: broyden_maxwt = 1.d5
     _REAL_ :: broyden_wtfact = 1.d-2

     _REAL_ :: diis_mixpar = 0.2d0
     integer :: diis_gen = 6
     logical :: diis_fromstart = .true.

     _REAL_ :: simple_mixpar = 0.05d0

     !logical :: d3disp = .false.
     _REAL_ :: d3a1 = -1.d0 ! 0.746d0
     _REAL_ :: d3a2 = -1.d0 ! 4.191d0
     _REAL_ :: d3s8 = -1.d0 ! 3.209d0

     ! logical :: d4disp = .false.
     ! _REAL_ :: d4a1 = 0.5467502d0
     ! _REAL_ :: d4a2 = 4.4955068d0
     ! _REAL_ :: d4s8 = 0.4727337d0

  end type dftbplus_nml_type

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

  
  
  type, extends(TQDepExtPotGen) :: DFTBPLUSCallbackT
     type(QMQMEwCommT) :: ewald
     _REAL_ :: ewald_qmqm
     _REAL_ :: ewald_qmmm
     _REAL_ :: ewald_ene
   contains
     procedure :: getExternalPot => dftbplus_sander_getExternalPot
     procedure :: getExternalPotGrad => dftbplus_sander_getExternalPotGrad
  end type DFTBPLUSCallbackT


  type DFTBPLUSIface
     type(dftbplus_nml_type) :: dftbplus_nml
     type(TDftbPlus) :: calc
  end type dftbplusiface

  
  type(DFTBPLUSIface),save :: DFTBPLUSstate
  type(DFTBPLUSCallbackT),save :: DFTBPLUScallback

#endif

  
contains

  
#ifdef DFTBPLUS


  subroutine dftbplus_sander_getExternalPot( &
       & this, chargePerAtom, &
       & chargePerShell, extPotAtom, extPotShell)
    
    implicit none
    
    !> Instance.
    class(DFTBPLUSCallbackT), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(wp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(wp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(wp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(wp), intent(out) :: extPotShell(:,:)

    real(wp),allocatable :: scf_mchg(:)

    integer :: i,n

    n = size(chargePerAtom)
    allocate( scf_mchg(n) )
    scf_mchg = - chargePerAtom

    extPotAtom = 0.d0
    extPotShell = 0.d0
    DFTBPLUScallback%ewald_qmmm = 0.d0
    DFTBPLUScallback%ewald_qmqm = 0.d0
    
    call QMQMEwComm_CalcEne( DFTBPLUScallback%ewald, scf_mchg, &
         & DFTBPLUScallback%ewald_qmmm, &
         & DFTBPLUScallback%ewald_qmqm, &
         & extPotAtom )
    
    !DFTBPLUScallback%ewald_ene = DFTBPLUScallback%ewald_qmmm + DFTBPLUScallback%ewald_qmqm
    DFTBPLUScallback%ewald_ene = - DFTBPLUScallback%ewald_qmqm
    
    extPotAtom = -extPotAtom
    
  end subroutine dftbplus_sander_getExternalPot

  
  subroutine dftbplus_sander_getExternalPotGrad( &
       & this, chargePerAtom, chargePerShell, extPotGrad )
    
    implicit none
    
    !> Class instance.
    class(DFTBPLUSCallbackT), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(wp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(wp), intent(in) :: chargePerShell(:,:)

    !> Gradient of the potential at each atom. Shape: [3, nAtom].
    real(wp), intent(out) :: extPotGrad(:,:)

    extPotGrad = 0.d0
    
  end subroutine dftbplus_sander_getExternalPotGrad



  subroutine newDFTBPLUSIface(self, nqmatoms, ntype, type_id, qm_coords)
    
    use qmmm_module, only: qmmm_nml
    use qmmm_module, only : qmmm_struct
    use ElementOrbitalIndex, only : ElementSymbol

    implicit none
    
    type(DFTBPLUSIface),intent(inout) :: self
    integer,intent(in) :: nqmatoms,ntype
    integer,intent(in) :: type_id(ntype)
    _REAL_,intent(in) :: qm_coords(3,nqmatoms) ! bohr

    type(TDftbPlusInput) :: input
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, &
                            pAnalysis, pFilling, pFermi, pHub, pHcorrection
    type(fnode), pointer :: pParserOpts, pOpts, pDamping, pMixer, pMixerOpt
    type(fnode), pointer :: pDisp, pD3 , pDispDamping, pZeroDamp, pH5
    character(10), allocatable :: speciesNames(:)
    character(1024), allocatable :: slako_files(:,:)
    integer :: i, j

    character(1), parameter :: maxAngNames(4) = ["s", "p", "d", "f"]
    integer :: maxAng3OB(118)
    integer, parameter :: iout = 6
    integer :: devNull

    logical :: isdftb
    
    call read_dftbplus_nml(self%dftbplus_nml)

    maxAng3OB = 2
    maxang3OB(1) = 1
    maxang3OB(15) = 3
    maxang3OB(16) = 3
    maxang3OB(17) = 3
    maxang3OB(30) = 3
    maxang3OB(35) = 3
    maxang3OB(53) = 3


    isdftb = .false.
    if ( self%dftbplus_nml%qm_level(1:4) == "DFTB" ) then
       isdftb = .true.
    end if
    
    allocate(speciesNames(ntype))
    do i = 1, ntype
       speciesNames(i) = elementSymbol(type_id(i))
    end do

    allocate(slako_files(ntype, ntype))
    call set_skf_filenames(speciesNames, slako_files, self%dftbplus_nml)

    ! Note: setting the global standard output to /dev/null 
    !        will also suppress run-time error messages
    devNull=6
    if ( qmmm_nml%verbosity < 1 ) then
       open(newunit=devNull, file="/dev/null", action="write")
    endif
    
    !if (self%dftbplus_nml%silent) then
    !   open(newunit=devNull, file="/dev/null", action="write")
    call TDftbPlus_init(self%calc, outputUnit=devNull)
    !else
    !   call TDftbPlus_init(self%calc)
    !end if

    ! You should provide the dftb_in.hsd and skfiles as found in the
    ! $DFTBPLUSHOME/test/prog/dftb+/non-scc/Si_2/ folder
    call self%calc%getEmptyInput(input)

    call input%getRootNode(pRoot)

    call setChild(pRoot, "Geometry", pGeo)

    call setChildValue(pGeo, "TypeNames", speciesNames)
    call setChildValue(pGeo, "TypesAndCoordinates", &
         & reshape(qmmm_struct%qm_atom_type, [1, nqmatoms]), &
         & qm_coords)

    call setChild(pRoot, "Hamiltonian", pHam)
    
    if ( isdftb ) then
    
       call setChild(pHam, "Dftb", pDftb)
       call setChildValue(pDftb, "Scc", .true.)
       call setChildValue(pDftb, "SccTolerance", self%dftbplus_nml%scftol)
       call setChildValue(pDftb, "MaxSCCIterations", self%dftbplus_nml%maxiter)
       call setChildValue(pDftb, "ConvergentSCCOnly",self%dftbplus_nml%reqconv)
       call setChildValue(pDftb, "Charge", real(qmmm_nml%qmcharge, wp))
       call setChildValue(pDftb, "ReadInitialCharges", .false.)
       call setChild(pDftb, "Filling", pFilling)
       call setChild(pFilling, "Fermi", pFermi)
       call setChildValue(pFermi, "Temperature", self%dftbplus_nml%tfermi, &
            & modifier="K")
       call setChildValue(pDftb, "ThirdOrderFull", self%dftbplus_nml%thirdorderfull)
       if (self%dftbplus_nml%thirdorderfull) then
          call setChild(pDftb, "HubbardDerivs", pHub)
          call set_uhubbder(type_id, speciesNames, pHub)
       end if
       call setChildValue(pDftb, "OrbitalResolvedSCC", &
            & self%dftbplus_nml%orbitalresolvedscc)
       call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
       do i = 1, ntype
          j = getMaxAngFromSlakoFile(slako_files(i, i)) + 1
          if ( self%dftbplus_nml%thirdorderfull ) then
             j = maxang3OB(type_id(i))
          end if
          call setChildValue(pMaxAng, speciesNames(i), &
               & maxAngNames(j))
       end do
       
       call setChild(pDftb, "SlaterKosterFiles", pSlakos)
       do i = 1, ntype
          do j = 1, ntype
             call setChildValue(pSlakos, &
                  & trim(speciesNames(i))//"-"//trim(speciesNames(j)), &
                  & trim(slako_files(i, j)))
          end do
       end do

       
       if (self%dftbplus_nml%hcorrection) then
          call setChild(pDftb, "HCorrection", pHcorrection)
          call setChild(pHcorrection, "Damping", pDamping)
          call setChildValue(pDamping, "Exponent", real(4.00, wp))
       end if
       
       if ( trim(self%dftbplus_nml%qm_level) == "DFTB3-D3H5" ) then
          call setChild(pDftb, "HCorrection", pHcorrection)
          call setChild(pHcorrection, "H5", pH5)
          call setChild(pDftb, "Dispersion", pDisp)
          call setChild(pDisp, "DftD3", pD3)
          call setChild(pD3, "Damping", pDispDamping)
          call setChild(pDispDamping, "ZeroDamping", pZeroDamp)
          call setChildValue(pZeroDamp, "sr6", 1.25_wp)
          call setChildValue(pZeroDamp, "alpha6", 29.61_wp)
          call setChildValue(pD3, "s6", 1.0_wp)
          call setChildValue(pD3, "s8", self%dftbplus_nml%d3s8)
          call setChildValue(pD3, "HHRepulsion", .true.)
       else if ( trim(self%dftbplus_nml%qm_level) == "DFTB3-D42B" ) then
          call setChild(pDftb, "Dispersion", pDisp)
          call setChild(pDisp, "DftD4", pD3)
          call setChildValue(pD3,"s6",1.d0)
          call setChildValue(pD3,"s9",0.d0)
          call setChildValue(pD3,"s8",self%dftbplus_nml%d3s8)
          call setChildValue(pD3,"a1",self%dftbplus_nml%d3a1)
          call setChildValue(pD3,"a2",self%dftbplus_nml%d3a2)
       else if ( trim(self%dftbplus_nml%qm_level) == "DFTB3-D43B" ) then
          call setChild(pDftb, "Dispersion", pDisp)
          call setChild(pDisp, "DftD4", pD3)
          call setChildValue(pD3,"s6",1.d0)
          call setChildValue(pD3,"s9",1.d0)
          call setChildValue(pD3,"s8",self%dftbplus_nml%d3s8)
          call setChildValue(pD3,"a1",self%dftbplus_nml%d3a1)
          call setChildValue(pD3,"a2",self%dftbplus_nml%d3a2)
       else if ( trim(self%dftbplus_nml%qm_level) == "DFTB3-D3" ) then
          ! Dispersion/SimpleDftD3/a1
          call setChild(pDftb, "Dispersion", pDisp)
          call setChild(pDisp, "SimpleDftD3", pD3)
          call setChildValue(pD3,"s6",1.d0)
          call setChildValue(pD3,"s8",self%dftbplus_nml%d3s8)
          call setChildValue(pD3,"a1",self%dftbplus_nml%d3a1)
          call setChildValue(pD3,"a2",self%dftbplus_nml%d3a2)
       end if
       
       ! end if
       
       if ( trim(self%dftbplus_nml%mixer) == "BROYDEN" ) then
          call setChild(pDftb, "Mixer", pMixer)
          call setChild(pMixer, "Broyden", pMixerOpt)
          call setChildValue(pMixerOpt, "MixingParameter",self%dftbplus_nml%broyden_mixpar)
          call setChildValue(pMixerOpt, "InverseJacobiWeight",self%dftbplus_nml%broyden_invjacwt)
          call setChildValue(pMixerOpt, "MinimalWeight",self%dftbplus_nml%broyden_minwt)
          call setChildValue(pMixerOpt, "MaximalWeight",self%dftbplus_nml%broyden_maxwt)
          call setChildValue(pMixerOpt, "WeightFactor",self%dftbplus_nml%broyden_wtfact)
       else if ( trim(self%dftbplus_nml%mixer) == "DIIS" ) then
          call setChild(pDftb, "Mixer", pMixer)
          call setChild(pMixer, "diis", pMixerOpt)
          call setChildValue(pMixerOpt, "InitMixingParameter",self%dftbplus_nml%diis_mixpar)
          call setChildValue(pMixerOpt, "Generations",self%dftbplus_nml%diis_gen)
          call setChildValue(pMixerOpt, "UseFromStart",self%dftbplus_nml%diis_fromstart)
       else if ( trim(self%dftbplus_nml%mixer) == "SIMPLE" ) then
          call setChild(pDftb, "Mixer", pMixer)
          call setChild(pMixer, "Simple", pMixerOpt)
          call setChildValue(pMixerOpt, "MixingParameter",self%dftbplus_nml%simple_mixpar)
       end if

          
    else

       call setChild(pHam, "xTB", pDftb)
       call setChildValue(pDftb, "Method", trim(self%dftbplus_nml%qm_level))
       call setChildValue(pDftb, "SccTolerance", self%dftbplus_nml%scftol)
       call setChildValue(pDftb, "MaxIterations", self%dftbplus_nml%maxiter)
       call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
       call setChildValue(pDftb, "Charge", real(qmmm_nml%qmcharge, wp))
       call setChild(pDftb, "Filling", pFilling)
       call setChild(pFilling, "Fermi", pFermi)
       call setChildValue(pFermi, "Temperature", self%dftbplus_nml%tfermi, &
            & modifier="K")
       
    end if
    
    call setChild(pRoot, "Analysis", pAnalysis)
    call setChildValue(pAnalysis, "CalculateForces", .true.)
    call setChildValue(pAnalysis, "WriteBandOut", .false.)

    call setChild(pRoot, "ParserOptions", pParserOpts)
    call setChildValue(pParserOpts, "ParserVersion", 5)


    call setChild(pRoot, "Options", pOpts)
    if ( .not. self%dftbplus_nml%silent ) then
       call setChildValue(pOpts, "WriteDetailedOut", .true.)
       call setChildValue(pParserOpts, "WriteHSDInput", .true.)
       call setChildValue(pOpts, "WriteCharges", .true.)
    else
       call setChildValue(pOpts, "WriteDetailedOut", .false.)
       call setChildValue(pParserOpts, "WriteHSDInput", .false.)
       call setChildValue(pOpts, "WriteCharges", .false.)
    end if

    
    if (self%dftbplus_nml%debug) then 
       write(output_unit,"(A)") 'Input tree in HSD format:'
       call dumpHsd(input%hsdTree, output_unit)
    end if

    call self%calc%setupCalculator(input)

    deallocate(slako_files)
    deallocate(speciesNames)

    if (self%dftbplus_nml%debug) then 
      call print_dftbplus_nml(self%dftbplus_nml)
      write(6,*) ''
      write(6,*) '<<< leaving subroutine dftbplus_input_setting'
    end if

  end subroutine newDFTBPLUSIface
  
#endif



  
  
  ! Calculate the DFTBPLUS qmmm force and energy
  subroutine get_dftbplus_qmmm_forces( &
       & natom, scaled_mm_charges, &
       & nqmatoms, ntype, type_id, qm_coords, &
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
    integer, intent(in) :: nqmatoms, ntype
    integer, intent(in) :: type_id(ntype)

    _REAL_, intent(in)  :: qm_coords(3,nqmatoms) 
    integer, intent(in) :: nclatoms
    _REAL_, intent(in)  :: clcoords(4,nclatoms)     ! MM atom coordinates
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    _REAL_, intent(out) :: scf_mchg(nqmatoms)

#ifdef DFTBPLUS
    ! local variables
    integer :: i,z
    integer :: ier = 0
    logical, save :: first_call = .true.
    _REAL_ :: qm_coords_bohr(3, nqmatoms)
    _REAL_ :: clcoords_bohr(3, nclatoms)
    _REAL_ :: clcharges(nclatoms)
    !_REAL_ :: extPot(nqmatoms)
    !_REAL_ :: extPotGrad(3, nqmatoms)
    !_REAL_ :: atomCharges(nqmatoms)
#endif
    
#ifdef MPI
    !_REAL_,allocatable :: buf(:)
    include 'mpif.h'
#endif


#ifndef DFTBPLUS
    
       call sander_bomb('get_dftbplus_qmmm_forces','DFTBPLUS is not enabled', &
            'Check your installation or reconfigure with the DFTBPLUS support.')
    
#else

    dxyzqm = 0.d0
    dxyzcl = 0.d0
    escf = 0.d0
    scf_mchg = 0.d0
    

    qm_coords_bohr(:,:) = qm_coords * AU_PER_ANG
    if (nclatoms > 0) then
       clcoords_bohr(1:3,:) = clcoords(1:3,:) * AU_PER_ANG
       clcharges(:) = clcoords(4,:)
       !do i=1,nclatoms
       !   mmatnums(i) = max(1,qmmm_struct%qm_mm_pair_atom_numbers(i))
       !end do
    end if

    !write(6,*)"nqmatoms",nqmatoms
    !write(6,*)"nclatoms",nclatoms
    !write(6,*)"natom",natom
    !write(6,*)clcharges
    !write(6,*)mmatnums

#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif
       if (qmmm_nml%verbosity > 4) then 
          write(6,*)
          write(6,*) ' DFTBPLUS input coordinates:'
          do i = 1, size(qm_coords_bohr, dim=2)
             write(6,'(3(f16.8,2x))') qm_coords_bohr(:,i)
          end do
          if ( nclatoms > 0 ) then
             write(6,*)
             write(6,*) ' DFTBPLUS external point coordinates and charges:'
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
          call newDFTBPLUSIFace(DFTBPLUSstate,nqmatoms,ntype,type_id,qm_coords_bohr)
          if ( qmmm_nml%qm_ewald > 0 ) then
             call DFTBPLUSstate%calc%setQDepExtPotGen(DFTBPLUScallback)
          end if
       end if
#ifdef MPI
    end if
#endif

    !
    ! Use all ranks to compute ewald quantities and then
    ! gather them on the first rank
    !

    call QMQMEwComm_Setup(DFTBPLUScallback%ewald,first_call,natom,scaled_mm_charges)
    DFTBPLUScallback%ewald_ene = 0.d0
    DFTBPLUScallback%ewald_qmqm = 0.d0
    DFTBPLUScallback%ewald_qmmm = 0.d0

    
#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif
       
       !
       ! Update the QM coordinates
       !

      call DFTBPLUSstate%calc%setGeometry(qm_coords_bohr)

       !
       ! Set the external MM point charges
       !

       if (nclatoms > 0) then
         call DFTBPLUSstate%calc%setExternalCharges(clcoords_bohr, clcharges)
       end if

       !
       ! Perform a single point calculation
       !

       dxyzqm = 0.d0
       escf = 0.d0
       call DFTBPLUSstate%calc%getEnergy(escf)
       call DFTBPLUSstate%calc%getGradients(dxyzqm)
       call DFTBPLUSstate%calc%getGrossCharges(scf_mchg)

       escf = (escf + DFTBPLUScallback%ewald_ene) / AU_PER_KCAL
       dxyzqm(:,:) = dxyzqm(:,:) * AU_PER_ANG / AU_PER_KCAL


       !
       ! Extract the MM gradients
       !

       if (nclatoms > 0 ) then         
          call DFTBPLUSstate%calc%getExtChargeGradients(dxyzcl)
          dxyzcl(:,:) = dxyzcl(:,:) * AU_PER_ANG / AU_PER_KCAL
       end if


       if (qmmm_nml%verbosity > 4) then 
          write(6,*)
          write(6,*) ' DFTBPLUS energy:', escf
          write(6,*) ' DFTBPLUS forces and mulliken charges'
          do i = 1, size(dxyzqm, dim=2)
             write(6,'(4es14.5)') dxyzqm(:,i), scf_mchg(i)
          end do
          write(6,*)
          write(6,*) ' DFTBPLUS MM forces'
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

#endif
    
  end subroutine get_dftbplus_qmmm_forces



  
#ifdef DFTBPLUS

  

  ! Read DFTBPLUS namelist from mdin file
  subroutine read_dftbplus_nml(dftbplus_nml)
#ifdef MPI
    use qmmm_module, only : qmmm_mpi
#endif
    implicit none
    type(dftbplus_nml_type) :: dftbplus_nml
    ! local variables
    integer, parameter :: iu_mdin = 5  ! assume mdin file connected to unit 5
    integer :: ierr
    logical :: is_open

    ! namelist variables
    character(len=80) :: qm_level,tmpstr
    _REAL_  :: scftol, tfermi
    integer :: hcorrection, maxiter
    logical :: silent, debug, orbitalresolvedscc
    character(len=80) :: mixer
    _REAL_ :: broyden_mixpar,broyden_invjacwt
    _REAL_ :: broyden_minwt,broyden_maxwt
    _REAL_ :: broyden_wtfact,diis_mixpar
    _REAL_ :: simple_mixpar
    integer :: diis_gen
    logical :: diis_fromstart
    logical :: reqconv
    _REAL_ :: d3s8,d3a1,d3a2
    
    namelist /dftbplus/ qm_level, scftol, tfermi, &
                        hcorrection,  &
                        silent, debug, maxiter, &
                        mixer, &
                        broyden_mixpar,broyden_invjacwt, &
                        broyden_minwt,broyden_maxwt, &
                        broyden_wtfact,diis_mixpar, &
                        simple_mixpar, &
                        diis_gen, &
                        diis_fromstart, &
                        reqconv, &
                        d3s8,d3a1,d3a2

    tmpstr=""
    
    ! Default namelist variable values
    qm_level = 'DFTB2'
    scftol = 1.0d-7
    tfermi = 0.0d0 ! Kelvin
    hcorrection = -1
    orbitalresolvedscc  = .false.
    silent = .true. ! suppress *all* output
    debug = .false.
    maxiter = 250
    reqconv = .true.
    
    mixer = "BROYDEN"
    broyden_mixpar = 0.2d0
    broyden_invjacwt = 0.01d0
    broyden_minwt = 1.0d0
    broyden_maxwt = 1.d5
    broyden_wtfact = 1.d-2
    diis_mixpar = 0.2d0
    diis_gen = 6
    diis_fromstart = .true.

    ! Recommended DFTB3/3ob-D3 parameters listed in Appendix G
    ! of the DFTB+ user documentation
    ! 
    ! _REAL_ :: d3a1 = 0.746d0
    ! _REAL_ :: d3a2 = 4.191d0
    ! _REAL_ :: d3s8 = 3.209d0

    ! logical :: d4disp = .false.
    ! _REAL_ :: d4a1 = 0.5467502d0
    ! _REAL_ :: d4a2 = 4.4955068d0
    ! _REAL_ :: d4s8 = 0.4727337d0

    ! values less than 0 will cause appropriate
    ! default values to be used, depending on qm_level
    d3a1 = -1.d0
    d3a2 = -1.d0
    d3s8 = -1.d0

    
#ifdef MPI
    if ( qmmm_mpi%commqmmm_master ) then
#endif
       
       ! Read namelist
       inquire(unit=iu_mdin, opened=is_open)
       if ( .not. is_open) then
          call sander_bomb('read_dftbplus_settings', &
               'mdin file not connected to unit 5', &
               'Stopping now.')
       end if
       rewind(unit=iu_mdin)
       ierr=0
       read(unit=iu_mdin, nml=dftbplus, iostat=ierr)
       if ( ierr /= 0 ) then
          call sander_bomb('read_dfbplus_settings', &
               '&dftbplus namelist read error', &
               'Please check your input.')
       end if

       tmpstr=""
       call str2upper(qm_level,tmpstr)
       qm_level = tmpstr

       dftbplus_nml%qm_level = qm_level
       dftbplus_nml%d3s8 = d3s8
       dftbplus_nml%d3a1 = d3a1
       dftbplus_nml%d3a2 = d3a2
       dftbplus_nml%scftol = scftol
       
       !write(6,'(3A)')"qm_level='",trim(qm_level),"'"
       ! Assign values to dftbplus_nml based on input namelist values
       if ( trim(qm_level) == 'DFTB3') then
          dftbplus_nml%thirdorderfull = .true.
          dftbplus_nml%hcorrection = .true.
       else if ( trim(qm_level) == 'DFTB3D3' .or. &
            & trim(qm_level) == 'DFTB3-D3' ) then
          dftbplus_nml%qm_level = "DFTB3-D3"
          dftbplus_nml%thirdorderfull = .true.
          dftbplus_nml%hcorrection = .true.
          if ( dftbplus_nml%d3a1 < 0 ) then
             dftbplus_nml%d3a1 = 0.746d0
          end if
          if ( dftbplus_nml%d3a2 < 0 ) then
             dftbplus_nml%d3a2 = 4.191d0
          end if
          if ( dftbplus_nml%d3s8 < 0 ) then
             dftbplus_nml%d3s8 = 3.209d0
          end if
       else if ( trim(qm_level) == 'DFTB3D42B' .or. &
            & trim(qm_level) == 'DFTB3-D42B' ) then
#ifndef WITH_SDFTD3
          write(6,*)"To use qm_level=",trim(qm_level),", you will need to"
          write(6,*)"reconfigure and recompile DFTB+ with -DWITH_SDFTB=TRUE."
          write(6,*)"You will then need to reconfigure and recompile AmberTools"
          write(6,*)"with -DUSE_DFTBPLUS=ON."
          call mexit(6,1)
#endif
          dftbplus_nml%qm_level = "DFTB3-D42B"
          dftbplus_nml%thirdorderfull = .true.
          dftbplus_nml%hcorrection = .true.

          if ( dftbplus_nml%d3a1 < 0 ) then
             dftbplus_nml%d3a1 = 0.5467502d0
          end if
          if ( dftbplus_nml%d3a2 < 0 ) then
             dftbplus_nml%d3a2 = 4.4955068d0
          end if
          if ( dftbplus_nml%d3s8 < 0 ) then
             dftbplus_nml%d3s8 = 0.4727337d0
          end if
       else if ( trim(qm_level) == 'DFTB3D43B' .or. &
            & trim(qm_level) == 'DFTB3-D43B' ) then
#ifndef WITH_SDFTD3
          write(6,*)"To use qm_level=",trim(qm_level),", you will need to"
          write(6,*)"reconfigure and recompile DFTB+ with -DWITH_SDFTB=TRUE."
          write(6,*)"You will then need to reconfigure and recompile AmberTools"
          write(6,*)"with -DUSE_DFTBPLUS=ON."
          call mexit(6,1)
#endif
          dftbplus_nml%qm_level = "DFTB3-D43B"
          dftbplus_nml%thirdorderfull = .true.
          dftbplus_nml%hcorrection = .true.

          if ( dftbplus_nml%d3a1 < 0 ) then
             dftbplus_nml%d3a1 = 0.5523240d0
          end if
          if ( dftbplus_nml%d3a2 < 0 ) then
             dftbplus_nml%d3a2 = 4.3537076d0
          end if
          if ( dftbplus_nml%d3s8 < 0 ) then
             dftbplus_nml%d3s8 = 0.6635015d0
          end if
       else if ( trim(qm_level) == 'DFTB3D3H5' .or. &
            & trim(qm_level) == 'DFTB3-D3H5' ) then
#ifndef WITH_SDFTD3
          write(6,*)"To use qm_level=",trim(qm_level),", you will need to"
          write(6,*)"reconfigure and recompile DFTB+ with -DWITH_SDFTB=TRUE."
          write(6,*)"You will then need to reconfigure and recompile AmberTools"
          write(6,*)"with -DUSE_DFTBPLUS=ON."
          call mexit(6,1)
#endif
          dftbplus_nml%qm_level = "DFTB3-D3H5"
          dftbplus_nml%thirdorderfull = .true.
          dftbplus_nml%hcorrection = .false.
          if ( dftbplus_nml%d3s8 < 0 ) then
             dftbplus_nml%d3s8 = 0.49d0
          end if
       else if ( trim(qm_level) == 'DFTB2') then
          dftbplus_nml%thirdorderfull = .false.
          dftbplus_nml%hcorrection = .false.
       else
          call sander_bomb('read_dfbplus_settings', &
               '&dftbplus namelist read error', &
               'Please check your qm_level.')
       end if


       if ( hcorrection == 0) then
          dftbplus_nml%hcorrection = .false. 
       else if ( hcorrection == 1 ) then
          dftbplus_nml%hcorrection = .true. 
       else if ( hcorrection == -1 ) then
          continue ! user has not set, default value has been applied already
       else
          call sander_bomb('read_dfbplus_settings', &
               '&dftbplus hcorrection read error', &
               'Please check your input.')
       end if

       dftbplus_nml%orbitalresolvedscc = orbitalresolvedscc
       dftbplus_nml%debug = debug
       dftbplus_nml%silent = silent
       dftbplus_nml%maxiter = max(1,maxiter)

       tmpstr=""
       call str2upper(mixer,tmpstr)
       dftbplus_nml%mixer=tmpstr
       dftbplus_nml%broyden_mixpar = broyden_mixpar
       dftbplus_nml%broyden_invjacwt = broyden_invjacwt
       dftbplus_nml%broyden_minwt = broyden_minwt
       dftbplus_nml%broyden_maxwt = broyden_maxwt
       dftbplus_nml%broyden_wtfact = broyden_wtfact
       dftbplus_nml%diis_mixpar = diis_mixpar
       dftbplus_nml%diis_gen = diis_gen
       dftbplus_nml%diis_fromstart = diis_fromstart
       dftbplus_nml%simple_mixpar = simple_mixpar
       dftbplus_nml%reqconv = reqconv
       
       if (dftbplus_nml%debug) then 
          write(6,*) '<<< leaving subroutine read_dftbplus_nml'
       end if

#ifdef MPI
    end if
#endif
       
  end subroutine read_dftbplus_nml


  
  ! Define Slater-Koster filenames based on atom types
  subroutine set_skf_filenames(speciesNames, slako_files, dftbplus_nml)
    
    implicit none
    character(*), intent(in) :: speciesNames(:)
    character(len=1024), intent(out) :: slako_files(:,:)
    type(dftbplus_nml_type), intent(in) :: dftbplus_nml

    ! local variables
    integer :: i, j, ntype
    character(len=1024) :: skroot
    character(*), parameter :: skf = ".skf"

    if (dftbplus_nml%debug) then 
      write(6,*) ''
      write(6,*) '>>> entered subroutine set_skf_filenames'
    end if

    ! Define the skroot: $AMBERHOME/dat/slko
    call getenv('AMBERHOME', skroot)
    
    ! For 3rd order, we use skf files from 3ob-3-1
    if (dftbplus_nml%thirdorderfull) then     
      if (dftbplus_nml%debug) then 
        write(6,*) 'qm_level is DFTB3/3ob'
      end if  
      
      skroot = TRIM(skroot) // '/dat/slko/3ob-3-1/'
      
      if (dftbplus_nml%debug) then 
        write(6,*) 'The skf file directory is: ', trim(skroot)
      end if
    ! Default setting will use skf files from mio-1-1 
    else
      if (dftbplus_nml%debug) then 
        write(6,*) 'Default qm_level is DFTB2/mio'
      end if

      skroot = TRIM(skroot) // '/dat/slko/mio-1-1/'
      
      if (dftbplus_nml%debug) then 
        write(6,*) 'The skf file directory is: ', trim(skroot)
      end if
    end if  

    ntype = size(speciesNames, dim=1)
    do i = 1, ntype
       do j = 1, ntype
          slako_files(i,j)=TRIM(skroot)//TRIM(speciesNames(i))// &
          & "-"//TRIM(speciesNames(j))//skf
       end do
    end do    
    
    if (dftbplus_nml%debug) then 
      write(6,*) '<<< leaving subroutine set_skf_filenames'
    end if

  end subroutine set_skf_filenames


  
  ! Set Hubbard derivative for each atom type
  subroutine set_uhubbder(type_id, speciesNames, pHub)
    implicit none
    integer, intent(in) :: type_id(:)
    character(*), intent(in) :: speciesNames(:)
    type(fNode), pointer, intent(in) :: pHub
    
    ! local variables
    _REAL_ :: hubbarDer(118)    ! All elements
    integer :: i

    ! Initialize all elements' Hubbard Derivs value
    ! >> The data is taken from www.dftb.org/parameters/download/3ob/3ob-3-1/
    hubbarDer(:)  = 0.0d0
    hubbarDer(1)  = -0.1857d0    ! H
    hubbarDer(6)  = -0.1492d0    ! C 
    hubbarDer(7)  = -0.1535d0    ! N 
    hubbarDer(8)  = -0.1575d0    ! O 
    hubbarDer(9)  = -0.1623d0    ! F 
    hubbarDer(11) = -0.0454d0    ! Na 
    hubbarDer(12) = -0.02d0      ! Mg
    hubbarDer(15) = -0.14d0      ! P 
    hubbarDer(16) = -0.11d0      ! S 
    hubbarDer(17) = -0.0697d0    ! Cl 
    hubbarDer(19) = -0.0339d0    ! K 
    hubbarDer(20) = -0.0340d0    ! Ca 
    hubbarDer(30) = -0.03d0      ! Zn 
    hubbarDer(35) = -0.0573d0    ! Br 
    hubbarDer(53) = -0.0433d0    ! I 

    do i = 1, size(speciesNames, dim=1)
       call setChildValue(pHub, trim(speciesNames(i)), hubbarDer(type_id(i)))
    end do

  end subroutine set_uhubbder


    ! maxang     = "p"
    ! maxang(1)  = "s" ! H
    ! maxang(15) = "d" ! P
    ! maxang(16) = "d" ! S
    ! maxang(17) = "d" ! Cl
    ! maxang(30) = "d" ! Zn
    ! maxang(35) = "d" ! Br
    ! maxang(53) = "d" ! I
    
  

  ! Print the DFTBPLUS namelist (Debug)
  subroutine print_dftbplus_nml(self)

    implicit none
    type(dftbplus_nml_type), intent(in) :: self
    
    write(6,'(/,a)')     '     ======== DFTBPLUS settings ======== '
    write(6,'(a,a)')     ' qm_level          : ', trim(self%qm_level)
    write(6,'(a,e10.2)') ' scftol            : ', self%scftol
    write(6,'(a,e10.2)') ' tfermi            : ', self%tfermi
    write(6,'(a,i6)')    ' maxiter           : ', self%maxiter
    write(6,'(a,l)')     ' hcorrection       : ', self%hcorrection
    write(6,'(a,l)')     ' thirdorderfull    : ', self%thirdorderfull
    write(6,'(a,l)')     ' debug             : ', self%debug
    write(6,'(a,l)')     ' silent            : ', self%silent
    write(6,'(a,l)')     ' reqconv           : ', self%reqconv
    write(6,'(a,a)')     ' mixer             : ', trim(self%mixer)
    write(6,'(a,e10.2)') ' broyden_mixpar    : ', self%broyden_mixpar
    write(6,'(a,e10.2)') ' broyden_invjacwt  : ', self%broyden_invjacwt
    write(6,'(a,e10.2)') ' broyden_minwt     : ', self%broyden_minwt
    write(6,'(a,e10.2)') ' broyden_maxwt     : ', self%broyden_maxwt
    write(6,'(a,e10.2)') ' broyden_wtfact    : ', self%broyden_wtfact
    write(6,'(a,e10.2)') ' diis_mixpar       : ', self%diis_mixpar
    write(6,'(a,i3)')    ' diis_gen          : ', self%diis_gen
    write(6,'(a,l)')     ' diis_fromstart    : ', self%diis_fromstart
    write(6,'(a,e10.2)') ' simple_mixpar     : ', self%simple_mixpar
    write(6,'(a,f9.5)')  ' d3s8              : ', self%d3s8
    write(6,'(a,f9.5)')  ' d3a1              : ', self%d3a1
    write(6,'(a,f9.5)')  ' d3a2              : ', self%d3a2
  end subroutine print_dftbplus_nml

  
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

#endif

end module dftbplus_module

