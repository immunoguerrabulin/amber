#include "../include/dprec.fh"

! Example Installation (15 Dec 2023)
! ----------------------------------
!
! 1. Download DeePMD-kit C library package from GitHub
!
! wget https://github.com/deepmodeling/deepmd-kit/releases/latest/download/libdeepmd_c.tar.gz
! tar xzf libdeepmd_c.tar.gz
!
! The package is compiled with specific versions of TensorFlow and CUDA. The versions may be
! updated over the time. The latest information will be updated in the documentation:
! https://docs.deepmodeling.com/projects/deepmd/en/master/install/install-from-c-library.html
! Go to https://github.com/deepmodeling/deepmd-kit/releases for other versions.
!
! 2. You should now see the following files:
!
! ${DEEPMDHOME}/lib/libdeepmd_c.so
! ${DEEPMDHOME}/include/deepmd.hpp
! ${DEEPMDHOME}/include/c_api.h
!
! 3. Change to the amber build directory
!
! cd /path/to/amber/build
!
! 4. Edit run_cmake to define USE_DEEPMDKIT and CMAKE_PREFIX_PATH; e.g.,
!
!   -DUSE_DEEPMDKIT=ON
!   -DCMAKE_PREFIX_PATH=${DEEPMDHOME}
!
! 5. Compile and install amber
!
! bash ./run_cmake
! make install
!
! 6. To enable GPU support, one need to load the specific version of CUDA Toolkit and
! cuDNN libraries during runtime by setting `LD_LIBRARY_PATH`.

module dprcmod

  use, intrinsic :: iso_c_binding

  implicit none

  public :: dprc_read_mdin
  public :: new_dprc_iface
  public :: cpt_dprc_iface

  private :: get_mdstep
  private :: set_mdstep
  private :: will_write_traj
  private :: calc_stats

  
  integer,parameter,private :: MAXMODELS = 4

  !  _REAL_,parameter,private :: AMBERELE = 18.2223d0
  !  _REAL_,parameter,private :: BOHRS_TO_A = 0.529177249d0
  !  _REAL_,parameter,private :: CODATA08_A_TO_BOHRS = 1. / BOHRS_TO_A;
  !  _REAL_,parameter,private :: CODATA08_AU_TO_KCAL = &
  !       & AMBERELE * AMBERELE * CODATA08_A_TO_BOHRS;

  ! _REAL_,parameter,private :: AU_PER_AMBER_CHARGE = 1. / AMBERELE;
  ! _REAL_,parameter,private :: AU_PER_AMBER_ANGSTROM = 1. / BOHRS_TO_A;
  ! _REAL_,parameter,private :: AU_PER_AMBER_KCAL_PER_MOL = &
  !      & AU_PER_AMBER_CHARGE * AU_PER_AMBER_CHARGE / AU_PER_AMBER_ANGSTROM;
  ! _REAL_,parameter,private :: AU_PER_AMBER_MMPOT = 1. / AU_PER_AMBER_ANGSTROM
  ! _REAL_,parameter,private :: AU_PER_AMBER_FORCE = &
  !      & AU_PER_AMBER_KCAL_PER_MOL / AU_PER_AMBER_ANGSTROM




  
  type dprc_type
     character(len=2560) :: mask = ""
     character(len=2560) :: interfile(MAXMODELS)
     character(len=2560) :: intrafile(MAXMODELS)
     _REAL_ :: rcut = 0.d0
     !logical :: cuda = .FALSE.
     integer :: idprc = 0
     integer,pointer :: tmask(:) => NULL()
     !_REAL_ :: althresh = 0.d0
     ! integer :: alskip = 50
     integer :: curstep = 0
     ! logical :: keepcrd = .false.
     logical :: al = .false.
     logical :: avg = .true.
     !logical :: usemmresnums = .true.
     !logical :: pairwise = .false.
  end type dprc_type

  type(dprc_type),public,save :: dprc_nml


#ifdef DPRC
  
  interface
     subroutine new_dprc( &
          & ninter, interfile, &
          & nintra, intrafile, &
          & ntb, avg &
          ! & pairwise &
#ifdef MPI
          & ,comm &
#endif
          & ) bind(c,name="new_dprc_")
       use, intrinsic :: iso_c_binding
       implicit none
       integer(kind=c_int),intent(in) :: ninter
       type(c_ptr) :: interfile
       integer(kind=c_int),intent(in) :: nintra
       type(c_ptr) :: intrafile
       integer(kind=c_int),intent(in) :: ntb
       integer(kind=c_int),intent(in) :: avg
       !integer(kind=c_int),intent(in) :: pairwise
#ifdef MPI
       integer(kind=c_int),intent(in) :: comm
#endif
       !integer(kind=c_int),intent(in) :: icuda
     end subroutine new_dprc
  end interface


  interface
     subroutine cpt_dprc( &
          & glb_crds, &
          & ml_ene, &
          & glb_frcs, &
          & maxstd, &
          & ntarget, &
          & nat, &
          & idxs, &
          & residxs, &
          & atomnames, &
          & nblist_ucell, &
          & calcstd, &
          & irank, isize ) bind(c,name="cpt_dprc_")
       use, intrinsic :: iso_c_binding
       implicit none
       real(kind=c_double),intent(in)    :: glb_crds
       real(kind=c_double),intent(out)   :: ml_ene
       real(kind=c_double),intent(inout) :: glb_frcs
       real(kind=c_double),intent(out)   :: maxstd
       integer(kind=c_int),intent(in)    :: ntarget
       integer(kind=c_int),intent(in)    :: nat
       integer(kind=c_int),intent(in)    :: idxs
       integer(kind=c_int),intent(in)    :: residxs
       type(c_ptr),dimension(nat) :: atomnames
       real(kind=c_double),intent(in)    :: nblist_ucell
       integer(kind=c_int),intent(in)    :: calcstd
       integer(kind=c_int),intent(in)    :: irank
       integer(kind=c_int),intent(in)    :: isize
     end subroutine cpt_dprc
  end interface



  interface
     subroutine img_dprc( &
          & nat, &
          & glb_crds, &
          & realcell, &
          & recipcell, &
          & ntarget, &
          & itarget, &
          & qmcut, &
          & imm ) bind(c,name="img_dprc_")
       use, intrinsic :: iso_c_binding
       implicit none
       integer(kind=c_int),intent(in)    :: nat
       real(kind=c_double),intent(inout) :: glb_crds
       real(kind=c_double),intent(in)    :: realcell
       real(kind=c_double),intent(in)    :: recipcell
       integer(kind=c_int),intent(in)    :: ntarget
       integer(kind=c_int),intent(in)    :: itarget
       real(kind=c_double),intent(in)    :: qmcut
       integer(kind=c_int),intent(out)   :: imm
     end subroutine img_dprc
  end interface


#endif

  
  
contains



  function calc_stats() result(res)
    implicit none
    logical :: res
    
    res = .false.
    if ( dprc_nml%al ) then
       if ( will_write_traj() ) then
          res = .true.
       end if
    end if
    
  end function calc_stats


  function will_write_traj() result(res)
    use file_io_dat, only : ntwx
    implicit none
#include "../include/md.h"
    logical :: res
    integer :: cstep

    integer :: tstep
    
    res = .false.

    tstep = 1
    if ( isynctraj > 0 ) then
       tstep = 0
    end if
    
    if ( nstlim == 0 .or. ntwx == 1 ) then
       res = .true.
    else if ( ntwx > 1 ) then
       cstep = get_mdstep()
       if ( mod(cstep,ntwx) == tstep .and. cstep > 0 ) then
          res = .true.
       end if
    end if
    
  end function will_write_traj


  function get_mdstep() result(curstep)
    implicit none
#include "../include/md.h"

    integer :: curstep

    curstep = dprc_nml%curstep

  end function get_mdstep



  subroutine set_mdstep()
    implicit none
#include "../include/md.h"
    
    dprc_nml%curstep = dprc_nml%curstep + 1
    if ( nstlim > 0 ) then
       if ( dprc_nml%curstep > nstlim ) then
          if ( irest > 0 ) then
             dprc_nml%curstep = 1
          else
             dprc_nml%curstep = 0
          end if
       end if
    else
       dprc_nml%curstep = 0
    end if
    
  end subroutine set_mdstep



  subroutine new_dprc_iface()
    use, intrinsic :: iso_c_binding

    use findmask, only : atommask
    use memory_module, only : x, ix, ih
    
    implicit none
#include "../include/memory.h"
#include "../include/md.h"
#include "box.h"
#ifdef MPI
#include "parallel.h"
    include 'mpif.h'
    integer :: ier
#endif
    integer :: i
    integer(kind=c_int) :: ninter,nintra
    character(kind=c_char,len=2561),target :: cinters(MAXMODELS)
    character(kind=c_char,len=2561),target :: cintras(MAXMODELS)
    type(c_ptr) :: cinterp(MAXMODELS)
    type(c_ptr) :: cintrap(MAXMODELS)
    integer :: doavg
    
#ifdef MPI
    ier=0
    call mpi_bcast(dprc_nml%idprc, 1, MPI_INTEGER, 0, commsander, ier)
#endif

    if ( dprc_nml%idprc < 1 ) then
       return
    endif

#ifdef MPI
    call mpi_bcast(dprc_nml%avg, 1, MPI_LOGICAL, 0, commsander, ier)
    call mpi_bcast(dprc_nml%al, 1, MPI_LOGICAL, 0, commsander, ier)
    !call mpi_bcast(dprc_nml%alskip, 1, MPI_INTEGER, 0, commsander, ier)
    call mpi_bcast(dprc_nml%mask, 2560, MPI_CHARACTER, 0, commsander, ier)
    !call mpi_bcast(dprc_nml%althresh, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
    do i=1,MAXMODELS
       call mpi_bcast(dprc_nml%interfile(i), 2560, MPI_CHARACTER, 0, commsander, ier)
    end do
    !call mpi_bcast(dprc_nml%intrafile, len(dprc_nml%intrafile), MPI_CHARACTER, 0, commsander, ier)
    call mpi_bcast(dprc_nml%rcut, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
    if ( sanderrank > 0 ) then
       dprc_nml%intrafile = ""
    end if
#endif

    doavg = 0
    if ( dprc_nml%avg ) then
       doavg = 1
    end if
    
    
    allocate( dprc_nml%tmask( natom ) )
    dprc_nml%tmask=0
    dprc_nml%curstep = 99999999

    !write(6,*)sanderrank,dprc_nml%mask
    !write(6,*)sanderrank,natom
    !write(6,*)sanderrank,nres
    !write(6,*)sanderrank,x(lcrd)


    call atommask( natom, nres, 0, ih(m04), ih(m06), &
         & ix(i02), ih(m02), x(lcrd), dprc_nml%mask, dprc_nml%tmask )

    !write(6,*)sanderrank,dprc_nml%tmask

    !write(6,*)"intrafile",dprc_nml%intrafile

    !icuda = 0
    !if ( dprc_nml%cuda ) then
    !   icuda = 1
    !end if

#if defined(MLCUDA) && defined(MPI)
    if ( sanderrank == 0 ) then
#endif

       ninter=0
       nintra=0
       cinters=""
       cintras=""
       do i=1,MAXMODELS
          if ( len_trim(adjustl(dprc_nml%interfile(i))) > 0 ) then
             ninter = ninter+1
          end if
          if ( len_trim(adjustl(dprc_nml%intrafile(i))) > 0 ) then
             nintra = nintra+1
          end if
          cinters(i) = trim(adjustl(dprc_nml%interfile(i))) // c_null_char
          cintras(i) = trim(adjustl(dprc_nml%intrafile(i))) // c_null_char
          cinterp(i) = c_loc( cinters(i) )
          cintrap(i) = c_loc( cintras(i) )
       end do

       !pairwise = 0
       !if ( dprc_nml%pairwise ) then
       !   pairwise = 1
       !end if

       
#ifdef DPRC

       call new_dprc( ninter, cinterp(1), nintra, cintrap(1), &
            & ntb, doavg &
#ifdef MPI
            & ,commsander &
#endif
            & )

#endif
 
#if defined(MLCUDA) && defined(MPI)
    end if
#endif



  end subroutine new_dprc_iface


  


  subroutine cpt_dprc_iface(crds,frcs,ecor)
    use, intrinsic :: iso_c_binding
    use nblist, only : ucell,recip
    
    use findmask, only : atommask
    use memory_module, only : x, ix, ih, residue_pointer
    use qmmm_module, only : get_atomic_number
    implicit none
#include "../include/memory.h"
#include "../include/md.h"
#include "box.h"
#ifdef MPI
#include "parallel.h"
    include 'mpif.h'
#endif
    
    _REAL_, intent(in) :: crds(*)
    _REAL_, intent(inout) :: frcs(*)
    _REAL_, intent(out) :: ecor



    integer,allocatable :: imask(:)
    integer :: ntmask,nimask,nat,i,j,k
    integer,allocatable :: idxs(:)
    integer,allocatable :: residxs(:)
    integer,allocatable :: atomic_nums(:)
    integer,allocatable :: allresnums(:)
    integer :: znum
    character(len=5),allocatable,target :: names(:)
    type(c_ptr),allocatable :: pnames(:)
    _REAL_, allocatable :: imgcrds(:)
    _REAL_, allocatable :: d3grd(:,:)
    _REAL_ :: d3ene
    integer :: irank,isize,ier
    _REAL_ :: maxstd
    integer :: calcstd


    if ( dprc_nml%idprc < 1 ) then
       return
    endif

    maxstd = 0.d0
    
    irank=0
    isize=1
#ifdef MPI
    irank = sanderrank
    isize = sandersize
    ier=0
#endif


    call set_mdstep()
    
    calcstd = 0
    if ( dprc_nml%al ) then
       if ( dprc_nml%avg ) then
          calcstd = 1
       else if ( calc_stats() ) then
          calcstd = 1
       end if
    end if

    
    ntmask = 0
    do i=1,natom
       if ( dprc_nml%tmask(i) > 0 ) then
          ntmask = ntmask + 1
       end if
    end do

    !write(6,*)"N target",ntmask
    allocate( idxs(ntmask) )
    j=0
    do i=1,natom
       if ( dprc_nml%tmask(i) > 0 ) then
          j = j + 1
          idxs(j) = i-1
          !write(6,*)"Target ",j,idxs(j)
       end if
    end do
    
    allocate(imgcrds(3*natom))
    do i=1,3*natom
       imgcrds(i) = crds(i)
    end do
    
    allocate( imask( natom ) )
    imask = 0

    
    allocate( allresnums(natom) )
    do i=1,nres
       allresnums( residue_pointer(i) : residue_pointer(i+1)-1 ) = i
    end do
    

#ifdef DPRC
    call img_dprc( natom, imgcrds(1), ucell(1,1), recip(1,1), ntmask, &
         & idxs(1), dprc_nml%rcut, imask(1) )

    !
    ! if any atom in a residue was marked, then mark all atoms within
    ! the residue
    !
    do i=1,natom
       if ( imask(i) > 0 ) then
          j = allresnums(i)
          !write(6,*)i,j,residue_pointer(j), residue_pointer(j+1)
          imask( residue_pointer(j) : residue_pointer(j+1)-1 ) = 1
       end if
    end do
    
#endif

!    if ( len_trim(dprc_nml%intermask) > 0 ) then
!       call atommask( natom, nres, 0, ih(m04), ih(m06), &
!            & ix(i02), ih(m02), imgcrds(1), dprc_nml%intermask, imask(1) )
!    end if

    do i=1,natom
       if ( dprc_nml%tmask(i) > 0 .and. imask(i) > 0 ) then
          imask(i) = 0
       end if
    end do

    nimask = 0
    do i=1,natom
       if ( imask(i) > 0 ) then
          nimask = nimask + 1
       end if
    end do

    !write(6,*)"nimask",nimask
    nat = ntmask + nimask

    deallocate( idxs )

    allocate( idxs(nat) )
    allocate( residxs(nat) )
    
    j=0
    do i=1,natom
       if ( dprc_nml%tmask(i) > 0 ) then
          j = j + 1
          idxs(j) = i-1
          residxs(j) = 0
       end if
    end do
    do i=1,natom
       if ( imask(i) > 0 ) then
          j = j + 1
          idxs(j) = i-1
          !if ( dprc_nml%usemmresnums ) then
          residxs(j) = allresnums(i)
          !else
          !   residxs(j) = 1
          !end if
       end if
    end do


    
    allocate( atomic_nums(ntmask) )
    
    ! ix(i100) is a flag that tells if the ATOMIC_NUMBERS section is in the
    ! prmtop. If it's set, get the atomic numbers from those read
    ! directly from the prmtop; otherwise use the name/mass to guess

    if (ix(i100).eq.0) then
       do j=1,ntmask
          i = idxs(j)+1
          call get_atomic_number(ih(m04+i-1),x(lmass+i-1),atomic_nums(j))
       end do
    else
       do j=1,ntmask
          i = idxs(j)+1
          atomic_nums(j) = ix(i100+i)
       end do
    endif

    allocate( names(nat) )
    allocate( pnames(nat) )
    do i=1,ntmask
       call get_atomic_symbol_from_number(atomic_nums(i),names(i))
    end do

    
    !do i=1,ntmask
    !   write(6,*)"target name ",i,idxs(i),names(i)
    !end do


    do i=1,nimask
       j=idxs(ntmask+i)+1
       names(ntmask+i) = ih(m06+j-1)
    end do

    if (ix(i100).eq.0) then
       do j=1,nimask
          i = idxs(ntmask+j)+1
          if ( trim(ih(m06+i-1)) .ne. "OW" .and. trim(ih(m06+i-1)) .ne. "HW" ) then
              call get_atomic_number(ih(m04+i-1),x(lmass+i-1),znum)
              call get_atomic_symbol_from_number(znum,names(ntmask+j))
              names(ntmask+j) = "m" // trim(names(ntmask+j))
          end if
       end do
    else
       do j=1,nimask
          i = idxs(ntmask+j)+1
          if ( trim(ih(m06+i-1)) .ne. "OW" .and. trim(ih(m06+i-1)) .ne. "HW" ) then
              znum = ix(i100+i)
              call get_atomic_symbol_from_number(znum,names(ntmask+j))
              names(ntmask+j) = "m" // trim(names(ntmask+j))
          end if
       end do
    endif


    ! do i=1,nat
    !   write(6,'(2i5,2xa4,i5)')i,idxs(i),names(i),residxs(i)
    ! end do
    

    do i=1,nat
       names(i) = trim(adjustl(names(i)))//char(0)
       pnames(i) = c_loc( names(i) )
    end do    

#ifdef DPRC
#if defined(MLCUDA) && defined(MPI)
    if ( irank == 0 ) then
       call cpt_dprc( imgcrds(1), ecor, frcs(1), maxstd, ntmask, nat, &
            & idxs(1), residxs(1), pnames(1), ucell(1,1), calcstd, 0, 1 )
    end if
#else
    call cpt_dprc( imgcrds(1), ecor, frcs(1), maxstd, ntmask, nat, &
         & idxs(1), residxs(1), pnames(1), ucell(1,1), calcstd, irank, isize )
#endif
#endif




    if ( calc_stats() ) then
       if ( irank == 0 ) then
          write(6,'(2A,F11.5,A)')"Active learning frame", &
               & " written with max. frc. std.: ", &
               & maxstd, " kcal/mol/A"
       end if
    end if
    


  end subroutine cpt_dprc_iface


  
  subroutine dprc_read_mdin(mdin_lun)
    implicit none
    integer, intent(in) :: mdin_lun
    call dprc_type_read_mdin(mdin_lun,dprc_nml)
  end subroutine dprc_read_mdin


  
  subroutine dprc_type_read_mdin(mdin_lun,data)
 
    use findmask, only : atommask
    use memory_module, only : x, ix, ih
    use file_io_dat, only : ntwx
    
    implicit none
#include "../include/memory.h"
#include "../include/md.h"
    
    integer, intent(in) :: mdin_lun
    type(dprc_type),intent(inout) :: data
    integer :: stat,i
    
    character(len=2560) :: mask
    character(len=2560) :: interfile(MAXMODELS)
    character(len=2560) :: intrafile(MAXMODELS)
    
    !logical :: cuda

    _REAL_ :: rcut !, althresh
    integer :: idprc !,alskip
    logical :: avg
    !logical :: usemmresnums
    !logical :: pairwise
    integer :: nmodel

    namelist /dprc/ &
         & idprc, &
         & mask, &
         & rcut, &
         & interfile, &
         & intrafile, &
         & avg

    idprc=0
    mask=""
    rcut=0.d0
    interfile=""
    intrafile=""
    !althresh = 0.d0
    !alskip=ntwx
    avg = .false.

    rewind(mdin_lun)
    read(unit=mdin_lun,nml=dprc,iostat=stat)

    if (stat /= 0) then
       write(6,'(A)') 'Error reading namelist &dprc.'
       call mexit(6,1)
    end if

#ifndef DPRC
    if ( idprc > 0 ) then
       write(6,*)"Cannot use idprc>0 within &dprc because amber was not"
       write(6,*)"configured to use deepmdkit. To use this feature, you"
       write(6,*)"need to download and install deepmdkit and reconfigure"
       write(6,*)"amber with the -DDPRC cmake option"
       call mexit(6,1)
    end if
#endif

    data%idprc=idprc
    data%mask=trim(adjustl(mask))
    data%rcut = rcut
    nmodel = 0
    do i=1,MAXMODELS
       if ( len_trim(adjustl(interfile(i))) > 0 .or. &
            & len_trim(adjustl(intrafile(i))) > 0 ) then
          nmodel = nmodel + 1
       end if
       data%interfile(i)=trim(adjustl(interfile(i)))
       data%intrafile(i)=trim(adjustl(intrafile(i)))
    end do


    if ( nmodel > 1 ) then
       data%al = .true.
       data%avg = avg       
    else
       data%al = .false.
       data%avg = .false.
    end if

    if ( data%idprc > 0 ) then
       if ( data%al .and. .not. data%avg ) then
          if ( isynctraj <= 0 ) then
             write(6,"(5a)")"When dprc uses multiple models, the standard", &
                  & " deviation of the atomic forces will be computed for each", &
                  & " saved frame (at a frequency of ntwx). This requires", &
                  & " setting isynctraj=1 in the &cntrl section. The default", &
                  & " (isynctraj=0) would save the coordinates of the NEXT frame."
             call mexit(6,1)
          end if
       end if
    end if
    
  end subroutine dprc_type_read_mdin


  subroutine get_atomic_symbol_from_number(z,sym)
    implicit none
    integer,intent(in) :: z
    character(len=*),intent(out) :: sym
  
    character(len=2),parameter :: pt(86) = (/ &
         & "H ", "He", & ! 2
         & "Li", "Be", & ! 2
         & "B ","C ","N ","O ","F ","Ne", & ! 6
         & "Na", "Mg", & ! 2
         & "Al","Si","P ","S ","Cl","Ar", & ! 6
         & "K ","Ca", & ! 2
         & "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", & ! 10
         & "Ga","Ge","As","Se","Br","Kr", & 
         & "Rb","Sr", &
         & "Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", &
         & "In","Sn","Sb","Te","I ","Xe", &
         & "Cs","Ba", &
         & "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd", &
         & "Tb","Dy","Ho","Er","Tm","Yb","Lu", &
         & "Hf","Ta","W ","Re","Os","It","Pt","Au","Hg",&
         & "Tl","Pb","Bi","Po","At","Rn" /)

    sym="XX"
    if ( z > 0 .and. z <= size(pt) ) then
       sym=pt(z)
    end if
    
  end subroutine get_atomic_symbol_from_number

  
end module dprcmod
