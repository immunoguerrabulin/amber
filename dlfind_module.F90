#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

module dlfind_module

  !use, intrinsic :: iso_c_binding
  use iso_c_binding
  use constants, only: CODATA08_A_TO_BOHRS, CODATA08_AU_TO_KCAL
  use state, only : state_rec

  implicit none
  private

  public :: run_dlfind

  _REAL_,parameter :: AU_PER_ANG = CODATA08_A_TO_BOHRS
  _REAL_,parameter :: AU_PER_KCAL = 1.d0 / CODATA08_AU_TO_KCAL
  _REAL_,parameter :: AU_PER_AMBER_FRC = AU_PER_KCAL / AU_PER_ANG
  _REAL_,parameter :: AU_PER_AMU = 1.82288848554095d+3

  integer,parameter :: MASKLEN = 2560
  
  type dlfind_nml_type
     !
     ! The restart file containing the coordinates used
     ! to select atoms from Amber-mask strings.
     ! If it is "", then the input coordinates (read with -c)
     ! are used.
     !
     ! string, default=""
     !
     character(len=MASKLEN) :: refcrd = ""
     !
     ! The Amber-mask which selects the movable atoms
     !
     ! string, default="@*"
     !
     character(len=MASKLEN) :: active = "@*"
     !
     !
     ! The coordinate representation. Must be either:
     !    CART : Cartesian
     !    MASS : Mass-weighted Cartesians
     !    DLC  : Delocalized coordinates (e.g., redundant internal coordinates)
     !    HDLC : Hybrid delocalized coordinates (inter-residues CART,
     !           and intra-residue DLC)
     !    TC   : Total connection (N*(N-1)/2 bonds)
     !    HDLC-TC : Hybrid delocalized/Total connection (inter-residue CART,
     !           and intra-residue TC)
     !
     !  string, default="HDLC"
     !
     character(len=20)   :: crdrep = "HDLC"
     !
     ! Optimization algorithm. Must be either:
     !    LBFGS : Minimum search with limited memory Broyden-Fletcher-
     !            Goldfarb-Shanno
     !    PRFO  : Minimum or TS search with partitioned rational function
     !            optimization
     !    CG    : Minimum search with conjugate gradient
     !    SD    : Minimum search with steepest descent
     !    NR    : Minimum or TS search using Newton-Raphson/quasi-Newton
     !            optimiser
     !    HESS  : Do not perform an optimization. Instead, calculate
     !            a Hessian and write it to hess_out
     !
     !  string, default="LBFGS"
     !
     character(len=20)   :: optalg = "LBFGS"
     !
     ! Trust radius approach. Must be either:
     !    CONST : Constant trust radius (specified with maxstep)
     !    ENE   : Trust radius based on the energy change
     !    GRAD  : Trust radius based on the projected gradient
     !
     ! string, default="CONST"
     !
     character(len=20)   :: trustrad = "CONST"
     !
     ! Initial Hessian. Must be either:
     !   ONEPOINT : Crude finite diff. approx. from single displacements
     !   TWOPOINT : Finite diff. approx. from symmetric displacements
     !   DIAGONAL : Only the diagonal elements are approximated from ONEPOINT
     !   IDENTITY : Identity matrix
     !
     ! string, default="IDENTITY"
     !
     character(len=20)   :: hessini = "IDENTITY"
     !
     ! Hessian update method. Must be either:
     !   NONE   : No update
     !   POWELL : Powell update
     !   BOFILL : Bofill update
     !   BFGS   : BFGS update
     !
     ! string, default="BOFILL"
     !
     character(len=20)   :: hessupd = "BOFILL"
     !
     ! Maximum number of Hessian updates.
     !
     ! integer, default = 5
     !
     integer :: maxupdate = 5
     !
     ! Finite difference displacement when calculating the Hessian.
     !
     ! float, atomic units, default=5.d-3
     !
     _REAL_ :: delta = 5.d-3
     !
     !
     ! The number of previous LBFGS Hessians to store
     ! when performing LBFGS optimizations.
     !
     ! integer, default=5
     !
     integer :: lbfgsmem = 5
     !
     ! Verbosity. Must be either:
     !   0 : Minimal output
     !   1 : Medium output
     !   2 : Full output
     !   3 : Debug output
     !
     ! integer, default=2
     !
     integer :: verbosity = 2
     !
     ! Convergence tolerance on max gradient component.
     !
     ! float, atomic units, default=4.5d-4
     !
     _REAL_ :: tol = 4.5d-4
     !
     ! Convergence criterion on max energy change.
     !
     ! float, atomic units, default=1.d-6
     !
     _REAL_ :: tole = 1.d-6
     !
     ! Maximum number of optimization cycles
     !
     ! integer, default=100
     !
     integer :: maxcycle = 1000
     !
     ! Maximum stepsize in internal coordinates
     !
     ! float, default=0.5
     !
     _REAL_ :: maxstep = 0.5d0
     !
     ! Do not update Hessian if predicted step size is less than minstep
     !
     ! float, default=1.d-7
     !
     _REAL_ :: minstep = 1.d-7
     !
     ! A constant factor use to scale the step size
     !
     ! float, default=1.d0
     !
     _REAL_ :: scalestep = 1.d0
     !
     ! Perform a TS search using the "dimer method"
     !
     ! logical, default=.false.
     !
     logical :: dimer = .false.
     !
     ! Dimer mode gives detailed control of the rotation line search
     !  0 : All rotations covered by the optimizer. Requires 2 energy
     !      evaluations per iteration.
     !  1 : Rotations by linesearch; No extrapolation of gradient.
     !      Two energy calculations are used per rotation. Requires
     !      at least 2 energy evaluations per iteration.
     !  2 : Rotations by linesearch; Extrapolation of gradient.
     !      One energy calculation is done per iteration, the other
     !      one is interpolated. Requires at least 2 energy evaluations
     !      per iteration.
     !
     ! integer, default=1
     !
     integer :: dimer_mode = 1
     !
     ! Input restart file containing the TS mode.
     ! The file must exist when dimer=.true.
     !
     ! string, default="tsmode_inp.rst7"
     !
     !character(len=MASKLEN) :: tsmode_inp = "tsmode_inp.rst7"
     !
     ! Output restart file containing the TS mode.
     !
     ! string, default="tsmode_out.rst7"
     !
     !character(len=MASKLEN) :: tsmode_out = "tsmode_out.rst7"
     !
     !_REAL_ :: dimer_delta
     !
     ! Angle tolerance in the dimer rotations
     !
     ! float, deg, default=2.0
     !
     _REAL_ :: dimer_tolrot = 2.0d0
     !
     ! Maximum number of rotations in each dimer step
     !
     ! integer, default=20
     !
     integer :: dimer_maxrot = 20
     !
     ! How/If nudge elastic band should be performed. Valid options are:
     !   NO      = Do not perform NEB (default)
     !   FREE    = Minimize end-states
     !   FROZEN  = Freeze end-states
     !   PERP    = Optimize end-states in the direction
     !             perpendicular to the path
     !
     ! string, default="NO"
     !
     character(LEN=20) :: neb = "NO"
     !
     ! This is a list of restart filenames.
     !
     ! When performing NEB, coords2 is a series of images (up to 12)
     ! that describe the path. When NEB is initialized it will
     ! interpolate neb_nimage images from this path.  The last
     ! supplied image should be the product state. For a linear
     ! interpolation, you only need to provides coords(1) = "product.rst7"
     ! The reactant state structure is the input coordinates, read
     ! on the command line with -c
     !
     ! When performing the dimer method, only coords2(1) is
     ! required: the initial guess at the transition state mode.
     ! That is, the approximate TS structure displaced by some
     ! amount in the direction of the product state.  The TS
     ! structure is read from the command line with -c
     !
     ! list of string, default=""
     !
     character(LEN=MASKLEN) :: coords2(12) = ""
     !
     ! The number of images to optimize. This should be >= the
     ! number of input images
     !
     ! integer, default = 1
     !
     integer :: neb_nimage = 1
     !
     ! Threshold scale factor for spawning climbing image
     ! (expressed as a multiple of the gradient tolerance for convergence)
     ! Setting neb_climb_test = neb_freeze_test = 0 may help
     ! if convergence issues are encountered during the later stage
     ! of the optimization.
     !
     ! float, default=3.0
     !
     _REAL_ :: neb_climb_test = 3.0d0

     !
     ! Threshold scale factor for freezing NEB images
     ! (expressed as a multiple of the gradient tolerance for convergence)
     !
     ! float, default=1.0
     !
     _REAL_ :: neb_freeze_test = 1.d0
     !
     !
     ! The NEB force constant
     !
     ! float, au, default=0.01
     !
     _REAL_ :: nebk = 0.01d0
     !
     ! If true, then interpolate the intermediate images
     ! using the specified coordinate system to create
     ! initial geometries, but perform the optimization
     ! in Cartesians.
     !
     ! logical, default=.true.
     !
     logical :: neb_cart = .true.
     !
     !
     ! Weight selection for Dimer and NEB methods.
     ! Only these atoms will contribute to the dimer rotations
     ! or NEB image distances.
     ! If the mask is empty, then all active atoms are selected.
     !
     ! string, default=""
     !
     character(len=MASKLEN) :: wtmask = ""
     !
     ! The weight values are 0 for all
     ! active atoms not in wtmask.
     ! The QM atoms within wtmask will be
     ! given a weight of 1.0.
     ! The MM atoms within wtmask will be
     ! given the weight specified by wtmm
     !
     ! float, default=1.0
     !
     _REAL_ :: wtmm = 1.0d0
     !
     ! The netcdf Hessian file used to store all (or part)
     ! of the Hessian. Only used when optalg="HESS".
     !
     ! character, default = "hess_out.nc"
     !
     character(len=MASKLEN) :: hess_out = "hess_out.nc"
     !
     ! If optalg="HESS", then a Hessian is calculated
     ! and written to hess_out using a special netcdf
     ! convention. hess_iseg and hess_nseg break the
     ! Hessian calculation into segments.
     ! The columns of the Hessian will is split into
     ! hess_nseg segments, and the current sander
     ! instance is only responsible for calculating
     ! segment hess_iseg.  If hess_nseg > 1, then
     ! each segment is written to a different netcdf
     ! output file: hess_out // "." // hess_iseg
     ! and these separate outputs can later be merged
     ! into a single hessian file.
     !
     ! One can concatenate netcdf files using:
     !    ncecat --rec_apn hess_out.nc.00* -O hess_out.nc
     ! where ncecat is part of the nco suite of utilities
     ! https://github.com/nco/nco.git
     !
     ! integer, default=1 (1 <= hess_iseg <= hess_nseg)
     !
     integer :: hess_iseg = 1
     !
     ! integer, default=1 (1 <= hess_nseg)
     !
     integer :: hess_nseg = 1
     !
     ! If true, then do not terminate the program
     ! when the system contains shaked bonds
     !
     ! logical, default=F
     !
     logical :: allow_shake = .false.
     !
     ! Disang file containing constraint coordinate definitions
     ! Only bond, angle, and dihedral constraints are supported.
     ! Biasing potentials are ignored.
     !
     character(LEN=MASKLEN) :: disang = ""
     !
     ! If a constraint crosses a residue boundary, merge the residues
     !
     ! logical, default=.false.
     !
     logical :: merge_residues = .false.
  end type dlfind_nml_type



  type active_data_type
     !
     ! Number of active atoms, determined from "active" mask
     ! This excludes 
     !
     integer :: nactive = 0
     !
     ! The 3N array of atomic coordinates solely used to
     ! select atoms from Amber-mask strings
     !
     _REAL_,pointer :: refcrd(:) => NULL()
     !
     ! A list of nactive integers. map(i) returns the atom
     ! index in the full list of coordinates
     !
     integer,pointer :: map(:) => null()
     !
     ! A list of nactive integers; the atomic numbers of
     ! the selected atoms
     !
     integer,pointer :: atnums(:) => null()
     !
     ! A list of nactive floats; the masses (au) of the
     ! selected atoms
     !
     _REAL_,pointer  :: atmass(:) => null()
     !
     ! A list of nactive integers; the residue numbers of
     ! the selected atoms
     !
     integer,pointer :: residx(:) => null()
     !
     ! A list of nactive weights; these are used in the
     ! dimer and neb methods to calculate distance
     ! metrics between images
     !
     _REAL_,pointer  :: wts(:) => null()
     !
     ! Scratch space for the nonbond list
     !
     integer,pointer :: ipairs(:) => null()
     !
     ! A flag indicating if the force routine needs to
     ! perform setup operations
     !
     logical :: qsetup = .false.
     !
     ! A list of 2*npairs integers; the atom pairs
     ! (indexed by the petite list of selected atoms)
     ! denoting shaked bonds
     !
     integer,pointer :: shaked_pairs(:) => null()
     !
     ! Constraint data size (5,ncon)
     ! The first element is 1, 2, or 3, denoting bond, angle, and torsion
     ! the remaining 4 elements are atom indexes
     !
     integer,pointer :: condata(:,:) => null()
     !
     !
     ! The number of extra coordinate sets given to
     ! dlfind. This is 0 for minimizations.
     ! This should be 1 for the dimer method (the extra
     ! coordinates are displaced from the input coordinates
     ! to describe the transition state mode).
     ! The should be 1 (or larger) for the neb method.
     ! When performing neb, the last set of extra coordinates
     ! should be the product state.
     !
     integer :: nframes = 0
     !
     ! A list of (3,nactive,nframes) floats; these are
     ! the extra coordinate sets given to dlfind
     !
     _REAL_,pointer :: frames(:,:,:) => null()
     !
     ! The method for calculating the Hessian.
     ! 2 = central finite difference
     ! 1 = forward finite difference
     ! 0 = error (in principle, this could be analytic)
     !
     integer :: hesstype = 2
     !
     ! Flag indicating if the optization is a search
     ! for a transition state. This lets us know
     ! how to write the output restart files when
     ! dlfind provides coordinates.
     !
     logical :: is_ts_search = .false.
     !
     ! Flag indicating if this is a hessian calculation
     ! rather than an optimization. This lets us know
     ! if we should even bother calling dlfind.
     !
     logical :: is_hess = .false.
     !
     ! Flag indicating if this is a neb calculation
     !
     logical :: is_neb = .false.
     !
     ! The dlfind interface packs floats and integers
     ! into huge arrays in a partiular order.
     ! We kind of have to determine this twice:
     ! once to figure out the storage requirements,
     ! and then a second time to pack the arrays.
     ! nspec and nvar2 are the storage requirements
     ! for the integers and floats.
     !
     integer :: nspec = 0
     integer :: nvar2 = 0
     !
     ! The number of energy evaluations. This is
     ! only used when printing the standard energy
     ! output to mdout.
     !
     integer :: neval = 0
     !
     ! The results from the latest energy evaluation
     !
     type(state_rec) :: ene
     !
     ! prefix used for writing restart files
     !
     character(LEN=MASKLEN) :: base
     !
     ! extension used for writing restart files
     !
     character(LEN=MASKLEN) :: ext
     
  end type active_data_type


  
  type(dlfind_nml_type),save :: dlfind_nml
  type(active_data_type),save :: active_data
  
  
contains




  subroutine read_dlfind_nml(options)
    
    implicit none
    
    type(dlfind_nml_type),intent(out) :: options
    
    integer, parameter :: iu_mdin = 5  ! assume mdin file connected to unit 5

    character(len=MASKLEN) :: refcrd
    character(len=MASKLEN) :: active
    character(len=20)   :: crdrep
    character(len=20)   :: optalg
    character(len=20)   :: trustrad
    character(len=20)   :: hessini
    character(len=20)   :: hessupd
    integer :: maxupdate
    _REAL_  :: delta
    integer :: lbfgsmem
    integer :: verbosity
    _REAL_  :: tol
    _REAL_  :: tole
    integer :: maxcycle
    _REAL_  :: maxstep
    _REAL_  :: minstep
    _REAL_  :: scalestep
    logical :: dimer
    integer :: dimer_mode
    _REAL_  :: dimer_tolrot
    integer :: dimer_maxrot
    !character(len=MASKLEN) :: tsmode_inp
    !character(len=MASKLEN) :: tsmode_out
    character(len=MASKLEN) :: wtmask
    _REAL_  :: wtmm

    character(len=MASKLEN) :: hess_out
    integer :: hess_iseg, hess_nseg
    logical :: allow_shake, merge_residues

    character(len=20) :: neb
    character(len=MASKLEN) :: coords2(12)
    integer :: neb_nimage
    _REAL_ :: neb_climb_test, neb_freeze_test
    _REAL_ :: nebk
    logical :: neb_cart
    character(len=MASKLEN) :: disang

    integer :: ierr

    namelist /dlfind/ &
         & refcrd, active, crdrep, optalg, trustrad, &
         & hessini, hessupd, maxupdate, delta, &
         & lbfgsmem, verbosity, &
         & tol, tole, maxcycle, &
         & maxstep, minstep, scalestep, &
         & dimer, dimer_mode, &
         & dimer_tolrot, dimer_maxrot, &
         !& tsmode_inp, &
         !& tsmode_out, &
         & wtmask, wtmm, &
         & hess_out, hess_iseg, hess_nseg, &
         & allow_shake, &
         & neb, coords2, nebk, &
         & neb_climb_test, neb_freeze_test, &
         & neb_nimage, neb_cart, &
         & disang,merge_residues

    refcrd=""
    active="@*"
    crdrep="HDLC"
    optalg="LBFGS"
    trustrad="CONST"
    hessini = "IDENTITY"
    hessupd = "BOFILL"
    maxupdate = 5
    delta = 5.d-3
    lbfgsmem = 5
    verbosity = 2
    tol = 4.5d-4
    tole = 1.d-6
    maxcycle = 1000
    maxstep = 0.5d0
    minstep = 1.d-7
    scalestep = 1.d0
    dimer = .false.
    dimer_mode = 1
    dimer_tolrot = 2.0d0
    dimer_maxrot = 20
    !tsmode_inp = "tsmode_inp.rst7"
    !tsmode_out = "tsmode_out.rst7"
    wtmask = ""
    wtmm = 1.0d0
    hess_out = "hess_out.nc"
    hess_iseg = 1
    hess_nseg = 1
    allow_shake = .false.

    neb = "NO"
    coords2 = ""
    neb_nimage = 1
    neb_climb_test = 3.d0
    neb_freeze_test = 1.d0
    nebk = 0.01d0
    neb_cart = .true.
    disang = ""
    merge_residues = .false.

    
    ierr = 0
    read(unit=iu_mdin, nml=dlfind, iostat=ierr)

    if ( ierr /= 0 ) then
       write(6,*)"Failed to read &dlfind namelist"
       call mexit(6,1)
    end if

    options%refcrd = adjustl(refcrd)
    options%active = active
    options%crdrep=""
    call str2upper(adjustl(crdrep),options%crdrep)
    options%optalg=""
    call str2upper(adjustl(optalg),options%optalg)
    options%trustrad=""
    call str2upper(adjustl(trustrad),options%trustrad)
    options%hessini=""
    call str2upper(adjustl(hessini),options%hessini)
    options%hessupd=""
    call str2upper(adjustl(hessupd),options%hessupd)
    options%neb=""
    call str2upper(adjustl(neb),options%neb)

    
    options%maxupdate = maxupdate
    options%delta = delta
    options%lbfgsmem = lbfgsmem
    options%verbosity = verbosity
    options%tol = tol
    options%tole = tole
    options%maxcycle = maxcycle
    options%maxstep = maxstep
    options%minstep = minstep
    options%scalestep = scalestep
    options%dimer = dimer
    options%dimer_mode = dimer_mode
    options%dimer_tolrot = dimer_tolrot
    options%dimer_maxrot = dimer_maxrot
    !options%tsmode_inp = tsmode_inp
    !options%tsmode_out = tsmode_out
    if ( len_trim(adjustl(wtmask)) == 0 ) then
       options%wtmask = active
    else
       options%wtmask = adjustl(wtmask)
    end if
    options%wtmm = wtmm

    options%hess_out = hess_out
    options%hess_iseg = hess_iseg
    options%hess_nseg = hess_nseg
    options%allow_shake = allow_shake

    options%coords2 = coords2
    options%neb_nimage = max(1,neb_nimage)
    options%neb_climb_test = max(0.d0,neb_climb_test)
    options%neb_freeze_test = max(0.d0,neb_freeze_test)
    options%nebk = max(1.d-8,nebk)
    options%neb_cart = neb_cart
    options%disang = adjustl(disang)
    options%merge_residues = merge_residues

  end subroutine read_dlfind_nml



  

  subroutine new_active_data(options,data,qsetup)

    use memory_module, only : atom_noshake
    use memory_module, only : lastpr
    use memory_module, only : natom,nres
    use memory_module, only : ix,ih,x
    use memory_module, only : m02,m04,m06,i02,lcrd,i100
    use memory_module, only : mass
    use memory_module, only : atom_name
    use memory_module, only : nbonh,nbona,nbper,iibh,ijbh
    use memory_module, only : residue_pointer
    use memory_module, only : coordinate
    use qmmm_module, only : qmmm_struct, get_atomic_number
    use findmask, only : atommask
    use state
    use file_io_dat, only : restrt,inpcrd
    
    implicit none

    type(dlfind_nml_type),intent(in) :: options
    type(active_data_type),intent(out) :: data
    logical,intent(in) :: qsetup

    integer, pointer :: amask(:) => null()
    integer, pointer :: residx(:) => null()
    integer, pointer :: wtmask(:) => null()
    integer, pointer :: qmmask(:) => null()
    integer, pointer :: revmap(:) => null()
    integer, pointer :: splitmask(:) => null()
    integer :: i,j,k,ia,ja,ires,jres
    logical :: err
    logical,parameter :: mergeres = .true.
    
    integer :: curidx,oldidx,ncons,nwts

    integer :: nbt,npairs,llind,ll,i3,j3,nsplit
    _REAL_,allocatable :: tcrd(:,:)

    character(len=MASKLEN) :: tmpstr

    logical :: has_constraint
    integer :: ppos

    integer,pointer :: discons(:,:) => NULL()
    integer :: ndiscons,nshake,naddshakemask
    integer,pointer :: addshakemask(:) => NULL()
    logical :: consmustbesameres
    
#include "../include/md.h"

    consmustbesameres = .false.
    if ( trim(options%crdrep) == "HDLC" ) then
       consmustbesameres = .true.
    else if ( trim(options%crdrep) == "HDLC-TC" ) then
       consmustbesameres = .true.
    end if
    
    has_constraint = .false.

    data%base=""
    data%ext=""
    ppos = scan(trim(restrt),".",BACK=.true.)
    if ( ppos > 0 ) then
       i = len_trim(restrt)
       data%base(1:ppos) = restrt(1:ppos)
       data%ext(1:i-ppos) = restrt(ppos+1:i)
    else
       data%base = trim(restrt)
       data%ext = "rst7"
    end if

    !write(6,'(3a)')"base='",trim(data%base),"'"
    !write(6,'(3a)')"ext ='",trim(data%ext),"'"

    if ( associated( data%condata ) ) then
       deallocate( data%condata )
    end if
    
    !
    ! Scratch space for the nonbond list
    !
    
    if ( associated(data%ipairs) ) then
       deallocate(data%ipairs)
    end if
    
    allocate( data%ipairs(lastpr) )
    
    data%ipairs(:) = 0

    !
    ! Flag indicating whether force() needs to
    ! do stuff related to setup
    !
    
    data%qsetup = qsetup

    !
    ! Energy initialization
    !
    
    data%ene = 0.d0


    data%is_neb = .false.
    if ( trim(options%neb) /= "NO" ) then
       data%is_neb = .true.
    end if
    
    data%is_ts_search = options%dimer
    if ( trim(options%optalg) == "PRFO" ) then
       data%is_ts_search = .true.
    end if

    data%is_hess = .false.
    if ( trim(options%optalg) == "HESS" ) then
       data%is_hess = .true.
       data%hesstype = 2
       if ( trim(options%hessini) == "ONEPOINT" ) then
          data%hesstype = 1
       end if
    end if
    

    if ( data%is_neb .and. options%dimer ) then
       write(6,'(2a)')"Cannot perform neb and dimer methods",&
            & "at the same time"
       call mexit(6,1)
    end if
    
    if ( data%is_neb .and. data%is_ts_search ) then
       write(6,'(2a)')"Cannot perform neb and optalg=PRFO ",&
            & "optimizations at the same time"
       call mexit(6,1)
    end if
    
    if ( data%is_neb .and. data%is_hess ) then
       write(6,'(2a)')"Cannot perform neb and optalg=HESS ",&
            & "calculations at the same time"
       call mexit(6,1)
    end if
    
    if ( data%is_hess .and. options%dimer ) then
       write(6,'(2a)')"Cannot perform dimer and optalg=HESS ",&
            & "calculations at the same time"
       call mexit(6,1)
    end if
    
    if ( trim(options%optalg) == "PRFO" .and. options%dimer ) then
       write(6,'(2a)')"Cannot perform dimer and optalg=PRFO ",&
            & "calculations at the same time"
       call mexit(6,1)
    end if

    !
    ! Full list of residue indexes
    !
    allocate( residx( natom ) )
    residx = 0
    do i=1,nres
       residx( residue_pointer(i) : residue_pointer(i+1) - 1 ) = i
    end do


    if ( associated(data%refcrd) ) then
       deallocate(data%refcrd)
    end if

    allocate( data%refcrd(3*natom) )
    data%refcrd(1:3*natom) = x(lcrd:lcrd+3*natom-1)
    
    if ( len_trim(options%refcrd) > 0 ) then
       call rdrest(natom,ntrx,options%refcrd,data%refcrd(1))
    end if
    
    !
    ! Find the active selection
    !
    
    allocate( amask( natom ) )

    tmpstr = options%active
    call atommask( natom, nres, 0, ih(m04), ih(m06), &
         & ix(i02), ih(m02), data%refcrd, tmpstr, &
         & amask )

    call PrintVMDSelection(natom,nres,amask,residue_pointer)


    !
    ! Do the UNSELECTED reference crds match the unselected
    ! current crds? If not, then print a warning.
    !
    if ( len_trim(options%refcrd) > 0 ) then
       ppos = CheckUnselectedCrds(natom,amask,data%refcrd)
       if ( ppos > 0 ) then
          write(6,'(/a,i6,4a/)')"NOTE: There are ",ppos, &
               & " unselected atoms in ",trim(options%refcrd), &
               & " whose coordinates differ from ",trim(inpcrd)
       end if
       ppos = 0
    end if
    
    !
    ! It's not worth trying to get the shake constraints
    ! to work as expected. If you uncomment this, shake
    ! constraints will be obeyed, but the dl-find optimized
    ! structures aren't actually stationary points, as
    ! can be verified by performing a frequency analysis.
    !
    ! It might have something to do with the way it builds
    ! up an approximate Hessian from BOFILL or BFGS.
    !
    ! The safer thing to do, is add all the active atoms
    ! to noshakemask.
    !
    ! Well, we can't do that here because this could be
    ! a MPI job with mpi_orig, so only sandermaster
    ! ever gets here.
    !
    ! Instead, we will mexit and print a message telling
    ! the used to use ntc=ntf=1
    !
    ! [[Can't call this here because only rank=1 gets here]]
    !call setnoshake(ix,amask,ntc,num_noshake,.true.)

    
    allocate( qmmask(natom) )
    qmmask = 0
    if ( qmmm_struct%nquant > 0 ) then
       do i = 1, qmmm_struct%nquant
          j = qmmm_struct%iqmatoms(i)
          qmmask(j) = 1
       end do
    else
       qmmask = 1
    end if
    
    if ( options%dimer .or. data%is_neb ) then
       
       allocate( wtmask(natom) )
       wtmask = qmmask
       if ( len_trim(options%wtmask) > 0 ) then
          wtmask = 0
          tmpstr = options%wtmask
          call atommask( natom, nres, 0, ih(m04), ih(m06), &
               & ix(i02), ih(m02), data%refcrd, tmpstr, &
               & wtmask )
       end if

    end if

    
    
    !
    ! Only particles with non-zero mass will be selected
    !
    
    data%nactive = 0
    do i=1,natom
       if ( amask(i) > 0 .and. mass(i) .gt. 0.1d0 ) then
          data%nactive = data%nactive + 1
       end if
    end do

    !
    ! Allocate space for the selected atoms
    !
    
    allocate( data%map(data%nactive) )
    data%map = 0
    allocate( data%atnums(data%nactive) )
    data%atnums = 0
    allocate( data%atmass(data%nactive) )
    data%atmass = 0
    allocate( data%residx(data%nactive) )
    data%residx = 0
    allocate( data%wts(data%nactive) )
    data%wts = 1.d0

    allocate( revmap(natom) )
    revmap = 0
    allocate( splitmask(natom) )
    splitmask = 0
    
    !
    ! Fill-in the selected atom properties
    !
    
    j = 0
    do i=1,natom
       if ( amask(i) > 0 .and. mass(i) .gt. 0.1d0 ) then
          j=j+1
          data%map(j) = i
          revmap(i) = j
          data%atmass(j) = mass(i) * AU_PER_AMU
          data%residx(j) = residx(i)
          
          ! If ix(i100) > 0, then parm7 has ATOMIC_NUMBERS section
          if ( ix(i100) == 0 ) then
             err = .false.
             call get_atomic_number(atom_name(i), &
                  & mass(i), data%atnums(j), err)
             if ( err ) then
                data%atnums(j) = 1
             end if
          else
             data%atnums(j) = ix(i100+i)
          end if

          if ( options%dimer .or. data%is_neb ) then
             if ( wtmask(i) > 0 ) then
                if ( qmmask(i) > 0 ) then
                   data%wts(j) = 1.d0
                else
                   data%wts(j) = options%wtmm
                end if
             else
                data%wts(j) = 0.d0
             end if
          end if
          
       end if
    end do

    !
    ! Reset residue indexes so they increment
    ! sequentially from 1 rather than the
    ! 'random'-like numbering from the extraction
    !

    oldidx = 0
    curidx = 0
    do i = 1, data%nactive
       if ( data%residx(i) /= oldidx ) then
          curidx = curidx + 1
          oldidx = data%residx(i)
       end if
       data%residx(i) = curidx
    end do

    !
    ! Create a list of shaked bonds
    !

    if ( associated( data%shaked_pairs ) ) then
       deallocate( data%shaked_pairs )
    end if

    npairs = 0
    if ( ntc > 1 ) then

       nbt = nbonh
       if ( ntc > 2 ) then
          nbt = nbonh+nbona+nbper
       end if

       do llind = 1,nbt
          if (atom_noshake(llind) == 1) cycle
          ll=llind
          i3 = ix( iibh + ll - 1 )
          j3 = ix( ijbh + ll - 1 )
          i  = i3/3+1
          j  = j3/3+1
          if ( mass(i) .le. 0.1d0 .or. mass(j) .le. 0.1d0 ) cycle
          if ( residx(i) /= residx(j) ) cycle
          if ( amask(i) + amask(j) == 0 ) cycle
          !
          ! Does the active region split a shaked bond?
          ! If so, then keep track of the atoms involved in the
          ! split so we can report to the user how they should
          ! modify noshakemask
          !
          if ( amask(i) + amask(j) == 1 ) then
             !write(6,*)i,j,amask(i),amask(j),residx(i),mass(i),mass(j)
             if ( amask(i) == 1 ) then
                splitmask(i) = 1
             else
                splitmask(j) = 1
             end if
             cycle
          end if
          npairs = npairs + 2
       end do

       has_constraint = .false.
       if ( npairs > 0 .or. sum(splitmask) > 0 ) then
          has_constraint = .true.
       end if
       
       if ( has_constraint ) then

          !if ( .not. data%is_hess ) then
          write(6,'(a)')"The DL-Find active atom selection includes shaked bonds."
          write(6,'(a)')"Either set ntc=1 and ntf=1 to turn off shake, or"
          write(6,'(a)')"include the active atoms in the noshakemask parameter."
          if ( .not. options%allow_shake ) then
             call mexit(6,1)
          end if
          !end if
          
          ! if ( trim(options%optalg) == "PRFO" ) then
          !    write(6,*)"This calculation requires a Hessian, but there are "
          !    write(6,*)"shake constraints involving the activated atoms."
          !    write(6,*)"You should add the &dlfind active mask to the"
          !    write(6,*)"noshakemask selection within the &cntrl section."
          !    call mexit(6,1)
          ! end if

          
          if ( ( trim(options%crdrep) == "CART" .or. &
               & trim(options%crdrep) == "MASS" ) .and. npairs > 0 ) then
             write(6,*)"The activated atoms include shaked bonds;"
             write(6,*)"however, this is only supported with:"
             write(6,*)"crdrep='DLC', 'HDLC', 'TC', or 'HDLC-TC'"
             write(6,*)"Either change crdrep in &dlfind or set noshakemask"
             write(6,*)"in &cntrl to include the active atoms"
             call mexit(6,1)
          end if
       end if

       
       allocate( data%shaked_pairs(npairs) )
       data%shaked_pairs = 0
       npairs = 0
       do llind = 1,nbt
          if (atom_noshake(llind) == 1) cycle
          ll=llind
          i3 = ix( iibh + ll - 1 )
          j3 = ix( ijbh + ll - 1 )
          i  = i3/3+1
          j  = j3/3+1
          
          if ( mass(i) .le. 0.1d0 .or. mass(j) .le. 0.1d0 ) cycle
          if ( residx(i) /= residx(j) ) cycle
          if ( amask(i) + amask(j) /= 2 ) cycle
          
          data%shaked_pairs(npairs+1) = revmap(i)
          data%shaked_pairs(npairs+2) = revmap(j)
          npairs = npairs + 2
       end do
       
    end if

    nsplit = sum(splitmask)
    if ( nsplit > 0 ) then
       call ReportNoshakemaskModification(natom,nres,splitmask,residue_pointer)
    end if

    
    naddshakemask = 0
    nshake = 0
    if ( associated(data%shaked_pairs) ) then
       nshake = size(data%shaked_pairs) / 2
       allocate( addshakemask(nshake) )
       addshakemask = 1
       naddshakemask = nshake
    end if
    
    ndiscons = 0
    if ( len_trim(options%disang) > 0 ) then
       call ReadDisang(options%disang,ndiscons,discons)
       
       do i=1,ndiscons
          if ( discons(1,i) == 1 ) then
             if ( amask( discons(2,i) ) + amask( discons(3,i) ) /= 2 ) then
                write(6,'(a,2i6,a)')"Bond constraint ", &
                     & discons(2,i),discons(3,i), &
                     & " is not fully contained in the active mask"
                call mexit(6,1)
             end if
          else if ( discons(1,i) == 2 ) then
             if ( amask( discons(2,i) ) + amask( discons(3,i) ) &
                  & + amask( discons(4,i) ) /= 3 ) then
                write(6,'(a,3i6,a)')"Angle constraint ", &
                     & discons(2,i),discons(3,i),discons(4,i), &
                     & " is not fully contained in the active mask"
                call mexit(6,1)
             end if
          else if ( discons(1,i) == 3 ) then
             if ( amask( discons(2,i) ) + amask( discons(3,i) ) &
                  & + amask( discons(4,i) ) + amask( discons(5,i) ) /= 4 ) then
                write(6,'(a,4i6,a)')"Torsion constraint ", &
                     & discons(2,i),discons(3,i),discons(4,i),discons(5,i), &
                     & " is not fully contained in the active mask"
                call mexit(6,1)
             end if
          end if
       end do


       if ( options%merge_residues ) then

          
          do i=1,ndiscons
             do k=3,2+discons(1,i)
                ires = data%residx( revmap( discons(2,i) ) )
                jres = data%residx( revmap( discons(k,i) ) )
                if ( ires < jres ) then
                   do j=1,data%nactive
                      if ( data%residx(j) == jres ) then
                         data%residx(j) = ires
                      end if
                   end do
                else if ( ires > jres ) then
                   do j=1,data%nactive
                      if ( data%residx(j) == ires ) then
                         data%residx(j) = jres
                      end if
                   end do
                end if
             end do
          end do
          
       else if ( consmustbesameres ) then
          
          do i=1,ndiscons
             if ( discons(1,i) == 1 ) then
                if ( residx( discons(2,i) ) /= residx( discons(3,i) ) ) then
                   write(6,'(a,2i6,2a)')"Bond constraint ", &
                        & discons(2,i),discons(3,i), &
                        & " crosses residue boundary; use crdrep='DLC' or 'TC'", &
                        & " or set merge_residues=T in &dlfind"
                   call mexit(6,1)
                end if
             else if ( discons(1,i) == 2 ) then
                if ( residx( discons(2,i) ) /= residx( discons(3,i) ) .or. &
                     &  residx( discons(2,i) ) /= residx( discons(4,i) ) ) then
                   write(6,'(a,3i6,2a)')"Angle constraint ", &
                        & discons(2,i),discons(3,i),discons(4,i), &
                        & " crosses residue boundary; use crdrep='DLC' or 'TC'", &
                        & " or set merge_residues=T in &dlfind"
                   call mexit(6,1)
                end if
             else if ( discons(1,i) == 3 ) then
                if ( residx( discons(2,i) ) /= residx( discons(3,i) ) .or. &
                     & residx( discons(2,i) ) /= residx( discons(4,i) ) .or. &
                     & residx( discons(2,i) ) /= residx( discons(5,i) ) ) then
                   write(6,'(a,4i6,2a)')"Torsion constraint ", &
                        & discons(2,i),discons(3,i),discons(4,i),discons(5,i), &
                        & " crosses residue boundary; use crdrep='DLC' or 'TC'", &
                        & " or set merge_residues=T in &dlfind"
                   call mexit(6,1)
                end if
             end if
          end do
       end if
       
       do i=1,ndiscons
          do j=2,2+discons(1,i)
             discons(j,i) = revmap( discons(j,i) )
          end do
       end do
       
       if ( associated(data%shaked_pairs) ) then
          do i=1,nshake
             do j=1,ndiscons
                if ( discons(1,j) == 1 ) then
                   if ( ( discons(2,j) == data%shaked_pairs(1+(i-1)*2) .and. &
                        & discons(3,j) == data%shaked_pairs(2+(i-1)*2) ) .or. &
                        ( discons(3,j) == data%shaked_pairs(1+(i-1)*2) .and. &
                        & discons(2,j) == data%shaked_pairs(2+(i-1)*2) ) ) then
                      addshakemask(i) = 0
                      naddshakemask = naddshakemask - 1
                      exit
                   end if
                end if
             end do
          end do
       end if
    end if

    if ( naddshakemask + ndiscons > 0 ) then
       
       write(6,'(2(a,i7),a/)')"Adding ",ndiscons," constraints from the disang file and ",&
            & naddshakemask," shake constraints"
       
       allocate( data%condata(5,naddshakemask + ndiscons) )
       data%condata = 0
       j=0
       do i=1,nshake
          if ( addshakemask(i) == 1 ) then
             j = j + 1
             data%condata(1:5,j) = (/ 1, data%shaked_pairs(1+(i-1)*2), &
                  & data%shaked_pairs(2+(i-1)*2), 0, 0 /)
          end if
       end do
       do i=1,ndiscons
          j=j+1
          data%condata(1:5,j) = discons(1:5,i)
       end do
    end if
    
    
    
    if ( trim(dlfind_nml%hessini) == "ONEPOINT" ) then
       data%hesstype = 1
    else if ( trim(dlfind_nml%hessini) == "TWOPOINT" ) then
       data%hesstype = 2
    else if ( trim(dlfind_nml%hessini) == "DIAGONAL" ) then
       data%hesstype = 0
    else if ( trim(dlfind_nml%hessini) == "IDENTITY" ) then
       data%hesstype = 0
    else
       data%hesstype = 0
    end if
    
    
    if ( options%dimer ) then

       if ( len_trim(options%coords2(1)) < 1 ) then
          write(6,'(3a)')"The dimer method requires you to provide ",&
               & "a guess at the transition state mode by setting ",&
               & "coords2(1) = 'tsmode.rst7' in &dlfind."
          write(6,'(5a)')"The TS mode can be made by performing a ",&
               & "minimization in which you harmonically restrain ",&
               & "the reaction coordinates slightly toward the ",&
               & "product state. Alternatively, the approximate TS ",&
               & "and TS mode will be printed during a NEB calculation"
          call mexit(6,1)
       end if
       
       data%nframes = 1
       allocate( data%frames(3,data%nactive,data%nframes) )
       data%frames = 0.d0

       allocate( tcrd(3,natom) )
       tcrd = 0.d0
       
       call rdrest(natom,ntrx,options%coords2(1),tcrd(1,1))
       do i = 1, data%nactive
          j = data%map(i)
          data%frames(1:3,i,1) = tcrd(1:3,j)
       end do

       ppos = CheckUnselectedCrds(natom,amask,tcrd(1,1))
       if ( ppos > 0 ) then
          write(6,'(/a,i6,4a/)')"NOTE: There are ",ppos, &
               & " unselected atoms in ",trim(options%refcrd), &
               & " whose coordinates differ from ",trim(inpcrd)
       end if
       ppos = 0
       
       deallocate(tcrd)
       
    else if ( data%is_neb ) then
       
       data%nframes = 0
       do i=1,size(options%coords2)
          if ( len_trim(options%coords2(i)) > 0 ) then
             if ( i > 1 ) then
                if ( options%coords2(i) == options%coords2(i-1) ) then
                   exit
                end if
             end if
             data%nframes = data%nframes + 1
          end if
       end do
       
       if ( data%nframes < 1 ) then
          write(6,'(2a)')"neb was activated, but a product ",&
               & "structure was not provided in coords2"
          call mexit(6,1)
       end if

       allocate( tcrd(3,natom) )
       tcrd = 0.d0
       
       allocate( data%frames(3,data%nactive,data%nframes) )
       data%frames = 0.d0
       
       j=0
       do i=1,size(options%coords2)
          if ( len_trim(options%coords2(i)) > 0 ) then
             
             if ( i > 1 ) then
                if ( options%coords2(i) == options%coords2(i-1) ) then
                   exit
                end if
             end if
             
             j=j+1
             
             call rdrest(natom,ntrx,options%coords2(i),tcrd(1,1))
             do ia = 1, data%nactive
                ja = data%map(ia)
                data%frames(1:3,ia,j) = tcrd(1:3,ja)
             end do

             ppos = CheckUnselectedCrds(natom,amask,tcrd(1,1))
             if ( ppos > 0 ) then
                write(6,'(/a,i6,4a/)')"NOTE: There are ",ppos, &
                     & " unselected atoms in ",trim(options%refcrd), &
                     & " whose coordinates differ from ",trim(inpcrd)
             end if
             ppos = 0
             
          end if
       end do
       
       deallocate(tcrd)

    end if

    nwts = 0
    if ( associated(data%wts) ) then
       nwts = data%nactive
    end if

    ncons = 0
    !if ( associated( data%shaked_pairs ) ) then
    !   ncons = size(data%shaked_pairs) / 2
    !end if
    if ( associated( data%condata ) ) then
       ncons = size(data%condata,2)
    end if
    
    data%nvar2 = data%nframes*(3*data%nactive) & ! frame coordinates
         & + nwts & ! weights
         & + data%nactive ! mass

    data%nspec = 2*data%nactive + 5*ncons + data%nactive


    if ( associated(amask) ) then
       deallocate(amask)
    end if

    if ( associated(residx) ) then
       deallocate(residx)
    end if

    if ( associated(wtmask) ) then
       deallocate(wtmask)
    end if

    if ( associated(qmmask) ) then
       deallocate(qmmask)
    end if

    if ( associated(revmap) ) then
       deallocate(revmap)
    end if

    if ( associated(splitmask) ) then
       deallocate(splitmask)
    end if

    if ( associated( discons ) ) then
       deallocate( discons )
    end if

    if ( associated( addshakemask ) ) then
       deallocate( addshakemask )
    end if

  end subroutine new_active_data


  

  subroutine enforce_shake()

    use memory_module, only : xx => x
    use memory_module, only : ix, ih
    use memory_module, only : natom,ibelly
    use memory_module, only : nbonh,nbona,iibh,ijbh,ibellygp
    use memory_module, only : conp, skip
    use memory_module, only : coordinate
    use memory_module, only : iifstwt,iifstwr
    use memory_module, only : noshake,nres,i02
    use fastwt, only : quick3

    implicit none
    
#include "../include/md.h"

    _REAL_ :: winv(nbonh+nbona)
    integer :: nitp
    logical :: belly
    logical qspatial
    
    winv(:) = 1.d0
    nitp = 0
    belly = ibelly.eq.1
    qspatial = .false.

    belly = ibelly.eq.1
    if( ntc .ne. 1 ) then
       call shake(nrp, nbonh, nbona, 0, ix(iibh), ix(ijbh), ix(ibellygp), &
            & winv, conp, skip, coordinate, coordinate, nitp, belly, &
            & ix(iifstwt), ix(noshake),qspatial)
       if (nitp == 0) then
          write(6,'(a)') 'SHAKE error encountered in dlfind_module.F90/enforce_shake'
          call mexit(6,1)
       end if
       call quick3(coordinate, coordinate, ix(iifstwr), natom, nres, ix(i02))
    end if
    
  end subroutine enforce_shake

  

  
  
  subroutine update_extrapts()
    
    use memory_module, only : x, ix, lcrd
    implicit none
    
#include "extra_pts.h"
    
    if( numextra > 0 ) call local_to_global(x(lcrd),x,ix)
    
  end subroutine update_extrapts


  
  
  
  subroutine set_crds_from_active_crds(acrds)
    
    use memory_module, only : coordinate
    implicit none
    
    real(kind=c_double),intent(in) :: acrds(3,active_data%nactive)
    integer :: i,j
    
    do i = 1 , active_data%nactive
       j = active_data%map(i)
       coordinate(1:3,j) = acrds(1:3,i) / AU_PER_ANG
    end do
    
    call update_extrapts()
    
  end subroutine set_crds_from_active_crds

  
  

  subroutine get_active_gradients(agrds)
    
    use memory_module, only : frc
    !use memory_module, only : natom
    implicit none

    real(kind=c_double),intent(out) :: agrds(3,active_data%nactive)
    integer :: i,j

    ! do i=natom-8,natom
    !    write(6,'(i7,3es13.4)')i,frc(1:3,i)
    ! end do
    
    do i = 1 , active_data%nactive
       j = active_data%map(i)
       agrds(1:3,i) =  - frc(1:3,j) * AU_PER_AMBER_FRC
    end do
    
  end subroutine get_active_gradients




  subroutine cpt_energy_and_forces(ener)

    use memory_module, only : xx => x
    use memory_module, only : ix, ih
    use memory_module, only : coordinate, frc
    use memory_module, only : l96,l97,l98,l99
    use memory_module, only : nbonh,nbona,iibh,ijbh,ibellygp,ibelly
    use memory_module, only : conp, skip, natom
    use memory_module, only : iifstwt,iifstwr,noshake,nres,i02
    use fastwt, only : quick3v
    use state
    
    implicit none

    type(state_rec),intent(inout) :: ener

#include "../include/md.h"
    
    logical :: do_list_update
    _REAL_ :: virials(4)
    _REAL_ :: winv(nbonh+nbona)
    integer :: nitp
    logical :: belly

    ener = 0.d0
    
    virials(:) = 0.d0
    winv(:) = 1.d0
    nitp = 0
    belly = ibelly.eq.1

    do_list_update = .true.
    ntnb = 1

    active_data%neval = active_data%neval + 1

    ! if( ntc .ne. 1 ) then
    !    call enforce_shake()
    ! end if

    frc=0.d0
    
    ! call force
    call force( xx,ix,ih,active_data%ipairs, &
         & coordinate,frc,active_data%ene,virials, &
         & xx(l96),xx(l97),xx(l98), xx(l99), active_data%qsetup, &
         & do_list_update, active_data%neval)

    active_data%ene%tot = active_data%ene%pot%tot
    
    ! if( ntc .ne. 1 ) then
    !    ! RATTLE-V, correct forces
    !    call rattlev(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
    !         & winv,conp,skip,coordinate,frc,nitp,belly,ix(iifstwt), &
    !         & ix(noshake), .false.)

    !    ! use SETTLE to deal with water model
    !    call quick3v(coordinate, frc, ix(iifstwr), natom, nres, ix(i02))
    ! end if

  end subroutine cpt_energy_and_forces
  
  

  
  subroutine get_energy_and_active_grds(acrds,ene,agrds)

    implicit none

    real(kind=c_double),intent(in)    :: acrds(3,active_data%nactive)
    real(kind=c_double),intent(out)   :: ene
    real(kind=c_double),intent(out)   :: agrds(3,active_data%nactive)

    call set_crds_from_active_crds(acrds)
    call cpt_energy_and_forces(active_data%ene)
    call get_active_gradients(agrds)
    ene = active_data%ene%pot%tot * AU_PER_KCAL

  end subroutine get_energy_and_active_grds




  subroutine print_energies()
    
    implicit none

    _REAL_ :: onefac(3)
    
    onefac(:) = 0.d0
    
    call prntmd(active_data%neval, &
         & 0.d0,active_data%ene,onefac,0,.false.)

  end subroutine print_energies
  

  
  
  subroutine write_rst_file()
    
    use memory_module, only : coordinate

    implicit none

#include "../include/md.h"
    
    !if ( master ) then
    call minrit(0,nrp,ntxo,coordinate(1,1))
    !end if

  end subroutine write_rst_file


  
  
  subroutine write_custom_rst_file(fname)

    use file_io_dat, only : owrite, title
    use binrestart, only : write_nc_restart
    use memory_module, only : x,lcrd,natom, coordinate
    use nblist, only: a,b,c,alpha,beta,gamma

    implicit none

    character(len=*),intent(in) :: fname
    
#include "../include/md.h"
#include "box.h"

    integer :: i
    integer :: fh

    fh = 50
    
    ! Open a file on filehandle 17, with name tsmode_rst
    if (ntxo == 2) then
       call write_nc_restart(trim(adjustl(fname)), &
            & title,owrite,nrp,ntb,.true.,x,x,0.0d0,.false.&
#ifdef MPI
            & , 0.0d0, 0, 0, (/ 0 /), (/ 0 /), (/ 0 /), 0, 0, 0, 0.d0, 0.d0 &
#endif
            & )
       return
    else
       call amopen(fh,trim(adjustl(fname)),'R','F','W')
       ! Write coordinates to filehandle 17
       !call minri2(fh,natom,ntxo,coordinate(1,1))
       !close(fh)
       write(fh,'(a4)')title
       if (natom>99999) then
          write(fh,'(i6)')natom
       else
          write(fh,'(i5)')natom
       end if
       write(fh,'(6f12.7)')(x(lcrd+i-1),i=1,3*natom)
       if ( ntb /= 0 ) then
          write(fh,'(6f12.7)') a,b,c,alpha,beta,gamma
       end if
       close(unit=fh)
    end if
    
    
  end subroutine write_custom_rst_file



  
  subroutine write_traj_frame()
    use memory_module, only : x,lcrd,natom
    use file_io_dat, only : ioutfm, MDCRD_UNIT
#ifdef BINTRAJ
   use netcdf
   use bintraj, only : end_binary_frame
#endif

    implicit none

#include "box.h"
    
    logical :: loutfm

    loutfm = (ioutfm <= 0)
    
    call corpac(x(lcrd),1,natom*3,MDCRD_UNIT,loutfm)
    !if (ntwf /= 0) call corpac(x(lforce),1,natom*3,MDFRC_UNIT,loutfm)
    if (ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)
#ifdef BINTRAJ
    if (ioutfm > 0) call end_binary_frame(MDCRD_UNIT)
#endif

  end subroutine write_traj_frame



  
  
  subroutine dlf_get_gradient_sander(nvar, coords, energy, &
       & gradient, iimage, kiter, status) bind(c)

    implicit none
    
    integer(kind=c_int), intent(in),value :: nvar ! number of xyz variables (3*nat)
    real(kind=c_double), intent(in) :: coords(nvar) ! coordinates
    real(kind=c_double), intent(out) :: energy ! energy
    real(kind=c_double), intent(out) :: gradient(nvar) ! gradient
    integer(kind=c_int), intent(in),value :: iimage ! current image (for NEB)
    integer(kind=c_int), intent(in),value :: kiter ! flag related to microiterations
    integer(kind=c_int), intent(out) :: status ! return code

    status = 0
    energy = 0.d0
    gradient(:) = 0.d0

    call get_energy_and_active_grds(coords(1),energy,gradient(1))
    call print_energies()
    
  end subroutine dlf_get_gradient_sander



  
  
  subroutine dlf_get_hessian_sander(nvar, coords, hessian, status) bind(c)
    
    implicit none
    
    integer(kind=c_int), intent(in), value :: nvar ! number of xyz variables (3*nat)
    real(kind=c_double), intent(in) :: coords(nvar) ! coordinates
    real(kind=c_double), intent(out) :: hessian(nvar, nvar) ! hessian
    integer(kind=c_int), intent(out) :: status ! return code

    _REAL_,allocatable :: gpos(:)
    _REAL_,allocatable :: gneg(:)
    _REAL_,allocatable :: tcrd(:)
    _REAL_ :: ene
    _REAL_ :: DELTA
    integer :: i,j
    
    hessian = 0.d0
    status = 0

    allocate(tcrd(nvar))
    tcrd = coords

    allocate(gpos(nvar))
    gpos = 0.d0

    allocate(gneg(nvar))
    gneg = 0.d0

    ene = 0.d0
    DELTA = dlfind_nml%delta
    
    if ( active_data%hesstype == 1 ) then

       write(6,'(a)')"Sander is calculating a one-point Hessian. This may take a while."
       
       call get_energy_and_active_grds(tcrd,ene,gneg)

       do i=1,nvar
          tcrd(i) = tcrd(i) + DELTA
          call get_energy_and_active_grds(tcrd,ene,gpos)
          tcrd(i) = tcrd(i) - DELTA
          hessian(:,i) = (gpos-gneg)/DELTA
       end do

    else if ( active_data%hesstype == 2 ) then
       
       write(6,'(a)')"Sander is calculating a two-point Hessian. This may take a while."

       do i=1,nvar
          tcrd(i) = tcrd(i) + DELTA
          call get_energy_and_active_grds(tcrd,ene,gpos)
          tcrd(i) = tcrd(i) - 2.d0*DELTA
          call get_energy_and_active_grds(tcrd,ene,gneg)
          tcrd(i) = tcrd(i) + DELTA
          hessian(:,i) = (gpos-gneg)/(2.d0*DELTA)
       end do

    else

       status = 1
       
    end if

    if ( status == 0 ) then
       do i=1,nvar-1
          do j=i+1,nvar
             ene = 0.5d0 * ( hessian(i,j) + hessian(j,i) )
             hessian(i,j) = ene
             hessian(j,i) = ene
          end do
       end do
    end if
       
  end subroutine dlf_get_hessian_sander


#ifdef BINTRAJ
  
  subroutine run_hessian_calc()
    
    use HessianNetCDF_mod, only : CheckHessianFileExists
    use HessianNetCDF_mod, only : OpenHessianNetCDF
    use HessianNetCDF_mod, only : CloseHessianNetCDF
    use HessianNetCDF_mod, only : WriteHessianColumn
    use HessianNetCDF_mod, only : WriteHessianGeometry
    use HessianNetCDF_mod, only : GetHessianStoredColIdxs
    use HessianNetCDF_mod, only : GetHessianNumStoredCols
    use HessianNetCDF_mod, only : GetHessianCrds

    use memory_module, only : coordinate

    implicit none
    
    integer :: nat,ncnat,nhess
    integer :: i,j,iseg,nseg
    logical :: exists
    _REAL_ :: DELTA
    _REAL_,allocatable :: acrds(:)
    _REAL_,allocatable :: tcrds(:)
    _REAL_,allocatable :: glo(:)
    _REAL_,allocatable :: ghi(:)
    character(len=MASKLEN) :: fname
    logical,allocatable :: colneeded(:)
    integer,allocatable :: segsizes(:)
    integer,allocatable :: segoffset(:)
    integer :: nstored
    integer,allocatable :: storedcols(:)
    _REAL_ :: etmp

    integer :: ncon
    integer,allocatable :: condata(:,:)
    character(len=20) :: existmsg
    
    nat = active_data%nactive
    DELTA = dlfind_nml%delta
    nhess = 3*nat
    iseg = dlfind_nml%hess_iseg
    nseg = dlfind_nml%hess_nseg

    
    if ( nseg < 1 .or. nseg > 999 ) then
       write(6,'(a,i0)')"hess_nseg must be >= 1 and <= 999, not ",nseg
       call mexit(6,1)
    end if
    if ( iseg < 1 .or. iseg > nseg ) then
       write(6,'(a,i0)')"hess_nseg must be >= 1 and <= hess_nseg, not ", &
            & iseg
       call mexit(6,1)
    end if
    
    if ( nseg == 1 ) then
       fname = dlfind_nml%hess_out
    else
       write(fname,'(a,a,i0.4)')trim(adjustl(dlfind_nml%hess_out)),".",iseg
    end if
    
    exists = CheckHessianFileExists(fname)
    
    if ( exists ) then
       existmsg = "Appending existing"
    else
       existmsg = "Writing new"
    end if
    
    ncon = 0
    if ( associated(active_data%condata) ) then
       ncon = size(active_data%condata,2)
       write(6,'(4a,i5)')trim(existmsg), &
            & " Hessian NetCDF3 to file: ",trim(fname), &
            & " with ncon=",ncon
       call OpenHessianNetCDF(fname,nat,ncon,active_data%condata)
    else
       write(6,'(4a,i5)')trim(existmsg), &
            & " Hessian NetCDF3 to file: ",trim(fname), &
            & " with ncon=",ncon
       allocate(condata(5,ncon))
       condata = 0
       call OpenHessianNetCDF(fname,nat,ncon,condata)
       deallocate(condata)
    end if
    
    allocate(acrds(nhess))
    acrds=0.d0
    allocate(glo(nhess))
    glo=0.d0
    allocate(ghi(nhess))
    ghi=0.d0

    allocate( segsizes(nseg) )
    segsizes = 0
    do i=1,nhess
       j = mod(i-1,nseg) + 1
       segsizes(j) = segsizes(j) + 1
    end do

    allocate( segoffset(nseg+1) )
    segoffset = 0
    do i=1,nseg
       segoffset(i+1) = segoffset(i) + segsizes(i)
    end do
    deallocate(segsizes)

    allocate( colneeded(nhess) )
    colneeded = .false.
    
    do i=segoffset(iseg)+1,segoffset(iseg+1)
       colneeded(i) = .true.
    end do
    
    do i=1,nat
       j = active_data%map(i)
       acrds(1+(i-1)*3) = coordinate(1,j) * AU_PER_ANG
       acrds(2+(i-1)*3) = coordinate(2,j) * AU_PER_ANG
       acrds(3+(i-1)*3) = coordinate(3,j) * AU_PER_ANG
    end do

    
    if ( exists ) then
       allocate(tcrds(nhess))
       tcrds=0.d0
       call GetHessianCrds(nat,tcrds(1))
       tcrds = abs(tcrds - acrds)
       do i=1,nat
          do j=1,3
             if ( tcrds(j+(i-1)*3) > 1.d-6 ) then
                write(6,'(a,i6,2a)')"Current coordinate atom=",i, &
                     & " differs from those stored in ", &
                     & trim(fname)
                call mexit(6,1)
             end if
          end do
       end do
       deallocate(tcrds)

       nstored = 0
       call GetHessianNumStoredCols(nstored)
       if ( nstored > 0 ) then
          allocate(storedcols(nstored))
          storedcols = 0
          call GetHessianStoredColIdxs(nstored,storedcols)
          do i=1,nstored
             j = storedcols(i)
             if ( j > 0 .and. j <= nhess ) then
                colneeded(j) = .false.
             end if
          end do
          deallocate(storedcols)
       end if
       
    else
       
       ! hessian file did not exist
       call WriteHessianGeometry(nat,acrds(1), &
            & active_data%map,active_data%atnums)
       
    end if

    if ( active_data%hesstype == 2 ) then

       write(6,'(a,2i7)')"Two-point Hessian calculation of columns ", &
            & segoffset(iseg)+1,segoffset(iseg+1)
       
       do i=segoffset(iseg)+1,segoffset(iseg+1)
          
          if ( .not. colneeded(i) ) then
             cycle
          end if

          acrds(i) = acrds(i) + DELTA
          call get_energy_and_active_grds(acrds(1),etmp,ghi(1))
          acrds(i) = acrds(i) - 2.d0*DELTA
          call get_energy_and_active_grds(acrds(1),etmp,glo(1))
          acrds(i) = acrds(i) + DELTA
          
          ghi = (ghi-glo)/(2.d0*DELTA)
          call WriteHessianColumn(nat,i,ghi)
       end do

    else

       write(6,'(a,2i7)')"One-point Hessian calculation of columns ", &
            & segoffset(iseg)+1,segoffset(iseg+1)
       
       call get_energy_and_active_grds(acrds(1),etmp,glo(1))
       
       do i=segoffset(iseg)+1,segoffset(iseg+1)
          
          if ( .not. colneeded(i) ) then
             cycle
          end if
          
          acrds(i) = acrds(i) + DELTA
          call get_energy_and_active_grds(acrds(1),etmp,ghi(1))
          acrds(i) = acrds(i) - DELTA
          
          ghi = (ghi-glo)/DELTA
          call WriteHessianColumn(nat,i,ghi)
       end do

       
    end if

    call CloseHessianNetCDF()
    
  end subroutine run_hessian_calc
  
#endif
  
  
  subroutine dlf_get_multistate_gradients_sander(nvar, coords, energy, &
       & gradient, coupling, needcoupling, iimage, status) bind(c)
    
    implicit none
    
    integer(kind=c_int), intent(in), value :: nvar ! number of xyz variables (3*nat)
    real(kind=c_double), intent(in) :: coords(nvar) ! coordinates
    real(kind=c_double), intent(out) :: energy(2) ! multistate energy
    real(kind=c_double), intent(out) :: gradient(3, nvar/3, 2) ! xyz multistate gradients
    real(kind=c_double), intent(out) :: coupling(3, nvar/3) ! xyz interstate coupling gradient
    integer(kind=c_int), intent(in), value :: needcoupling ! true if interstate coupling gradients should be calculated
    integer(kind=c_int), intent(in), value :: iimage ! current image (for NEB)
    integer(kind=c_int), intent(out) :: status ! return code
   
    energy(:) = 0.d0
    gradient(:,:,:) = 0.d0
    coupling(:,:) = 0.d0
    status = 0
    
    write(6,*)"ERROR: dlf_get_multistate_gradients_sander is unimplemented"
    status = 1
 
  end subroutine dlf_get_multistate_gradients_sander


  
  subroutine dlf_get_params_sander( &
       & nvar, nvar2, nspec, coords, coords2, spec, ierr, tolerance, &
       & printl, maxcycle, maxene, tatoms, icoord, iopt, iline, maxstep, &
       & scalestep, lbfgs_mem, nimage, nebk, dump, restart, nz, ncons, nconn, &
       & update, maxupd, delta, soft, inithessian, carthessian, &
       & tsrel, maxrot, tolrot, nframe, nmass, nweight, timestep, fric0, &
       & fricfac, fricp, imultistate, state_i, state_j, pf_c1, pf_c2, &
       & gp_c3, gp_c4, ln_t1, ln_t2, printf, tolerance_e, distort, massweight, &
       & minstep, maxdump, task, temperature, po_pop_size, &
       & po_radius, po_contraction, po_tolerance_r, po_tolerance_g, &
       & po_distribution, po_maxcycle, po_init_pop_size, po_reset, &
       & po_mutation_rate, po_death_rate, po_scalefac, po_nsave, ntasks, &
       & tdlf_farm, n_po_scaling, neb_climb_test, neb_freeze_test, &
       & nzero, coupled_states, qtsflag, imicroiter, maxmicrocycle, &
       & micro_esp_fit ) bind(c)
    
    use memory_module, only : coordinate
    
    implicit none
    
    integer(kind=c_int), intent(in), value :: nvar ! number of xyz variables (3*nat)
    integer(kind=c_int), intent(in), value :: nvar2 ! number of variables to read in the second array (coords2)
    integer(kind=c_int), intent(in), value :: nspec ! number of values in the integer array spec
    real(kind=c_double), intent(inout) :: coords(nvar) ! start coordinates
    real(kind=c_double), intent(inout) :: coords2(nvar2) ! a real array that can be used depending on the calculation e.g. a second set of coordinates
    integer(kind=c_int), intent(inout) :: spec(nspec) ! (nat) fragment number, or -1: frozen, -2: x frozen, see dlf_coords.f90
    integer(kind=c_int), intent(out) :: ierr ! error code
    real(kind=c_double), intent(inout) :: tolerance ! main convergence criterion (Max grad comp.)
    integer(kind=c_int), intent(inout) :: printl ! how verbosely to write info to stdout
    integer(kind=c_int), intent(inout) :: maxcycle ! maximum number of cycles
    integer(kind=c_int), intent(inout) :: maxene ! maximum number of E&G evaluations
    integer(kind=c_int), intent(inout) :: tatoms ! atoms or arbitrary DOF
    integer(kind=c_int), intent(inout) :: icoord ! type of internal coordinates
    integer(kind=c_int), intent(inout) :: iopt ! type of optimisation algorithm
    integer(kind=c_int), intent(inout) :: iline ! type of line search or trust radius
    real(kind=c_double), intent(inout) :: maxstep ! maximum length of the step in internals
    real(kind=c_double), intent(inout) :: scalestep ! constant factor with which to scale the step
    integer(kind=c_int), intent(inout) :: lbfgs_mem ! number of steps in LBFGS memory
    integer(kind=c_int), intent(inout) :: nimage ! Number of images (e.g. in NEB)
    real(kind=c_double), intent(inout) :: nebk ! force constant for NEB
    integer(kind=c_int), intent(inout) :: dump ! after how many E&G calculations to dump a checkpoint file?
    integer(kind=c_int), intent(inout) :: restart ! restart mode: 0 new, 1 read dump file ...
    integer(kind=c_int), intent(inout) :: nz ! entries of nuclear charges (same order as coords)
    integer(kind=c_int), intent(inout) :: ncons ! number of constraints
    integer(kind=c_int), intent(inout) :: nconn ! number of user provided connections
    integer(kind=c_int), intent(inout) :: update ! Hessian update scheme
    integer(kind=c_int), intent(inout) :: maxupd ! Maximum number of Hessian updates
    real(kind=c_double), intent(inout) :: delta ! Delta-x in finite-difference Hessian
    real(kind=c_double), intent(inout) :: soft ! Abs(eigval(hess)) < soft -> ignored in P-RFO
    integer(kind=c_int), intent(inout) :: inithessian ! Option for method of calculating the initial Hessian
    integer(kind=c_int), intent(inout) :: carthessian ! Hessian update in cartesians?
    integer(kind=c_int), intent(inout) :: tsrel ! Transition vector I/O absolute or relative?
    integer(kind=c_int), intent(inout) :: maxrot ! maximum number of rotations in each DIMER step
    real(kind=c_double), intent(inout) :: tolrot ! angle tolerance for rotation (deg) in DIMER
    integer(kind=c_int), intent(inout) :: nframe ! number of structures
    integer(kind=c_int), intent(inout) :: nmass ! entries of atomic masses (nat or 0)
    integer(kind=c_int), intent(inout) :: nweight ! entries of weights (nat or 0)
    real(kind=c_double), intent(inout) :: timestep ! time step
    real(kind=c_double), intent(inout) :: fric0 ! start friction
    real(kind=c_double), intent(inout) :: fricfac ! factor to reduce friction (<1) whenever the energy is decreasing
    real(kind=c_double), intent(inout) :: fricp ! friction to use whenever energy increasing
    integer(kind=c_int), intent(inout) :: imultistate ! type of multistate calculation (0 = none)
    integer(kind=c_int), intent(inout) :: state_i ! lower state
    integer(kind=c_int), intent(inout) :: state_j ! upper state
    real(kind=c_double), intent(inout) :: pf_c1 ! penalty function parameter (aka alpha)
    real(kind=c_double), intent(inout) :: pf_c2 ! penalty function parameter (aka beta)
    real(kind=c_double), intent(inout) :: gp_c3 ! gradient projection parameter (aka alpha0)
    real(kind=c_double), intent(inout) :: gp_c4 ! gradient projection parameter (aka alpha1)
    real(kind=c_double), intent(inout) :: ln_t1 ! Lagrange-Newton orthogonalisation on threshold
    real(kind=c_double), intent(inout) :: ln_t2 ! Lagrange-Newton orthogonalisation off threshold
    integer(kind=c_int), intent(inout) :: printf ! how verbosely files should be written
    real(kind=c_double), intent(inout) :: tolerance_e ! convergence criterion on energy change
    real(kind=c_double), intent(inout) :: distort ! shift start structure along coords2 (+ or -)
    integer(kind=c_int), intent(inout) :: massweight ! use mass-weighted coordinates
    real(kind=c_double), intent(inout) :: minstep ! Hessian is not updated if step < minstep
    integer(kind=c_int), intent(inout) :: maxdump ! do only dump restart file after the at most maxdump E&G evaluations
    integer(kind=c_int), intent(inout) :: task ! number of taks for the task manager
    real(kind=c_double), intent(inout) :: temperature ! temperature for thermal analysis
    integer(kind=c_int), intent(inout) :: po_pop_size ! sample population size
    real(kind=c_double), intent(inout) :: po_radius  ! per-atom search radii (the prevailing units)
    real(kind=c_double), intent(inout) :: po_contraction ! factor by which the search radius decreases between search cycles. Cycle = 1 energy eval for each member of the sample population.
    real(kind=c_double), intent(inout) :: po_tolerance_r ! tolerances on po_radius: stop if any component of po_radius shrinks to less than the corresponding component of po_tolerance_r
    real(kind=c_double), intent(inout) :: po_tolerance_g ! convergence criterion: max abs component of g
    integer(kind=c_int), intent(inout) :: po_distribution ! type of distribn of sample points in space
    integer(kind=c_int), intent(inout) :: po_maxcycle ! maximum number of cycles
    integer(kind=c_int), intent(inout) :: po_init_pop_size ! size of initial population
    integer(kind=c_int), intent(inout) :: po_reset ! number of cycles before population resetting
    real(kind=c_double), intent(inout) :: po_mutation_rate ! Fraction of the total number of coordinates in the population to be mutated (randomly shifted) per cycle
    real(kind=c_double), intent(inout) :: po_death_rate ! Fraction of the population to be replaced by offspring per cycle
    real(kind=c_double), intent(inout) :: po_scalefac ! Multiplying factor for the absolute gradient vector in the force_bias stoch. search scheme
    integer(kind=c_int), intent(inout) :: po_nsave ! number of low-energy minima to store
    integer(kind=c_int), intent(inout) :: ntasks ! number of taskfarms (workgroups)
    integer(kind=c_int), intent(inout) :: tdlf_farm ! Some flag for the task farm
    integer(kind=c_int), intent(inout) :: n_po_scaling ! entries of radii scaling factors in the parallel optimization (0 [meaning all radii set to the base value], or a pre-known nivar) i.e. nvarin2= nframe*nat*3 + nweight + nmass + n_po_scaling
    real(kind=c_double), intent(inout) :: neb_climb_test ! threshold scale factor for spawning climbing image
    real(kind=c_double), intent(inout) :: neb_freeze_test ! threshold scale factor for freezing NEB images
    integer(kind=c_int), intent(inout) :: nzero ! number of zero vibrational modes in system
    integer(kind=c_int), intent(inout) :: coupled_states   ! Do we need to calculate the interstate coupling gradient? If coupled_states is false, coupling = zero.
    integer(kind=c_int), intent(inout) :: qtsflag ! additional info, like if tunnelig splittings are to be calculated (see dlf_qts.f90)
    integer(kind=c_int), intent(inout) :: imicroiter ! flag for microiterative calculations =0 : standard, non-microiterative calculation >0 : microiterative calculation [=1 : inside macroiterative loop [=2 : inside microiterative loop]
    integer(kind=c_int), intent(inout) :: maxmicrocycle ! max number of microiterative cycles before switching back to macro
    integer(kind=c_int), intent(inout) :: micro_esp_fit ! fit ESP charges to inner region during microiterations

    integer :: i,j,iframe
    integer :: icrdtype,idimer,ineb,ineb_cart
    

    do i=1,active_data%nactive
       j = active_data%map(i)
       coords(1+(i-1)*3) = coordinate(1,j) * AU_PER_ANG
       coords(2+(i-1)*3) = coordinate(2,j) * AU_PER_ANG
       coords(3+(i-1)*3) = coordinate(3,j) * AU_PER_ANG
       do iframe=1,active_data%nframes
          coords2(1+(i-1)*3+(iframe-1)*nvar) = active_data%frames(1,i,iframe) * AU_PER_ANG
          coords2(2+(i-1)*3+(iframe-1)*nvar) = active_data%frames(2,i,iframe) * AU_PER_ANG
          coords2(3+(i-1)*3+(iframe-1)*nvar) = active_data%frames(3,i,iframe) * AU_PER_ANG
       end do
    end do

    ! write(6,*)"nframes=",active_data%nframes
    ! do i=1,active_data%nactive
    !    write(6,'(i7,4f12.5)')i,active_data%wts(i), active_data%frames(1:3,i,1)
    ! end do

    j = nvar*active_data%nframes
    if ( associated( active_data%wts ) ) then
       do i=1,active_data%nactive
          coords2(i+j) = active_data%wts(i)
       end do    
       j = j + active_data%nactive
    end if
    
    do i=1,active_data%nactive
       coords2(i+j) = active_data%atmass(i)
    end do

    do i=1,active_data%nactive
       spec(i) = active_data%residx(i)
       spec(i+active_data%nactive) = active_data%atnums(i)
    end do



    ncons = 0
    ! if ( associated( active_data%shaked_pairs ) ) then
    !    ncons = size( active_data%shaked_pairs ) / 2
    !    j = 2*active_data%nactive
    !    do i=1,ncons
    !       spec(j+(i-1)*5 + 1) = 1
    !       spec(j+(i-1)*5 + 2) = active_data%shaked_pairs(1+(i-1)*2)
    !       spec(j+(i-1)*5 + 3) = active_data%shaked_pairs(2+(i-1)*2)
    !       spec(j+(i-1)*5 + 4) = 0
    !       spec(j+(i-1)*5 + 5) = 0
    !    end do
    ! end if
    
    if ( associated( active_data%condata ) ) then
       ncons = size(active_data%condata,2)
       j = 2*active_data%nactive
       do i=1,ncons
          spec(j+(i-1)*5 + 1) = active_data%condata(1,i)
          spec(j+(i-1)*5 + 2) = active_data%condata(2,i)
          spec(j+(i-1)*5 + 3) = active_data%condata(3,i)
          spec(j+(i-1)*5 + 4) = active_data%condata(4,i)
          spec(j+(i-1)*5 + 5) = active_data%condata(5,i)
       end do
    end if


    j = 2*active_data%nactive + 5*ncons
    do i=1,active_data%nactive
       spec(i+j) = 0
    end do

    
    ierr = 0
    tolerance = dlfind_nml%tol
    maxcycle = dlfind_nml%maxcycle
    maxene = 7*dlfind_nml%maxcycle
    tatoms = 1 ! this must be 1

    
    icoord = -1
    icrdtype = 0
    massweight = 0
    if ( trim(dlfind_nml%crdrep) == "CART" ) then
       icrdtype = 0
    else if ( trim(dlfind_nml%crdrep) == "MASS" ) then
       icrdtype = 0
       massweight = 1
    else if ( trim(dlfind_nml%crdrep) == "DLC" ) then
       icrdtype = 3
    else if ( trim(dlfind_nml%crdrep) == "HDLC" ) then
       icrdtype = 1
    else if ( trim(dlfind_nml%crdrep) == "TC" ) then
       icrdtype = 4
    else if ( trim(dlfind_nml%crdrep) == "HDLC-TC" ) then
       icrdtype = 2
    else
       write(6,'(3a)')"ERROR: Invalid crdrep='",trim(dlfind_nml%crdrep),"'"
    endif

    if ( icrdtype == 3 .or. icrdtype == 4 ) then
       do i=1,active_data%nactive
          spec(i) = 1
       end do
    end if
    
    iopt = -1
    if ( trim(dlfind_nml%optalg) == "LBFGS" ) then
       iopt = 3
    else if ( trim(dlfind_nml%optalg) == "PRFO" ) then
       iopt = 10
    else if ( trim(dlfind_nml%optalg) == "CG" ) then
       iopt = 1
    else if ( trim(dlfind_nml%optalg) == "SD" ) then
       iopt = 0
    else if ( trim(dlfind_nml%optalg) == "NR" ) then
       iopt = 20
    else
       write(6,'(3a)')"ERROR: Invalid optalg='",trim(dlfind_nml%optalg),"'"
    end if

    
    iline = -1
    if ( trim(dlfind_nml%trustrad) == "CONST" ) then
       iline = 0
    else if ( trim(dlfind_nml%trustrad) == "ENE" ) then
       iline = 1
    else if ( trim(dlfind_nml%trustrad) == "GRAD" ) then
       iline = 2
    else
       write(6,'(3a)')"ERROR: Invalid trustrad='",trim(dlfind_nml%trustrad),"'"
    end if

    inithessian = -1
    if ( trim(dlfind_nml%hessini) == "ONEPOINT" ) then
       inithessian = 1
    else if ( trim(dlfind_nml%hessini) == "TWOPOINT" ) then
       inithessian = 2
    else if ( trim(dlfind_nml%hessini) == "DIAGONAL" ) then
       inithessian = 3
    else if ( trim(dlfind_nml%hessini) == "IDENTITY" ) then
       inithessian = 4
    else
       write(6,'(3a)')"ERROR: Invalid hessini='",trim(dlfind_nml%hessini),"'"
    end if
    
    update = -1
    if ( trim(dlfind_nml%hessupd) == "NONE" ) then
       update = 0
    else if ( trim(dlfind_nml%hessupd) == "POWELL" ) then
       update = 1
    else if ( trim(dlfind_nml%hessupd) == "BOFILL" ) then
       update = 2
    else if ( trim(dlfind_nml%hessupd) == "BFGS" ) then
       update = 3
    else
       write(6,'(3a)')"ERROR: Invalid hessupd='",trim(dlfind_nml%hessupd),"'"
    end if

    if ( dlfind_nml%verbosity < 1 ) then
       printl = 0
       printf = 0
    else if ( dlfind_nml%verbosity == 1 ) then
       printl = 2
       printf = 2
    else if ( dlfind_nml%verbosity == 2 ) then
       printl = 4
       printf = 4
    else
       printl = 6
       printf = 6
    end if

    idimer = 0
    if ( dlfind_nml%dimer ) then
       if ( dlfind_nml%dimer_mode >= 0 .and. dlfind_nml%dimer_mode <= 2 ) then
          idimer = 200 + 10 * dlfind_nml%dimer_mode
       else
          write(6,'(a,i4)')"ERROR: Invalid dimer_mode=",dlfind_nml%dimer_mode
       end if
    end if

    ineb = 0
    if ( active_data%is_neb ) then
       if ( trim(dlfind_nml%neb) == "FREE" ) then
          ineb = 100
       else if ( trim(dlfind_nml%neb) == "PERP" ) then
          ineb = 110
       else if ( trim(dlfind_nml%neb) == "FROZEN" ) then
          ineb = 120
       else
          write(6,'(2a)')"Invalid neb option ",trim(dlfind_nml%neb)
          call mexit(6,1)
       end if
    end if
    ! ineb = 100 if FREE
    ! ineb = 110 if PERP
    ! ineb = 120 if FROZEN
    
    ineb_cart = 0
    if ( active_data%is_neb .and. dlfind_nml%neb_cart ) then
       ineb_cart = 1
    end if
    ! ineb_cart = 1 if neb_cart = .true.

    icoord = icrdtype + ineb + ineb_cart*30 + idimer
       
    maxstep = dlfind_nml%maxstep
    scalestep = dlfind_nml%scalestep
    lbfgs_mem = dlfind_nml%lbfgsmem
    nimage = 0
    if ( active_data%is_neb ) then
       nimage = dlfind_nml%neb_nimage
    end if
    nebk = dlfind_nml%nebk
    dump = -1
    restart = -1
    nz = active_data%nactive
    nconn = 0
    maxupd = dlfind_nml%maxupdate
    delta = dlfind_nml%delta
    soft = 1.d+20
    carthessian = -1
    tsrel = 0
    maxrot = dlfind_nml%dimer_maxrot
    tolrot = dlfind_nml%dimer_tolrot
    nframe = active_data%nframes
    nmass = active_data%nactive
    
    nweight = 0
    if ( associated(active_data%wts) ) then
       nweight = size(active_data%wts)
    end if
    
    timestep = -1.d0
    fric0 = -1.d0
    fricfac = -1.d0
    fricp = -1.d0
    imultistate = 0
    state_i = 1
    state_j = 2
    pf_c1 = 5.d0
    pf_c2 = 5.d0
    gp_c3 = 1.d0
    gp_c4 = 0.9d0
    ln_t1 = 1.d-4
    ln_t2 = 1.d0
    tolerance_e = dlfind_nml%tole
    distort = 0.d0
    minstep = dlfind_nml%minstep
    maxdump = -1
    task = -1
    temperature = -1.d0
    po_pop_size = -1
    po_radius = -1.d0
    po_contraction = -1.d0
    po_tolerance_r = -1.d0
    po_tolerance_g = -1.d0
    po_distribution = -1
    po_maxcycle = -1
    po_init_pop_size = -1
    po_reset = -1
    po_mutation_rate = -1.d0
    po_death_rate = -1.d0
    po_scalefac = -1.d0
    po_nsave = -1
    ntasks = 1
    tdlf_farm = 1
    n_po_scaling = 0
    
    neb_climb_test = dlfind_nml%neb_climb_test
    neb_freeze_test = dlfind_nml%neb_freeze_test
    nzero = -1
    coupled_states = 1
    qtsflag = -1
    imicroiter = 0
    maxmicrocycle = -1
    micro_esp_fit = -1
    
  end subroutine dlf_get_params_sander



  
  
  subroutine dlf_put_coords_sander(nvar, switch, energy, coords, iam) bind(c)
    
    use file_io_dat, only : restrt,mdcrd,ntwx
    
    implicit none
    
    integer(kind=c_int), intent(in), value :: nvar ! number of xyz variables (3*nat)
    integer(kind=c_int), intent(in), value :: switch ! 1: coords contains actual geometry, 2: coords contains transition mode
    real(kind=c_double), intent(in), value :: energy ! energy
    real(kind=c_double), intent(in) :: coords(nvar) ! coordinates
    integer(kind=c_int), intent(in), value :: iam ! flag for MPI runs
    integer, save :: num_steps = 0
    character(len=MASKLEN) :: fname
    
    if ( active_data%is_ts_search ) then
       if ( switch == 1 ) then
          ! ignore
          return
       else if ( switch == 2 ) then
          fname = trim(active_data%base) // "tsmode." // trim(active_data%ext)
          write(6,'(2a)')"Saving tsmode to restart: ",trim(fname)
          call set_crds_from_active_crds(coords(1))
          call write_custom_rst_file(fname)
       else if ( switch == 3 ) then
          num_steps = num_steps + 1
          write(6,'(2a)')"Saving ts to restart: ",trim(restrt)
          call set_crds_from_active_crds(coords(1))
          call write_rst_file()
          if ( ntwx > 0 ) then
             if ( mod(num_steps,ntwx) == 0 ) then
                write(6,'(2a)')"Saving trajectory frame to: ",trim(mdcrd)
                call write_traj_frame()
             end if
          end if
       else if ( switch == 4 ) then
          fname = trim(active_data%base) // "product." // trim(active_data%ext)
          write(6,'(2a)')"Saving product-side minimum to restart: ",trim(fname)
          call set_crds_from_active_crds(coords(1))
          call write_custom_rst_file(fname)
       else if ( switch == 5 ) then
          fname = trim(active_data%base) // "reactant." // trim(active_data%ext)
          write(6,'(2a)')"Saving reactant-side minimum to restart: ",trim(fname)
          call set_crds_from_active_crds(coords(1))
          call write_custom_rst_file(fname)
       else
          write(6,'(a,i5)')"dlf_put_coords: unexpected switch with dimer=T:",&
               & switch
       end if
       return
    end if

    if ( active_data%is_neb ) then
       if ( switch == 1 ) then
          ! ignore
          return
       else if ( switch == 2 ) then
          fname = trim(active_data%base) // "tsmode." // trim(active_data%ext)
          write(6,'(2a)')"Saving tsmode to restart: ",trim(fname)
          call set_crds_from_active_crds(coords(1))
          call write_custom_rst_file(fname)
       else if ( switch == 3 ) then
          fname = trim(active_data%base) // "ts." // trim(active_data%ext)
          write(6,'(2a)')"Saving ts to restart: ",trim(fname)
          call set_crds_from_active_crds(coords(1))
          call write_custom_rst_file(fname)
       else if ( switch < 0 ) then
          write(fname,'(i0.4)')-switch
          fname = trim(active_data%base) // trim(fname) // "." // trim(active_data%ext)
          write(6,'(a,i3,2a)')"Saving img ",-switch," to restart: ",trim(fname)
          call set_crds_from_active_crds(coords(1))
          call write_custom_rst_file(fname)
       else
          write(6,'(a,i5)')"dlf_put_coords: unexpected switch with neb=T:",switch
       end if
       return
    end if
    
    if ( switch == 1 ) then
       num_steps = num_steps + 1
       write(6,'(2a)')"Saving coordinates to restart: ",trim(restrt)
       call set_crds_from_active_crds(coords(1))
       call write_rst_file()
       if ( ntwx > 0 ) then
          if ( mod(num_steps,ntwx) == 0 ) then
             write(6,'(2a)')"Saving trajectory frame to: ",trim(mdcrd)
             call write_traj_frame()
          end if
       end if
    else
       write(6,'(a,i5)')"dlf_put_coords: unexpected switch with dimer=F, neb=F:", &
            & switch
    end if
    
  end subroutine dlf_put_coords_sander



  
  
  subroutine dlf_update_sander() bind(c)

    implicit none

    write(6,*)"ERROR: dlfind_module/dlf_update_sander is unimplemented"
    
  end subroutine dlf_update_sander

  
  
  
  subroutine dlf_error_sander() bind(c)

    implicit none

    write(6,*)"DL-Find encountered an error. Aborting optimization."
    call mexit(6,1)
    
  end subroutine dlf_error_sander

  
  

  subroutine run_dlfind(crds,frcs,ene,qsetup)
    !use iso_c_binding, only: c_funloc
    !use iso_c_binding, only: c_int
    use iso_c_binding
    
    implicit none
    
    !_REAL_,  intent(inout) :: xx(*)          ! real dynamic memory
    !integer, intent(inout) :: ix(*)          ! integer dynamic memory
    !character(len=4), intent(inout) :: ih(*) ! hollerith dynamic memory
    !integer, intent(inout) :: ipairs(*)      ! nonbond pair list dynamic memory
    _REAL_,  intent(inout) :: crds(*)
    _REAL_,  intent(inout) :: frcs(*)
    type(state_rec),  intent(inout) :: ene
    logical,  intent(inout) :: qsetup

    integer(kind=c_int) :: nvarin
    integer(kind=c_int) :: nvarin2
    integer(kind=c_int) :: nspec
    integer(kind=c_int) :: imaster

    interface
       
       subroutine api_dl_find(nvarin, nvarin2, nspec, master, &
            & dlf_error_c, dlf_get_gradient_c, dlf_get_hessian_c, &
            dlf_get_multistate_gradients_c, dlf_get_params_c, &
            & dlf_put_coords_c, dlf_update_c) bind(c)
         use iso_c_binding, only: c_int, c_double, c_funptr, c_f_procpointer

         implicit none
         integer(kind=c_int), intent(in), value :: nvarin ! number of variables to read in 3*nat
         integer(kind=c_int), intent(in), value :: nvarin2 ! number of variables to read in in the second array (coords2)
         integer(kind=c_int), intent(in), value :: nspec ! number of values in the integer array spec
         integer(kind=c_int), intent(in), value :: master ! 1 if this task is the master of a parallel run, 0 otherwise
         type(c_funptr), intent(in), value :: dlf_error_c, dlf_get_gradient_c, &
              & dlf_get_hessian_c, dlf_get_multistate_gradients_c, &
              & dlf_get_params_c, dlf_put_coords_c, dlf_update_c ! Functions received from C side

       end subroutine api_dl_find

    end interface

    
    call read_dlfind_nml(dlfind_nml)

    call new_active_data(dlfind_nml,active_data,qsetup)

    nvarin = 3 * active_data%nactive
    nvarin2 = active_data%nvar2
    nspec = active_data%nspec
    imaster = 1

    if ( active_data%is_hess ) then

#ifdef BINTRAJ
       call run_hessian_calc()
#else
       write(6,'(2a)')"Cannot use optalg='HESS' because amber was ", &
            & "not compiled with netcdf support"
       call mexit(6,1)
#endif

    else
       
       call enforce_shake()
    
       call api_dl_find( nvarin, nvarin2, nspec, imaster, &
            & c_funloc(dlf_error_sander), &
            & c_funloc(dlf_get_gradient_sander), &
            & c_funloc(dlf_get_hessian_sander), &
            & c_funloc(dlf_get_multistate_gradients_sander), &
            & c_funloc(dlf_get_params_sander), &
            & c_funloc(dlf_put_coords_sander), &
            & c_funloc(dlf_update_sander) )
       
       ene = active_data%ene

       if ( dlfind_nml%verbosity < 1 .and. .not. active_data%is_neb ) then
          call write_rst_file()
       end if
       
    end if
    
    
  end subroutine run_dlfind


  


  
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



  subroutine PrintVMDSelection(natom,nres,amask,residue_pointer)
    implicit none
    integer, intent(in) :: natom
    integer, intent(in) :: nres
    integer, intent(in) :: amask(natom)
    integer, intent(in) :: residue_pointer(nres+1)


    integer :: i,j
    logical :: full_match, part_match
    
    integer,pointer :: anyresmask(:) => null()
    integer,pointer :: fullresmask(:) => null()
    integer,pointer :: patommask(:) => null()
    
    integer,pointer :: res_from(:) => null()
    integer,pointer :: res_to(:) => null()
    integer,pointer :: atom_from(:) => null()
    integer,pointer :: atom_to(:) => null()

    integer :: resrange
    integer :: atomrange

    character(len=8) :: s1,s2
    
    ! vmd selection
    allocate( fullresmask(nres) )
    fullresmask = 0

    allocate( anyresmask(nres) )
    anyresmask = 0

    allocate( patommask(natom) )
    patommask = 0

    allocate( res_from(nres) )
    res_from = 0

    allocate( res_to(nres) )
    res_to = 0

    allocate( atom_from(natom) )
    atom_from = 0

    allocate( atom_to(natom) )
    atom_to = 0


    do i=1,nres
       
       full_match = .true.
       part_match = .false.
       do j=residue_pointer(i),residue_pointer(i+1)-1
          if ( amask(j) == 0 ) then
             full_match = .false.
          else
             part_match = .true.
             anyresmask(i) = 1
          end if
       end do

       if ( full_match ) then
          fullresmask(i) = 1
       else if ( part_match ) then
          do j=residue_pointer(i),residue_pointer(i+1)-1
             patommask(j) = amask(j)
          end do
       end if
       
    end do

    resrange = 0
    do i=1,nres
       if ( fullresmask(i) > 0 ) then
          resrange = resrange + 1
          res_from(resrange) = i
          res_to(resrange) = i
          exit
       end if
    end do
    if ( resrange > 0 ) then
       do i=res_from(resrange)+1,nres
          if ( fullresmask(i) > 0 ) then
             if ( fullresmask(i-1) > 0 ) then
                res_to(resrange) = i
             else
                resrange = resrange + 1
                res_from(resrange) = i
                res_to(resrange) = i
             end if
          end if
       end do
    end if



    atomrange = 0
    do i=1,natom
       if ( patommask(i) > 0 ) then
          atomrange = atomrange + 1
          atom_from(atomrange) = i
          atom_to(atomrange) = i
          exit
       end if
    end do
    if ( atomrange > 0 ) then
       do i=atom_from(atomrange)+1,natom
          if ( patommask(i) > 0 ) then
             if ( patommask(i-1) > 0 ) then
                atom_to(atomrange) = i
             else
                atomrange = atomrange + 1
                atom_from(atomrange) = i
                atom_to(atomrange) = i
             end if
          end if
       end do
    end if
    
    
    write(6,'(/a)')"&dlfind active selection in VMD mask format:"

    if ( resrange > 0 ) then
       write(6,'(a)',advance="NO")"resid"
       do i=1,resrange
          write(s1,'(i8)')res_from(i)
          write(s2,'(i8)')res_to(i)
          if ( res_from(i) == res_to(i) ) then
             write(6,'(1x,a)',advance="NO")trim(adjustl(s1))
          else
             write(6,'(1x,3a)',advance="NO")trim(adjustl(s1))," to ",trim(adjustl(s2))
          end if
       end do
    end if

    if ( atomrange > 0 ) then
       if ( resrange > 0 ) then
          write(6,'(a)',advance="NO")" or "
       end if
       write(6,'(a)',advance="NO")"serial"
       do i=1,atomrange
          write(s1,'(i8)')atom_from(i)
          write(s2,'(i8)')atom_to(i)
          if ( atom_from(i) == atom_to(i) ) then
             write(6,'(1x,a)',advance="NO")trim(adjustl(s1))
          else
             write(6,'(1x,3a)',advance="NO")trim(adjustl(s1))," to ",trim(adjustl(s2))
          end if
       end do
    end if
    
    write(6,'(//a)')"&dlfind active selection in Amber mask format:"
    
    if ( resrange > 0 ) then
       write(6,'(a)',advance="NO")":"
       do i=1,resrange
          write(s1,'(i8)')res_from(i)
          write(s2,'(i8)')res_to(i)
          if ( i > 1 ) then
             write(6,'(a)',advance="NO")","
          end if
          if ( res_from(i) == res_to(i) ) then
             write(6,'(a)',advance="NO")trim(adjustl(s1))
          else
             write(6,'(3a)',advance="NO")trim(adjustl(s1)),"-",trim(adjustl(s2))
          end if
       end do
    end if

    if ( atomrange > 0 ) then
       if ( resrange > 0 ) then
          write(6,'(a)',advance="NO")"|"
       end if
       write(6,'(a)',advance="NO")"@"
       do i=1,atomrange
          write(s1,'(i8)')atom_from(i)
          write(s2,'(i8)')atom_to(i)
          if ( i > 1 ) then
             write(6,'(a)',advance="NO")","
          end if
          if ( atom_from(i) == atom_to(i) ) then
             write(6,'(a)',advance="NO")trim(adjustl(s1))
          else
             write(6,'(3a)',advance="NO")trim(adjustl(s1)),"-",trim(adjustl(s2))
          end if
       end do
    end if
    write(6,'(a)')""



    resrange = 0
    do i=1,nres
       if ( anyresmask(i) > 0 ) then
          resrange = resrange + 1
          res_from(resrange) = i
          res_to(resrange) = i
          exit
       end if
    end do
    if ( resrange > 0 ) then
       do i=res_from(resrange)+1,nres
          if ( anyresmask(i) > 0 ) then
             if ( anyresmask(i-1) > 0 ) then
                res_to(resrange) = i
             else
                resrange = resrange + 1
                res_from(resrange) = i
                res_to(resrange) = i
             end if
          end if
       end do
    end if

    
    write(6,'(/a)')"Residues with at least 1 active atom in VMD mask format:"

    if ( resrange > 0 ) then
       write(6,'(a)',advance="NO")"resid"
       do i=1,resrange
          write(s1,'(i8)')res_from(i)
          write(s2,'(i8)')res_to(i)
          if ( res_from(i) == res_to(i) ) then
             write(6,'(1x,a)',advance="NO")trim(adjustl(s1))
          else
             write(6,'(1x,3a)',advance="NO")trim(adjustl(s1))," to ",trim(adjustl(s2))
          end if
       end do
    end if


    write(6,'(//a)')"Residues with at least 1 active atom in Amber mask format:"
    
    if ( resrange > 0 ) then
       write(6,'(a)',advance="NO")":"
       do i=1,resrange
          write(s1,'(i8)')res_from(i)
          write(s2,'(i8)')res_to(i)
          if ( i > 1 ) then
             write(6,'(a)',advance="NO")","
          end if
          if ( res_from(i) == res_to(i) ) then
             write(6,'(a)',advance="NO")trim(adjustl(s1))
          else
             write(6,'(3a)',advance="NO")trim(adjustl(s1)),"-",trim(adjustl(s2))
          end if
       end do
    end if
    write(6,'(/)')

    
    if ( associated(anyresmask) ) then
       deallocate(anyresmask)
    end if

    if ( associated(fullresmask) ) then
       deallocate(fullresmask)
    end if

    if ( associated(patommask) ) then
       deallocate(patommask)
    end if

    if ( associated(res_from) ) then
       deallocate(res_from)
    end if

    if ( associated(res_to) ) then
       deallocate(res_to)
    end if

    if ( associated(atom_from) ) then
       deallocate(atom_from)
    end if

    if ( associated(atom_to) ) then
       deallocate(atom_to)
    end if
    
  end subroutine PrintVMDSelection


  subroutine ReportNoshakeMaskModification(natom,nres,amask,residue_pointer)
    implicit none
    integer,intent(in) :: natom
    integer,intent(in) :: nres
    integer,intent(in) :: amask(natom)
    integer,intent(in) :: residue_pointer(nres+1)

    integer :: i,j,resrange
    
    integer,pointer :: fullresmask(:) => null()    
    integer,pointer :: res_from(:) => null()
    integer,pointer :: res_to(:) => null()


    character(len=8) :: s1,s2
    
    ! vmd selection
    allocate( fullresmask(nres) )
    fullresmask = 0

    allocate( res_from(nres) )
    res_from = 0

    allocate( res_to(nres) )
    res_to = 0


    do i=1,nres
       do j=residue_pointer(i),residue_pointer(i+1)-1
          if ( amask(j) > 0 ) then
             fullresmask(i) = 1
             exit
          end if
       end do
    end do

    resrange = 0
    do i=1,nres
       if ( fullresmask(i) > 0 ) then
          resrange = resrange + 1
          res_from(resrange) = i
          res_to(resrange) = i
          exit
       end if
    end do
    if ( resrange > 0 ) then
       do i=res_from(resrange)+1,nres
          if ( fullresmask(i) > 0 ) then
             if ( fullresmask(i-1) > 0 ) then
                res_to(resrange) = i
             else
                resrange = resrange + 1
                res_from(resrange) = i
                res_to(resrange) = i
             end if
          end if
       end do
    end if


    write(6,'(/a)')"The &dlfind active mask severs shaked bonds."
    write(6,'(a)')"You should take one of the following actions."
    write(6,'(a)')"   Set ntc=1 in &cntrl"
    write(6,'(a)')"   -or- set ntc>=2 and include the following mask within &cntrl noshakemask"
    write(6,'(a/)')"   -or- set ntc=2 and include the following mask within &dlfind active"

    write(6,'(a)',advance="NO")":"
    do i=1,resrange
       write(s1,'(i8)')res_from(i)
       write(s2,'(i8)')res_to(i)
       if ( i > 1 ) then
          write(6,'(a)',advance="NO")","
       end if
       if ( res_from(i) == res_to(i) ) then
          write(6,'(a)',advance="NO")trim(adjustl(s1))
       else
          write(6,'(3a)',advance="NO")trim(adjustl(s1)),"-",trim(adjustl(s2))
       end if
    end do
    write(6,'(//)')
    
    
    if ( associated(fullresmask) ) then
       deallocate(fullresmask)
    end if

    if ( associated(res_from) ) then
       deallocate(res_from)
    end if

    if ( associated(res_to) ) then
       deallocate(res_to)
    end if

    call mexit(6,1)

  end subroutine ReportNoshakeMaskModification

  integer function CheckUnselectedCrds(natom,amask,crds)
    use memory_module, only : coordinate
    implicit none
    integer,intent(in) :: natom
    integer,intent(in) :: amask(natom)
    _REAL_,intent(in) :: crds(3,natom)

    integer :: i,j
    _REAL_ :: dx,dr2
    _REAL_,parameter :: TOL=1.d-12

    CheckUnselectedCrds = 0
    do i=1,natom
       if ( amask(i) == 0 ) then
          dr2 = 0.d0
          do j=1,3
             dx = crds(j,i) - coordinate(j,i)
             dr2 = dr2 + dx*dx
          end do
          if ( dr2 > TOL ) then
            CheckUnselectedCrds = CheckUnselectedCrds + 1
          end if
       end if
    end do
    
  end function CheckUnselectedCrds



  subroutine ReadDisang(fname,ncon,condata)
    implicit none

    character(len=*),intent(in) :: fname

    integer,intent(out) :: ncon
    integer,pointer :: condata(:,:)
    
    !------------------------------------------------------------------------------

    integer :: nmrnum
    integer :: i,ii
    logical :: isvalid

    integer,parameter :: maxcdata = 500
    integer,allocatable :: cdata(:,:)
    
    !------------------------------------------------------------------------------
    integer :: iin,ifind
    
    integer,parameter :: maxigr = 200
    _REAL_, parameter :: rstwttol = 1.d-7
    _REAL_, parameter :: ZERO = 0.d0
    
    integer :: iat(8)
    integer :: iat1,iat2,iat3,iat4,iat5,iat6,iat7,iat8
    character(256) :: restraint
    _REAL_  :: rstwt(4)
    integer :: nstep1,nstep2,irstyp,ninc
    integer :: iresid,imult
    
    character(len=4) :: atnam(8),grnam1(maxigr),grnam2(maxigr), &
         & grnam3(maxigr),grnam4(maxigr),grnam5(maxigr), &
         & grnam6(maxigr),grnam7(maxigr),grnam8(maxigr) !,grnamarr,grnampass

    integer :: igr1(maxigr),igr2(maxigr),igr3(maxigr),igr4(maxigr), &
         & igr5(maxigr),igr6(maxigr),igr7(maxigr),igr8(maxigr)

    _REAL_  :: r0,r1,r2,r3,r4,k0,k2,k3,r0a,r1a,r2a,r3a,r4a,k0a,rk2a,rk3a,rk2,rk3
    integer :: ir6,ifntyp,ifvari
    _REAL_  :: rjcoef(3)
    integer :: ixpk,nxpk,ialtd,iconstr,fxyz(5),outxyz
    
    
    equivalence (iat(1),iat1),(iat(2),iat2)
    equivalence (iat(3),iat3),(iat(4),iat4)
    equivalence (iat(5),iat5),(iat(6),iat6)
    equivalence (iat(7),iat7),(iat(8),iat8)
    
    namelist /rst/ iat,restraint,rstwt,nstep1,nstep2,irstyp,ninc, &
         iresid,imult,atnam,igr1,igr2,igr3,igr4,igr5,igr6,igr7,igr8, &
         grnam1,grnam2,grnam3,grnam4,grnam5,grnam6,grnam7,grnam8,&
         r0,r1,r2,r3,r4,k0,rk2,rk3,r0a,r1a,r2a,r3a,r4a,k0a,rk2a,rk3a,ir6,ifntyp, &
         ifvari,rjcoef,ixpk,nxpk,ialtd,iconstr,fxyz,outxyz
    !------------------------------------------------------------------------------

    write(6,'(/2a)')"Reading constraints from ",trim(fname)
    
    iin = 62
    call amopen(iin,fname,"O","F","R")
    

    allocate( cdata(5,maxcdata) )
    cdata = 0
    
    nmrnum=0
    
    do i=1,99999
       do ii=1,8
          iat(ii) = 0
          atnam(ii) = '    '
       end do
       do ii = 1,4
          rstwt(ii) = 0.0
       end do
       restraint = ' '
       do ii = 1,3
          rjcoef(ii) = ZERO
       end do
       do ii=1,maxigr
          igr1(ii) = 0
          igr2(ii) = 0
          igr3(ii) = 0
          igr4(ii) = 0
          igr5(ii) = 0
          igr6(ii) = 0
          igr7(ii) = 0
          igr8(ii) = 0
          grnam1(ii) = '    '
          grnam2(ii) = '    '
          grnam3(ii) = '    '
          grnam4(ii) = '    '
          grnam5(ii) = '    '
          grnam6(ii) = '    '
          grnam7(ii) = '    '
          grnam8(ii) = '    '
       end do
       iconstr = 0
       ixpk = 0
       nxpk = 0
       nstep1 = 0
       nstep2 = 0
       irstyp = 0
       ifvari = 0
       iresid = 0
       ifntyp = 0
       r0 = ZERO
       r0a = ZERO
       k0 = ZERO
       k0a = ZERO
       fxyz = 0
       
       ! Look for "rst" namelist. If not found, assume we are done reading them
       
       isvalid = .true.
       call nmlsrc('rst',iin,ifind)
       if (ifind == 0) goto 227
       read (iin,nml=rst,end=227)
      
227    if (restraint /= ' ' .and. iat1 ==0) then
          write(6,'(2a)')"WARNING: Ignoring natural language constraint. Use &rst iat variable",&
               & " to define constrained coordinates"
          isvalid = .false.
       else if (restraint /= ' ' .and. iat1 /=0) then
          write(6,'(2a)')"WARNING: Natural language constraint and iat variables used in &rst",&
               & "The natural language definition will be ignored."
          isvalid = .true.
       end if

       
       !If IAT1=0, assume the last restraint has been read.
       
       if (iat1 == 0) then
          isvalid = .false.
          goto 100
       end if
       
       if (fxyz(1) .ne. 0 .or. fxyz(2) .ne. 0 .or. fxyz(3) .ne. 0 &
            & .or. fxyz(4) .ne. 0 .or. fxyz(5) .ne. 0 ) then
          write(6,'(a)')"Skipping constraint; fxyz(:)/=0 is unsupported"
          isvalid = .false.
          cycle
       end if
       
       do ii=5,8
          if ( iat(ii) /= 0 ) then
             write(6,'(a,i1,a)')"Skipping constraint; iat(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
       end do
       do ii=1,4
          if ( abs(rstwt(ii)) > rstwttol ) then
             write(6,'(a,i1,a)')"Skipping constraint; abs(rstwt(",ii,"))>0 is unsupported"
             isvalid = .false.
          end if
       end do
       do ii=1,maxigr
          if ( igr1(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr1(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr2(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr2(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr3(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr3(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr4(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr4(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr5(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr5(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr6(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr6(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr7(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr7(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
          if ( igr8(ii) /= 0 ) then
             write(6,'(a,i6,a)')"Skipping constraint; igr8(",ii,")/=0 is unsupported"
             isvalid = .false.
          end if
       end do

       if ( iconstr /= 0 ) then
          write(6,'(a)')"Skipping constraint; icnstr/=0 is unsupported"
          isvalid=.false.
       end if

       if ( ixpk /= 0 ) then
          write(6,'(a)')"Skipping constraint; ixpk/=0 is unsupported"
          isvalid=.false.
       end if

       if ( nxpk /= 0 ) then
          write(6,'(a)')"Skipping constraint; nxpk/=0 is unsupported"
          isvalid=.false.
       end if

       if ( nstep1 /= 0 ) then
          write(6,'(a)')"Skipping constraint; nstep1/=0 is unsupported"
          isvalid=.false.
       end if

       if ( nstep2 /= 0 ) then
          write(6,'(a)')"Skipping constraint; nstep2/=0 is unsupported"
          isvalid=.false.
       end if

       if ( irstyp /= 0 ) then
          write(6,'(a)')"Skipping constraint; irstyp/=0 is unsupported"
          isvalid=.false.
       end if
       
       if ( ifvari /= 0 ) then
          write(6,'(a)')"Skipping constraint; ifvari/=0 is unsupported"
          isvalid=.false.
       end if

       if ( iresid /= 0 ) then
          write(6,'(a)')"Skipping constraint; iresid/=0 is unsupported"
          isvalid=.false.
       end if

       if ( ifntyp /= 0 ) then
          write(6,'(a)')"Skipping constraint; ifntyp/=0 is unsupported"
          isvalid=.false.
       end if

       if ( iat(2) < 1 ) then
          write(6,'()')"Skipping constraint; iat(2) should be > 0"
          isvalid=.false.
       end if
       
       if ( isvalid ) then
          nmrnum = nmrnum + 1

          if ( nmrnum > maxcdata ) then
             write(6,'(a,i5,a)')"ERROR: A maximum of ",maxcdata," constraints can be read"
             call mexit(6,1)
          end if
          
          if ( iat3 > 0 ) then
             if ( iat4 > 0 ) then
                cdata(1:5,nmrnum) = (/3,iat1,iat2,iat3,iat4/)
             else
                cdata(1:5,nmrnum) = (/2,iat1,iat2,iat3,0/)
             end if
          else
             cdata(1:5,nmrnum) = (/1,iat1,iat2,0,0/)
          end if
          
       end if
    end do
    
100 continue


    if ( associated(condata) ) then
       deallocate(condata)
    end if

    ncon = nmrnum
    allocate(condata(5,ncon))
    condata(1:5,1:ncon) = cdata(1:5,1:ncon)
    deallocate( cdata )

    close(unit=iin)

    write(6,'(a,i5,2a/)')"Read ",ncon," constraint definitions from ",trim(fname)
    
  end subroutine ReadDisang

  
  
end module dlfind_module
