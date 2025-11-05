! <compile=optimized>
#include "../include/dprec.fh"
module qmmm_fires_module
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Implementation of the Flexible Inner Region Ensemble Separator (FIRES)
!
!References: 
!Dmitrii Beglov, Benoît Roux; Finite representation of an infinite bulk system:
!    Solvent boundary potential for computer simulations. 
!    J. Chem. Phys. 15 June 1994; 100 (12): 9050–9063.
!Rowley CN, Roux B. The Solvation Structure of Na(+) and K(+) in Liquid
!    Water Determined from High Level ab Initio Molecular Dynamics Simulations. 
!    J Chem Theory Comput. 2012
!Implemented by:
!    Date: 8/15/2024
!    Benoît Roux (UChicago) roux@uchicago.edu
!    Jose Guerra (UChicago) jlguerra@uchicago.edu
! Inpute File Information:
!    fires = 0 (FIRES off), fires = 1 (FIRES on)
!    fires_inner = (atom indices of the inner region)
!    fires_outer = (atom indices of the outer region)
!    fires_k = (force constant for the boundary potential)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    use findmask
    implicit none
#ifdef MPI
#  include "parallel.h"
  include 'mpif.h'
#endif
    private
    public :: FireMasks, calculate_fires_force, fires_force, setup_fires, fires_set_local_bounds, fires_restraint_enabled, fires_prepare_if_needed, fires_mask_signature, fires_set_step
    ! Debug/diagnostics: last computed FIRES force components
    _REAL_, public, save :: last_f_pot(3)           = 0.0d0
    _REAL_, public, save :: last_f_inner_water(3)   = 0.0d0
    _REAL_, public, save :: last_f_pen_water(3)     = 0.0d0
    integer, public, save :: fires = 0
    integer, public, save :: last_n_pen = 0              ! diagnostics: # of penetrating outer atoms this eval
    integer, public, save :: last_n_inw_out = 0          ! diagnostics: # of inner-solvent atoms outside rmax
    _REAL_,  public, save :: last_rmax = 0.0d0           ! diagnostics: boundary radius used this eval
    integer, public, save :: last_anchor_idx = 0         ! diagnostics: anchor atom index used this eval
    _REAL_,  public, save :: last_rmax_inw = 0.0d0       ! diagnostics: farthest inner-solvent distance
    integer, public, save :: last_inw_max_idx = 0         ! diagnostics: atom index for farthest inner-solvent
    integer, public, save :: last_pen_atom_idx = 0        ! diagnostics: most inward penetrating outer atom index
    _REAL_,  public, save :: last_pen_r = 0.0d0           ! diagnostics: its distance to COM
    _REAL_,  public, save :: last_pen_depth = 0.0d0       ! diagnostics: penetration depth (rmax - r)
    !JOSE MOD: MTS controls
    !_REAL_ small-step size in ps; when > 0 and mts_n>1, subcycle FIRES with dt_fast=mts_fires
    _REAL_,  public, save :: mts_fires = 0.0d0
    ! number of FIRES microsteps per outer step
    integer, public, save :: mts_n = 1
    ! Optional band filter: when enabled, only atoms with |r - rmax| <= halfwidth are considered
    logical, public, save :: fires_band_filter_enabled = .false.
    _REAL_,  public, save :: fires_band_halfwidth = 6.0d0
    ! Freezing controls
    logical, public, save :: fires_static_mask     = .false.
    integer, public, save :: fires_freeze_on_nstep = 1
    integer, save :: fires_current_step            = 0  ! cached md step for freeze gating
    ! Hysteresis control to stabilize mask membership near the boundary
    logical, public, save :: fires_hyst_enabled    = .false.
    _REAL_,  public, save :: fires_hyst_rin  = 3.9d0
    _REAL_,  public, save :: fires_hyst_rout = 4.1d0
    ! Track whether the last call to force() already injected FIRES into the total force
    logical, public, save :: fires_in_force = .false.
    integer, dimension(:), allocatable :: fire_inner_mask, fire_outer_mask, fire_inner_solvent_mask
    integer, dimension(:), allocatable, save :: fire_inner_mask_prev, fire_inner_solvent_mask_prev
    ! Removed unused module-level natom/nres to avoid name collisions with COMMON from memory.h
    character(len=256), public, save :: fire_inner = ""
    character(len=256), public, save :: fire_inner_solvent = ""
    character(len=256), public, save :: fire_outer = ""
    ! (deleted unused ipres/igraph/isymbl/lbres/crd/nsc/inner_atoms/outer_atoms)
    ! Control whether to auto-refresh FIRES masks at runtime (each step). Default: off
    logical, public, save :: fires_refresh_runtime = .false.

    ! Boundary radius control: 0 = derive from inner solute extent (default)
    !                           1 = use constant radius (fires_rmax_const) from inner COM
    integer, public, save :: fires_rmax_mode = 0
    _REAL_,  public, save :: fires_rmax_const = 0.0d0
    _REAL_,  public, save :: fires_rmax_pad   = 0.0d0
    _REAL_,  public, save :: fires_rmax_min   = 0.0d0
    ! Optional smoothing of rmax to avoid "buzz" when the farthest atom flips (1.0 = off)
    _REAL_,  public, save :: fires_rmax_ema_alpha = 1.0d0
    _REAL_,  save        :: rmax_ema = -1.0d0

    ! FIRES parameters
    _REAL_, public, save :: fires_k = 0.0D0
    _REAL_, parameter :: EPS = 1.0D-12
    _REAL_, parameter :: MASS_EPS = 1.0D-12

    type :: FireMasks
      integer, dimension(:), allocatable :: inner_mask
      integer, dimension(:), allocatable :: inner_solvent_mask
      integer, dimension(:), allocatable :: outer_mask
    end type FireMasks
    type(FireMasks) :: masks

    !JOSE MOD: guard to avoid double-initialization/allocations in setup_fires
    logical, save :: fires_initialized = .false.
    ! Cache sizes so we can refresh masks later if needed
    integer, save :: fires_natom = 0
    integer, save :: fires_nres  = 0
    ! Ensure distance-based masks are refreshed once after coords are available
    logical, save :: fires_need_first_refresh = .true.
    logical, save :: fires_notice_printed = .false.
    ! Local slice (this-rank) bounds for FIRES force writes (3*indices)
    integer, save :: g_loc_istart3 = 1
    integer, save :: g_loc_iend3   = -1
    ! Center control: 0 = use COM of inner mask; >0 = use this atom index as center
    integer, public, save :: fires_center_atom = 0

contains

    subroutine setup_fires(natom_in, nres_in, ix_in, ih_in, x_in, nbonh_in, nbona_in)
        ! Use global topology/coordinates so atommask has real data
    use memory_module, only: atom_name, amber_atom_type, residue_label, residue_pointer, x
        implicit none
#include "../include/memory.h"
#include "../include/md.h"
        integer, intent(in) :: natom_in, nres_in, nbonh_in, nbona_in
        character(len=4) ih_in(*)
        integer, intent(in) :: ix_in(*)
        _REAL_, intent(in) :: x_in(*)
    integer :: num_atoms, i, inner_count, outer_count,counter, n_overlap
        character(len=256) :: s_inner, s_inw, s_outer
        integer :: prnlev
#ifdef MPI
        integer :: ierr
        integer :: mask_sizes(3)
#endif
        
        !JOSE MOD: if already initialized, skip to avoid re-allocation crash
        if (fires_initialized) then
            print *, 'FIRES: setup_fires already initialized; skipping reinit.'
            return
        end if

    num_atoms = natom_in  ! Use natom as the number of atoms
    fires_natom = natom_in
    fires_nres  = nres_in

    ! Note: The mask engine (findmask) requires valid topology arrays:
    ! atom_name (igraph), amber_atom_type (isymbl), residue_label (lbres)
    ! residue_pointer (ipres with size nres+1)
    ! coordinates as a flat 3*natom array (x)
    ! These are provided globally via memory_module after memory_init.

    !JOSE MOD: allocate mask buffers only if not already allocated
        if (.not. allocated(fire_inner_mask))          allocate(fire_inner_mask(natom_in))
        if (.not. allocated(fire_inner_solvent_mask))  allocate(fire_inner_solvent_mask(natom_in))
        if (.not. allocated(fire_outer_mask))          allocate(fire_outer_mask(natom_in))
    ! Prepare local copies (avoid side-effects of appending ';' to module strings)
    s_inner = trim(fire_inner)
    s_inw   = trim(fire_inner_solvent)
    s_outer = trim(fire_outer)
    prnlev = 0
    ! Evaluate masks using global topology/coordinates
    call atommask(natom_in, nres_in, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, x(lcrd), s_inner,           fire_inner_mask)
    call atommask(natom_in, nres_in, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, x(lcrd), s_inw,             fire_inner_solvent_mask)
    call atommask(natom_in, nres_in, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, x(lcrd), s_outer,           fire_outer_mask)
    !JOSE MOD: ensure mask arrays are (re)allocated with correct sizes
    if (allocated(masks%inner_mask))         deallocate(masks%inner_mask)
    if (allocated(masks%inner_solvent_mask)) deallocate(masks%inner_solvent_mask)
    if (allocated(masks%outer_mask))         deallocate(masks%outer_mask)
    allocate(masks%inner_mask( sum(fire_inner_mask(1:natom_in))))
    allocate(masks%inner_solvent_mask( sum(fire_inner_solvent_mask(1:natom_in))))
    allocate(masks%outer_mask( sum(fire_outer_mask(1:natom_in))))
        counter=1
        do i = 1, natom_in
            if (fire_inner_mask(i) == 1) then
                masks%inner_mask(counter) = i
                counter = counter + 1
            end if
        end do
        counter=1
        do i = 1, natom_in
            if (fire_inner_solvent_mask(i) == 1) then
                masks%inner_solvent_mask(counter) = i
                counter = counter + 1
            end if
        end do
        counter=1
        do i = 1, natom_in
            if (fire_outer_mask(i) == 1) then
                masks%outer_mask(counter) = i
                counter = counter + 1
            end if
        end do

#ifdef MPI
        if (numtasks > 1) then
            call mpi_bcast(fire_inner_mask, natom_in, MPI_INTEGER, 0, commsander, ierr)
            call mpi_bcast(fire_inner_solvent_mask, natom_in, MPI_INTEGER, 0, commsander, ierr)
            call mpi_bcast(fire_outer_mask, natom_in, MPI_INTEGER, 0, commsander, ierr)

            mask_sizes(1) = size(masks%inner_mask)
            mask_sizes(2) = size(masks%inner_solvent_mask)
            mask_sizes(3) = size(masks%outer_mask)
            call mpi_bcast(mask_sizes, 3, MPI_INTEGER, 0, commsander, ierr)

            if (allocated(masks%inner_mask)) then
                if (size(masks%inner_mask) /= mask_sizes(1)) then
                    deallocate(masks%inner_mask)
                    allocate(masks%inner_mask(mask_sizes(1)))
                end if
            else
                allocate(masks%inner_mask(mask_sizes(1)))
            end if

            if (allocated(masks%inner_solvent_mask)) then
                if (size(masks%inner_solvent_mask) /= mask_sizes(2)) then
                    deallocate(masks%inner_solvent_mask)
                    allocate(masks%inner_solvent_mask(mask_sizes(2)))
                end if
            else
                allocate(masks%inner_solvent_mask(mask_sizes(2)))
            end if

            if (allocated(masks%outer_mask)) then
                if (size(masks%outer_mask) /= mask_sizes(3)) then
                    deallocate(masks%outer_mask)
                    allocate(masks%outer_mask(mask_sizes(3)))
                end if
            else
                allocate(masks%outer_mask(mask_sizes(3)))
            end if

            if (mask_sizes(1) > 0) then
                call mpi_bcast(masks%inner_mask, mask_sizes(1), MPI_INTEGER, 0, commsander, ierr)
            end if
            if (mask_sizes(2) > 0) then
                call mpi_bcast(masks%inner_solvent_mask, mask_sizes(2), MPI_INTEGER, 0, commsander, ierr)
            end if
            if (mask_sizes(3) > 0) then
                call mpi_bcast(masks%outer_mask, mask_sizes(3), MPI_INTEGER, 0, commsander, ierr)
            end if
        end if
#endif
        write(6,'(a)')    'FIRES: Using masks:'
        write(6,'(a,1x,a)') '  fire_inner         =', trim(s_inner)
        write(6,'(a,1x,a)') '  fire_inner_solvent =', trim(s_inw)
        write(6,'(a,1x,a)') '  fire_outer         =', trim(s_outer)
        write(6,'(a,i0)') 'FIRES: Inner region size: ', size(masks%inner_mask)
        write(6,'(a,i0)') 'FIRES: Inner solvent size: ', size(masks%inner_solvent_mask)
        write(6,'(a,i0)') 'FIRES: Outer region size: ', size(masks%outer_mask)
        write(6,'(a,1pe12.4)') 'FIRES: fires_k = ', fires_k
        write(6,'(a,1pe12.4,a,i0)') 'FIRES: mts_fires = ', mts_fires, ', mts_n = ', mts_n
        if (size(masks%inner_mask)==0 .or. size(masks%outer_mask)==0) then
            write(6,*) 'FIRES WARNING: One or more masks are empty; FIRES forces will be zero.'
        end if
        if (fires_k <= 0.0d0) then
            write(6,*) 'FIRES WARNING: fires_k <= 0; FIRES forces will be zero.'
        end if

        ! NOTE: Mass printing removed here to avoid dependency on undefined symbol 'amass'.
        ! Masses are now printed (when needed) inside calculate_fires_force() where the
        ! per-atom mass array is available via the 'mass' dummy argument.

    ! One-time mask overlap diagnostic: warn if inner and outer selections overlap
        n_overlap = 0
        do i = 1, natom_in
            if (fire_inner_mask(i) == 1 .and. fire_outer_mask(i) == 1) n_overlap = n_overlap + 1
        end do
        if (n_overlap > 0) then
            write(6,'(a,i0)') 'FIRES WARNING: inner ∩ outer overlap atoms = ', n_overlap
        end if

        ! If selection looks suspicious, re-run with prnlev=1 to print parser state
        if ( size(masks%outer_mask) == 0 .or. size(masks%inner_solvent_mask) > natom_in/2 ) then
            write(6,'(a)') 'FIRES NOTICE: Re-evaluating masks with prnlev=1 for additional details.'
            prnlev = 1
            call atommask(natom_in, nres_in, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, x(lcrd), s_inner,           fire_inner_mask)
            call atommask(natom_in, nres_in, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, x(lcrd), s_inw,             fire_inner_solvent_mask)
            call atommask(natom_in, nres_in, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, x(lcrd), s_outer,           fire_outer_mask)
            if (allocated(masks%inner_mask))         deallocate(masks%inner_mask)
            if (allocated(masks%inner_solvent_mask)) deallocate(masks%inner_solvent_mask)
            if (allocated(masks%outer_mask))         deallocate(masks%outer_mask)
            allocate(masks%inner_mask( sum(fire_inner_mask(1:natom_in))))
            allocate(masks%inner_solvent_mask( sum(fire_inner_solvent_mask(1:natom_in))))
            allocate(masks%outer_mask( sum(fire_outer_mask(1:natom_in))))
            counter=1
            do i = 1, natom_in
                if (fire_inner_mask(i) == 1) then
                    masks%inner_mask(counter) = i
                    counter = counter + 1
                end if
            end do
            counter=1
            do i = 1, natom_in
                if (fire_inner_solvent_mask(i) == 1) then
                    masks%inner_solvent_mask(counter) = i
                    counter = counter + 1
                end if
            end do
        counter=1
        do i = 1, natom_in
            if (fire_outer_mask(i) == 1) then
                masks%outer_mask(counter) = i
                counter = counter + 1
            end if
        end do
#ifdef MPI
        if (numtasks > 1) then
            call mpi_bcast(fire_inner_mask, natom_in, MPI_INTEGER, 0, commsander, ierr)
            call mpi_bcast(fire_inner_solvent_mask, natom_in, MPI_INTEGER, 0, commsander, ierr)
            call mpi_bcast(fire_outer_mask, natom_in, MPI_INTEGER, 0, commsander, ierr)

            mask_sizes(1) = size(masks%inner_mask)
            mask_sizes(2) = size(masks%inner_solvent_mask)
            mask_sizes(3) = size(masks%outer_mask)
            call mpi_bcast(mask_sizes, 3, MPI_INTEGER, 0, commsander, ierr)

            if (allocated(masks%inner_mask)) then
                if (size(masks%inner_mask) /= mask_sizes(1)) then
                    deallocate(masks%inner_mask)
                    allocate(masks%inner_mask(mask_sizes(1)))
                end if
            else
                allocate(masks%inner_mask(mask_sizes(1)))
            end if

            if (allocated(masks%inner_solvent_mask)) then
                if (size(masks%inner_solvent_mask) /= mask_sizes(2)) then
                    deallocate(masks%inner_solvent_mask)
                    allocate(masks%inner_solvent_mask(mask_sizes(2)))
                end if
            else
                allocate(masks%inner_solvent_mask(mask_sizes(2)))
            end if

            if (allocated(masks%outer_mask)) then
                if (size(masks%outer_mask) /= mask_sizes(3)) then
                    deallocate(masks%outer_mask)
                    allocate(masks%outer_mask(mask_sizes(3)))
                end if
            else
                allocate(masks%outer_mask(mask_sizes(3)))
            end if

            if (mask_sizes(1) > 0) then
                call mpi_bcast(masks%inner_mask, mask_sizes(1), MPI_INTEGER, 0, commsander, ierr)
            end if
            if (mask_sizes(2) > 0) then
                call mpi_bcast(masks%inner_solvent_mask, mask_sizes(2), MPI_INTEGER, 0, commsander, ierr)
            end if
            if (mask_sizes(3) > 0) then
                call mpi_bcast(masks%outer_mask, mask_sizes(3), MPI_INTEGER, 0, commsander, ierr)
            end if
        end if
#endif
        end if

    call ensure_prev_mask(fire_inner_mask, fire_inner_mask_prev)
    call ensure_prev_mask(fire_inner_solvent_mask, fire_inner_solvent_mask_prev)

    !JOSE Comment: mark as initialized to prevent double allocation on later calls
    fires_initialized = .true.
    fires_need_first_refresh = .true.
    end subroutine setup_fires

    ! Re-evaluate FIRES masks on demand using current coordinates
    subroutine refresh_fires_masks(xcur)
        use memory_module, only: atom_name, amber_atom_type, residue_label, residue_pointer, mass
        implicit none
#include "../include/memory.h"
#include "../include/md.h"
        _REAL_, intent(in) :: xcur(*)
        _REAL_ :: center(3)
        logical :: have_center
        integer :: i, counter
        character(len=256) :: s_inner, s_inw, s_outer
        integer :: prnlev
#ifdef MPI
        integer :: ierr
        integer :: mask_sizes(3)
#endif

        if (fires_static_mask) return

        if (fires_natom <= 0 .or. fires_nres <= 0) return

        have_center = .false.
        center = 0.0d0

        if (allocated(masks%inner_mask)) then
            if (size(masks%inner_mask) > 0) then
                call get_center(xcur, mass, center)
                have_center = .true.
            end if
        end if

        if (.not. allocated(fire_inner_mask))          allocate(fire_inner_mask(fires_natom))
        if (.not. allocated(fire_inner_solvent_mask))  allocate(fire_inner_solvent_mask(fires_natom))
        if (.not. allocated(fire_outer_mask))          allocate(fire_outer_mask(fires_natom))
        s_inner = trim(fire_inner)
        s_inw   = trim(fire_inner_solvent)
        s_outer = trim(fire_outer)
        prnlev = 0
        call atommask(fires_natom, fires_nres, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, xcur, s_inner,          fire_inner_mask)
        call atommask(fires_natom, fires_nres, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, xcur, s_inw,            fire_inner_solvent_mask)
        call atommask(fires_natom, fires_nres, prnlev, atom_name, amber_atom_type, residue_pointer, residue_label, xcur, s_outer,          fire_outer_mask)

        if (fires_hyst_enabled) then
            call apply_mask_hysteresis(fire_inner_mask, fire_inner_mask_prev, xcur, center, have_center)
            call apply_mask_hysteresis(fire_inner_solvent_mask, fire_inner_solvent_mask_prev, xcur, center, have_center)
        else
            call ensure_prev_mask(fire_inner_mask, fire_inner_mask_prev)
            call ensure_prev_mask(fire_inner_solvent_mask, fire_inner_solvent_mask_prev)
        end if

        if (allocated(masks%inner_mask))         deallocate(masks%inner_mask)
        if (allocated(masks%inner_solvent_mask)) deallocate(masks%inner_solvent_mask)
        if (allocated(masks%outer_mask))         deallocate(masks%outer_mask)
        allocate(masks%inner_mask( sum(fire_inner_mask(1:fires_natom))))
        allocate(masks%inner_solvent_mask( sum(fire_inner_solvent_mask(1:fires_natom))))
        allocate(masks%outer_mask( sum(fire_outer_mask(1:fires_natom))))

        counter=1
        do i = 1, fires_natom
            if (fire_inner_mask(i) == 1) then
                masks%inner_mask(counter) = i
                counter = counter + 1
            end if
        end do
        counter=1
        do i = 1, fires_natom
            if (fire_inner_solvent_mask(i) == 1) then
                masks%inner_solvent_mask(counter) = i
                counter = counter + 1
            end if
        end do
        counter=1
        do i = 1, fires_natom
            if (fire_outer_mask(i) == 1) then
                masks%outer_mask(counter) = i
                counter = counter + 1
            end if
        end do
#ifdef MPI
        if (numtasks > 1) then
            call mpi_bcast(fire_inner_mask, fires_natom, MPI_INTEGER, 0, commsander, ierr)
            call mpi_bcast(fire_inner_solvent_mask, fires_natom, MPI_INTEGER, 0, commsander, ierr)
            call mpi_bcast(fire_outer_mask, fires_natom, MPI_INTEGER, 0, commsander, ierr)

            mask_sizes(1) = size(masks%inner_mask)
            mask_sizes(2) = size(masks%inner_solvent_mask)
            mask_sizes(3) = size(masks%outer_mask)
            call mpi_bcast(mask_sizes, 3, MPI_INTEGER, 0, commsander, ierr)

            if (allocated(masks%inner_mask)) then
                if (size(masks%inner_mask) /= mask_sizes(1)) then
                    deallocate(masks%inner_mask)
                    allocate(masks%inner_mask(mask_sizes(1)))
                end if
            else
                allocate(masks%inner_mask(mask_sizes(1)))
            end if

            if (allocated(masks%inner_solvent_mask)) then
                if (size(masks%inner_solvent_mask) /= mask_sizes(2)) then
                    deallocate(masks%inner_solvent_mask)
                    allocate(masks%inner_solvent_mask(mask_sizes(2)))
                end if
            else
                allocate(masks%inner_solvent_mask(mask_sizes(2)))
            end if

            if (allocated(masks%outer_mask)) then
                if (size(masks%outer_mask) /= mask_sizes(3)) then
                    deallocate(masks%outer_mask)
                    allocate(masks%outer_mask(mask_sizes(3)))
                end if
            else
                allocate(masks%outer_mask(mask_sizes(3)))
            end if

            if (mask_sizes(1) > 0) then
                call mpi_bcast(masks%inner_mask, mask_sizes(1), MPI_INTEGER, 0, commsander, ierr)
            end if
            if (mask_sizes(2) > 0) then
                call mpi_bcast(masks%inner_solvent_mask, mask_sizes(2), MPI_INTEGER, 0, commsander, ierr)
            end if
            if (mask_sizes(3) > 0) then
                call mpi_bcast(masks%outer_mask, mask_sizes(3), MPI_INTEGER, 0, commsander, ierr)
            end if
        end if
#endif
        write(6,'(a)') 'FIRES: Refreshed masks using current coordinates.'
        write(6,'(a,i0)') 'FIRES: Inner region size: ', size(masks%inner_mask)
        write(6,'(a,i0)') 'FIRES: Inner solvent size: ', size(masks%inner_solvent_mask)
        write(6,'(a,i0)') 'FIRES: Outer region size: ', size(masks%outer_mask)
    end subroutine refresh_fires_masks

    subroutine calculate_fires_force_slow(x,natom,f,efires)
        implicit none
        _REAL_, intent(in) :: x(*)
        integer, intent(in) :: natom
        _REAL_, intent(inout) :: f(*)
        _REAL_ :: efires
    integer  i
    i = 1
        do while(i <= size(masks%inner_mask))
        !    print *, 'FIRES Inner Region:', (masks%inner_mask(i))
        !    print *, 'FIRES Inner Region Mass:',
            i = i + 1
        end do
    end subroutine calculate_fires_force_slow

    subroutine center_of_mass(x, mass, com)
        implicit none
        _REAL_, intent(in) :: x(*), mass(*)
        _REAL_, intent(out) :: com(3)
        integer :: i, a, atom_index, ref_idx, aref
        _REAL_ :: msum
        _REAL_ :: rref(3), rref_local(3), dr(3), rimag(3)
#include "../include/md.h"

        com = 0.0D0
        msum = 0.0D0

        if (size(masks%inner_mask) <= 0) return

        ref_idx = masks%inner_mask(1)
        if (ref_idx <= 0) return
        aref = 3*ref_idx - 2
        rref = x(aref:aref+2)

        i = 1
        do while (i <= size(masks%inner_mask))
            atom_index = masks%inner_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                dr = x(a:a+2) - rref
                call min_image_vec(dr)
                rimag = rref + dr
                com = com + mass(atom_index) * rimag
                msum = msum + mass(atom_index)
            else
                exit
            end if
            i = i + 1
        end do
        if (msum > 0.0D0) com = com / msum
    end subroutine center_of_mass

    ! Return the FIRES center: either a designated core atom or the COM of inner mask
    subroutine get_center(x, mass, center)
        implicit none
        _REAL_, intent(in)  :: x(*), mass(*)
        _REAL_, intent(out) :: center(3)
        integer :: a, i, atom_index, ref_idx, aref
        integer :: have_ref_local, have_ref_global
        _REAL_ :: msum_local, msum_global
        _REAL_ :: rref(3), rref_local(3), dr(3), rimag(3)
        _REAL_ :: center_local(3), center_global(3)
#ifdef MPI
        integer :: ierr
#endif
#include "../include/md.h"

        center = 0.0d0
        rref = 0.0d0
        dr = 0.0d0
        rimag = 0.0d0
        center_local = 0.0d0
        center_global = 0.0d0
        msum_local = 0.0d0
        msum_global = 0.0d0
        have_ref_local = 0
        have_ref_global = 0

        if (fires_center_atom > 0) then
            ref_idx = fires_center_atom
        else if (size(masks%inner_mask) > 0) then
            ref_idx = masks%inner_mask(1)
        else
            return
        end if

        if (ref_idx <= 0) return
        aref = 3*ref_idx - 2
        if (aref >= g_loc_istart3 .and. aref+2 <= g_loc_iend3) then
            rref = x(aref:aref+2)
            have_ref_local = 1
        end if
#ifdef MPI
        call mpi_allreduce(have_ref_local, have_ref_global, 1, MPI_INTEGER, mpi_sum, commsander, ierr)
        rref_local = rref
        call mpi_allreduce(rref_local, rref, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
#else
        have_ref_global = have_ref_local
#endif

        if (have_ref_global <= 0) return

        if (fires_center_atom > 0) then
            center = rref
            return
        end if

        if (size(masks%inner_mask) <= 0) return

        i = 1
        do while (i <= size(masks%inner_mask))
            atom_index = masks%inner_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = x(a:a+2) - rref
                    call min_image_vec(dr)
                    rimag = rref + dr
                    center_local = center_local + mass(atom_index) * rimag
                    msum_local = msum_local + mass(atom_index)
                end if
            else
                exit
            end if
            i = i + 1
        end do

#ifdef MPI
        call mpi_allreduce(center_local, center_global, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
        call mpi_allreduce(msum_local, msum_global, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
#else
        center_global = center_local
        msum_global = msum_local
#endif

        if (msum_global > 0.0D0) then
            center = center_global / msum_global
        else
            ! Fallback: if there is exactly one inner atom, use its coordinate even if massless
            if (size(masks%inner_mask) == 1) then
                atom_index = masks%inner_mask(1)
                if (atom_index > 0) then
                    center = rref
                end if
            end if
            if (.not. fires_notice_printed) then
                write(6,'(a)') 'FIRES WARNING: inner mask total mass is zero; using fallback center.'
                fires_notice_printed = .true.
            end if
        end if
    end subroutine get_center

    ! Apply minimum-image convention to a 3-vector using the current box (orthorhombic)
    subroutine min_image_vec(dr)
        implicit none
        _REAL_, intent(inout) :: dr(3)
#include "../include/md.h"
#include "box.h"
        integer :: k
        if (ifbox /= 0) then
            do k = 1, 3
                if (box(k) > 0.0d0) dr(k) = dr(k) - anint(dr(k)/box(k))*box(k)
            end do
        end if
    end subroutine min_image_vec

    logical function fires_in_band(d, was_in)
        implicit none
        _REAL_, intent(in) :: d
        logical, intent(in) :: was_in
        _REAL_ :: rin_local, rout_local

        rin_local  = max(0.0d0, fires_hyst_rin)
        rout_local = max(rin_local, fires_hyst_rout)

        if (.not. fires_hyst_enabled) then
            fires_in_band = (d <= rout_local)
            return
        end if

        if (was_in) then
            fires_in_band = (d <= rout_local)
        else
            fires_in_band = (d <= rin_local)
        end if
    end function fires_in_band

    subroutine ensure_prev_mask(curr, prev)
        implicit none
        integer, intent(in) :: curr(:)
        integer, allocatable, intent(inout) :: prev(:)

        if (allocated(prev)) then
            if (size(prev) /= size(curr)) then
                deallocate(prev)
                allocate(prev(size(curr)))
                prev = curr
            end if
        else
            allocate(prev(size(curr)))
            prev = curr
        end if
    end subroutine ensure_prev_mask

    subroutine apply_mask_hysteresis(mask_curr, mask_prev, xcur, center, have_center)
        implicit none
        integer, intent(inout) :: mask_curr(:)
        integer, allocatable, intent(inout) :: mask_prev(:)
        _REAL_, intent(in) :: xcur(*), center(3)
        logical, intent(in) :: have_center
        integer :: i, a
        _REAL_ :: dr(3), dist
        logical :: was_in
        integer, allocatable :: prev_snapshot(:), curr_local(:)
#ifdef MPI
        integer :: ierr
#endif
#include "../include/md.h"

        call ensure_prev_mask(mask_curr, mask_prev)

        if (.not. have_center) then
            mask_curr = mask_prev
#ifdef MPI
            call mpi_allreduce(MPI_IN_PLACE, mask_curr, size(mask_curr), MPI_INTEGER, mpi_max, commsander, ierr)
#endif
            mask_prev = mask_curr
            return
        end if

        allocate(prev_snapshot(size(mask_prev)))
        prev_snapshot = mask_prev

        allocate(curr_local(size(mask_curr)))
        curr_local = mask_curr

        do i = 1, size(mask_curr)
            if (curr_local(i) == 0) then
                was_in = (prev_snapshot(i) == 1)
                a = 3*i - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = xcur(a:a+2) - center
                    call min_image_vec(dr)
                    dist = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                    if (fires_in_band(dist, was_in)) curr_local(i) = 1
                else
                    if (was_in) curr_local(i) = 1
                end if
            end if
        end do

#ifdef MPI
        call mpi_allreduce(curr_local, mask_curr, size(mask_curr), MPI_INTEGER, mpi_max, commsander, ierr)
#else
        mask_curr = curr_local
#endif

        mask_prev = mask_curr

        deallocate(prev_snapshot)
        deallocate(curr_local)
    end subroutine apply_mask_hysteresis

    subroutine calculate_fires_force(x, mass, natom, f, efires)
        implicit none
        integer, intent(in) :: natom
        _REAL_, intent(in) :: x(3*natom), mass(natom)
        _REAL_, intent(inout) :: f(3*natom)
        _REAL_ :: efires

! Bring in MD counters (nstep, etc.)
#include "../include/md.h"

        ! locals
#ifdef MPI
        integer :: ierr
#endif
    _REAL_ :: com(3)
    _REAL_ :: r, rmax, delta, ui(3), dr(3)
        _REAL_ :: r_inw_max
        integer :: inw_max_idx
        _REAL_ :: eacc
    _REAL_ :: f_react_com(3)
    _REAL_ :: m_tot_inner
        _REAL_ :: f_sum_pen(3)          ! sum of forces on penetrating (outer) water
        _REAL_ :: f_sum_inw(3)          ! sum of forces on inner water (inner solvent)
        _REAL_ :: fa(3)                 ! per-atom applied force
    integer :: i, atom_index, a, j_anchor, idx
        _REAL_ :: min_outer_r
        integer :: min_outer_idx
    logical :: did_refresh
        _REAL_ :: rmax_minus, rmax_plus, rmaxm2, rmaxp2, r2
        integer :: n_in, n_inw, n_out
    ! (debug variables removed)
    ! Prepass extremes for single-atom inner fallback boundary
    _REAL_ :: pre_inw_max_r, pre_out_min_r
    integer :: pre_inw_max_idx, pre_out_min_idx
    ! Effective inner list skipping zero-mass (e.g. TIP4P extra points)
    integer, allocatable :: inner_eff(:)
    integer :: n_inner_eff
    _REAL_ :: pre_inw_max_r_glob, pre_out_min_r_glob
    _REAL_ :: f_react_com_glob(3), m_tot_inner_glob
    _REAL_ :: rmax_local
    integer :: j_anchor_local, anchor_candidate
    _REAL_ :: anchor_tol

        ! Initialize accumulators (do not zero caller's force buffer here;
        ! callers are responsible for clearing when needed)
        eacc = 0.0D0
        f_react_com = 0.0D0
    rmax = -1.0D0
        r_inw_max = 0.0D0
        inw_max_idx = 0
        did_refresh = .false.
        j_anchor = -1
    ! Ensure module global bounds are sane for this call (fallback)
    if (g_loc_istart3 < 1) g_loc_istart3 = 1
    if (g_loc_iend3 < 1) g_loc_iend3 = 3*natom
    f_sum_pen = 0.0D0
    f_sum_inw = 0.0D0
    last_n_pen = 0
    last_n_inw_out = 0
    min_outer_r = 1.0d99
    min_outer_idx = 0
    m_tot_inner = 0.0d0
    fa=0.d0

    ! Debug printing removed to reduce runtime verbosity.

    ! If any selection looks empty, recompute masks with current coordinates (only if enabled)
    associate(nstep => fires_current_step)
        if (fires_need_first_refresh) then
            call refresh_fires_masks(x)
            did_refresh = .true.
            fires_need_first_refresh = .false.

            if (nstep >= fires_freeze_on_nstep) then
                fires_static_mask     = .true.
                fires_refresh_runtime = .false.
                fires_hyst_enabled    = .false.
            end if

        else
            if (.not. fires_static_mask .and. nstep >= fires_freeze_on_nstep) then
                fires_static_mask     = .true.
                fires_refresh_runtime = .false.
                fires_hyst_enabled    = .false.
            end if

            if (fires_refresh_runtime) then
                if (size(masks%outer_mask) == 0 .or. size(masks%inner_mask) == 0) then
                    call refresh_fires_masks(x)
                    did_refresh = .true.
                end if
            end if
        end if
    end associate

    call get_center(x, mass, com)

    ! Build effective inner list (skip zero/near-zero mass atoms)
    if (allocated(inner_eff)) deallocate(inner_eff)
    if (size(masks%inner_mask) > 0) then
        allocate(inner_eff(size(masks%inner_mask)))
        n_inner_eff = 0
        do i = 1, size(masks%inner_mask)
            atom_index = masks%inner_mask(i)
            if (atom_index > 0) then
                if (mass(atom_index) > MASS_EPS) then
                    n_inner_eff = n_inner_eff + 1
                    inner_eff(n_inner_eff) = atom_index
                end if
            end if
        end do
        if (n_inner_eff == 0) then
            if (.not. fires_notice_printed) then
                write(6,'(a)') 'FIRES NOTICE: Inner mask has no massive atoms (likely virtual sites). Will derive rmax from inner-solvent cloud.'
                write(6,'(a)') 'FIRES NOTICE: Consider selecting a heavy atom for fire_inner or set fires_center_atom.'
                fires_notice_printed = .true.
            end if
            ! Continue; reaction redistribution to inner is skipped when m_tot_inner==0.
        end if
    else
        n_inner_eff = 0
    end if

    ! Short early-call debug prints removed.

        ! Prepass: farthest inner-solvent and nearest outer distances from COM
        pre_inw_max_r = 0.0d0
        pre_inw_max_idx = 0
        i = 1
        do while (i <= size(masks%inner_solvent_mask))
            atom_index = masks%inner_solvent_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                    if (r > pre_inw_max_r) then
                        pre_inw_max_r = r
                        pre_inw_max_idx = atom_index
                    end if
                end if
            else
                exit
            end if
            i = i + 1
        end do
        pre_out_min_r = 1.0d99
        pre_out_min_idx = 0
        i = 1
        do while (i <= size(masks%outer_mask))
            atom_index = masks%outer_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                    if (r < pre_out_min_r) then
                        pre_out_min_r = r
                        pre_out_min_idx = atom_index
                    end if
                end if
            else
                exit
            end if
            i = i + 1
        end do
#ifdef MPI
        pre_inw_max_r_glob = pre_inw_max_r
        pre_out_min_r_glob = pre_out_min_r
        call mpi_allreduce(pre_inw_max_r_glob, pre_inw_max_r, 1, MPI_DOUBLE_PRECISION, mpi_max, commsander, ierr)
        call mpi_allreduce(pre_out_min_r_glob, pre_out_min_r, 1, MPI_DOUBLE_PRECISION, mpi_min, commsander, ierr)
#endif

        ! First pass (JOSE): define boundary
        ! Mode 0: rmax from farthest INNER (solute) atom relative to COM
        ! Mode 1: rmax is constant (fires_rmax_const) relative to COM of INNER
        if (fires_rmax_mode == 0) then
            if (n_inner_eff > 0) then
                do i = 1, n_inner_eff
                    atom_index = inner_eff(i)
                    a = 3*atom_index - 2
                    if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                        dr = x(a:a+2) - com
                        call min_image_vec(dr)
                        r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                        if (r > rmax) then
                            rmax = r
                            j_anchor = atom_index
                        end if
                        m_tot_inner = m_tot_inner + mass(atom_index)
                    end if
                end do
            else
                ! No massive inner atoms: use farthest inner-solvent distance as boundary
                rmax = pre_inw_max_r
                if (size(masks%inner_mask) > 0) then
                    j_anchor = masks%inner_mask(1)
                else
                    j_anchor = 0
                end if
            end if
            ! If inner is a single atom (degenerate extent), infer a reasonable boundary
            if (rmax < 1.0d-6 .and. size(masks%inner_mask) <= 1) then
                if (pre_inw_max_r > 0.0d0 .and. pre_out_min_r < 1.0d90) then
                    rmax = 0.5d0*(pre_inw_max_r + pre_out_min_r)
                else if (pre_inw_max_r > 0.0d0) then
                    rmax = pre_inw_max_r
                end if
            end if
            rmax = rmax + max(0.0d0, fires_rmax_pad)
            if (fires_rmax_min > 0.0d0) rmax = max(rmax, fires_rmax_min)
            if (rmax < 1.0d-6) rmax = 1.0d-6
        else
            ! Constant-radius mode: choose an anchor (first INNER atom) and set rmax
            if (size(masks%inner_mask) > 0) then
                j_anchor = masks%inner_mask(1)
                ! accumulate total inner mass for reaction distribution
                do i = 1, n_inner_eff
                    atom_index = inner_eff(i)
                    a = 3*atom_index - 2
                    if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                        m_tot_inner = m_tot_inner + mass(atom_index)
                    end if
                end do
            end if
            rmax = fires_rmax_const
            rmax = max(rmax, fires_rmax_min)
            rmax = rmax + max(0.0d0, fires_rmax_pad)
        end if

        ! Optional smoothing for auto mode; reset when using constant mode
        if (fires_rmax_mode == 0) then
            if (rmax_ema < 0.0d0 .or. did_refresh) rmax_ema = rmax
            rmax_ema = fires_rmax_ema_alpha*rmax + (1.0d0 - fires_rmax_ema_alpha)*rmax_ema
            rmax = rmax_ema
        else
            rmax_ema = -1.0d0
        end if

    rmax_local = rmax
    j_anchor_local = j_anchor
    anchor_candidate = 0
    anchor_tol = 1.0d-12*max(1.0d0, abs(rmax_local))
#ifdef MPI
    call mpi_allreduce(rmax_local, rmax, 1, MPI_DOUBLE_PRECISION, mpi_max, commsander, ierr)
    if (rmax > 0.0d0 .and. j_anchor_local > 0) then
        if (abs(rmax - rmax_local) <= max(1.0d-12, anchor_tol)) anchor_candidate = j_anchor_local
    end if
    call mpi_allreduce(anchor_candidate, j_anchor, 1, MPI_INTEGER, mpi_max, commsander, ierr)
#else
    j_anchor = j_anchor_local
#endif

    ! Reset notice throttle if FIRES is active this step
    if (j_anchor > 0 .and. rmax >= EPS) fires_notice_printed = .false.

    if (j_anchor <= 0 .or. rmax < EPS) then
            if (fires_refresh_runtime .and. .not. did_refresh) then
                call refresh_fires_masks(x)
                did_refresh = .true.
                ! Recompute COM and boundary once after refresh
                call get_center(x, mass, com)
                rmax = -1.0D0
                j_anchor = -1
                if (fires_rmax_mode == 0) then
                    i = 1
                    do while (i <= size(masks%inner_mask))
                        atom_index = masks%inner_mask(i)
                        if (atom_index .gt. 0) then
                            a = 3*atom_index - 2
                            if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                                dr = x(a:a+2) - com
                                call min_image_vec(dr)
                                r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                                if (r > rmax) then
                                    rmax = r
                                    j_anchor = atom_index
                                end if
                                m_tot_inner = m_tot_inner + mass(atom_index)
                            end if
                        else
                            exit
                        end if
                        i = i + 1
                    end do
                    rmax = rmax + max(0.0d0, fires_rmax_pad)
                    if (fires_rmax_min > 0.0d0) rmax = max(rmax, fires_rmax_min)
                    if (rmax < 1.0d-6) rmax = 1.0d-6
                else
                    if (size(masks%inner_mask) > 0) then
                        j_anchor = masks%inner_mask(1)
                        ! accumulate total inner mass for reaction distribution
                        i = 1
                        do while (i <= size(masks%inner_mask))
                            atom_index = masks%inner_mask(i)
                            if (atom_index .gt. 0) then
                                a = 3*atom_index - 2
                                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                                    m_tot_inner = m_tot_inner + mass(atom_index)
                                end if
                            else
                                exit
                            end if
                            i = i + 1
                        end do
                    end if
                    rmax = fires_rmax_const
                    rmax = max(rmax, fires_rmax_min)
                    rmax = rmax + max(0.0d0, fires_rmax_pad)
                end if
                if (fires_rmax_mode == 0) then
                    rmax_ema = rmax
                else
                    rmax_ema = -1.0d0
                end if
            end if
        end if

    rmax_local = rmax
    j_anchor_local = j_anchor
    anchor_candidate = 0
    anchor_tol = 1.0d-12*max(1.0d0, abs(rmax_local))
#ifdef MPI
    call mpi_allreduce(rmax_local, rmax, 1, MPI_DOUBLE_PRECISION, mpi_max, commsander, ierr)
    if (rmax > 0.0d0 .and. j_anchor_local > 0) then
        if (abs(rmax - rmax_local) <= max(1.0d-12, anchor_tol)) anchor_candidate = j_anchor_local
    end if
    call mpi_allreduce(anchor_candidate, j_anchor, 1, MPI_INTEGER, mpi_max, commsander, ierr)
#else
    j_anchor = j_anchor_local
#endif

        if (j_anchor <= 0 .or. rmax < EPS) then
            ! Even if FIRES is inactive, compute farthest inner-solvent distance for diagnostics
            r_inw_max = 0.0d0
            i = 1
            do while (i <= size(masks%inner_solvent_mask))
                atom_index = masks%inner_solvent_mask(i)
                if (atom_index .gt. 0) then
                    a = 3*atom_index - 2
                    if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                        dr = x(a:a+2) - com
                        call min_image_vec(dr)
                        r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                        if (r > r_inw_max) then
                            r_inw_max = r
                            inw_max_idx = atom_index
                        end if
                    end if
                else
                    exit
                end if
                i = i + 1
            end do
            efires = 0.0D0
            last_f_pot(1:3) = 0.0D0
            last_f_inner_water(1:3) = 0.0D0
            last_f_pen_water(1:3) = 0.0D0
            last_rmax = 0.0d0
            last_anchor_idx = 0
            last_rmax_inw = r_inw_max
            last_inw_max_idx = inw_max_idx
            last_pen_atom_idx = 0
            last_pen_r = 0.0d0
            last_pen_depth = 0.0d0
            ! Print a one-time hint if FIRES inactive due to masks
            if (.not. fires_notice_printed) then
                n_in  = size(masks%inner_mask); n_inw = size(masks%inner_solvent_mask); n_out = size(masks%outer_mask)
                write(6,'(a)') 'FIRES NOTICE: zero boundary radius or no anchor; FIRES inactive this step.'
                write(6,'(a,i0,a,i0,a,i0)') 'FIRES NOTICE: selection sizes: inner=', n_in, &
                     ', inner_solvent=', n_inw, ', outer=', n_out
                write(6,'(a,1x,a)') 'FIRES NOTICE: fire_inner         =', trim(fire_inner)
                write(6,'(a,1x,a)') 'FIRES NOTICE: fire_inner_solvent =', trim(fire_inner_solvent)
                write(6,'(a,1x,a)') 'FIRES NOTICE: fire_outer         =', trim(fire_outer)
                fires_notice_printed = .true.
            end if
            return
        end if

        ! Second pass (JOSE):
        ! Repel OUTER atoms if they penetrate inside (r < rmax)
        ! Keep INNER atoms inside by pulling back if they try to leave (r > rmax)

        ! Precompute squared thresholds to reduce sqrt calls
        rmax_minus = max(0.0d0, rmax - EPS)
        rmax_plus  = rmax + EPS
        rmaxm2 = rmax_minus*rmax_minus
        rmaxp2 = rmax_plus*rmax_plus

        ! OUTER mask: push outward if inside boundary
        i = 1
        do while (i <= size(masks%outer_mask))
            atom_index = masks%outer_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
                    if (r2 < rmaxm2) then
                        r = sqrt(r2)
                        if (fires_band_filter_enabled) then
                            ! For OUTER (inside), keep only those within halfwidth of boundary
                            if ((rmax - r) > fires_band_halfwidth) then
                                i = i + 1
                                cycle
                            end if
                        end if
                        if (r < min_outer_r) then
                            min_outer_r = r
                            min_outer_idx = atom_index
                        end if
                        ui = dr / max(r, EPS)
                        delta = r - rmax  ! negative
                        ! Applied force on penetrating outer atom
                        fa = -fires_k * delta * ui
                        f(a:a+2) = f(a:a+2) + fa
                        eacc = eacc + 0.5D0 * fires_k * (rmax - r)*(rmax - r)
                        ! Reaction on inner COM is negative of applied forces
                        f_react_com = f_react_com - fa
                        f_sum_pen = f_sum_pen + fa
                        last_n_pen = last_n_pen + 1
                    end if
                end if
            else
                exit
            end if
            i = i + 1
        end do

        ! INNER mask: pull inward if outside boundary
        i = 1
        do while (i <= size(masks%inner_mask))
            atom_index = masks%inner_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
                    if (r2 > rmaxp2) then
                        r = sqrt(r2)
                        if (fires_band_filter_enabled) then
                            ! For INNER (outside), keep only those within halfwidth of boundary
                            if ((r - rmax) > fires_band_halfwidth) then
                                i = i + 1
                                cycle
                            end if
                        end if
                        ui = dr / max(r, EPS)
                        delta = r - rmax  ! positive
                        ! Applied force on inner (non-water) atom: pull inward
                        fa = -fires_k * delta * ui
                        f(a:a+2) = f(a:a+2) + fa
                        eacc = eacc + 0.5D0 * fires_k * (r - rmax)*(r - rmax)
                        ! Reaction applied to inner COM to conserve momentum
                        f_react_com = f_react_com - fa
                    end if
                end if
            else
                exit
            end if
            i = i + 1
        end do

        ! INNER solvent mask: treat same as INNER (kept inside)
    i = 1
        do while (i <= size(masks%inner_solvent_mask))
            atom_index = masks%inner_solvent_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
                    r = sqrt(r2)
                    if (r > r_inw_max) then
                        r_inw_max = r
                        inw_max_idx = atom_index
                    end if
                    if (r2 > rmaxp2) then
                        if (fires_band_filter_enabled) then
                            if ((r - rmax) > fires_band_halfwidth) then
                                i = i + 1
                                cycle
                            end if
                        end if
                        ui = dr / max(r, EPS)
                        delta = r - rmax
                        fa = -fires_k * delta * ui
                        f(a:a+2) = f(a:a+2) + fa
                        eacc = eacc + 0.5D0 * fires_k * (r - rmax)*(r - rmax)
                        f_react_com = f_react_com - fa
                        f_sum_inw = f_sum_inw + fa
                        last_n_inw_out = last_n_inw_out + 1
                    end if
                end if
            else
                exit
            end if
            i = i + 1
        end do

        ! Apply the force to the INNER group mass-weighted to conserve momentum and avoid torque
#ifdef MPI
        m_tot_inner_glob = m_tot_inner
        f_react_com_glob = f_react_com
        call mpi_allreduce(m_tot_inner_glob, m_tot_inner, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
        call mpi_allreduce(f_react_com_glob, f_react_com, 3, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
#endif
        if (m_tot_inner > EPS) then
            do i = 1, n_inner_eff
                atom_index = inner_eff(i)
                a = 3*atom_index - 2
                if (a >= g_loc_istart3 .and. a+2 <= g_loc_iend3) then
                    f(a:a+2) = f(a:a+2) + (mass(atom_index)/m_tot_inner) * f_react_com
                end if
            end do
        end if

    efires = eacc

    ! Export for printing at the driver level
    last_f_pot(1:3)         = f_react_com(1:3)
    last_f_inner_water(1:3) = f_sum_inw(1:3)
    last_f_pen_water(1:3)   = f_sum_pen(1:3)
    last_rmax               = rmax
    last_anchor_idx         = j_anchor
    if (r_inw_max > 0.0d0) then
        last_rmax_inw    = r_inw_max
        last_inw_max_idx = inw_max_idx
    else
        last_rmax_inw    = pre_inw_max_r
        last_inw_max_idx = pre_inw_max_idx
    end if
    if (min_outer_idx > 0) then
        last_pen_atom_idx = min_outer_idx
        last_pen_r        = min_outer_r
        last_pen_depth    = max(0.0d0, rmax - min_outer_r)
    else
        last_pen_atom_idx = 0
        last_pen_r        = 0.0d0
        last_pen_depth    = 0.0d0
    end if

        if (allocated(inner_eff)) deallocate(inner_eff)
    end subroutine calculate_fires_force

            ! --- FIRES helper: set local 3*index bounds for this rank ---
    subroutine fires_set_local_bounds(istart3_in, iend3_in)
        implicit none
        integer, intent(in) :: istart3_in, iend3_in
        g_loc_istart3 = istart3_in
        g_loc_iend3   = iend3_in
end subroutine fires_set_local_bounds

    subroutine fires_set_step(step)
        implicit none
        integer, intent(in) :: step
        fires_current_step = step
    end subroutine fires_set_step

    subroutine fires_force(x,mass,natom, f, epot)
        implicit none
        integer, intent(in) :: natom
        _REAL_, intent(in)  :: x(3*natom),mass(natom)
        _REAL_, intent(inout) :: f(3*natom)
        _REAL_,intent(out) ::  epot
        _REAL_ :: efires
        
        epot = 0.d0
        ! Early gate if FIRES disabled or zero strength
        if (fires == 0 .or. fires_k <= 0.0d0) then
            epot = 0.0d0
            return
        end if
        if (.not. allocated(masks%inner_mask)) then
            epot = 0.0d0
            return
        end if
        if (size(masks%inner_mask) <= 0) then
            epot = 0.0d0
            return
        end if
        call calculate_fires_force(x,mass,natom, f, efires)
        epot = efires
    end subroutine fires_force

    subroutine fires_mask_signature(sig_inner, sig_inner_solvent, sig_outer, count_inner, count_inner_solvent, count_outer)
        use iso_fortran_env, only: int64
        implicit none
        integer(int64), intent(out) :: sig_inner, sig_inner_solvent, sig_outer
        integer, intent(out) :: count_inner, count_inner_solvent, count_outer
        integer(int64), parameter :: fnv_offset = int(Z'CBF29CE484222325', int64)
        integer(int64), parameter :: fnv_prime  = int(Z'100000001B3', int64)
        integer(int64) :: acc
        integer :: i

        sig_inner         = fnv_offset
        sig_inner_solvent = fnv_offset
        sig_outer         = fnv_offset

        if (allocated(masks%inner_mask)) then
            count_inner = size(masks%inner_mask)
            do i = 1, count_inner
                acc = int(masks%inner_mask(i), int64) + int(i, int64)
                sig_inner = ieor(sig_inner, acc)
                sig_inner = sig_inner * fnv_prime
            end do
        else
            count_inner = 0
        end if

        if (allocated(masks%inner_solvent_mask)) then
            count_inner_solvent = size(masks%inner_solvent_mask)
            do i = 1, count_inner_solvent
                acc = int(masks%inner_solvent_mask(i), int64) + int(i, int64)
                sig_inner_solvent = ieor(sig_inner_solvent, acc)
                sig_inner_solvent = sig_inner_solvent * fnv_prime
            end do
        else
            count_inner_solvent = 0
        end if

        if (allocated(masks%outer_mask)) then
            count_outer = size(masks%outer_mask)
            do i = 1, count_outer
                acc = int(masks%outer_mask(i), int64) + int(i, int64)
                sig_outer = ieor(sig_outer, acc)
                sig_outer = sig_outer * fnv_prime
            end do
        else
            count_outer = 0
        end if
    end subroutine fires_mask_signature

    subroutine fires_prepare_if_needed(x, mass, natom)
        implicit none
        integer, intent(in) :: natom
        _REAL_, intent(in)  :: x(3*natom), mass(natom)
        _REAL_, allocatable :: f_dummy(:)
        _REAL_ :: efires_dummy

        if (fires == 0) return
        if (.not. fires_need_first_refresh) return
        if (.not. allocated(masks%inner_mask)) return
        if (size(masks%inner_mask) <= 0) return

        allocate(f_dummy(3*natom))
        f_dummy = 0.0d0
        call calculate_fires_force(x, mass, natom, f_dummy, efires_dummy)
        deallocate(f_dummy)
    end subroutine fires_prepare_if_needed

    logical function fires_restraint_enabled()
        implicit none
        fires_restraint_enabled = .false.
        if (fires == 0) return
        if (fires_k <= EPS) return
        if (.not. allocated(masks%inner_mask)) return
        if (size(masks%inner_mask) <= 0) return
        fires_restraint_enabled = .true.
    end function fires_restraint_enabled

end module qmmm_fires_module
