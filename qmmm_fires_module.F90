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
    private
    public :: FireMasks, calculate_fires_force, fires_force, setup_fires, &
              fires_prepare_if_needed, fires_set_local_bounds, fires_set_step
    ! Debug/diagnostics: last computed FIRES force components
    _REAL_, public, save :: last_f_pot(3)           = 0.0d0
    _REAL_, public, save :: last_f_inner_water(3)   = 0.0d0
    _REAL_, public, save :: last_f_pen_water(3)     = 0.0d0
    integer, public, save :: fires = 0
    integer, public, save :: last_n_pen = 0              ! diagnostics: # of penetrating outer atoms this eval
    integer, public, save :: last_n_inw_out = 0          ! diagnostics: # of inner-solvent atoms outside rmax
    _REAL_,  public, save :: last_rmax = 0.0d0           ! diagnostics: boundary radius used this eval
    integer, public, save :: last_anchor_idx = 0         ! diagnostics: anchor atom index used this eval
    _REAL_,  public, save :: last_rmax_inw = 0.0d0       ! diagnostics: farthest inner-solvent atom distance from COM
    integer, public, save :: last_inw_max_idx = 0         ! diagnostics: atom index for farthest inner-solvent
    integer, public, save :: last_pen_atom_idx = 0        ! diagnostics: most inward penetrating outer atom index
    _REAL_,  public, save :: last_pen_r = 0.0d0           ! diagnostics: its distance to COM
    _REAL_,  public, save :: last_pen_depth = 0.0d0       ! diagnostics: penetration depth (rmax - r)
    !JOSE MOD: MTS controls
    !_REAL_ small-step size in ps; when > 0 and mts_n>1, subcycle FIRES with dt_fast=mts_fires
    _REAL_,  public, save :: mts_fires = 0.0d0
    ! number of FIRES microsteps per outer step
    integer, public, save :: mts_n = 1
    integer, dimension(:), allocatable :: fire_inner_mask, fire_outer_mask, fire_inner_solvent_mask
    ! Removed unused module-level natom/nres to avoid name collisions with COMMON from memory.h
    character(len=256), public, save :: fire_inner = ""
    character(len=256), public, save :: fire_inner_solvent = ""
    character(len=256), public, save :: fire_outer = ""
    ! (deleted unused ipres/igraph/isymbl/lbres/crd/nsc/inner_atoms/outer_atoms)
    ! Control whether to auto-refresh FIRES masks at runtime (each step). Default: off
    logical, public, save :: fires_refresh_runtime = .false.
    logical, public, save :: fires_restraint_enabled = .true.
    logical, public, save :: fires_static_mask = .false.
    logical, public, save :: fires_band_filter_enabled = .false.
    _REAL_,  public, save :: fires_band_halfwidth = 0.0d0
    logical, public, save :: fires_in_force = .false.
    integer, public, save :: fires_freeze_on_nstep = 0
    integer, parameter :: fires_hash_kind = selected_int_kind(18)
    integer(fires_hash_kind), public, save :: fires_mask_signature = 0_fires_hash_kind

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
    ! Center control: 0 = use COM of inner mask; >0 = use this atom index as center
    integer, public, save :: fires_center_atom = 0
    integer, save :: g_loc_istart3 = 1
    integer, save :: g_loc_iend3   = 0

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
        write(6,'(a)')    'FIRES: Using masks:'
        write(6,'(a,1x,a)') '  fire_inner         =', trim(s_inner)
        write(6,'(a,1x,a)') '  fire_inner_solvent =', trim(s_inw)
        write(6,'(a,1x,a)') '  fire_outer         =', trim(s_outer)
        write(6,'(a,i0)') 'FIRES: Inner region size: ', size(masks%inner_mask)
        write(6,'(a,i0)') 'FIRES: Inner solvent size: ', size(masks%inner_solvent_mask)
        write(6,'(a,i0)') 'FIRES: Outer region size: ', size(masks%outer_mask)
        write(6,'(a,1pe12.4)') 'FIRES: fires_k = ', fires_k
        write(6,'(a,1pe12.4,a,i0)') 'FIRES: mts_fires = ', mts_fires, ', mts_n = ', mts_n
    if (fires_rmax_mode == 0) then
            write(6,'(a,1x,1pe12.4)') 'FIRES: rmax mode = inner-extent (auto), pad =', fires_rmax_pad
        else
            write(6,'(a,1x,1pe12.4)') 'FIRES: rmax mode = constant, rmax =', fires_rmax_const
        end if
    write(6,'(a,1x,1pe12.4)') 'FIRES: rmax_min =', fires_rmax_min
        write(6,'(a,1x,1pe12.4)') 'FIRES: rmax_ema_alpha =', fires_rmax_ema_alpha
        if (size(masks%inner_mask)==0 .or. size(masks%outer_mask)==0) then
            write(6,*) 'FIRES WARNING: One or more masks are empty; FIRES forces will be zero.'
        end if
        if (fires_k <= 0.0d0) then
            write(6,*) 'FIRES WARNING: fires_k <= 0; FIRES forces will be zero.'
        end if

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
            write(6,*) 'FIRES DEBUG: Re-evaluating masks with prnlev=1 to print parse/eval details.'
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
            write(6,'(a)') 'FIRES DEBUG: Post-debug selection sizes:'
            write(6,'(a,i0)') '  Inner region size: ', size(masks%inner_mask)
            write(6,'(a,i0)') '  Inner solvent size: ', size(masks%inner_solvent_mask)
            write(6,'(a,i0)') '  Outer region size: ', size(masks%outer_mask)
        end if

    !JOSE Comment: mark as initialized to prevent double allocation on later calls
    call update_fires_mask_signature()

    fires_initialized = .true.
    fires_need_first_refresh = .true.
    end subroutine setup_fires

    subroutine fires_set_local_bounds(istart3_in, iend3_in)
        implicit none
        integer, intent(in) :: istart3_in, iend3_in

        g_loc_istart3 = max(1, istart3_in)
        g_loc_iend3   = iend3_in
    end subroutine fires_set_local_bounds

    subroutine fires_prepare_if_needed(x, mass, natom)
        implicit none
        integer, intent(in) :: natom
        _REAL_, intent(in)  :: x(3*natom)
        _REAL_, intent(in)  :: mass(natom)

        if (fires == 0) return
        if (fires_k <= 0.0d0) return

        if (fires_need_first_refresh) then
            call refresh_fires_masks(x)
            fires_need_first_refresh = .false.
        end if

        call update_fires_mask_signature()
    end subroutine fires_prepare_if_needed

    subroutine fires_set_step(nstep)
        use memory_module, only : x
        implicit none
        integer, intent(in) :: nstep
#include "../include/memory.h"

        if (fires == 0) return

        if (fires_freeze_on_nstep > 0 .and. nstep == fires_freeze_on_nstep) then
            if (.not. fires_static_mask) then
                call refresh_fires_masks(x(lcrd))
                fires_static_mask = .true.
            end if
        else if (.not. fires_static_mask .and. fires_refresh_runtime) then
            call refresh_fires_masks(x(lcrd))
        end if

        call update_fires_mask_signature()
    end subroutine fires_set_step

    ! Re-evaluate FIRES masks on demand using current coordinates
    subroutine update_fires_mask_signature()
        implicit none
        integer :: i
        integer(fires_hash_kind) :: sig
        sig = 0_fires_hash_kind

        if (allocated(masks%inner_mask)) then
            do i = 1, size(masks%inner_mask)
                sig = ieor(sig, int(masks%inner_mask(i), fires_hash_kind))
            end do
        end if
        if (allocated(masks%inner_solvent_mask)) then
            do i = 1, size(masks%inner_solvent_mask)
                sig = ieor(sig, ishft(int(masks%inner_solvent_mask(i), fires_hash_kind), 1))
            end do
        end if
        if (allocated(masks%outer_mask)) then
            do i = 1, size(masks%outer_mask)
                sig = ieor(sig, ishft(int(masks%outer_mask(i), fires_hash_kind), 2))
            end do
        end if
        fires_mask_signature = sig
    end subroutine update_fires_mask_signature

    logical function index_in_local_range(idx3)
        implicit none
        integer, intent(in) :: idx3

        if (g_loc_iend3 >= g_loc_istart3 .and. g_loc_iend3 > 0) then
            index_in_local_range = (idx3 >= g_loc_istart3 .and. idx3 <= g_loc_iend3)
        else
            index_in_local_range = .true.
        end if
    end function index_in_local_range

    logical function atom_is_local(atom_index)
        implicit none
        integer, intent(in) :: atom_index
        integer :: idx3

        idx3 = 3*atom_index - 2
        atom_is_local = index_in_local_range(idx3) .and. index_in_local_range(idx3+2)
    end function atom_is_local

    subroutine refresh_fires_masks(xcur)
        use memory_module, only: atom_name, amber_atom_type, residue_label, residue_pointer
        implicit none
#include "../include/memory.h"
#include "../include/md.h"
        integer :: i, counter
        _REAL_, intent(in) :: xcur(*)
        character(len=256) :: s_inner, s_inw, s_outer
        integer :: prnlev

        if (fires_natom <= 0 .or. fires_nres <= 0) return

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
        write(6,'(a)') 'FIRES: Refreshed masks using current coordinates.'
        write(6,'(a,i0)') 'FIRES: Inner region size: ', size(masks%inner_mask)
        write(6,'(a,i0)') 'FIRES: Inner solvent size: ', size(masks%inner_solvent_mask)
        write(6,'(a,i0)') 'FIRES: Outer region size: ', size(masks%outer_mask)
        call update_fires_mask_signature()
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
        _REAL_ :: rref(3), dr(3), rimag(3)
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
      use memory_module, only : natom
        implicit none
        _REAL_, intent(in)  :: x(*), mass(*)
        _REAL_, intent(out) :: center(3)
        integer :: a, i, atom_index, ref_idx, aref
        _REAL_ :: msum, rref(3), dr(3), rimag(3)
#include "../include/md.h"

        msum = 0.d0
        rref = 0.d0
        dr = 0.d0
        rimag = 0.d0
        a = 0
        i = 0
        atom_index = 0
        ref_idx = 0
        aref = 0
        center = 0.d0
        
        if (fires_center_atom > 0) then
            a = 3*fires_center_atom - 2
            center = x(a:a+2)
            return
        end if

        ! Fallback: COM of inner mask (same logic as center_of_mass)
        center = 0.0D0
        msum   = 0.0D0
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
                center = center + mass(atom_index) * rimag
                msum   = msum   + mass(atom_index)
            else
                exit
            end if
            i = i + 1
        end do
        if (msum > 0.0D0) center = center / msum
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

    subroutine calculate_fires_force(x, mass, natom, f, efires)
        implicit none
        integer, intent(in) :: natom
        _REAL_, intent(in) :: x(3*natom), mass(natom)
        _REAL_, intent(inout) :: f(3*natom)
        _REAL_ :: efires

        ! locals
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
    ! Prepass extremes for single-atom inner fallback boundary
    _REAL_ :: pre_inw_max_r, pre_out_min_r
    integer :: pre_inw_max_idx, pre_out_min_idx

        ! Initialize accumulators (do not zero caller's force buffer here;
    ! callers are responsible for clearing when needed)
        eacc = 0.0D0
        f_react_com = 0.0D0
    rmax = -1.0D0
        r_inw_max = 0.0D0
        inw_max_idx = 0
        did_refresh = .false.
        j_anchor = -1
    f_sum_pen = 0.0D0
    f_sum_inw = 0.0D0
    last_n_pen = 0
    last_n_inw_out = 0
    min_outer_r = 1.0d99
    min_outer_idx = 0
    m_tot_inner = 0.0d0
    fa=0.d0

    ! If any selection looks empty, recompute masks with current coordinates (only if enabled)
    if (fires_need_first_refresh) then
        call refresh_fires_masks(x)
        did_refresh = .true.
        fires_need_first_refresh = .false.
    else if (fires_refresh_runtime) then
        if (size(masks%outer_mask) == 0 .or. size(masks%inner_mask) == 0) then
            call refresh_fires_masks(x)
            did_refresh = .true.
        end if
    end if

    call get_center(x, mass, com)

        ! Prepass: farthest inner-solvent and nearest outer distances from COM
        pre_inw_max_r = 0.0d0
        pre_inw_max_idx = 0
        i = 1
        do while (i <= size(masks%inner_solvent_mask))
            atom_index = masks%inner_solvent_mask(i)
            if (atom_index .gt. 0) then
                a = 3*atom_index - 2
                dr = x(a:a+2) - com
                call min_image_vec(dr)
                r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                if (r > pre_inw_max_r) then
                    pre_inw_max_r = r
                    pre_inw_max_idx = atom_index
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
                dr = x(a:a+2) - com
                call min_image_vec(dr)
                r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                if (r < pre_out_min_r) then
                    pre_out_min_r = r
                    pre_out_min_idx = atom_index
                end if
            else
                exit
            end if
            i = i + 1
        end do

        ! First pass (JOSE): define boundary
        ! Mode 0: rmax from farthest INNER (solute) atom relative to COM
        ! Mode 1: rmax is constant (fires_rmax_const) relative to COM of INNER
        if (fires_rmax_mode == 0) then
            i = 1
            do while (i <= size(masks%inner_mask))
                atom_index = masks%inner_mask(i)
                if (atom_index .gt. 0) then
                    a = 3*atom_index - 2
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                    if (r > rmax) then
                        rmax = r
                        j_anchor = atom_index
                    end if
                    m_tot_inner = m_tot_inner + mass(atom_index)
                else
                    exit
                end if
                i = i + 1
            end do
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
                i = 1
                do while (i <= size(masks%inner_mask))
                    atom_index = masks%inner_mask(i)
                    if (atom_index .gt. 0) then
                        m_tot_inner = m_tot_inner + mass(atom_index)
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

        ! Optional smoothing for auto mode; reset when using constant mode
        if (fires_rmax_mode == 0) then
            if (rmax_ema < 0.0d0 .or. did_refresh) rmax_ema = rmax
            rmax_ema = fires_rmax_ema_alpha*rmax + (1.0d0 - fires_rmax_ema_alpha)*rmax_ema
            rmax = rmax_ema
        else
            rmax_ema = -1.0d0
        end if

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
                            dr = x(a:a+2) - com
                            call min_image_vec(dr)
                            r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                            if (r > rmax) then
                                rmax = r
                                j_anchor = atom_index
                            end if
                            m_tot_inner = m_tot_inner + mass(atom_index)
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
                                m_tot_inner = m_tot_inner + mass(atom_index)
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

        if (j_anchor <= 0 .or. rmax < EPS) then
            ! Even if FIRES is inactive, compute farthest inner-solvent distance for diagnostics
            r_inw_max = 0.0d0
            i = 1
            do while (i <= size(masks%inner_solvent_mask))
                atom_index = masks%inner_solvent_mask(i)
                if (atom_index .gt. 0) then
                    a = 3*atom_index - 2
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                    r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
                    if (r > r_inw_max) then
                        r_inw_max = r
                        inw_max_idx = atom_index
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
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
                if (r2 < rmaxm2) then
                    r = sqrt(r2)
                    if (r < min_outer_r) then
                        min_outer_r = r
                        min_outer_idx = atom_index
                    end if
                    ui = dr / max(r, EPS)
                    delta = r - rmax  ! negative
                    ! Applied force on penetrating outer atom
                    fa = -fires_k * delta * ui
                    if (atom_is_local(atom_index)) then
                        f(a:a+2) = f(a:a+2) + fa
                    end if
                    ! Reaction on inner COM is negative of applied forces
                    f_react_com = f_react_com - fa
                    f_sum_pen = f_sum_pen + fa
                    last_n_pen = last_n_pen + 1
                    eacc = eacc + 0.5D0 * fires_k * (rmax - r)*(rmax - r)
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
                    dr = x(a:a+2) - com
                    call min_image_vec(dr)
                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
                if (r2 > rmaxp2) then
                    r = sqrt(r2)
                    ui = dr / max(r, EPS)
                    delta = r - rmax  ! positive
                    ! Applied force on inner (non-water) atom: pull inward
                    fa = -fires_k * delta * ui
                    if (atom_is_local(atom_index)) then
                        f(a:a+2) = f(a:a+2) + fa
                    end if
                    ! Reaction applied to inner COM to conserve momentum
                    f_react_com = f_react_com - fa
                    eacc = eacc + 0.5D0 * fires_k * (r - rmax)*(r - rmax)
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
                dr = x(a:a+2) - com
                call min_image_vec(dr)
                r2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
                r = sqrt(r2)
        if (r > r_inw_max) then
            r_inw_max = r
            inw_max_idx = atom_index
        end if
                if (r2 > rmaxp2) then
                    ui = dr / max(r, EPS)
                    delta = r - rmax
                    fa = -fires_k * delta * ui
                    if (atom_is_local(atom_index)) then
                        f(a:a+2) = f(a:a+2) + fa
                    end if
                    f_react_com = f_react_com - fa
                    f_sum_inw = f_sum_inw + fa
                    last_n_inw_out = last_n_inw_out + 1
                    eacc = eacc + 0.5D0 * fires_k * (r - rmax)*(r - rmax)
                end if
            else
                exit
            end if
            i = i + 1
        end do

        ! Apply the force to the INNER group mass-weighted to conserve momentum and avoid torque
        if (m_tot_inner > EPS) then
            i = 1
            do while (i <= size(masks%inner_mask))
                atom_index = masks%inner_mask(i)
                if (atom_index .gt. 0) then
                    a = 3*atom_index - 2
                    if (atom_is_local(atom_index)) then
                        f(a:a+2) = f(a:a+2) + (mass(atom_index)/m_tot_inner) * f_react_com
                    end if
                else
                    exit
                end if
                i = i + 1
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

    end subroutine calculate_fires_force

    subroutine fires_force(x,mass,natom, f, epot)
        implicit none
        integer, intent(in) :: natom
        _REAL_, intent(in)  :: x(3*natom),mass(natom)
        _REAL_, intent(inout) :: f(3*natom)
        _REAL_,intent(out) ::  epot
        _REAL_ :: efires
        
        integer :: lstart, lend

        epot = 0.d0

        ! Determine local bounds (defaults to the whole array in serial runs)
        lstart = max(1, g_loc_istart3)
        lend   = g_loc_iend3
        if (lend <= 0 .or. lend < lstart) then
            lstart = 1
            lend   = 3*natom
        else
            lend   = min(lend, 3*natom)
            lstart = min(max(lstart, 1), lend)
        end if

        if (lend >= lstart) f(lstart:lend) = 0.0d0

        ! Early gate if FIRES disabled or zero strength or globally disabled
        if (fires == 0 .or. fires_k <= 0.0d0 .or. .not. fires_restraint_enabled) then
            epot = 0.0d0
            return
        end if

        call calculate_fires_force(x,mass,natom, f, efires)
        epot = efires
    end subroutine fires_force

end module qmmm_fires_module
