#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

module HessianNetCDF_mod
  
#ifdef BINTRAJ
  use AmberNetcdf_mod, only : checkNCerror
  use AmberNetcdf_mod, only : NetcdfFileExists, NC_openRead
  use AmberNetcdf_mod, only : NC_close, GetDimInfo
  use AmberNetcdf_mod, only : NC_error
#endif
  
  implicit none
  private
  
#ifdef BINTRAJ

  integer, parameter :: mdout = 6

  !
  ! Atomic coordinates and unit cell information
  ! This follows the same variable conventions as the
  ! amber netcdf trajectory format
  !
  
  character(5), parameter  :: NCLABEL = "label"
  integer, parameter       :: NCLABELLEN = 5
  character(12), parameter :: NCSPATIAL = "spatial"
  character(12), parameter :: NCCELL_SPATIAL = "cell_spatial"
  character(12), parameter :: NCCELL_LENGTHS = "cell_lengths"
  character(12), parameter :: NCCELL_ANGULAR = "cell_angular"
  character(11), parameter :: NCCELL_ANGLES = "cell_angles"
  character(11), parameter :: NCCOORDS = "coordinates"
  character(6),  parameter :: NCNATOMS = "natoms"

  !
  ! The 1-based atom indexes of the atoms saved within the netcdf file
  ! This is usually only a very small portion of the full system.
  !
  character(6),  parameter :: NCATIDXS = "atidxs"
  !
  ! The atomic numbers of the stored atoms
  !
  character(6),  parameter :: NCATNUMS = "atnums"
  !
  ! The length of the conmat array (5,ncon)
  !
  character(6), parameter :: NCCONMAT = "conmat"
  !
  ! The list of constraints
  !
  character(4),  parameter :: NCNCON = "ncon"
  !
  ! A label for a dimension of length 5. Each constraint has 5 integers.
  !
  character(4),  parameter :: NCCON5 = "con5"
  !
  ! The hessian dimension (row or column) is 3*natoms
  !
  character(8),  parameter :: NCNATOMS3 = "natoms3"
  !
  ! A dimension of length 1, so we can store a scalar value
  !
  character(8),  parameter :: NCCNT = "cnt"
  !
  ! The number of hessian columns saved in within this netcdf file
  !
  character(10),  parameter :: NCNSAVEDCOLS = "nsavedcols"
  !
  ! The 1-based column indexes of the hessian columns saved within
  ! this file
  !
  character(6),  parameter :: NCCOLIDX = "colidx"
  !
  ! The column vectors of hessian values
  !
  character(7),  parameter :: NCCOLVALS = "colvals"

  !
  ! This is space to save isotopologue analysis in the
  ! future
  !
  character(8),  parameter :: NCNISO    = "niso"
  character(8),  parameter :: NCNVIB    = "nvib"
  character(8),  parameter :: NCISOMASK = "isomask"
  character(8),  parameter :: NCATMASS  = "atmass"
  character(8),  parameter :: NCREDMASS = "redmass"
  character(8),  parameter :: NCXMAT    = "xmat"
  character(8),  parameter :: NCFREQS   = "freqs"
  character(8),  parameter :: NCDISPMAT = "dispmat"


  integer :: hess_ncid = -1
  integer :: coordVarID = -1
  integer :: atidxsVarID = -1
  integer :: atnumsVarID = -1
  integer :: cell_lengthVarID = -1
  integer :: cell_angleVarID = -1
  integer :: colidxVarID = -1
  integer :: colvalsVarID = -1
  integer :: nsavedcolsVarID = -1
  integer :: cntDimID = -1
  integer :: hess_entry = -1
  integer :: HESSIAN_UNIT = 79

  
  public :: CheckHessianFileExists
  public :: OpenHessianNetCDF
  public :: CloseHessianNetCDF
  public :: WriteHessianGeometry
  public :: WriteHessianColumn
  public :: GetHessianNumStoredCols
  public :: GetHessianStoredColIdxs
  public :: GetHessianCrds
  
#endif
  
contains

#ifdef BINTRAJ

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_OPENWRITE()
!> @brief Open netcdf file for writing.
!> @param fname Netcdf file to open
!> @param ncid Netcdf ID of opened file.
!> @return false if successful, true if error occurs.
logical function NC_openWrite(fname, ncid)
  use netcdf
  implicit none
  character(*), intent(in) :: fname
  integer, intent(out)     :: ncid
  NC_openWrite=.true.
  !
  ! The NF90_WRITE opens in "classic" mode, which only allows
  ! a single unlimited dimension. We need to use the "extended"
  ! mode by specifying NF90_NETCDF4
  !
  if (NC_error( nf90_open( fname, NF90_WRITE, ncid ))) return
  !
  ! Actually, if Amber was not linked to an external version of
  ! netcdf, then netcdf4 is not available.
  !
  !if (NC_error( nf90_open( fname, NF90_NETCDF4, ncid ))) return
  NC_openWrite=.false.
end function NC_openWrite


  
!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION NC_DEFINE_VAR
!> @brief Wrapper around nf90_def_var for multiple dimensions.
!> @param vcid NetCDF ID of file to define variable for.
!> @param varname Name of variable to define.
!> @param dtype Data type of variable.
!> @param ndim Number of dimensions in variable.
!> @param dimID Array holding size of each dimension.
  integer function NC_define_var(ncid, varname, dtype, ndim, dimID)
    use netcdf
    implicit none
    ! Passed vars
    integer, intent(in)      :: ncid
    character(*), intent(in) :: varname
    integer, intent(in)      :: dtype
    integer, intent(in)      :: ndim
    integer, dimension(NF90_MAX_DIMS), intent(in) :: dimID
    
    if (ndim .eq. 0) then
       if (NC_error( nf90_def_var(ncid, varname, dtype, NC_define_var) )) &
            & NC_define_var=-1
    else if (ndim .eq. 1) then
       if (NC_error( nf90_def_var(ncid, varname, dtype, (/ dimID(1) /), &
            & NC_define_var) )) NC_define_var=-1
    else if (ndim .eq. 2) then
       if (NC_error( nf90_def_var(ncid, varname, dtype, (/ dimID(1), dimID(2) /), &
            & NC_define_var) )) NC_define_var=-1
    else ! ndim .eq. 3
       if (NC_error( nf90_def_var(ncid, varname, dtype, (/ dimID(1), dimID(2), dimID(3) /),&
            & NC_define_var) )) NC_define_var=-1
    endif
  end function NC_define_var

!--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION GETVARINFO
!> @return variable ID of specified variable if it exists.
!> @return -1 if the variable does not exist.
  integer function GetVarInfo(ncid, attribute)
    use netcdf
    implicit none
    integer, intent(in)      :: ncid
    character(*), intent(in) :: attribute
    integer VID
    GetVarInfo=-1
    if ( nf90_inq_varid(ncid, attribute, VID) .eq. NF90_NOERR ) &
         & GetVarInfo=VID
  end function GetVarInfo


  !--------------------------------------------------------------------
!> MODULE AMBERNETCDF FUNCTION CHECKATTRTEXT
!> @brief Check that text associated with attribute of vid matches textIn
  logical function CheckAttrText(ncid, vid, attribute, textIn)
    use netcdf
    implicit none
    integer, intent(in)      :: ncid, vid
    character(*), intent(in) :: attribute, textIn
    ! Local vars
    character(80) :: text
    CheckAttrText=.true.
    if (NC_error( nf90_get_att(ncid, vid, attribute, text),&
         & attribute)) return
    if ( text .ne. textIn ) &
         & write(mdout,'(6a)') 'Warning: In netcdf file, expected ', &
         & textIn, ' for attribute ', attribute, ', got ', trim(text)
    CheckAttrText=.false.
  end function CheckAttrText



!--------------------------------------------------------------------
!> MODULE AMBERNETCDF SUBROUTINE NC_SETUPBOX
!> @brief Check box info, setup cellLengthVID and cellAngleVID.
!> @param cellLengthVID Cell lengths variable, set to -1 if no box or error.
!> @param cellAngleVID Cell angles variable, set to -1 if no box or error.
  subroutine NC_setupBox(ncid, cellLengthVID, cellAngleVID)
    use netcdf
    implicit none
    integer, intent(in)  :: ncid
    integer, intent(out) :: cellLengthVID, cellAngleVID
    
    cellAngleVID = -1
    cellLengthVID = GetVarInfo( ncid, NCCELL_LENGTHS )
    if ( cellLengthVID.ne.-1 ) then
       if (NC_error(nf90_inq_varid(ncid,NCCELL_ANGLES,cellAngleVID),&
            & 'Getting cell angles')) then
          cellLengthVID = -1
          return
       endif
    endif
  end subroutine NC_setupBox



  
  logical function NC_checkHessianConvention(ncid)
    use netcdf
    implicit none
    
    integer,intent(in) :: ncid
    
    character(80) :: attribute

    NC_checkHessianConvention = .true.
    
    if (NC_error( nf90_get_att(ncid, NF90_GLOBAL, 'Conventions', attribute), &
         & 'Getting Conventions attribute')) return

    if (attribute.ne."HESSIAN") then
       write(mdout,'(a)') 'ERROR: NetCDF traj has Conventions that are not HESSIAN.'
       return
    endif

    NC_checkHessianConvention = .false.

  end function NC_checkHessianConvention
  
    
    
  logical function NC_setupHessCols(ncid, ncols, nsavedcolsVID, colidxVID, colvalsVID)
    use netcdf
    implicit none
    integer, intent(in)  :: ncid
    integer, intent(out) :: ncols, nsavedcolsVID, colidxVID, colvalsVID
    integer :: t1(1)
    
    NC_setupHessCols=.true.
    !colDID = GetDimInfo( ncid, NCHESSCOLS, ncols )
    !if (colDID.eq.-1) return
    

    nsavedcolsVID = GetVarInfo( ncid, NCNSAVEDCOLS )
    if ( nsavedcolsVID .eq. -1 ) return

    t1 = 0
    if (NC_error(nf90_get_var(ncid, nsavedcolsVID, t1, &
         start = (/ 1 /),&
         count = (/ 1 /)),&
         'Getting nsavedcols')) return
    ncols = t1(1)

    colidxVID = GetVarInfo( ncid, NCCOLIDX )
    if ( colidxVID .eq. -1 ) return

    colvalsVID = GetVarInfo( ncid, NCCOLVALS )
    if ( colvalsVID .eq. -1 ) return
    
    NC_setupHessCols=.false.
    
  end function NC_setupHessCols

  
  
  logical function NC_setupCoords(ncid, natoms, coordVID, atidxsVID, atnumsVID)
    
    use netcdf
    implicit none
    
    integer, intent(in)  :: ncid
    integer, intent(out) :: natoms
    integer, intent(out) :: coordVID, atidxsVID, atnumsVID

    integer :: spatial
    integer :: natomsDID,spatialDID,spatialVID
    
    NC_setupCoords = .true.

    natomsDID = GetDimInfo( ncid, NCNATOMS, natoms )
    if (natomsDID .eq. -1) return

    coordVID = GetVarInfo( ncid, NCCOORDS )
    if ( coordVID .eq. -1 ) return
    
    if (CheckAttrText(ncid, coordVID, "units", "bohr")) return

    spatialDID = GetDimInfo( ncid, NCSPATIAL, spatial )
    if (spatialDID .eq. -1) return
    
    if (spatial .ne. 3) then
       write(mdout,'(a,i6)') 'Error: Netcdf: Expected 3 spatial dimensions, got ', &
            & spatial
       return
    endif
    
    if ( NC_error(nf90_inq_varid(ncid, NCSPATIAL, spatialVID),&
         & 'Getting spatial VID') ) return

    atidxsVID = GetVarInfo( ncid, NCATIDXS )
    if ( atidxsVID .eq. -1 ) return

    atnumsVID = GetVarInfo( ncid, NCATNUMS )
    if ( atnumsVID .eq. -1 ) return
    
    NC_setupCoords = .false.
    
  end function NC_setupCoords

  


    
  logical function NC_createHessian( &
       & filename, owrite, natomIn, hasBox, &
       & ncon, conmat, &
       & ncid, coordVID, atidxsVID, atnumsVID, &
       & cellLengthVID, cellAngleVID, &
       & nsavedcolsVID, colidxVID, colvalsVID )
    
    use netcdf
    implicit none

    character(*), intent(in)  :: filename
    character,    intent(in)  :: owrite
    integer,      intent(in)  :: natomIn
    logical,      intent(in)  :: hasBox
    integer,      intent(in)  :: ncon
    integer,      intent(in)  :: conmat(5,ncon)
    
    integer,      intent(out) :: ncid
    integer,      intent(out) :: coordVID
    integer,      intent(out) :: atidxsVID, atnumsVID
    integer,      intent(out) :: cellLengthVID, cellAngleVID
    integer,      intent(out) :: nsavedcolsVID, colidxVID, colvalsVID
    !integer,      intent(out) :: hesscolsDID
    
    
    ! Local variables
    integer, dimension(NF90_MAX_DIMS) :: dimensionID
    integer  :: NDIM, dataType, cmode, ierr
    integer  :: spatialDID, spatialVID
    integer  :: cell_spatialDID, labelDID, cell_angularDID
    integer  :: cell_spatialVID, cell_angularVID
    integer  :: natoms3DID, natomsDID

    integer :: nisoDID
    integer :: nvibVID
    integer :: isomaskVID
    integer :: atmassVID
    integer :: redmassVID
    integer :: xmatVID
    integer :: freqsVID
    integer :: dispmatVID
    
    integer :: nconDID
    integer :: con5DID
    integer :: conmatVID
    integer :: cntDID
    
    NC_createHessian=.true. ! Assume error until success
    
    dataType = NF90_DOUBLE
    cmode = nf90_64bit_offset
    !cmode = NF90_NETCDF4
    
    if (owrite == 'N') cmode = ior(cmode, nf90_noclobber)
    
    if ( len(filename) .eq. 0 ) return
    
    ierr = nf90_create(path=filename, cmode=cmode, ncid=ncid)
    if (ierr == nf90_eexist) then
       write(mdout,'(a,a)') 'Error: NetCDF file exists: ',filename
       return
    end if
    
    if (NC_error( ierr )) return
    
    !!write(6,*)"nisoDID"
    ! hesscolsDID, dim=VARIABLE
    if (NC_error( nf90_def_dim( ncid, NCNISO, NF90_UNLIMITED, nisoDID ), &
         & 'Defining niso dimension.' )) return

    !!write(6,*)"hesscolsDID"
    ! hesscolsDID, dim=VARIABLE
    if (NC_error( nf90_def_dim( ncid, NCCNT, 1, cntDID ), &
         & 'Defining cnt dimension.' )) return

    
    !write(6,*)"spatialDID"
    ! spatialDID, dim=3
    if ( NC_error( nf90_def_dim( ncid, NCSPATIAL, 3, spatialDID), &
         & 'Defining spatial dimension')) return
    
    !write(6,*)"natomsDID"
    ! natomsDID, dim=natomIn
    if (NC_error( nf90_def_dim( ncid, NCNATOMS, natomIn, natomsDID ), &
         & 'Defining natoms dimension.' )) return
    
    !write(6,*)"natoms3DID"
    ! natoms3DID, dim=3*natomIn
    if (NC_error( nf90_def_dim( ncid, NCNATOMS3, 3*natomIn, natoms3DID ), &
         & 'Defining natoms3 dimension.' )) return

    if ( ncon > 0 ) then
       !write(6,*)"nbondpairsDID"
       ! nbondpairsDID, dim=nbondpairs
       if (NC_error( nf90_def_dim( ncid, NCNCON, ncon, nconDID ), &
            & 'Defining ncon dimension.' )) return
       if (NC_error( nf90_def_dim( ncid, NCCON5, 5, con5DID ), &
            & 'Defining con5 dimension.' )) return
    end if


    
    !write(6,*)"spatialVID"
    ! spatialVID is size (3), type=char
    dimensionID(1) = spatialDID
    spatialVID = NC_define_var(ncid, NCSPATIAL, NF90_CHAR, 1, dimensionID)
    if (spatialVID .eq. -1) return ! 'Defining spatial variable'
    
    !write(6,*)"coordVID"
    ! coordVID is size (3,natom), type=double
    dimensionID(1) = spatialDID
    dimensionID(2) = natomsDID
    coordVID = NC_define_var(ncid, NCCOORDS, dataType, 2, dimensionID)
    if (coordVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, coordVID, "units", "bohr"), &
         & 'Writing coordinates variable units.')) return
    
    !write(6,*)"atidxsVID"
    ! atidxsVID is size (natom), type=int
    dimensionID(1) = natomsDID
    atidxsVID = NC_define_var(ncid, NCATIDXS, NF90_INT, 1, dimensionID)
    if (atidxsVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, atidxsVID, "units", "1-based index"), &
         & 'Writing atidxs variable units.')) return
    
    !write(6,*)"atnumsVID"
    ! atnumsVID is size (natom), type=int
    dimensionID(1) = natomsDID
    atnumsVID = NC_define_var(ncid, NCATNUMS, NF90_INT, 1, dimensionID)
    if (atnumsVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, atnumsVID, "units", "unitless"), &
         & 'Writing atnums variable units.')) return

    if ( ncon > 0 ) then
       !write(6,*)"bondpairsVID"
       ! bondpairsVID is size (nbondpairs), type=int
       dimensionID(1) = con5DID
       dimensionID(2) = nconDID
       conmatVID = NC_define_var(ncid, NCCONMAT, NF90_INT, 2, dimensionID)
       if (conmatVID .eq. -1) return ! 'Defining bondpairs variable.')) return
       if (NC_error( nf90_put_att( ncid, conmatVID, "units", "unitless"), &
            & 'Writing conmat variable units.')) return
    end if
   
    
    
    !write(6,*)"colidxVID"
    ! colidxVID is size (VARIABLE), type=double
    dimensionID(1) = cntDID
    nsavedcolsVID = NC_define_var(ncid, NCNSAVEDCOLS, NF90_INT, 1, dimensionID)
    if (nsavedcolsVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, nsavedcolsVID, "units", "unitless"), &
         & 'Writing nsavedcols variable units.')) return
    
    !write(6,*)"colidxVID"
    ! colidxVID is size (VARIABLE), type=double
    dimensionID(1) = natoms3DID
    colidxVID = NC_define_var(ncid, NCCOLIDX, NF90_INT, 1, dimensionID)
    if (colidxVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, colidxVID, "units", "1-based column index"), &
         & 'Writing colidx variable units.')) return
    
    !write(6,*)"colvalsVID"
    ! colvalsVID is size (3*natom)
    dimensionID(1) = natoms3DID
    dimensionID(2) = natoms3DID
    colvalsVID = NC_define_var(ncid, NCCOLVALS, datatype, 2, dimensionID)
    if (colvalsVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, colvalsVID, "units", "hartree/bohr/bohr"), &
         & 'Writing colvals variable units.')) return
    
    
    !write(6,*)"nvibVID"
    ! nvibVID is size (3*natom,niso)
    !dimensionID(1) = natoms3DID
    dimensionID(1) = nisoDID
    nvibVID = NC_define_var(ncid, NCNVIB, NF90_INT, 1, dimensionID)
    if (nvibVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, nvibVID, "units", "unitless"), &
         & 'Writing nvib variable units.')) return

    !write(6,*)"isomaskVID"
    ! isomaskVID is size (natom,niso)
    dimensionID(1) = natomsDID
    dimensionID(2) = nisoDID
    isomaskVID = NC_define_var(ncid, NCISOMASK, NF90_INT, 2, dimensionID)
    if (isomaskVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, isomaskVID, "units", "1 if isosub"), &
         & 'Writing isomask variable units.')) return

    !write(6,*)"atmassVID"
    ! atmassVID is size (natom,niso)
    dimensionID(1) = natomsDID
    dimensionID(2) = nisoDID
    atmassVID = NC_define_var(ncid, NCATMASS, datatype, 2, dimensionID)
    if (atmassVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, atmassVID, "units", "au mass"), &
         & 'Writing atmass variable units.')) return

    !write(6,*)"redmassVID"
    ! redmassVID is size (3*natom,niso)
    dimensionID(1) = natoms3DID
    dimensionID(2) = nisoDID
    redmassVID = NC_define_var(ncid, NCREDMASS, datatype, 2, dimensionID)
    if (redmassVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, redmassVID, "units", "au mass"), &
         & 'Writing redmass variable units.')) return

    !write(6,*)"xmatVID"
    ! xmatVID is size (nvib,3*natom,niso)
    dimensionID(1) = natoms3DID
    dimensionID(2) = natoms3DID
    dimensionID(3) = nisoDID
    xmatVID = NC_define_var(ncid, NCXMAT, datatype, 3, dimensionID)
    if (xmatVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, xmatVID, "units", "au"), &
         & 'Writing xmat variable units.')) return

    !write(6,*)"evecsVID"
    ! evecsVID is size (nvib,3*natom,niso)
    dimensionID(1) = natoms3DID
    dimensionID(2) = natoms3DID
    dimensionID(3) = nisoDID
    dispmatVID = NC_define_var(ncid, NCDISPMAT, datatype, 3, dimensionID)
    if (xmatVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, dispmatVID, "units", "au"), &
         & 'Writing dispmat variable units.')) return

    !write(6,*)"evalsVID"
    ! evalsVID is size (nvib,niso)
    dimensionID(1) = natoms3DID
    dimensionID(2) = nisoDID
    freqsVID = NC_define_var(ncid, NCFREQS, datatype, 2, dimensionID)
    if (freqsVID .eq. -1) return ! 'Defining coordinates variable.')) return
    if (NC_error( nf90_put_att( ncid, freqsVID, "units", "au"), &
         & 'Writing freqs variable units.')) return


    
    !write(6,*)"hasbox",hasbox

    if (hasBox) then
       
       !write(6,*)"cell_spatialDID"
       ! cell_spatialDID, dim=3
       if (NC_error( nf90_def_dim( ncid, NCCELL_SPATIAL, 3, cell_spatialDID), &
            & 'Defining cell spatial dimension.' )) return
       
       !write(6,*)"labelDID"
       ! labelDID, dim=NCLABELLEN
       if (NC_error( nf90_def_dim( ncid, NCLABEL, NCLABELLEN, labelDID), &
            & 'Defining label dimension.' )) return
       
       !write(6,*)"cell_angularDID"
       ! cell_angularDID, dim=3
       if (NC_error( nf90_def_dim( ncid, NCCELL_ANGULAR, 3, cell_angularDID), &
            & 'Defining cell angular dimension' )) return
       
       !write(6,*)"cell_spatialVID"
       ! cell_spatialVID is size (3), type=char
       dimensionID(1) = cell_spatialDID
       cell_spatialVID = NC_define_var(ncid, NCCELL_SPATIAL, NF90_CHAR, 1, dimensionID)
       if (cell_spatialVID .eq. -1) return ! 'Defining cell spatial variable.'
       
       !write(6,*)"cell_angularVID"
       ! cell_angularVID is size (NCLABELLEN,3), type=char
       dimensionID(1) = labelDID;
       dimensionID(2) = cell_angularDID;
       cell_angularVID = NC_define_var(ncid, NCCELL_ANGULAR, NF90_CHAR, 2, dimensionID)
       if (cell_angularVID .eq. -1) return ! 'Defining cell angular variable.'
       
       !write(6,*)"cell_lengthVID"
       ! cellLengthVID is size (3), type=double
       dimensionID(1) = cell_spatialDID
       cellLengthVID = NC_define_var(ncid, NCCELL_LENGTHS, NF90_DOUBLE, 1, &
            & dimensionID)
       if (cellLengthVID .eq. -1) return !  'Defining cell length variable.'
       if (NC_error( nf90_put_att(ncid, cellLengthVID, "units", "bohr"),&
            & 'Writing cell length variable units' )) return

       
       !write(6,*)"cell_angleVID"
       ! cellAngleVID is size (3), type=double
       dimensionID(1) = cell_angularDID;
       cellAngleVID = NC_define_var(ncid, NCCELL_ANGLES, NF90_DOUBLE, 1, &
            dimensionID)
       if (cellAngleVID .eq. -1) return ! 'Defining cell angle variable.'
       if (NC_error( nf90_put_att(ncid, cellAngleVID, "units", "degree"), &
            & 'Writing cell angle variable units.' )) return
       
    endif
    
    !write(6,*)"Conventions"
    if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","HESSIAN"),&
         'Writing hessian conventions')) return
    !write(6,*)"Version"
    if (NC_error(nf90_put_att(ncid,NF90_GLOBAL,"ConventionVersion","1.0"),&
         'Writing conventions version')) return
    
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set fill mode
    !write(6,*)"fill"
    if (NC_error(nf90_set_fill(ncid, NF90_NOFILL, ierr), &
         & 'setting fill value.')) return
    ! End netcdf definitions
    !write(6,*)"enddef"
    if (NC_error(nf90_enddef(ncid), 'ending definitions')) return
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! We can populate the values of the labels
    !write(6,*)"put spatialVID"
    if (NC_error(nf90_put_var(ncid, spatialVID, &
         & (/ 'x', 'y', 'z' /), start = (/ 1 /), count = (/ 3 /)), &
         & 'write spatial variable')) return
    
    if ( hasBox ) then
       
       !write(6,*)"put cell_spatialVID"
       if (NC_error(nf90_put_var(ncid, cell_spatialVID, &
            & (/ 'a','b','c' /), start=(/ 1 /), count=(/ 3 /)), &
            & 'write cell spatial variable')) return
       
       !write(6,*)"put cell_angularVID"
       if (NC_error(nf90_put_var(ncid, cell_angularVID, &
            & (/ 'alpha','beta ','gamma' /), &
            & start=(/ 1, 1 /), count=(/ NCLABELLEN, 3 /)), &
            & 'write spatial variable')) return
       
    end if

    if ( ncon > 0 ) then
       if (NC_error(nf90_put_var(ncid, conmatVID, &
            & conmat, &
            & start = (/ 1, 1 /), &
            & count = (/ 5, ncon /)), &
            & 'write conmat variable')) return
    end if


    if (NC_error(nf90_put_var(ncid, nsavedcolsVID, &
         & (/ 0 /), start=(/ 1 /), count=(/ 1 /)), &
         & 'write nsavedcols')) return
    
    ! All is well, exit false.
    NC_createHessian=.false.
  end function NC_createHessian


  
  logical function NC_setupHessian( &
       & ncid, ncols, ncatom, &
       & coordVID, atidxsVID, atnumsVID, &
       & cellLengthVID, cellAngleVID, &
       & nsavedcolsVID, colidxVID, colvalsVID )
    
    use netcdf
    implicit none

    integer, intent(in)   :: ncid
    integer, intent(out)  :: ncols
    integer, intent(out)  :: ncatom
    integer, intent(out)  :: coordVID, atidxsVID, atnumsVID
    integer, intent(out)  :: cellLengthVID, cellAngleVID
    integer, intent(out)  :: nsavedcolsVID,colidxVID, colvalsVID

    character(80) :: attribute
    
    NC_setupHessian = .true.

    if (NC_checkHessianConvention(ncid)) return
    if (NC_setupHessCols(ncid, ncols, nsavedcolsVID, colidxVID, colvalsVID )) return
    if (NC_setupCoords( ncid, ncatom, coordVID, atidxsVID, atnumsVID)) return
    if ( coordVID .eq. -1 ) then
       write(mdout,'(a)') 'Error: NetCDF file has no coordinates.'
       return
    endif
    if ( atidxsVID .eq. -1 ) then
       write(mdout,'(a)') 'Error: NetCDF file has no atom indexes.'
       return
    endif
    if ( atnumsVID .eq. -1 ) then
       write(mdout,'(a)') 'Error: NetCDF file has no atomic numbers.'
       return
    endif
    
    call NC_setupBox(ncid, cellLengthVID, cellAngleVID)
    
    NC_setupHessian = .false.
  end function NC_setupHessian


  logical function NC_readHessianBox(ncid,a,b,c,alpha,beta,gamma)
    use netcdf
    implicit none
    
    integer, intent(in)           :: ncid 
    _REAL_, intent(out) :: a,b,c,alpha,beta,gamma
    double precision box(3)
     integer  :: cellLengthVID, cellAngleVID
   
    NC_readHessianBox = .true.

    a = 0.d0
    b = 0.d0
    c = 0.d0
    alpha = 0.d0
    beta = 0.d0
    gamma = 0.d0
    box = 0.d0

    cellLengthVID=-1
    cellAngleVID=-1

    
    call NC_setupBox(ncid, cellLengthVID, cellAngleVID)
    if (cellLengthVID.eq.-1 .or. cellAngleVID.eq.-1) return
    
    if (NC_error(nf90_get_var(ncid, cellLengthVID, box(1:3), &
         & start = (/ 1 /), count = (/ 3 /)),&
         & 'Getting box lengths')) return
    a = box(1)
    b = box(2)
    c = box(3)
    
    ! Cell angles
    if (NC_error(nf90_get_var(ncid, cellAngleVID, box(1:3), &
         & start = (/ 1 /), count = (/ 3 /)),&
         & 'Getting box angles')) return
    alpha = box(1)
    beta  = box(2)
    gamma = box(3)
    
    NC_readHessianBox = .true.
    
  end function NC_readHessianBox


  logical function NC_checkHessian(filename)
    !
    ! returns true if filename exists and is a hessian netcdf file
    !
    use netcdf
    implicit none
    character(*), intent(in) :: filename
    ! local
    integer ncid, ierr
    NC_checkHessian=.false.
    ierr = nf90_open( filename, NF90_NOWRITE, ncid )
    if (ierr .eq. NF90_NOERR) then
       if (.not.NC_checkHessianConvention(ncid)) NC_checkHessian=.true.
       call NC_close(ncid)
    endif
  end function NC_checkHessian

  

  subroutine OpenHessianNetCDF( fname, nactive, ncon, conmat )
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nactive
    integer,intent(in) :: ncon
    integer,intent(in) :: conmat(5,ncon)
    
    logical :: exists
    integer :: ncatom

#include "box.h"
    
    exists = NetcdfFileExists(fname)

    hess_ncid    = -1
    coordVarID   = -1
    atidxsVarID  = -1
    atnumsVarID  = -1
    cell_lengthVarID = -1
    cell_angleVarID  = -1
    nsavedcolsVarID = -1
    colidxVarID  = -1
    colvalsVarID = -1
    !hesscolsDimID = -1

    if ( .not. exists ) then
       
       if ( NC_createHessian( &
            & fname, 'N', nactive, ntb.gt.0, &
            & ncon, conmat, &
            & hess_ncid, coordVarID, atidxsVarID, atnumsVarID, &
            & cell_lengthVarId, cell_angleVarID, &
            & nsavedcolsVarID, colidxVarID, colvalsVarID ) ) then
          
          write(mdout,'(a)')"Creation of Hessian NetCDF trajectory file failed."
          call mexit(6,1)
          
       end if
       
       hess_entry = 1
       
    else
       
       if (NC_openWrite( fname, hess_ncid )) then
          
          write(mdout,'(2a)')"Failed to open Hessian NetCDF file: ",trim(adjustl(fname))
          call mexit(6,1)
          
       end if
       
       if (NC_setupHessian( &
            & hess_ncid, hess_entry, ncatom, &
            & coordVarID, atidxsVarID, atnumsVarID, &
            & cell_lengthVarID, cell_angleVarID, &
            & nsavedcolsVarID, colidxVarID, colvalsVarID )) then
          
          write(mdout,'(2a)')"Failed to read Hessian NetCDF file: ",trim(adjustl(fname))
   
       end if
       
       hess_entry = hess_entry + 1

       if ( (ntb.gt.0).neqv.(Cell_lengthVarID.ne.-1) ) then
          
          write(mdout,'(a)') 'Error: Cannot append to Hessian NetCDF, box info mismatch.'
          call mexit(6,1)
          
       endif
    
       if (ncatom .ne. nactive) then
          
          write(mdout,'(2(a,i8))') 'Error: NetCDF Hessian file has ', ncatom, &
               & ' atoms, but current selection has ', nactive
          call mexit(6,1)
          
       endif
       
    end if
    
  end subroutine OpenHessianNetCDF

  
  subroutine CloseHessianNetCDF()
    implicit none
    call NC_close(hess_ncid)
  end subroutine CloseHessianNetCDF

  
  subroutine WriteHessianGeometry(nactive,acrds,atidxs,atnums)
    use netcdf
    use constants, only: CODATA08_A_TO_BOHRS
    use nblist, only : a,b,c,alpha,beta,gamma
    implicit none

    integer,intent(in) :: nactive
    _REAL_,intent(in) :: acrds(3,nactive)
    integer,intent(in) :: atidxs(nactive)
    integer,intent(in) :: atnums(nactive)


    _REAL_,parameter :: AU_PER_ANG = CODATA08_A_TO_BOHRS
    _REAL_ :: mya, myb, myc

#include "box.h"
    
    if ( ntb>0 ) then
       
       mya = a * AU_PER_ANG
       myb = b * AU_PER_ANG
       myc = c * AU_PER_ANG
       
       call checkNCerror(nf90_put_var(hess_ncid, cell_lengthVarID, &
            & (/ mya,myb,myc/), start = (/ 1 /), &
            & count = (/ 3 /)), &
            & 'write cell lengths')
       
       call checkNCerror(nf90_put_var(hess_ncid, cell_angleVarID, &
            & (/ alpha,beta,gamma /), &
            & start = (/ 1 /), &
            & count = (/ 3 /)), &
            & 'write cell angles')
       
    end if

    call checkNCerror(nf90_put_var(hess_ncid,coordVarID, acrds, &
         & start = (/ 1, 1 /), &
         & count = (/ 3, nactive /)), &
         & 'write atom coords')

    call checkNCerror(nf90_put_var(hess_ncid,atidxsVarID, atidxs, &
         & start = (/ 1 /), &
         & count = (/ nactive /)), &
         & 'write atom indexes')

    call checkNCerror(nf90_put_var(hess_ncid,atnumsVarID, atnums, &
         & start = (/ 1 /), &
         & count = (/ nactive /)), &
         & 'write atomic numbers')
    
  end subroutine WriteHessianGeometry


  
  subroutine WriteHessianColumn(nactive,icol,data)
    use netcdf
    implicit none

    integer,intent(in) :: nactive
    integer,intent(in) :: icol
    _REAL_,intent(in) :: data(3*nactive)

    
    call checkNCerror(nf90_put_var(hess_ncid,colvalsVarID, data, &
         & start = (/ 1, hess_entry /), &
         & count = (/ 3*nactive, 1 /)), &
         & 'write hessian column values')

    call checkNCerror(nf90_put_var(hess_ncid,colidxVarID, (/icol/), &
         & start = (/ hess_entry /), &
         & count = (/ 1 /)), &
         & 'write hessian column index')
    
    call checkNCerror(nf90_sync(hess_ncid))
    
    call checkNCerror(nf90_put_var(hess_ncid,nsavedcolsVarID, (/hess_entry/), &
         & start = (/ 1 /), &
         & count = (/ 1 /)), &
         & 'write nsasvedcols')
    
    call checkNCerror(nf90_sync(hess_ncid))
    hess_entry = hess_entry + 1

  end subroutine WriteHessianColumn


  subroutine GetHessianNumStoredCols(ncols)
    use netcdf
    implicit none
    integer,intent(out) :: ncols
    
    integer :: nsavedcolsVID
    integer :: t1(1)

    t1 = 0
    ncols = 0

    !colDID = GetDimInfo( hess_ncid, NCHESSCOLS, ncols )
    !if ( colDID == -1 ) then
    !   ncols = 0
    !end if

    
    nsavedcolsVID = GetVarInfo( hess_ncid, NCNSAVEDCOLS )
    if ( nsavedcolsVID .eq. -1 ) return
    
    if (NC_error(nf90_get_var(hess_ncid, nsavedcolsVID, t1, &
         start = (/ 1 /),&
         count = (/ 1 /)),&
         'Getting nsavedcols')) return

    ncols = t1(1)

  end subroutine GetHessianNumStoredCols


  subroutine GetHessianStoredColIdxs(ncols,idxs)
    use netcdf
    implicit none
    integer,intent(in) :: ncols
    integer,intent(out) :: idxs(ncols)

    if ( ncols > 0 ) then
       idxs = 0
       
       if ( NC_error(nf90_get_var( &
            & hess_ncid,colidxVarID,idxs, &
            & start = (/ 1 /), &
            & count = (/ ncols /)), &
            & 'reading hessian column indexes') ) then
          return
       end if
    end if
    
  end subroutine GetHessianStoredColIdxs


  
  subroutine GetHessianCrds(nactive,acrds)
    use netcdf
    implicit none
    integer,intent(in) :: nactive
    _REAL_,intent(out) :: acrds(3,nactive)

    
    if ( nactive > 0 ) then
       acrds = 0.d0

       if ( NC_error(nf90_get_var( &
            & hess_ncid,coordVarID,acrds, &
            & start = (/ 1,1 /), &
            & count = (/ 3,nactive /)), &
            & 'reading hessian coords') ) then
          return
       end if
    end if
    
  end subroutine GetHessianCrds


  
  logical function CheckHessianFileExists(fname)
    implicit none
    character(len=*) fname
    CheckHessianFileExists = NetcdfFileExists(fname)
  end function CheckHessianFileExists

  
#endif
  
end module HessianNetCDF_mod
  
  
