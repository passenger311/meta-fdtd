!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd / meta3
!
!  fdtd core algorithm.
!
!  subs:
!
!  InitializeFdtd
!    AllocateFields
!    ReadEpsilonField
!  FinalizeFdtd
!  StepE
!  StepH
!
!----------------------------------------------------------------------


!======================================================================
!
! m4 macro-preprocessor runs over this file and replaces
! M4 FTYPE -> real(kind=8) or complex(kind=8)
!

module fdtd
 
  use constant
  use grid
 
  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'FDTD'
  logical, private :: modconfigured

  ! --- Public Methods

  public :: ReadConfigFdtd
  public :: InitializeFdtd
  public :: FinalizeFdtd
  public :: StepE
  public :: StepH

  ! --- Public Data

  public :: Ex, Ey, Ez, Hx, Hy, Hz, EPSINV
  public :: pfxepsilon

  ! --- Constants

  character(len=20), parameter :: pfxepsilon = 'epsilon'

  ! --- Data

  M4_FTYPE, allocatable, dimension(:, :, :) :: Ex, Ey, Ez
  M4_FTYPE, allocatable, dimension(:, :, :) :: Hx, Hy, Hz
  real(kind=8), allocatable, dimension(:, :, :) :: EPSINV

 
contains

!----------------------------------------------------------------------

 subroutine ReadConfigFdtd(funit,string)

   integer :: funit
   character(len=*) :: string

   integer :: ios

   M4_WRITE_DBG({". enter ReadConfigFdtd/fdtd"})
      
   if ( string .ne. "(FDTD" ) then
      M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigFdtd/fdtd"})
   endif
   
! TODO: read field output structures here

   read(funit,*,iostat=ios) string
   M4_WRITE_DBG({"read terminator: ", TRIM(string)})
   
   if ( string .ne. ")" ) then
      M4_FATAL_ERROR({"BAD SECTION TERMINATOR: ReadConfigFdtd/fdtd"})
   endif
   
   modconfigured = .true.

   M4_WRITE_DBG({". exit ReadConfigFdtd/fdtd"})
   
 end subroutine ReadConfigFdtd
 

!----------------------------------------------------------------------

 subroutine InitializeFdtd
   
   integer :: err
    
    M4_WRITE_DBG({". enter InitializeFdtd/fdtd"})

   if ( .not. modconfigured ) then
      M4_FATAL_ERROR({"NOT CONFIGURED: InitializeFdtd/fdtd"})
   endif
 
   call AllocateFields
   call ReadEpsilonField

    M4_WRITE_DBG({". exit InitializeFdtd/fdtd"})
   
 contains
      
   subroutine AllocateFields

     M4_WRITE_DBG({". enter InitializeFdtd/fdtd . AllocateFields"})
     
     allocate(Ex(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Ey(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Ez(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Hx(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")

     allocate(Hy(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Hz(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(EPSINV(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     Ex = 0.0
     Ey = 0.0
     Ez = 0.0
     
     Hx = 0.0
     Hy = 0.0
     Hz = 0.0
     
     EPSINV = 1.0

     M4_WRITE_DBG({". exit InitializeFdtd/fdtd . AllocateFields"})

   end subroutine AllocateFields
   
   
   subroutine ReadEpsilonField

     integer :: ios, i, j, k
     real(kind=8) :: val
     
     character(len=STRLNG) :: file
     
     M4_WRITE_DBG({". enter InitializeFdtd/fdtd . ReadEpsilonField"})

     file = cat2(pfxepsilon,mpi_sfxin)
     
     M4_WRITE_DBG({"trying to open ", file})
     
     open(UNITTMP, FILE=file, STATUS='old', IOSTAT=ios)
     M4_OPEN_ERROR(ios,file)
     
     do k=KBEG,KEND+1
        do j=JBEG,JEND+1
           do i=IBEG, IEND+1
              read(UNITTMP,*) val
              EPSINV(i,j,k)=1.0/val
           end do
        end do
     end do
     
     close(UNITTMP)

     M4_WRITE_DBG({". exit InitializeFdtd/fdtd . ReadEpsilonField"})
     
   end subroutine ReadEpsilonField
   
 end subroutine InitializeFdtd

!----------------------------------------------------------------------

 subroutine FinalizeFdtd
   
   M4_WRITE_DBG({". enter FinalizeFdtd/fdtd"})
   
   deallocate(Hz)
   deallocate(Hy)
   deallocate(Hx)
   deallocate(Ez)
   deallocate(Ey)
   deallocate(Ex)
   deallocate(EPSINV)

   M4_WRITE_DBG({". exit FinalizeFdtd/fdtd"})
   
 end subroutine FinalizeFdtd
 
 !----------------------------------------------------------------------

  subroutine StepH

    implicit none
    
    real(kind=8) :: dtx, dty, dtz
    M4_FTYPE :: Exh,Eyh,Ezh
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz


    ! H in GT+1/2DT

    do k=KSIG, KEIG-1
       do j=JSIG, JEIG-1
          do i=ISIG, IEIG-1

             Exh=Ex(i,j,k)
             Eyh=Ey(i,j,k)
             Ezh=Ez(i,j,k)

             Hx(i,j,k) =  Hx(i,j,k) &
                  - dty*( Ez(i,j+1,k) - Ezh ) &
                  + dtz*( Ey(i,j,k+1) - Eyh )
             Hy(i,j,k) = Hy(i,j,k) &
                  - dtz*( Ex(i,j,k+1) - Exh ) &
                  + dtx*( Ez(i+1,j,k) - Ezh )
             Hz(i,j,k) = Hz(i,j,k) &
                  - dtx*( Ey(i+1,j,k) - Eyh ) &
                  + dty*( Ex(i,j+1,k) - Exh )

          enddo
       enddo
    enddo
        
  end subroutine StepH

!----------------------------------------------------------------------

  subroutine StepE

    implicit none
    
    real(kind=8) :: dtx, dty, dtz
    M4_FTYPE :: Hxh, Hyh, Hzh
    real(kind=8) :: epsinvx, epsinvy, epsinvz
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz
    
    ! E in GT + DT
    do k=KSIG, KEIG-1
       do j=JSIG, JEIG-1
          do i=ISIG, IEIG-1

             Hxh=Hx(i,j,k)
             Hyh=Hy(i,j,k)
             Hzh=Hz(i,j,k) 

             epsinvx = 0.5*(EPSINV(i,j,k) +  EPSINV(i+1,j,k))
             epsinvy = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j+1,k))
             epsinvz = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j,k+1))

             Ex(i,j,k) =  Ex(i,j,k) +  epsinvx* &
                  ( dty*( Hzh - Hz(i,j-1,k) ) &
                  - dtz*( Hyh - Hy(i,j,k-1) ))
             Ey(i,j,k) =  Ey(i,j,k) +  epsinvy* &
                  ( dtz*( Hxh - Hx(i,j,k-1) ) &
                  - dtx*( Hzh - Hz(i-1,j,k) ))
             Ez(i,j,k) =  Ez(i,j,k) +  epsinvz* &
                  ( dtx*( Hyh - Hy(i-1,j,k) )  &
                  - dty*( Hxh - Hx(i,j-1,k) ))

          enddo
       enddo
    enddo
    
  end subroutine StepE

!----------------------------------------------------------------------

end module fdtd

!
! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
