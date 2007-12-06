!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd / max3d
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
! M4_FT -> real(kind=8) or complex(kind=8)
!

module fdtd
 
  use constant
  use grid
 
  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'fdtd'

  ! --- Public Methods

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

  M4_FT, allocatable, dimension(:, :, :) :: Ex, Ey, Ez
  M4_FT, allocatable, dimension(:, :, :) :: Hx, Hy, Hz
  real(kind=8), allocatable, dimension(:, :, :) :: EPSINV

 
contains

!----------------------------------------------------------------------

  subroutine InitializeFdtd(sfx)

    implicit none

    character(len=*) :: sfx

    integer :: err

    call AllocateFields
    call ReadEpsilonField(sfx)


    contains

      subroutine AllocateFields

        implicit none

        allocate(Ex(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif
        allocate(Ey(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif
        allocate(Ez(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif

        allocate(Hx(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif
        allocate(Hy(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif
        allocate(Hz(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif

        allocate(EPSINV(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
        if(err .ne. 0) then
           write(STDERR,*) "!ERROR OUT OF MEMORY: AllocateFields/fdtd"
           stop
        endif

        Ex = 0.0
        Ey = 0.0
        Ez = 0.0

        Hx = 0.0
        Hy = 0.0
        Hz = 0.0

        EPSINV = 1.0

      end subroutine AllocateFields


      subroutine ReadEpsilonField(sfx)

        implicit none

        character(len=*) :: sfx

        integer :: ios, i, j, k
        real(kind=8) :: val

        character(len=STRLNG) :: file

        file = cat2(pfxepsilon,sfx)

        open(UNITTMP, FILE=file, STATUS='unknown')
        do k=KBEG,KEND+1
           do j=JBEG,JEND+1
              do i=IBEG, IEND+1
                 read(UNITTMP,*) val
                 EPSINV(i,j,k)=1.0/val
              end do
           end do
        end do
        close(UNITTMP)

      end subroutine ReadEpsilonField


  end subroutine InitializeFdtd

!----------------------------------------------------------------------

  subroutine FinalizeFdtd

    implicit none

    deallocate(Hz)
    deallocate(Hy)
    deallocate(Hx)
    deallocate(Ez)
    deallocate(Ey)
    deallocate(Ex)
    deallocate(EPSINV)

  end subroutine FinalizeFdtd

!----------------------------------------------------------------------

  subroutine StepH

    implicit none
    
    real(kind=8) :: dtx, dty, dtz
    M4_FT :: Exh,Eyh,Ezh
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
    M4_FT :: Hxh, Hyh, Hzh
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
