!-*- F90 -*------------------------------------------------------------
!
!  module: grid / max3d
!
!  grid definition.
!
!  subs:
!
!    InitializeGrid
!      ReadConfig
!      Initialize
!    FinalizeGrid
!    EchoGrid
!
!----------------------------------------------------------------------

!======================================================================
!
!
!

module grid

  use constant
  use strings
  use mpiworld

  implicit none
  save

  ! --- Constants

  character(len=STRLNG), parameter :: pfxgrid = 'grid'

  ! --- Variables  

  integer :: IBEG, IEND
  integer :: JBEG, JEND
  integer :: KBEG, KEND
  integer :: GRIDSIZE

  integer :: IMIN = -1, IMAX = 2
  integer :: JMIN = -1, JMAX = 2   
  integer :: KMIN = -1, KMAX = 2
  integer :: ISIG = 0
  integer :: IEIG = 2 
  integer :: JSIG = 0
  integer :: JEIG = 2
  integer :: KSIG = 0
  integer :: KEIG = 2

  real(8) :: SX  = 1.0    
  real(8) :: SY  = 1.0
  real(8) :: SZ  = 1.0

  integer :: NCYCMAX = 0
  real(8) :: GT  = 0.0
  real(8) :: DT  = 0.577

  integer :: PARTITIONS

contains

!----------------------------------------------------------------------

  subroutine InitializeGrid(sfx)

    implicit none

    character(len=*) :: sfx

    call ReadConfig(sfx)
    call Initialize()

    contains

      subroutine ReadConfig(sfx)

        implicit none

        character(len=*) :: sfx

        character(len=STRLNG) :: file 
        integer :: err, i

        file = cat2(pfxgrid,sfx)

        open(UNITTMP,FILE=file,STATUS='unknown')
        read(UNITTMP,*) PARTITIONS
        read(UNITTMP,*) NCYCMAX
        read(UNITTMP,*) DT
        read(UNITTMP,*) IBEG, IEND    ! from ... to ranges
        read(UNITTMP,*) JBEG, JEND
        read(UNITTMP,*) KBEG, KEND

        close(UNITTMP) 


      end subroutine ReadConfig

      subroutine Initialize

        implicit none

        if ( PARTITIONS .ne. numproc .and. mpi_started .ne. 0 ) then
           write(STDERR,*) "!ERROR GRID PARTITIONS MISMATCH: InitializeGrid/grid"
           stop
        end if

        ! inner ranges are [IBEG,IEND] etc., outer ranges [IMIN,IMAX]

        IMIN = IBEG-1          
        IMAX = IEND+1
        JMIN = JBEG-1
        JMAX = JEND+1
        KMIN = KBEG-1
        KMAX = KEND+1

       ! ranges for PML sheets [ISIG,IEIG]

        ISIG=IBEG
        IEIG=IMAX
        JSIG=IBEG
        JEIG=JMAX
        KSIG=KBEG
        KEIG=KMAX

        GRIDSIZE = (IEND - IBEG)*(JEND - JBEG)*(KEND - KBEG)

      end subroutine Initialize

  end subroutine InitializeGrid

!----------------------------------------------------------------------


  subroutine FinalizeGrid

  end subroutine FinalizeGrid

!----------------------------------------------------------------------

  subroutine EchoGrid

    implicit none
    
    write(6,*) 
    write(6,*) '------------ EchoGrid() -----------'
    write(6,*) 
    write(6,*) 'IMAX, JMAX, KMAX: ', IMAX, JMAX, KMAX
    write(6,'(A20, 1E12.4)') 'DT: ', DT
    write(6,'(A20, 3E12.4)') 'SX, SY, SZ: ', SX, SY, SZ
    write(6,*) 'NCYCMAX:          ', NCYCMAX
    write(6,*) 'ISIG, IEIx/G:       ', ISIG, IEIG
    write(6,*) 'JSIG, JEIG:       ', JSIG, JEIG
    write(6,*) 'KSIG, KEIG:       ', KSIG, KEIG

  end subroutine EchoGrid

!----------------------------------------------------------------------

end module grid

!
! Authors:  J.Hamm, A.Klaedtke, S.Scholz, 
! Modified: 4/12/2007
!
!======================================================================
