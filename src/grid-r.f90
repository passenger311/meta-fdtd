!----------------------------------------------------------------------
!
!  module: grid(-r) / max3d
!
!  grid definition and field allocation.
!
!  subs:
!
!    CreateGrid
!      ReadConfig
!      Initialize
!    DestroyGrid
!    EchoGrid
!
!----------------------------------------------------------------------


module grid

  use constant
  use strings
  use mpiworld

  implicit none
  save

  ! --- Constants

  character(len=255), parameter :: pfxgrid = 'grid'

  ! --- Variables  

  integer :: IBEG, IEND
  integer :: JBEG, JEND
  integer :: KBEG, KEND

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

  subroutine CreateGrid(sfx)

    implicit none

    character(len=*) :: sfx

    call ReadConfig(sfx)
    call Initialize()

    contains

      subroutine ReadConfig(sfx)

        implicit none

        character(len=*) :: sfx

        character(len=255) :: file 
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

        if ( PARTITIONS .ne. numproc .and. mpi_started .ne. 0 ) then
           write(6,*) "config: number of read in parititons does not match mpi numproc"
           stop
        end if

        ! Inner ranges are [IBEG,IEND] etc., outer ranges [IMIN,IMAX]

        IMIN = IBEG-1          
        IMAX = IEND+1
        JMIN = JBEG-1
        JMAX = JEND+1
        KMIN = KBEG-1
        KMAX = KEND+1

       ! Ranges for PML sheaths [ISIG,IEIG]

        ISIG=IBEG
        IEIG=IMAX
        JSIG=IBEG
        JEIG=JMAX
        KSIG=KBEG
        KEIG=KMAX

      end subroutine ReadConfig

      subroutine Initialize() 

        implicit none


      end subroutine Initialize

  end subroutine CreateGrid



  subroutine DestroyGrid()

    implicit none

  end subroutine DestroyGrid



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


end module grid







