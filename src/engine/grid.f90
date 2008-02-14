!-*- F90 -*------------------------------------------------------------
!
!  module: grid / meta3
!
!  grid definition.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


!======================================================================
!
!
!

module grid

  use constant
  use strings
  use parse
  use mpiworld
  use reglist

  implicit none
  public
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'GRID'
  logical, private :: modconfigured = .false.

  ! --- Public Variables  

  integer :: IBEG = 0, IEND = 1
  integer :: JBEG = 0, JEND = 1
  integer :: KBEG = 0, KEND = 1
  integer :: GRIDSIZE, FIELDSIZE

  integer :: IMIN = -1, IMAX = 2
  integer :: JMIN = -1, JMAX = 2   
  integer :: KMIN = -1, KMAX = 2

  integer :: IBIG = 0, IEIG = 1 
  integer :: JBIG = 0, JEIG = 1
  integer :: KBIG = 0, KEIG = 1

  real(8) :: SX  = 1.0    
  real(8) :: SY  = 1.0
  real(8) :: SZ  = 1.0

  integer :: NCYCMAX = 0
  real(8) :: GT  = 0.0
  real(8) :: DT  = 0.577

  integer :: DIM

  integer :: PARTITIONS

  type(T_REG) :: fdtdreg

contains

!----------------------------------------------------------------------

  subroutine ReadConfigGrid(funit,lcount,string,sdim)

    integer :: funit, lcount
    character(len=*) :: string
    integer :: sdim
    
    character(len=LINELNG) :: line
    integer :: v(2)

    M4_WRITE_DBG({". enter ReadConfigGrid"})

    M4_WRITE_DBG({"received token: ", TRIM(string)})
    if ( string .ne. "(GRID" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigGrid"})
    endif

    JBEG = 0
    JEND = 0
    KBEG = 0
    KEND = 0
    
    call readint(funit, lcount, DIM)
    M4_WRITE_DBG({"read DIM: ", DIM})
    if ( DIM .ne. sdim ) then
       M4_FATAL_ERROR({"ENGINE HAS NOT BEEN COMPILED FOR DIM = ",TRIM(i2str(DIM))})
    end if
    call readint(funit, lcount, PARTITIONS)
    M4_WRITE_DBG({"read PARTITIONS: ", PARTITIONS})
    call readint(funit, lcount, NCYCMAX)
    M4_WRITE_DBG({"read NCYCMAX: ", NCYCMAX})
    call readfloat(funit, lcount, DT)
    M4_WRITE_DBG({"read DT: ", DT})
    call readints(funit, lcount, v ,2)
    IBEG=v(1)
    IEND=v(2)
    M4_WRITE_DBG({"read IBEG/IEND: ", IBEG, IEND})
    call readints(funit, lcount, v ,2)
    JBEG=v(1)
    JEND=v(2)
    if ( DIM .lt. 2 ) JEND = JBEG
    M4_WRITE_DBG({"read JBEG/JEND: ", JBEG, JEND})
    call readints(funit, lcount, v ,2)
    KBEG=v(1)
    KEND=v(2)
    if ( DIM .lt. 3 ) KEND = KBEG
    M4_WRITE_DBG({"read KBEG/KEND: ", KBEG, KEND})
    
    call readtoken(funit, lcount, ")GRID")

    modconfigured = .true.

    M4_WRITE_DBG({"initializing grid early!"})

    call InitializeGrid

    fdtdreg = CreateBoxRegObj(IBEG, IEND, 1, JBEG, JEND, 1, KBEG, KEND, 1)

    M4_WRITE_DBG({". exit ReadConfigGrid"})

  end subroutine ReadConfigGrid

!----------------------------------------------------------------------

  subroutine InitializeGrid
    
    M4_WRITE_DBG({". enter InitializeGrid"})

    if ( .not. modconfigured ) then
       M4_FATAL_ERROR({"NOT CONFIGURED: InitializeGrid"})
    endif

    if ( PARTITIONS .ne. numproc .and. mpi_started .ne. 0 ) then
       M4_FATAL_ERROR({"GRID PARTITIONS MISMATCH: InitializeGrid"})
    end if
    
    ! inner ranges are [IBEG,IEND] etc., outer ranges [IMIN,IMAX]

    IMIN = IBEG-1          
    IMAX = IEND+1
    JMIN = JBEG-1
    JMAX = JEND+1
    KMIN = KBEG-1
    KMAX = KEND+1

    M4_WRITE_DBG({"init IMIN/IMAX: ", IMIN, IMAX})
    M4_WRITE_DBG({"init JMIN/JMAX: ", JMIN, JMAX})
    M4_WRITE_DBG({"init KMIN/KMAX: ", KMIN, KMAX})

    IBIG=IBEG
    IEIG=IEND
    JBIG=JBEG
    JEIG=JEND
    KBIG=KBEG
    KEIG=KEND

    M4_WRITE_DBG({"set IBIG/IEIG: ", IBIG, IEIG})
    M4_WRITE_DBG({"set JBIG/JEIG: ", JBIG, JEIG})
    M4_WRITE_DBG({"set KBIG/KEIG: ", KBIG, KEIG})

    FIELDSIZE = (IMAX - IMIN + 1)*(JMAX - JMIN + 1)*(KMAX - KMIN + 1)
    GRIDSIZE = (IEND - IBEG + 1)*(JEND - JBEG + 1)*(KEND - KBEG + 1)
    M4_WRITE_DBG({"init GRIDSIZE: ", GRIDSIZE})

    M4_WRITE_DBG({". exit InitializeGrid"})
    
  end subroutine InitializeGrid

!----------------------------------------------------------------------

  subroutine FinalizeGrid

    M4_WRITE_DBG({". FinalizeGrid"})

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
    write(6,*) 'DIM:              ', DIM
    write(6,*) 'NCYCMAX:          ', NCYCMAX
    write(6,*) 'IBIG, IEIx/G:     ', IBIG, IEIG
    write(6,*) 'JBIG, JEIG:       ', JBIG, JEIG
    write(6,*) 'KBIG, KEIG:       ', KBIG, KEIG

  end subroutine EchoGrid

!----------------------------------------------------------------------

end module grid

!
! Authors:  J.Hamm, A.Klaedtke, S.Scholz, 
! Modified: 4/12/2007
!
!======================================================================
