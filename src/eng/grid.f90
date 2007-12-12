!-*- F90 -*------------------------------------------------------------
!
!  module: grid / meta3
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
  use reglist

  implicit none
  public
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'GRID'
  logical, private :: modconfigured = .false.

  ! --- Public Variables  

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

  type(T_REG) :: fdtdreg

contains

!----------------------------------------------------------------------

  subroutine ReadConfigGrid(funit,string)

    integer :: funit
    character(len=*) :: string
    
    integer :: ios

    M4_WRITE_DBG({". enter ReadConfigGrid"})

    M4_WRITE_DBG({"received token: ", TRIM(string)})
    if ( string .ne. "(GRID" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigGrid"})
    endif

    read(funit,*) PARTITIONS
    M4_WRITE_DBG({"read PARTITIONS: ", PARTITIONS})
    read(funit,*) NCYCMAX
    M4_WRITE_DBG({"read NCYCMAX: ", NCYCMAX})
    read(funit,*) DT
    M4_WRITE_DBG({"read DT: ", DT})
    read(funit,*) IBEG, IEND    ! from ... to ranges
    M4_WRITE_DBG({"read IBEG/IEND: ", IBEG, IEND})
    read(funit,*) JBEG, JEND
    M4_WRITE_DBG({"read JBEG/JEND: ", JBEG, JEND})
    read(funit,*) KBEG, KEND
    M4_WRITE_DBG({"read KBEG/KEND: ", KBEG, KEND})
    read(funit,*,iostat=ios) string
    M4_WRITE_DBG({"read terminator: ", TRIM(string)})

    ! TODO: add some checks on numerical values

    if ( string(1:1) .ne. ")" ) then
       M4_FATAL_ERROR({"BAD SECTION TERMINATOR: ReadConfigGrid"})
    endif

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
    M4_WRITE_DBG({"init KMIN/KMAX: ", JMIN, JMAX})
    M4_WRITE_DBG({"init KMIN/KMAX: ", KMIN, KMAX})
    
    ! ranges for PML sheets [ISIG,IEIG]
    
    ISIG=IBEG
    IEIG=IMAX
    JSIG=IBEG
    JEIG=JMAX
    KSIG=KBEG
    KEIG=KMAX
    M4_WRITE_DBG({"init ISIG/IEIG: ", ISIG, IEIG})
    M4_WRITE_DBG({"init KSIG/KEIG: ", JSIG, JEIG})
    M4_WRITE_DBG({"init KSIG/KEIG: ", KSIG, KEIG})

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
