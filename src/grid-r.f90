!----------------------------------------------------------------------
!
!  module: grid(-r) / max3d
!
!  grid definition and field allocation.
!
!----------------------------------------------------------------------


module grid

  use mpistart
  use constant
  implicit none
  save

  ! Variables, Defaults  

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

  real(8), allocatable, dimension(:, :, :) :: Ex, Ey, Ez
  real(8), allocatable, dimension(:, :, :) :: Hx, Hy, Hz
  real(8), allocatable, dimension(:, :, :) :: EPSINV

  integer :: NCYCMAX = 0
  real(8) :: GT  = 0.0
  real(8) :: DT  = 0.577

  integer :: PARTITIONS

contains

  subroutine InitGrid()

    implicit none

    ! Data is being read from file FNCFG (defined in constant)
    ! Dataflag: #GRID

    ! Variables
    integer :: error = 0  
    integer :: ios
    character(len=STRLNG) :: str
    integer :: err

    ! Inner ranges are [IBEG,IEND] etc., outer ranges [IMIN,IMAX]
    ! these are calculated in config.f90/ReadGrid

    allocate(Ex(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of Ex failed."
       stop
    endif
    allocate(Ey(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of Ey failed."
       stop
    endif
    allocate(Ez(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of Ez failed." 
       stop
    endif

    allocate(Hx(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of Hx failed."
       stop
    endif
    allocate(Hy(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of Hy failed."
       stop
    endif
    allocate(Hz(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of Hz failed." 
       stop
    endif

    allocate(EPSINV(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
    if(err .ne. 0) then
       write(STDERR,*) "Allocation of EPSINV failed." 
       stop
    endif

    Ex = 0.0
    Ey = 0.0
    Ez = 0.0
 
    Hx = 0.0
    Hy = 0.0
    Hz = 0.0

    EPSINV = 1.0

  end subroutine InitGrid

  subroutine EchoGrid
    
    ! Bildschirmausgabe der Grid Daten
    
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


  subroutine FinaliseGrid ()
    implicit none

    deallocate(Hz)
    deallocate(Hy)
    deallocate(Hx)
    deallocate(Ez)
    deallocate(Ey)
    deallocate(Ex)
    deallocate(EPSINV)

  end subroutine FinaliseGrid


  subroutine ReadEpsilon

    implicit none

    integer :: ios, i, j, k
    character(len=255) :: file
    real(8) :: val

    file = cat2(pfxepsilon,mpi_sfxin)

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
  

  end subroutine ReadEpsilon

end module grid







