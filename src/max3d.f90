!----------------------------------------------------------------------
!
!  program: max3d / max3d
!
!  main program.
! 
!
!
!----------------------------------------------------------------------


program max3d

  use mpistart
  use config
  use constant
  use grid
  use fdtd3d
  use upml
  use bound
  use output
  use tfsf

  implicit none

  integer :: ncyc, error, l, ides, isrc
  character(len=12) :: str

  ! --- time measurment

  double precision :: tt
  double precision :: time_waited (2)
  double precision :: mysecond
  external mysecond

  time_waited = 0
  
  ! --- startup MPI

  call MPIinit

  ! --- init grid

  if (myrank .eq. 0) then

     write(6,*) '* --------- max3d start'

  end if

  do l=0, numproc-1  ! sort output

     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

     if (myrank .eq. l) then

        write(6,*) '* initialising for ',ranklbl
        write(6,*) '* ... ReadConfig'
        call ReadConfig(mpi_sfxin)        ! read config files
        write(6,*) '* ... InitGrid'
        call InitGrid()                   ! initialize grid, allocate fields
        write(6,*) '* ... ReadEpsilon'
        call ReadEpsilon
        write(6,*) '* ... InitPML'
        call InitPML()
        write(6,*) '* ... InitTFSF '
        call InitTFSF()
        write(6,*) '* ... InitOutput '
        call InitOutput()
        write(6,*) '* ... InitPzDFT '
        call InitPzDFT()

     end if

  end do

  mpisize=(IMAX-IMIN+1)*(JMAX-JMIN+1)  ! set mpi mpiet size
 
  do l=0, numproc-1  ! sort output

     if (myrank .eq. l) then
        
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
        
        str = i2str(NCYCMAX)
        write(6,*) '* entering time loop ',ranklbl 
        write(6,*) '* ... steps = ', str 
        
     end if
     
  end do
  
!  call MSGexchangeE
!  call MSGexchangeH
  

  ncyc = 0

  do ncyc = 1, NCYCMAX

     if (myrank .eq. 0 .and. mod(ncyc*10000/NCYCMAX,100) .eq. 0  ) then

        str = cat2(i2str(ncyc*100/NCYCMAX),'%')
        write(6,*) 'running ... ', str

     end if

     call StepH()
     call PMLH()
!     call NewTFSF_H()

! ** send tangential H fields to right

     if ( myrank .ne. numproc-1 ) then
        ides= myrank+1
        call MPI_ISEND( Hx(IMIN,JMIN,KEND),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Hy(IMIN,JMIN,KEND),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )
     endif

! receive tangential H fields from left

     if ( myrank .ne. 0 ) then
        isrc= myrank-1
        call MPI_IRECV( Hx(IMIN,JMIN,KMIN),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Hy(IMIN,JMIN,KMIN),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )
     endif

     tt = mysecond()

! and wait ...

     if ( myrank .ne. 0 ) then
        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )
     endif

     if ( myrank .ne. numproc-1 ) then
        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )
     endif

     tt = mysecond() - tt
     time_waited(1) = time_waited(1) + tt

! ** finished tangential H comm


!     call DataPrepOutput(ncyc)
     call StepE()
     call PMLE() 
!     call NewTFSF_E()
     call PEC()


! ** send tangential E fields to left

     if ( myrank .ne. 0 ) then
        ides = myrank-1
        call MPI_ISEND( Ex(IMIN,JMIN,KBEG),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Ey(IMIN,JMIN,KBEG),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )
     endif

! receive tangential E fields from right

     if ( myrank .ne. numproc-1 ) then
        isrc = myrank+1
        call MPI_IRECV( Ex(IMIN,JMIN,KMAX),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Ey(IMIN,JMIN,KMAX),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )
     endif

     tt = mysecond()

! and wait ...

     if ( myrank .ne. numproc-1 ) then
        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )
     endif


     if ( myrank .ne. 0 ) then
        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )
     endif

     tt = mysecond() - tt
     time_waited(2) = time_waited(2) + tt

! ** finished tangential H comm
    

!     call DataOutGPL(ncyc)
     call PzDFT(ncyc)

     GT = GT + DT

  end do

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  do l=0, numproc-1  ! sort output

     if (myrank .eq. l) then
 
        write(6,*) '* finalising ', ranklbl

        call OutPzDFT()
        call FinaliseOutput()
        call FinalisePML()
        call FinaliseGrid()

        write(6,*) ranklbl, "Time spent in H-field communication", time_waited(1)
        write(6,*) ranklbl, "Time spent in E-field communication", time_waited(2)
     end if

  end do

! --- terminate MPI

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  call MPIfinalise


  if (myrank .eq. 0) then

     write(6,*) '* --------- max3d end'

  end if


contains
  ! From the stream benchmark as the second_wall.c function mysecond

  ! A semi-portable way to determine the clock granularity
  ! Adapted from a code by John Henning of Digital Equipment Corporation

  INTEGER FUNCTION checktick()
    IMPLICIT NONE
                                                                               
    !     .. Parameters ..
    integer, parameter :: n = 20

    !     .. Local Scalars ..
    DOUBLE PRECISION :: t1, t2
    INTEGER :: i,j,jmin

    !     .. Local Arrays ..
    DOUBLE PRECISION :: timesfound(n)

    !     .. External Functions ..
    DOUBLE PRECISION :: mysecond
    EXTERNAL mysecond

    !     .. Intrinsic Functions ..
    INTRINSIC :: max,min,nint

    do i = 0, n-1
       t2 = mysecond()
       
       do while (t2.ne.t1)
          t2 = mysecond()
       end do
       
       t1 = t2
       timesfound(i) = t1
       
    end do
    
    jmin = 1000000
    do i = 2,n
       j = nint((timesfound(i)-timesfound(i-1))*1d6)
       jmin = min(jmin,max(j,0))
    end do
       
    if (jmin.GT.0) then
       checktick = jmin
    else
       write(6,*) 'Your clock granularity appears to be less ', \
       'than one microsecond'
       checktick = 1
    end if
    return
    
  end function checktick

  
end program max3d
