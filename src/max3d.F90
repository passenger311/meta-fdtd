!----------------------------------------------------------------------
!
!  program: max3d / max3d
!
!  main program.
!
!----------------------------------------------------------------------

#ifndef MPI
#undef MPE_LOG
#endif


program max3d

#if MPI
  use mpistart
#endif /* MPI */
  use config
  use constant
  use grid
  use fdtd3d
  use upml
  use bound
  use output
  use tfsf
  use pointsource

  implicit none

  integer :: ncyc, error, l, ides, isrc
  character(len=12) :: str

  ! --- comms time measurement

#if MPI
  double precision :: tt
  double precision :: time_waited (2)
  time_waited = 0
#endif /* MPI */

! =========================== INITIALIZATION ==========================

  ! --- startup MPI

#if MPI
     call MPIinit
#endif /* MPI */

  ! --- init grid

#if MPE_LOG

  if (myrank .eq. 0) then

     write(6,*) '* --------- max3d start'

  ! --- init logging

     call MPE_INIT_LOG() !*

     call MPE_Describe_state(1,2,"ISEND","red:gray0")
     call MPE_Describe_state(3,4,"IRECV","blue:gray1")
     call MPE_Describe_state(5,6,"WAIT","green:gray2")
     call MPE_Describe_state(7,8,"StepH","yellow:gray3")
     call MPE_Describe_state(9,10,"PMLH","violet:gray4")
     call MPE_Describe_state(11,12,"StepE","orange1:gray5")
     call MPE_Describe_state(13,14,"PMLE","snow:gray6")
     call MPE_Describe_state(15,16,"PEC","pink1:gray7")
     call MPE_Describe_state(17,18,"PzDFT","chocolate:gray8")

  end if
#endif /* MPE_LOG */

#if MPI
  do l=0, numproc-1  ! sort output

     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

     if (myrank .eq. l) then

        write(6,*) '* initialising for process = ',ranklbl
#endif /* MPI */

        write(6,*) '* -> ReadConfig'
        call ReadConfig(mpi_sfxin)        ! read config files
        write(6,*) '* -> InitGrid'
        call InitGrid()                   ! initialize grid, allocate fields
        write(6,*) '* -> ReadEpsilon'
        call ReadEpsilon
        write(6,*) '* -> InitPML'
        call InitPML()
        write(6,*) '* -> InitTFSF '
        call InitTFSF()
        write(6,*) '* -> InitOutput '
        call InitOutput()
        write(6,*) '* -> InitPzDFT '
        call InitPzDFT()
        write(6,*) '* -> InitSource '
        call InitSource()

#if MPI
     end if

  end do
#endif /* MPI */

#if MPI

  mpisize=(IMAX-IMIN+1)*(JMAX-JMIN+1)  ! set mpi mpiet size
 
  do l=0, numproc-1  ! sort output

     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

     if (myrank .eq. l) then
                
        write(6,*) '* process entering loop, rank = ',ranklbl 

#endif /* MPI */

        str = i2str(NCYCMAX)
        write(6,*) '* time steps = ', str 

#if MPI
        
     end if
     
  end do

#endif /* MPI */
  

! ============================ TIMELOOP ===============================

  ncyc = 0
  do ncyc = 1, NCYCMAX

     if (myrank .eq. 0 .and. mod(ncyc*10000/NCYCMAX,100) .eq. 0  ) then

        str = cat2(i2str(ncyc*100/NCYCMAX),'%')
        write(6,*) 'running ... ', str

     end if

! ------------------------------ StepH --------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(7,"StepH",mpierr)
#endif /* MPE_LOG */

     call StepH()

#if MPE_LOG
     call MPE_LOG_EVENT(8,"StepH",mpierr)
#endif /* MPE_LOG */

! ----------------------------- SourceHy ------------------------------

!     call SourceHy(ncyc)

! ------------------------------ PMLH ---------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(9,"PMLH",mpierr) 
#endif /* MPE_LOG */

     call PMLH()

#if MPE_LOG
     call MPE_LOG_EVENT(10,"PMLH",mpierr)
#endif /* MPE_LOG */

! ---------------------------- TFSFH -------------------------------

!     call NewTFSF_H()

! -------------------------- COMMS H -------------------------------

! 1. send tangential H fields to right

#if MPI

     if ( myrank .ne. numproc-1 ) then
        ides= myrank+1

#if MPE_LOG
        call MPE_LOG_EVENT(1,"ISEND",mpierr)
#endif /* MPE_LOG */

        call MPI_ISEND( Hx(IMIN,JMIN,KEND),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Hy(IMIN,JMIN,KEND),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )
        
#if MPE_LOG
        call MPE_LOG_EVENT(2,"ISEND",mpierr)
#endif /* MPE_LOG */

     endif

! 2. receive tangential H fields from left

     if ( myrank .ne. 0 ) then
        isrc= myrank-1

#if MPE_LOG
        call MPE_LOG_EVENT(3,"IRECV",mpierr) 
#endif /* MPE_LOG */

        call MPI_IRECV( Hx(IMIN,JMIN,KMIN),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Hy(IMIN,JMIN,KMIN),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )

#if MPE_LOG
        call MPE_LOG_EVENT(4,"IRECV",mpierr) 
#endif /* MPE_LOG */

     endif

     tt = MPI_WTIME()

#endif /* MPI */

! ---------------------------- WAIT H ---------------------------------

! wait for pending comms to finish

#if MPI

     if ( myrank .ne. 0 ) then

#if MPE_LOG
        call MPE_LOG_EVENT(5,"WAIT",mpierr)
#endif /* MPE_LOG */

        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )

#if MPE_LOG
        call MPE_LOG_EVENT(6,"WAIT",mpierr)
#endif /* MPE_LOG */

     endif

     if ( myrank .ne. numproc-1 ) then

#if MPE_LOG
        call MPE_LOG_EVENT(5,"WAIT",mpierr)
#endif /* MPE_LOG */

        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )
 
#if MPE_LOG
        call MPE_LOG_EVENT(6,"WAIT",mpierr)
#endif /* MPE_LOG */

     endif

     tt = MPI_WTIME() - tt
     time_waited(1) = time_waited(1) + tt

#endif /* MPI */

! ---------------------------- OUTPUT ---------------------------------

     call DataPrepOutput(ncyc)

! ---------------------------- StepE ----------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(11,"StepE",mpierr)
#endif /* MPE_LOG */

     call StepE()

#if MPE_LOG
     call MPE_LOG_EVENT(12,"StepE",mpierr)
#endif /* MPE_LOG */

! -------------------------- SourceEy ---------------------------------

     call SourceEy(ncyc)

! ---------------------------- PMLE -----------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(13,"PMLE",mpierr)
#endif /* MPE_LOG */

     call PMLE() 

#if MPE_LOG
     call MPE_LOG_EVENT(14,"PMLE",mpierr)
#endif /* MPE_LOG */

! --------------------------- TFSFE -----------------------------------

!     call NewTFSF_E()

! ---------------------------- PEC ------------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(15,"PEC",mpierr) 
#endif /* MPE_LOG */

     call PEC()

#if MPE_LOG
     call MPE_LOG_EVENT(16,"PEC",mpierr) 
#endif /* MPE_LOG */

! -------------------------- COMMS E ----------------------------------

#if MPI 

! 1. send tangential E fields to the left

     if ( myrank .ne. 0 ) then
        ides = myrank-1

#if MPE_LOG
        call MPE_LOG_EVENT(1,"ISEND",mpierr)
#endif /* MPE_LOG */

        call MPI_ISEND( Ex(IMIN,JMIN,KBEG),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Ey(IMIN,JMIN,KBEG),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )

#if MPE_LOG
        call MPE_LOG_EVENT(2,"ISEND",mpierr)
#endif /* MPE_LOG */

     endif

! 2. receive tangential E fields from right

     if ( myrank .ne. numproc-1 ) then
        isrc = myrank+1

#if MPE_LOG
        call MPE_LOG_EVENT(3,"IRECV",mpierr)
#endif /* MPE_LOG */

        call MPI_IRECV( Ex(IMIN,JMIN,KMAX),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Ey(IMIN,JMIN,KMAX),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )

#if MPE_LOG
        call MPE_LOG_EVENT(4,"IRECV",mpierr)
#endif /* MPE_LOG */

     endif

     tt = MPI_WTIME()

#endif /* MPI */

! --------------------------- WAIT E ----------------------------------

#if MPI

     if ( myrank .ne. numproc-1 ) then

#if MPE_LOG
        call MPE_LOG_EVENT(5,"WAIT",mpierr) 
#endif /* MPE_LOG */

        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )

#if MPE_LOG
        call MPE_LOG_EVENT(6,"WAIT",mpierr) 
#endif /* MPE_LOG */

     endif


     if ( myrank .ne. 0 ) then

#if MPE_LOG
        call MPE_LOG_EVENT(5,"WAIT",mpierr) 
#endif /* MPE_LOG */

        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )

#if MPE_LOG
        call MPE_LOG_EVENT(6,"WAIT",mpierr) 
#endif /* MPE_LOG */

     endif

     tt = MPI_WTIME() - tt
     time_waited(2) = time_waited(2) + tt

#endif /* MPI */
    
! ----------------------------- OUTPUT --------------------------------

     call DataOutGPL(ncyc)

! ------------------------------ PzDFT --------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(17,"PzDFT",mpierr)  
#endif /* MPE_LOG */

     call PzDFT(ncyc)

#if MPE_LOG
     call MPE_LOG_EVENT(18,"PzDFT",mpierr)  
#endif /* MPE_LOG */

! ---------------------------------------------------------------------

     GT = GT + DT

  end do

! =========================== FINALIZING ==============================

! end of time loop.

#if MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

#if MPI
  do l=0, numproc-1  ! sort output

     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

     if (myrank .eq. l) then
 
        write(6,*) '* process finalising, rank = ', ranklbl
#endif /* MPI */

        write(6,*) '* -> OutPzDFTn'
        call OutPzDFT()
        write(6,*) '* -> FinaliseOutput'
        call FinaliseOutput()
        write(6,*) '* -> FinalisePML'
        call FinalisePML()
        write(6,*) '* -> FinaliseGrid'
        call FinaliseGrid()

#if MPI
        write(6,*) ranklbl, "* time for Comms H = ", time_waited(1)
        write(6,*) ranklbl, "* time for Comms E = ", time_waited(2)
     end if

  end do
#endif /* MPI */

! --- terminate logging

#if MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

#if MPE_LOG
  call MPE_Finish_log("max3d.log") !
#endif /* MPE_LOG */

! --- terminate MPI

#if MPI
  call MPIfinalise

  if (myrank .eq. 0) then
#endif /* MPI */

     write(6,*) '* --------- max3d end'

#if MPI
  end if
#endif /* MPI */
  
end program max3d


!----------------------------------------------------------------------
!----------------------------------------------------------------------
