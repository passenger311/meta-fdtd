!-*- F90 -*------------------------------------------------------------
!
!  program: max3d
!
!  main program
!
!----------------------------------------------------------------------

#ifndef MPI
#undef MPE_LOG
#endif


program max3d

  use constant
  use strings
  use mpiworld
  use grid
  use output
  use fdtd
  use pec
  use upml
  use tfsf
  use fdtd_output
  use src_point

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

  call InitializeMPIWorld

  ! --- init grid

#if MPE_LOG

  if (myrank .eq. 0) then

     write(6,*) '* --------- max3d start'

  ! --- init logging

     call MPE_INIT_LOG() 

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

  do l=0, numproc-1  ! sort output

#if MPI
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

     if (myrank .eq. l) then

        write(6,*) '* initialising modules: ',ranklbl

        write(6,*) '* -> InitializeGrid'
        call InitializeGrid(mpi_sfxin)
        write(6,*) '* -> InitializeFdtd'
        call InitializeFdtd(mpi_sfxin)  
        write(6,*) '* -> InitializeUPML'
        call InitializeUPML
        write(6,*) '* -> InitializeTFSF'
        call InitializeTFSF
        write(6,*) '* -> InitializeOutAsc'
        call InitializeOutAsc
        write(6,*) '* -> InitializeFdtdOutput'
        call InitializeFdtdOutput
        write(6,*) '* -> InitializePzDFT '
        call InitializePzDFT
        write(6,*) '* -> InitializeSource '
        call InitializeSource

     end if

  end do

  mpisize=(IMAX-IMIN+1)*(JMAX-JMIN+1)  ! set mpi mpiet size
 
  do l=0, numproc-1  ! sort output

#if MPI
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

     if (myrank .eq. l) then
                
        write(6,*) '* process entering loop, rank = ',ranklbl 


        str = i2str(NCYCMAX)
        write(6,*) '* time steps = ', str 

     end if
     
  end do
  

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
     call MPE_LOG_EVENT(9,"StepUPMLH",mpierr) 
#endif /* MPE_LOG */

     call StepUPMLH()

#if MPE_LOG
     call MPE_LOG_EVENT(10,"StepUPMLH",mpierr)
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
     call MPE_LOG_EVENT(13,"StepUPMLE",mpierr)
#endif /* MPE_LOG */

     call StepUPMLE() 

#if MPE_LOG
     call MPE_LOG_EVENT(14,"StepUPMLE",mpierr)
#endif /* MPE_LOG */

! --------------------------- TFSFE -----------------------------------

!     call NewTFSF_E()

! ---------------------------- PEC ------------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(15,"PEC",mpierr) 
#endif /* MPE_LOG */

     call SetPEC()

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

  do l=0, numproc-1  ! sort output

#if MPI
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

     if (myrank .eq. l) then
 
        write(6,*) '* process finalising, rank = ', ranklbl

        write(6,*) '* -> OutPzDFTn'
        call OutPzDFT
        write(6,*) '* -> FinalizeUPML'
        call FinalizeUPML
        write(6,*) '* -> FinalizeFdtd'
        call FinalizeFdtd
        write(6,*) '* -> FinalizeOutput'
        call FinalizeOutAsc
        write(6,*) '* -> FinalizeFdtdOutput'
        call FinalizeFdtdOutput
        write(6,*) '* -> FinalizeGrid'
        call FinalizeGrid

#if MPI
        write(6,*) ranklbl, "* time for Comms H = ", time_waited(1)
        write(6,*) ranklbl, "* time for Comms E = ", time_waited(2)
#endif /* MPI */

     end if

  end do

! --- terminate logging

#if MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

#if MPE_LOG
  call MPE_FINISH_LOG("max3d.log") !
#endif /* MPE_LOG */

! --- terminate MPI

  call FinalizeMPIWorld

  if (myrank .eq. 0) then

     write(6,*) '* --------- max3d end'

  end if
  
end program max3d


!----------------------------------------------------------------------
!----------------------------------------------------------------------
