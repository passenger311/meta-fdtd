!-*- F90 -*------------------------------------------------------------
!
!  program: meta3 / meta3
!
!  main program
!
!----------------------------------------------------------------------

!======================================================================
!
! Compile with preprocessor 
! -DMPI     : to activate parallel MPI code
! -DMPE_LOG : to activate MPE logging
! 

#ifndef MPI
#undef MPE_LOG
#endif

program meta3

  use constant
  use strings
  use list
  use mpiworld
  use grid
  use outgpl
  use fdtd
  use pec
  use pml
  use mat
  use diag
  use tfsf
  use outgpl_fdtd
  use matsource

  implicit none

  integer :: ncyc, error, l, ides, isrc
  character(len=12) :: str

  ! --- comms time measurement

#if MPI
  real(kind=8) :: tt
  real(kind=8) :: time_waited (2)
  time_waited = 0
#endif /* MPI */

! =========================== INITIALIZATION ==========================

  ! --- startup MPI

  call InitializeMPIWorld

  ! --- init grid

#if MPE_LOG

  if (myrank .eq. 0) then

     write(6,*) '* --------- meta3 start'

  ! --- init logging

     call MPE_INIT_LOG() 

     call MPE_Describe_state(1,2,"ISEND","red:gray0")
     call MPE_Describe_state(3,4,"IRECV","blue:gray1")
     call MPE_Describe_state(5,6,"WAIT","green:gray2")
     call MPE_Describe_state(7,8,"StepH","yellow:gray3")
     call MPE_Describe_state(9,10,"PmlH","violet:gray4")
     call MPE_Describe_state(11,12,"StepE","orange1:gray5")
     call MPE_Describe_state(13,14,"PmlE","snow:gray6")
     call MPE_Describe_state(15,16,"StepHMat","pink1:gray7")
     call MPE_Describe_state(17,18,"StepEMat","chocolate:gray8")
     call MPE_Describe_state(19,20,"StepHDiag","pink1:gray7")
     call MPE_Describe_state(21,22,"StepEDiag","chocolate:gray8")

  end if
#endif /* MPE_LOG */

  do l=0, numproc-1

#if MPI
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

     if (myrank .eq. l) then

        write(6,*) '* initialising modules: ',ranklbl

        write(6,*) '* -> InitializeList'
        call InitializeList
        write(6,*) '* -> InitializeGrid'
        call InitializeGrid(mpi_sfxin)
        write(6,*) '* -> InitializeFdtd'
        call InitializeFdtd(mpi_sfxin)  
        write(6,*) '* -> InitializePml'
        call InitializePml
        write(6,*) '* -> InitializeMat'
        call InitializeMat
        write(6,*) '* -> InitializeOut'        !
        call InitializeOut
        write(6,*) '* -> WriteHeaderOut'
        call WriteHeaderOut

     end if

  end do

  mpisize=(IMAX-IMIN+1)*(JMAX-JMIN+1)  ! set mpisize
 
  do l=0, numproc-1

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


! ------------------------------ PMLH ---------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(9,"StepHPml",mpierr) 
#endif /* MPE_LOG */

     call StepHPml()

#if MPE_LOG
     call MPE_LOG_EVENT(10,"StepHPml",mpierr)
#endif /* MPE_LOG */


! -------------------------- StepHMat ---------------------------------


#if MPE_LOG
     call MPE_LOG_EVENT(15,"StepHMat",mpierr) 
#endif /* MPE_LOG */

     call StepHMat(ncyc)

#if MPE_LOG
     call MPE_LOG_EVENT(16,"StepHMat",mpierr) 
#endif /* MPE_LOG */


! -------------------------- COMMS H -------------------------------

! H field calculation if finished -> initiate comms

#if MPI

! 1. send tangential H fields to right

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


! -------------------------- StepHDiag --------------------------------

! use the time while comms are pending to do some useful stuff

#if MPE_LOG
     call MPE_LOG_EVENT(19,"StepHDiag",mpierr) 
#endif /* MPE_LOG */

     call StepHDiag(ncyc)

#if MPE_LOG
     call MPE_LOG_EVENT(20,"StepHDiag",mpierr) 
#endif /* MPE_LOG */


! ------------------------ LoadDataOut --------------------------------

! use the time while comms are pending to do some useful stuff


     call LoadDataOut(ncyc) ! buffer output data from this half-step


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


! ---------------------------- StepE ----------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(11,"StepE",mpierr)
#endif /* MPE_LOG */

     call StepE()

#if MPE_LOG
     call MPE_LOG_EVENT(12,"StepE",mpierr)
#endif /* MPE_LOG */


! ---------------------------- PMLE -----------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(13,"StepEPml",mpierr)
#endif /* MPE_LOG */

     call StepEPml() 

     call SetPec()

#if MPE_LOG
     call MPE_LOG_EVENT(14,"StepEPml",mpierr)
#endif /* MPE_LOG */


! -------------------------- StepEMat ---------------------------------

#if MPE_LOG
     call MPE_LOG_EVENT(17,"StepEMat",mpierr) 
#endif /* MPE_LOG */

     call StepEMat(ncyc)

#if MPE_LOG
     call MPE_LOG_EVENT(18,"StepEMat",mpierr) 
#endif /* MPE_LOG */

! -------------------------- COMMS E ----------------------------------

! E field calculation is finished -> initiate comms

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


! -------------------------- StepEDiag --------------------------------

! use the time while comms are pending to do some useful stuff

#if MPE_LOG
     call MPE_LOG_EVENT(21,"StepEDiag",mpierr) 
#endif /* MPE_LOG */

     call StepEDiag(ncyc)

#if MPE_LOG
     call MPE_LOG_EVENT(22,"StepEDiag",mpierr) 
#endif /* MPE_LOG */

! ------------------------ WriteDataOut -------------------------------

! use the time while comms are pending to do some useful stuff

     call WriteDataOut(ncyc) ! write output data after full-step

! --------------------------- WAIT E ----------------------------------

! complete comms

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

! ---------------------------------------------------------------------

     GT = GT + DT

  end do

! =========================== FINALIZING ==============================

! end of time loop.

#if MPI
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

  do l=0, numproc-1

#if MPI
     call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif /* MPI */

     if (myrank .eq. l) then
 
        write(6,*) '* process finalising, rank = ', ranklbl

        write(6,*) '* -> FinalizeOut'
        call FinalizeOut
        write(6,*) '* -> FinalizeMat'
        call FinalizeMat
        write(6,*) '* -> FinalizePml'
        call FinalizePml
        write(6,*) '* -> FinalizeFdtd'
        call FinalizeFdtd
        write(6,*) '* -> FinalizeGrid'
        call FinalizeGrid
        write(6,*) '* -> FinalizeList'
        call FinalizeList

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
  call MPE_FINISH_LOG("meta3.log") !
#endif /* MPE_LOG */

! --- terminate MPI

  call FinalizeMPIWorld

  if (myrank .eq. 0) then

     write(6,*) '* --------- meta3 end'

  end if
  
end program meta3

! ---------------------------------------------------------------------

!
! Authors:  J.Hamm, A.Klaedtke, S.Scholz
! Modified: 4/12/2007
!
!======================================================================
