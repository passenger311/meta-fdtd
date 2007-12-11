!-*- F90 -*------------------------------------------------------------
!
!  program: meta3 / meta3
!
!  main program
!
!----------------------------------------------------------------------

!======================================================================
!
!

program meta3

  use constant
  use strings
  use list
  use mpiworld
  use grid
  use fdtd
  use pec
  use pml
  use mat
  use out
  use diag
  use config

  implicit none


  character(len=STRLNG), parameter :: modname = "META3"
  integer :: ncyc, error, l, ides, isrc
  character(len=12) :: str

  ! --- comms time measurement

M4_IFELSE_MPI({
  real(kind=8) :: tt
  real(kind=8) :: time_waited (2)
  time_waited = 0
},{})

! =========================== INITIALIZATION ==========================


  ! --- startup MPI

  call InitializeMPIWorld

  ! --- init grid

M4_IFELSE_MPELOG({

  if (myrank .eq. 0) then

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
})

  if (myrank .eq. 0) then

     write(6,*) "* ------------------------ META-3 ENGINE STARTS  ------------------------ "
     write(6,*) "*"
     call DisplayVersionLine
     write(6,*) "*"
  end if

  do l=0, numproc-1

M4_IFELSE_MPI({call MPI_BARRIER(MPI_COMM_WORLD,mpierr)},{})

     if (myrank .eq. l) then

        write(6,*) "* initialising modules: myrank = ", ranklbl

        write(6,*) "* -> InitializeList"
        call InitializeList
        write(6,*) "* -> ReadConfig"
        call ReadConfig
        write(6,*) "* -> InitializeFdtd"
        call InitializeFdtd
        write(6,*) "* -> InitializePml"
        call InitializePml
        write(6,*) "* -> InitializeMat"
        call InitializeMat
        write(6,*) "* -> InitializeDiag"
        call InitializeDiag
        write(6,*) "* -> InitializeOut"
        call InitializeOut
       
     end if

  end do

  mpisize=(IMAX-IMIN+1)*(JMAX-JMIN+1)  ! set mpisize
 
  do l=0, numproc-1

M4_IFELSE_MPI({call MPI_BARRIER(MPI_COMM_WORLD,mpierr)})

     if (myrank .eq. l) then
                
        write(6,*) "* process entering loop: myrank = ",ranklbl 


        str = i2str(NCYCMAX)
        write(6,*) "* time steps = ", str 

     end if
     
  end do
  

! ============================ TIMELOOP ===============================

  ncyc = 0
  do ncyc = 0, NCYCMAX

     if (myrank .eq. 0 .and. mod(ncyc*10000/NCYCMAX,100) .eq. 0  ) then

        str = cat2(i2str(ncyc*100/NCYCMAX),"%")
        write(6,*) "running ... ", str

     end if

! ------------------------------ StepH --------------------------------

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(7,"StepH",mpierr)})

     if ( ncyc .gt. 0 ) call StepH()

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(8,"StepH",mpierr)})

! ------------------------------ PMLH ---------------------------------

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(9,"StepHPml",mpierr)})

     if ( ncyc .gt. 0 )  call StepHPml()

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(10,"StepHPml",mpierr)})

! -------------------------- StepHMat ---------------------------------

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(15,"StepHMat",mpierr)})

     if ( ncyc .gt. 0 ) call StepHMat(ncyc)

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(16,"StepHMat",mpierr)})

! -------------------------- COMMS H -------------------------------

! H field calculation if finished -> initiate comms

M4_IFELSE_MPI({

! 1. send tangential H fields to right

     if ( myrank .ne. numproc-1 ) then
        ides= myrank+1

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(1,"ISEND",mpierr)})

        call MPI_ISEND( Hx(IMIN,JMIN,KEND),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Hy(IMIN,JMIN,KEND),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )
        
M4_IFELSE_MPELOG({call MPE_LOG_EVENT(2,"ISEND",mpierr)})

     endif

! 2. receive tangential H fields from left

     if ( myrank .ne. 0 ) then
        isrc= myrank-1

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(3,"IRECV",mpierr)})

        call MPI_IRECV( Hx(IMIN,JMIN,KMIN),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Hy(IMIN,JMIN,KMIN),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(4,"IRECV",mpierr)})

     endif

     tt = MPI_WTIME()

})


! -------------------------- StepHDiag --------------------------------

! use the time while comms are pending to do some useful stuff

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(19,"StepHDiag",mpierr)})

     call StepHDiag(ncyc)

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(20,"StepHDiag",mpierr)})

! ------------------------ LoadDataOut --------------------------------

! use the time while comms are pending to do some useful stuff

     call WriteDataOut(ncyc,.false.) 

! ---------------------------- WAIT H ---------------------------------

! wait for pending comms to finish

M4_IFELSE_MPI({

     if ( myrank .ne. 0 ) then

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(5,"WAIT",mpierr)})

        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(6,"WAIT",mpierr)})

     endif

     if ( myrank .ne. numproc-1 ) then

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(5,"WAIT",mpierr)})

        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )
 
M4_IFELSE_MPELOG({call MPE_LOG_EVENT(6,"WAIT",mpierr)})

     endif

     tt = MPI_WTIME() - tt
     time_waited(1) = time_waited(1) + tt

})


! ---------------------------- StepE ----------------------------------

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(11,"StepE",mpierr)})

      if ( ncyc .gt. 0 ) call StepE()

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(12,"StepE",mpierr)})

! ---------------------------- PMLE -----------------------------------

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(13,"StepEPml",mpierr)})

      if ( ncyc .gt. 0 ) call StepEPml() 
      if ( ncyc .gt. 0 ) call SetPec()

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(14,"StepEPml",mpierr)})

! -------------------------- StepEMat ---------------------------------

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(17,"StepEMat",mpierr)})

      if ( ncyc .gt. 0 ) call StepEMat(ncyc)

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(18,"StepEMat",mpierr)})

! -------------------------- COMMS E ----------------------------------

! E field calculation is finished -> initiate comms

M4_IFELSE_MPI({ 

! 1. send tangential E fields to the left

     if ( myrank .ne. 0 ) then
        ides = myrank-1

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(1,"ISEND",mpierr)})

        call MPI_ISEND( Ex(IMIN,JMIN,KBEG),mpisize,mpitype,ides, &
             xtag,mpicomm,ireqs(1),mpierr )
        call MPI_ISEND( Ey(IMIN,JMIN,KBEG),mpisize,mpitype,ides, &
             ytag,mpicomm,ireqs(2),mpierr )

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(2,"ISEND",mpierr)})

     endif

! 2. receive tangential E fields from right

     if ( myrank .ne. numproc-1 ) then
        isrc = myrank+1

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(3,"IRECV",mpierr)})

        call MPI_IRECV( Ex(IMIN,JMIN,KMAX),mpisize,mpitype,isrc, &
             xtag,mpicomm,ireqr(1),mpierr )
        call MPI_IRECV( Ey(IMIN,JMIN,KMAX),mpisize,mpitype,isrc, &
             ytag,mpicomm,ireqr(2),mpierr )

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(4,"IRECV",mpierr)})

     endif

     tt = MPI_WTIME()

})


! -------------------------- StepEDiag --------------------------------

! use the time while comms are pending to do some useful stuff

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(21,"StepEDiag",mpierr)})

     call StepEDiag(ncyc)

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(22,"StepEDiag",mpierr)})

! ------------------------ WriteDataOut -------------------------------

! use the time while comms are pending to do some useful stuff

     call WriteDataOut(ncyc, .true.) 

! --------------------------- WAIT E ----------------------------------

! complete comms

M4_IFELSE_MPI({

     if ( myrank .ne. numproc-1 ) then

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(5,"WAIT",mpierr)})

        call MPI_WAIT(ireqr(2),mpistatus, mpierr )
        call MPI_WAIT(ireqr(1),mpistatus, mpierr )

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(6,"WAIT",mpierr)})

     endif

     if ( myrank .ne. 0 ) then

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(5,"WAIT",mpierr)})

        call MPI_WAIT(ireqs(2),mpistatus, mpierr )
        call MPI_WAIT(ireqs(1),mpistatus, mpierr )

M4_IFELSE_MPELOG({call MPE_LOG_EVENT(6,"WAIT",mpierr)})

     endif

     tt = MPI_WTIME() - tt
     time_waited(2) = time_waited(2) + tt

})

! ---------------------------------------------------------------------

     GT = GT + DT

  end do

! =========================== FINALIZING ==============================

! end of time loop.

M4_IFELSE_MPI({call MPI_BARRIER(MPI_COMM_WORLD,mpierr)})

  do l=0, numproc-1

M4_IFELSE_MPI({call MPI_BARRIER(MPI_COMM_WORLD,mpierr)})

     if (myrank .eq. l) then
 
        write(6,*) "* process finalising, rank = ", ranklbl

        write(6,*) "* -> FinalizeOut"
        call FinalizeOut
        write(6,*) "* -> FinalizeMat"
        call FinalizeMat
        write(6,*) "* -> FinalizePml"
        call FinalizePml
        write(6,*) "* -> FinalizeFdtd"
        call FinalizeFdtd
        write(6,*) "* -> FinalizeGrid"
        call FinalizeGrid
        write(6,*) "* -> FinalizeList"
        call FinalizeList

M4_IFELSE_MPI({
        write(6,*) ranklbl, "* time for Comms H = ", time_waited(1)
        write(6,*) ranklbl, "* time for Comms E = ", time_waited(2)
})

     end if

  end do

! --- terminate logging

M4_IFELSE_MPI({call MPI_BARRIER(MPI_COMM_WORLD,mpierr)})

M4_IFELSE_MPELOG({call MPE_FINISH_LOG("meta3.log")})

! --- terminate MPI

  call FinalizeMPIWorld

  if (myrank .eq. 0) then

     write(6,*) "* ------------------------ META-3 ENGINE STOPPED ------------------------ "

  end if
  
end program meta3

! ---------------------------------------------------------------------

!
! Authors:  J.Hamm, A.Klaedtke, S.Scholz
! Modified: 4/12/2007
!
!======================================================================
