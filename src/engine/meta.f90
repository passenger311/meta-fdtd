!-*- F90 -*------------------------------------------------------------
!
!  program: meta / meta
!
!  main program
!
!----------------------------------------------------------------------

!======================================================================
!
!

program meta

  use strings
  use list
  use mpiworld
  use grid
  use fdtd
  use fdtd_calc
  use bound
  use mat
  use out
  use diag
  use config
  use mpicomms
  use manifest
  use timer

  implicit none

  character(len=STRLNG), parameter :: modname = "META"
  integer :: ncyc, l , prog0, prog1
  real(kind=8) :: cells
  character(len=12) :: str
  character(len=1) :: busychar ='|'


! =========================== INITIALIZATION ==========================


  ! --- startup MPI

  call InitializeMPIWorld

  if (myrank .eq. 0) then

     write(6,*) "* ------------------------- META-ENGINE STARTS -------------------------- "
     write(6,*) "*"
     call DisplayManifest
     write(6,*) "*"
     write(6,*) "* ----------------------------------------------------------------------- "
  	 M4_WRITE_INFO("START")

  end if

  do l=0, numproc-1

M4_IFELSE_MPI(call SynchronizeMPIWorld)
  
     if (myrank .eq. l) then

        write(6,*) "*"
        write(6,*) "* initialising modules: myrank = ", ranklbl
        call CheckTimer
        write(6,*) "* -> InitializeList"
        call InitializeList
        write(6,*) "* -> ReadConfig(M4_SDIM{})"
        call ReadConfig(M4_SDIM{})
        write(6,*) "* -> InitializeFdtd"
        call InitializeFdtd
        write(6,*) "* -> InitializeFdtdCalc"
        call InitializeFdtdCalc
        write(6,*) "* -> InitializeBound"
        call InitializeBound
        write(6,*) "* -> InitializeMat"
        call InitializeMat
        write(6,*) "* -> InitializeDiag"
        call InitializeDiag
        write(6,*) "* -> InitializeOut"
        call InitializeOut
       
     end if

  end do
       
  if (myrank .eq. 0) then

     write(6,*) "*"
     write(6,*) "* ----------------------------------------------------------------------- "
  	 M4_WRITE_INFO("LOOP")

  end if

 
  do l=0, numproc-1

M4_IFELSE_MPI(call SynchronizeMPIWorld)

     if (myrank .eq. l) then
        write(6,*) "* "
        write(6,*) "* process entering loop: myrank = ",ranklbl 

        str = i2str(NCYCMAX)
        write(6,*) "* time steps = ", str 

     end if
     
  end do

M4_IFELSE_MPI(call InitializeMPIComms)

! ============================ TIMELOOP ===============================

  do l = 1, MAXTIMER
     call ResetTimer(l) ! start timer
  end do

  ncyc = 0
  prog0 = 0
  prog1 = -1

  do ncyc = 0, NCYCMAX

! do some pretty print, a cursor and a 100.0% progress indicator
     
     !  progress 0...1000
     prog0 = int(ncyc*1000./NCYCMAX) 

     if (myrank .eq. 0 .and. prog0 .gt. prog1 ) then

        prog1 = prog0

        str = cat4(i2str(prog0/10),".",i2str(mod(prog0,10)),"%")
        select case ( busychar )
        case ( '|' ) 
           busychar = '/'
        case ( '/' ) 
           busychar = '-'
        case ( '-' ) 
M4_IFELSE_ARCH({PGI}, {
           busychar = '\\'
        case ( '\\' ) 
}, {
           busychar = '\'
        case ( '\' ) 
})
           busychar = '|'
        end select
        write(6,*) busychar," running ... ", str

     end if

! ------------------------------ StepH --------------------------------

call StartTimer(2)
M4_MPE_SECTION(4,{
	if ( ncyc .gt. 0 ) call StepH
})
call StopTimer(2)


! -------------------------- StepHMat ---------------------------------

call StartTimer(3)
M4_MPE_SECTION(8,{
	if ( ncyc .gt. 0 ) call StepHMat(ncyc)
})
call StopTimer(3)

! -------------------------- StepHBound -------------------------------

call StartTimer(4)
M4_MPE_SECTION(7,{
	if ( ncyc .gt. 0 )  call StepHBound
})
call StopTimer(4)

! ---------------------------- COMMS H --------------------------------

! H field calculation if finished -> initiate comms

call StartTimer(10)
M4_IFELSE_MPI({
	call InitiateHMPIComms
})
call StartTimer(10)

! -------------------------- StepHDiag --------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(5)
M4_MPE_SECTION(10,{
	call StepHDiag(ncyc)
})
call StopTimer(5)

! ------------------------ LoadDataOut --------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(9)
     call WriteDataOut(ncyc,.false.) 
call StopTimer(9)

! ---------------------------- WAIT H ---------------------------------

! complete pending H comms

call StartTimer(10)
M4_IFELSE_MPI({
	call CompleteHMPIComms
})
call StopTimer(10)

! ---------------------------- StepE ----------------------------------

call StartTimer(2)
M4_MPE_SECTION(6,{ 
	if ( ncyc .gt. 0 ) call StepE 
})
call StopTimer(2)

! -------------------------- StepEMat ---------------------------------

call StartTimer(3)
M4_MPE_SECTION(9,{
	if ( ncyc .gt. 0 ) call StepEMat(ncyc)
})
call StopTimer(3)

! -------------------------- StepEBound -------------------------------

call StartTimer(4)
M4_MPE_SECTION(7,{
	if ( ncyc .gt. 0 ) call StepEBound
})
call StopTimer(4)


! -------------------------- COMMS E ----------------------------------

! E field calculation is finished -> initiate comms

call StartTimer(10)
M4_IFELSE_MPI({
	call InitiateEMPIComms
})
call StopTimer(10)
	
! -------------------------- StepEDiag --------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(5)
M4_MPE_SECTION(11,{
	call StepEDiag(ncyc)
})
call StopTimer(5)

! ------------------------ WriteDataOut -------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(9)
     call WriteDataOut(ncyc, .true.) 
call StopTimer(9)

! --------------------------- WAIT E ----------------------------------

! complete pending E comms

call StartTimer(10)
M4_IFELSE_MPI({
	call CompleteEMPIComms
})
call StopTimer(10)

! ----------------------- END OF TIME LOOP ----------------------------

     GT = GT + DT

  end do

  call StopTimer(1)

! ---------------------------------------------------------------------

M4_IFELSE_MPI(call FinalizeMPIComms)

  
! =========================== FINALIZING ==============================

! end of time loop.

  if (myrank .eq. 0) then

     write(6,*) "*"
     write(6,*) "* ----------------------------------------------------------------------- "
     M4_WRITE_INFO("TIMING")
     write(6,*) "*"
     call DisplayTotalTimer("Fdtd   : ",2)
     call DisplayTotalTimer("Mat    : ",3)
     call DisplayTotalTimer("Bound  : ",4)
     call DisplayTotalTimer("Diag   : ",5)
     call DisplayTotalTimer("Output : ",9)
     call DisplayTotalTimer("Comms  : ",10)
     call DisplayTotalTimer("Total  : ",1)
     write(6,*) "*"
     cells = real(NCYCMAX) * real(IEIG-IBIG+1) * real(JEIG-JBIG+1) * real(KEIG-KBIG+1)
     call DisplayMcpsTimer( "Fdtd   : ",2, cells)
     cells = real(NCYCMAX) * real(IBIG-IBEG+1) * real( JBIG-JBEG+1) * real(KBIG-KBEG+1)
     call DisplayMcpsTimer( "Bound  : ",4, cells)

     write(6,*) "*"
     write(6,*) "* ----------------------------------------------------------------------- "
     M4_WRITE_INFO("FINALIZE")
     write(6,*) "*"

  end if

M4_IFELSE_MPI(call SynchronizeMPIWorld)

  do l=0, numproc-1

M4_IFELSE_MPI(call SynchronizeMPIWorld)

     if (myrank .eq. l) then
 
        write(6,*) "* process finalising, rank = ", ranklbl

        write(6,*) "* -> FinalizeOut"
        call FinalizeOut
        write(6,*) "* -> FinalizeMat"
        call FinalizeMat
        write(6,*) "* -> FinalizeDiag"
        call FinalizeDiag
        write(6,*) "* -> FinalizeBound"
        call FinalizeBound
        write(6,*) "* -> FinalizeFdtd"
        call FinalizeFdtd
        write(6,*) "* -> FinalizeGrid"
        call FinalizeGrid
        write(6,*) "* -> FinalizeList"
        call FinalizeList

        write(6,*) "*"

     end if

  end do

! --- terminate MPI

  call FinalizeMPIWorld

  if (myrank .eq. 0) then

  	 M4_WRITE_INFO("END")
     write(6,*) "* ------------------------- META-ENGINE STOPPED ------------------------- "

  end if
  
end program meta

! ---------------------------------------------------------------------

!
! Authors:  J.Hamm, A.Klaedtke, S.Scholz
! Modified: 6/1/2008
!
!======================================================================
