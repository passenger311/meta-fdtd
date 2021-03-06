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
  use checkpoint
  use list
  use mpiworld
  use grid
  use fdtd
  use fdtd_calc
  use bound
  use out
  use src
  use mat
  use diag
  use config
  use mpicomms
  use manifest
  USE lumped
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

        if ( load_state ) then
           write(6,*) "* loading checkpoint file"
           if ( load_state ) then
              write(6,*) "* checkpoint.in found (loading)"
              open(unit=UNITCHK, file=checkpoint_fn, status = "old", form = "unformatted"  )
              write(6,*) "!WRN Ensure in config file that even first timestep follows on odd last timestep &
		and vice versa (Beware that timestep ncyc=0 is not performed)."
              write(6,*) "!WRN The diagebal module evaluates mod(2*ncyc,3) and enforcing odd following on &
                even timesteps only is not sufficient for this module."
           else
              M4_FATAL_ERROR({"COULD NOT LOAD 'checkpoint.in'"})
          !    write(6,*) "* checkpoint.in not found (skipping)"
           end if
        end if

        write(6,*) "* -> InitializeFdtd"
        call InitializeFdtd
        write(6,*) "* -> InitializeFdtdCalc"
        call InitializeFdtdCalc
        write(6,*) "* -> InitializeLumped"
        call InitializeLumped
        write(6,*) "* -> InitializeBound"
        call InitializeBound
        write(6,*) "* -> InitializeSrc"
        call InitializeSrc
        write(6,*) "* -> InitializeMat"
        call InitializeMat
        write(6,*) "* -> InitializeDiag"
        call InitializeDiag
        write(6,*) "* -> InitializeOut"
        call InitializeOut

        if ( load_state ) then
           close(UNITCHK) 
        end if

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

        str = i2str(NCYCMAX-NCYCMIN)
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

  do ncyc = NCYCMIN, NCYCMAX

! do some pretty print, a cursor and a 100.0% progress indicator
     
     !  progress 0...1000
     prog0 = int((ncyc-NCYCMIN)*1000./(NCYCMAX-NCYCMIN)) 

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

!        call flush()
      end if

! ------------------------------ StepH --------------------------------

call StartTimer(2)
M4_MPE_SECTION(2,{
	if ( ncyc .gt. 0 ) call StepH
})
call StopTimer(2)

! ------------------------------ StepHLumped --------------------------------

call StartTimer(3)
M4_MPE_SECTION(3,{
	if ( ncyc .gt. 0 ) call StepHLumped(ncyc)
})
call StopTimer(3)


! -------------------------- StepHSrc ---------------------------------

call StartTimer(4)
M4_MPE_SECTION(4,{
	if ( ncyc .gt. 0 ) call StepHSrc(ncyc)
})
call StopTimer(4)

! -------------------------- StepHMat ---------------------------------

call StartTimer(5)
M4_MPE_SECTION(5,{
	if ( ncyc .gt. 0 ) call StepHMat(ncyc)
})
call StopTimer(5)

! -------------------------- StepHBound -------------------------------

call StartTimer(6)
M4_MPE_SECTION(6,{
	if ( ncyc .gt. 0 )  call StepHBound
})
call StopTimer(6)

! ---------------------------- COMMS H --------------------------------

! H field calculation if finished -> initiate comms

call StartTimer(10)
M4_IFELSE_MPI({
M4_MPE_SECTION(10,{
	call InitiateHMPIComms
})
})
call StartTimer(10)

! -------------------------- StepHDiag --------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(7)
M4_MPE_SECTION(7,{
	call StepHDiag(ncyc)
})
call StopTimer(7)

! ------------------------ LoadDataOut --------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(9)
M4_MPE_SECTION(9,{
     call WriteDataOut(ncyc,.false.) 
})
call StopTimer(9)

! ---------------------------- WAIT H ---------------------------------

! complete pending H comms

call StartTimer(10)
M4_IFELSE_MPI({
M4_MPE_SECTION(10,{
	call CompleteHMPIComms
})
})
call StopTimer(10)

! ---------------------------- StepE ----------------------------------

call StartTimer(2)
M4_MPE_SECTION(6,{ 
	if ( ncyc .gt. 0 ) call StepE 
})
call StopTimer(2)

! -------------------------- StepELumped ---------------------------------

call StartTimer(3)
M4_MPE_SECTION(3,{
	if ( ncyc .gt. 0 ) call StepELumped(ncyc)
})
call StopTimer(3)

! -------------------------- StepESrc ---------------------------------

call StartTimer(4)
M4_MPE_SECTION(4,{
	if ( ncyc .gt. 0 ) call StepESrc(ncyc)
})
call StopTimer(4)

! -------------------------- StepEMat ---------------------------------

call StartTimer(5)
M4_MPE_SECTION(5,{
	if ( ncyc .gt. 0 ) call StepEMat(ncyc)
})
call StopTimer(5)

! -------------------------- StepEBound -------------------------------

call StartTimer(6)
M4_MPE_SECTION(6,{
	if ( ncyc .gt. 0 ) call StepEBound
})
call StopTimer(6)


! -------------------------- COMMS E ----------------------------------

! E field calculation is finished -> initiate comms

call StartTimer(10)
M4_IFELSE_MPI({
M4_MPE_SECTION(10,{
	call InitiateEMPIComms
})
})
call StopTimer(10)
	
! -------------------------- StepEDiag --------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(7)
M4_MPE_SECTION(7,{
	call StepEDiag(ncyc)
})
call StopTimer(7)

! ------------------------ WriteDataOut -------------------------------

! use the time while comms are pending to do some useful stuff

call StartTimer(9)
M4_MPE_SECTION(9,{
     call WriteDataOut(ncyc, .true.) 
})
call StopTimer(9)

! --------------------------- WAIT E ----------------------------------

! complete pending E comms

call StartTimer(10)
M4_IFELSE_MPI({
M4_MPE_SECTION(10,{
	call CompleteEMPIComms
})
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
     CALL DisplayTotalTimer("Lumped : ",3)
     call DisplayTotalTimer("Src    : ",4)
     call DisplayTotalTimer("Mat    : ",5)
     call DisplayTotalTimer("Bound  : ",6)
     call DisplayTotalTimer("Diag   : ",7)
     call DisplayTotalTimer("Output : ",9)
     call DisplayTotalTimer("Comms  : ",10)
     call DisplayTotalTimer("Total  : ",1)
     write(6,*) "*"
     cells = real(NCYCMAX) * real(IEIG-IBIG+1) * real(JEIG-JBIG+1) * real(KEIG-KBIG+1)
     call DisplayMcpsTimer( "Fdtd   : ",2, cells)
     cells = real(NCYCMAX) * real(IBIG-IBEG+1) * real( JBIG-JBEG+1) * real(KBIG-KBEG+1)
     call DisplayMcpsTimer( "Bound  : ",6, cells)

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

        if ( save_state ) then
           write(6,*) "* writing checkpoint file"
           open(unit=UNITCHK, file=checkpoint_fn, status = "unknown", form = "unformatted" ) 
        end if

        write(6,*) "* -> FinalizeList"
        call FinalizeList
        write(6,*) "* -> FinalizeGrid"
        call FinalizeGrid
        write(6,*) "* -> FinalizeFdtd"
        call FinalizeFdtd
        write(6,*) "* -> FinalizeLumped"
        call FinalizeLumped
        write(6,*) "* -> FinalizeBound"
        call FinalizeBound
        write(6,*) "* -> FinalizeSrc"
        call FinalizeSrc
        write(6,*) "* -> FinalizeMat"
        call FinalizeMat
        write(6,*) "* -> FinalizeDiag"
        call FinalizeDiag
        write(6,*) "* -> FinalizeOut"
        call FinalizeOut

        if ( save_state ) then
           close(UNITCHK) 
        end if

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
! Authors:  J.Hamm, A.Klaedtke, S.Scholz, R.Crowter
! Modified: 02/07/2009
!
!======================================================================
