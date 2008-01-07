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

  use constant
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

  implicit none

  character(len=STRLNG), parameter :: modname = "META"
  integer :: ncyc, l 
  character(len=12) :: str
  character(len=1) :: busychar ='|'


! =========================== INITIALIZATION ==========================


  ! --- startup MPI

  call InitializeMPIWorld

  if (myrank .eq. 0) then

     write(6,*) "* ------------------------- META-ENGINE STARTS -------------------------- "
     write(6,*) "*"
     call DisplayVersion
     write(6,*) "*"
     write(6,*) "* ----------------------------------------------------------------------- "
  	 M4_WRITE_INFO("START")

  end if

  do l=0, numproc-1

M4_IFELSE_MPI(call SynchronizeMPIWorld)
  
     if (myrank .eq. l) then

        write(6,*) "*"
        write(6,*) "* initialising modules: myrank = ", ranklbl

        write(6,*) "* -> InitializeList"
        call InitializeList
        write(6,*) "* -> ReadConfig"
        call ReadConfig
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

  ncyc = 0
  do ncyc = 0, NCYCMAX

! do some pretty print, a cursor and a 100.0% progress indicator
     if (myrank .eq. 0 .and. mod(ncyc*10000/NCYCMAX,10) .eq. 0  ) then

        str = cat4(i2str(ncyc*100/NCYCMAX),".",i2str(mod(ncyc*1000/NCYCMAX,10)),"%")
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

M4_MPE_SECTION(4,{
	if ( ncyc .gt. 0 ) call StepH
})

! -------------------------- StepHBound -------------------------------

M4_MPE_SECTION(7,{
	if ( ncyc .gt. 0 )  call StepHBound
})

! -------------------------- StepHMat ---------------------------------

M4_MPE_SECTION(8,{
	if ( ncyc .gt. 0 ) call StepHMat(ncyc)
})

! ---------------------------- COMMS H --------------------------------

! H field calculation if finished -> initiate comms
M4_IFELSE_MPI({
	call InitiateHMPIComms
})

! -------------------------- StepHDiag --------------------------------

! use the time while comms are pending to do some useful stuff
M4_MPE_SECTION(10,{
	call StepHDiag(ncyc)
})

! ------------------------ LoadDataOut --------------------------------

! use the time while comms are pending to do some useful stuff
     call WriteDataOut(ncyc,.false.) 

! ---------------------------- WAIT H ---------------------------------

! complete pending H comms
M4_IFELSE_MPI({
	call CompleteHMPIComms
})

! ---------------------------- StepE ----------------------------------

M4_MPE_SECTION(6,{ 
	if ( ncyc .gt. 0 ) call StepE 
})

! -------------------------- StepEBound -------------------------------

M4_MPE_SECTION(7,{
	if ( ncyc .gt. 0 ) call StepEBound
})

! -------------------------- StepEMat ---------------------------------

M4_MPE_SECTION(9,{
	if ( ncyc .gt. 0 ) call StepEMat(ncyc)
})

! -------------------------- COMMS E ----------------------------------

! E field calculation is finished -> initiate comms
M4_IFELSE_MPI({
	call InitiateEMPIComms
})
	
! -------------------------- StepEDiag --------------------------------

! use the time while comms are pending to do some useful stuff
M4_MPE_SECTION(11,{
	call StepEDiag(ncyc)
})

! ------------------------ WriteDataOut -------------------------------

! use the time while comms are pending to do some useful stuff
     call WriteDataOut(ncyc, .true.) 

! --------------------------- WAIT E ----------------------------------

! complete pending E comms
M4_IFELSE_MPI({
	call CompleteEMPIComms
})

! ----------------------- END OF TIME LOOP ----------------------------

     GT = GT + DT

  end do

! ---------------------------------------------------------------------

M4_IFELSE_MPI(call FinalizeMPIComms)

  
! =========================== FINALIZING ==============================

! end of time loop.

  if (myrank .eq. 0) then

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
