!-*- F90 -*------------------------------------------------------------
!
!  module: mpiworld / meta
!
!  create mpi world 
!
!----------------------------------------------------------------------

!======================================================================
!
!


module mpiworld

  use strings
  use constant
  implicit none
  save 

M4_IFELSE_MPI({include 'mpif.h'})

  integer :: numproc = 1
  integer :: myrank = 0
  integer :: mpi_started = 0
  character(len=STRLNG) :: mpi_sfxin, mpi_sfxout, mpi_sfx
  character(len=20)   :: ranklbl
  logical :: mpibound(6)

M4_IFELSE_MPI({
  integer :: mpicomm
  integer :: mpierr
  integer :: mpistatus(MPI_STATUS_SIZE) 
})

M4_IFELSE_MPELOG({
  character(len=STRLNG),dimension(40) :: mpestrings
})

contains

!----------------------------------------------------------------------

  subroutine InitializeMPIWorld

    integer :: i

M4_IFELSE_MPI({

    mpicomm = MPI_COMM_WORLD

    call MPI_INIT(mpierr)
    
    if ( mpierr .eq. 0 ) call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, mpierr)
    if ( mpierr .eq. 0 ) call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpierr)

})

    mpi_sfxin = cat3(".",i2str(myrank),".in")
    mpi_sfxout = cat3(".",i2str(myrank),".out")
    mpi_sfx = cat2(".",i2str(myrank))
    ranklbl = cat5("(",i2str(myrank),"/",i2str(numproc),")")
    mpibound = .false. ! this will be setup properly in mpicomms

M4_IFELSE_MPI({

    if ( mpierr .ne. 0 ) then
       M4_FATAL_ERROR({"MPISTART FAILED!"})
    end if
    
    mpi_started = 1

  M4_IFELSE_MPELOG({

  if (myrank .eq. 0) then

  ! --- init logging

     call MPE_INIT_LOG 

     call RegMPEState(1,"ISEND","red:gray0")
     call RegMPEState(2,"IRECV","blue:gray1")
     call RegMPEState(3,"WAIT","green:gray2")
     call RegMPEState(4,"StepH","yellow:gray3")
     call RegMPEState(5,"StepHBound","violet:gray4")
     call RegMPEState(6,"StepE","orange1:gray5")
     call RegMPEState(7,"StepEBound","snow:gray6")
     call RegMPEState(8,"StepHMat","pink1:gray7")
     call RegMPEState(9,"StepEMat","chocolate:gray8")
     call RegMPEState(10,"StepHDiag","pink1:gray7")
     call RegMPEState(11,"StepEDiag","chocolate:gray8")

  end if
  })
})


  end subroutine InitializeMPIWorld

!----------------------------------------------------------------------

  subroutine RegMPEState(num,string, color)
    
    integer :: num
    character(len=*) :: string, color
    
M4_IFELSE_MPELOG({
  
    call MPE_Describe_state(num,40+num,string,color)    
    mpestrings(num) = string

})        
  end subroutine RegMPEState

!----------------------------------------------------------------------

  subroutine StartMPELog(num)
    
    integer :: num

M4_IFELSE_MPELOG({
        
    character(len=STRLNG) :: string

        string = mpestrings(num)
    call MPE_LOG_EVENT(num,string,mpierr)
    
})  
    
  end subroutine StartMPELog

  !----------------------------------------------------------------------

  subroutine StopMPELog(num)
    
    integer :: num

M4_IFELSE_MPELOG({

    character(len=STRLNG) :: string
    string = mpestrings(num)
    call MPE_LOG_EVENT(40+num,string,mpierr)    

})  
 
  end subroutine StopMPELog

!----------------------------------------------------------------------

  subroutine SynchronizeMPIWorld
  
M4_IFELSE_MPI({call MPI_BARRIER(MPI_COMM_WORLD,mpierr)},{})
  
  end subroutine SynchronizeMPIWorld
   
!----------------------------------------------------------------------

  subroutine FinalizeMPIWorld

M4_IFELSE_MPI({

   call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

   M4_IFELSE_MPELOG({call MPE_FINISH_LOG("meta-mpe.log")})

   call MPI_FINALIZE(mpierr)

})

  end subroutine FinalizeMPIWorld

!----------------------------------------------------------------------

end module mpiworld

!
! Authors:  J.Hamm, A.Klaedtke
! Modified: 4/12/2007
!
!======================================================================









