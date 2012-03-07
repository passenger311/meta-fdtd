!-*- F90 -*------------------------------------------------------------
!
!  module: bound / meta
!
!  generic boundary conditions interface.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------

!======================================================================
!

module bound

  use constant
  use strings
  use mpiworld
  use grid  
  use fdtd
  use pec 
  use pmc
  use pbc
  use pml
  use cpml

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'BOUND'
  logical, private :: modconfigured = .false.

  ! --- Public Methods

  public :: ReadConfigBound
  public :: InitializeBound
  public :: FinalizeBound
  public :: StepEBound
  public :: StepHBound

  ! --- Public Data

  ! --- Constants

  ! --- Data

  integer, dimension(6) :: planebound=(/ 0,0,0,0,0,0 /)

contains

!----------------------------------------------------------------------

  subroutine ReadConfigBound(funit,lcount,string)

    integer :: funit, lcount
    character(len=*) :: string
    character(len=LINELNG) :: line
      
    logical :: err, eof
    integer :: ios, i

    M4_WRITE_DBG({". enter ReadConfigBound"})

    if ( string .ne. "(BOUND" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigBound"})
    endif

    call readintvec(funit, lcount, planebound, 6)
    M4_WRITE_DBG({"read planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    err = .false.

    do

      call readline(funit,lcount,eof,line)
      call getstring(line,string,err)
      
      M4_SYNTAX_ERROR(err,lcount,{"[STRING]"})
      M4_WRITE_DBG({"got token ",TRIM(string)})
 
      select case (string)
      case( "(PML" )
         M4_WRITE_INFO({"--> ReadConfigPml"})
         call ReadConfigPml(funit,lcount,string)
      case( "(CPML" )
         M4_WRITE_INFO({"--> ReadConfigCPml"})
         call ReadConfigCPml(funit,lcount,string)
      case( "(MPIBC" )
          M4_WRITE_DBG({"-> invoking ReadConfigMpibc"})
!          call ReadConfigMpibc(UNITTMP,lcount,string) 
      case( "(PBC" )
         M4_WRITE_INFO({"--> ReadConfigPbc"})
         call ReadConfigPbc(funit,lcount,string) 
      case( ")BOUND" )
         exit
      case default	
         M4_BADTOKEN_ERROR(.true.,lcount,{string})
      end select

    enddo
    
    modconfigured = .true.
    
    M4_WRITE_DBG({". exit ReadConfigBound"})

  end subroutine ReadConfigBound

!----------------------------------------------------------------------

  subroutine InitializeBound

    integer :: i
    integer :: planeset(6)

    M4_WRITE_DBG({". enter InitializeBound"})

    if ( .not. modconfigured ) then
       M4_WRITE_WARN({"NO BOUNDARY GIVEN -> DEFAULTING TO PEC"})
       return
    end if  

    do i = 1,M4_SDIM*2
       planebound(i) = Max(planebound(i),0)
    end do

    do i = M4_SDIM*2+1, 6
       planebound(i) = 0
    end do
    
    call InitializePec(planebound,0)
    call InitializePml(planebound,1)
    call InitializeCPml(planebound,4)
    call InitializePmc(planebound,2)
    call InitializePbc(planebound,3)

    M4_WRITE_DBG({". exit InitializeBound"})

  end subroutine InitializeBound

!----------------------------------------------------------------------

  subroutine FinalizeBound

    M4_WRITE_DBG({". enter FinalizeBound"})

    call FinalizePml
    call FinalizeCPml
    call FinalizePec
    call FinalizePmc
    call FinalizePbc

    M4_WRITE_DBG({". exit FinalizeBound"})

  end subroutine FinalizeBound

!----------------------------------------------------------------------

  ! update the H-fields of all boundary layers
  
  subroutine StepHBound

    integer :: i


    ! do pmls first!
    do i = 1, M4_SDIM*2

       if ( mpibound(i) ) cycle
       if ( planebound(i) .eq. 1 ) call StepHBoundPml(i)
       if ( planebound(i) .eq. 4 ) call StepHBoundCPml(i)

    end do

    do i = 1, M4_SDIM*2

       if ( mpibound(i) ) cycle

       select case ( planebound(i) )

          case ( 0 )
             call StepHBoundPec(i)
          case ( 1 )
             ! call StepHBoundPml(i)
          case ( 2 )
             call StepHBoundPmc(i)
          case ( 3 )
             call StepHBoundPbc(i)
          case ( 4 )
             ! call StepHBoundCPml(i)

          case default
             M4_FATAL_ERROR({"BOUNDARY CONDITION # ",TRIM(i2str(planebound(i))),&
                  " NOT IMPLEMENTED!"})

       end select

    end do

  end subroutine StepHBound

!----------------------------------------------------------------------

  ! update the E-fields of all boundary layers
  
  subroutine StepEBound

    integer :: i

    ! do pmls first!
    do i = 1, M4_SDIM*2

       if ( mpibound(i) ) cycle
       if ( planebound(i) .eq. 1 ) call StepEBoundPml(i)
       if ( planebound(i) .eq. 4 ) call StepEBoundCPml(i)

    end do

    do i = 1, M4_SDIM*2

       if ( mpibound(i) ) cycle

       select case ( planebound(i) )
          
       case ( 0 )
          call StepEBoundPec(i)
       case ( 1 )
          ! call StepEBoundPml(i)
       case ( 2 )
          call StepEBoundPmc(i)
       case ( 3 )
          call StepEBoundPbc(i)
       case ( 4 )
          ! call StepEBoundCPml(i)

       case default
          M4_FATAL_ERROR({"BOUNDARY CONDITION # ",TRIM(i2str(planebound(i))),&
               " NOT IMPLEMENTED!"})
          
       end select

    end do
    
  end subroutine StepEBound
  
!----------------------------------------------------------------------

end module bound

!
! Authors:  J.Hamm
! Modified: 4/12/2007
!
!======================================================================
