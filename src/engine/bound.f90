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
  use sbc
  use pbc
  use pml

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

  subroutine ReadConfigBound(funit,string)

    integer :: funit
    character(len=*) :: string
    character(len=STRLNG) :: skiptill, line
      
    integer :: ios, i

    M4_WRITE_DBG({". enter ReadConfigBound"})

    M4_WRITE_DBG({"received token: ", TRIM(string)})
    if ( string .ne. "(BOUND" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigBound"})
    endif

    read(UNITTMP,*) (planebound(i),i=1, M4_SDIM*2)
    M4_WRITE_DBG({"read planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    skiptill = ""
    do	
       read(funit,*) line
       string = TRIM(ADJUSTL(line))
       
       if ( skiptill .ne. "" ) then 
          M4_WRITE_DBG({"skipping line ",TRIM(string)})
          if ( string .eq. skiptill ) skiptill = ""  
          cycle              
       endif
       
       select case (string)
       case( "(PML" )
          M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigPml"})
          call ReadConfigPml(UNITTMP,string) 
   
        case( "(PBC" )
          M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigPbc"})
          call ReadConfigPbc(UNITTMP,string) 
       
          ! add optional configs for other boundaries here

       case default	
          if ( string(1:2) .eq. "(!" ) then
             skiptill = cat2(")",string(3:))
             M4_WRITE_DBG({"got token (! -> skiptill = ", TRIM(skiptill)})  
             cycle
          end if
          M4_WRITE_DBG({"read terminator: ", TRIM(string)})
          if ( string(1:1) .ne. ")" ) then
             M4_FATAL_ERROR({"BAD TERMINATOR: ReadConfigBound"})
          end if
          exit
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
    
    call InitializePec(planebound,0)
    call InitializePml(planebound,1)
    call InitializePbc(planebound,2)

    M4_WRITE_DBG({". exit InitializeBound"})

  end subroutine InitializeBound

!----------------------------------------------------------------------

  subroutine FinalizeBound

    M4_WRITE_DBG({". enter FinalizeBound"})

    call FinalizePml
    call FinalizePec
    call FinalizePbc

    M4_WRITE_DBG({". exit FinalizeBound"})

  end subroutine FinalizeBound

!----------------------------------------------------------------------

  ! update the H-fields of all boundary layers
  
  subroutine StepHBound

    integer :: i

    do i = 1, M4_SDIM*2

!       M4_WRITE_DBG({"step h bound",i})

       if ( mpibound(i) ) cycle

       select case ( planebound(i) )

          case ( 0 )
             call StepHBoundPec(i)
          case ( 1 )
             call StepHBoundPml(i)
          case ( 2 )
             call StepHBoundSbc(i)
          case ( 3 )
             call StepHBoundPbc(i)

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

    do i = 1, M4_SDIM*2

 !      M4_WRITE_DBG({"step e bound",i})

       if ( mpibound(i) ) cycle

       select case ( planebound(i) )
          
       case ( 0 )
          call StepEBoundPec(i)
       case ( 1 )
          call StepEBoundPml(i)
       case ( 2 )
          call StepEBoundPbc(i)

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
