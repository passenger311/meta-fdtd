!-*- F90 -*------------------------------------------------------------
!
!  module: config / meta
!
!  config file reader. 
!
!----------------------------------------------------------------------


! =====================================================================
!


module config

  use strings
  use parse
  use constant
  use reglist

  use grid
  use checkpoint
  use fdtd
  use bound
  use src
  use mat 
  use diag
  use lumped

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'CONFIG'

  ! --- Public Methods

  public :: ReadConfig

  ! --- Constants

  character(len=STRLNG), parameter :: pfxconfig = 'config'


contains

!----------------------------------------------------------------------


  subroutine ReadConfig(dim)
      
      integer :: dim
      character(len=STRLNG) :: file
      character(len=LINELNG) :: string, line, skiptill
      logical :: gotgrid = .false.
      logical :: gotfdtd = .false.
      logical :: err, eof
      integer :: ios, lcount = 1
      
!OMP SERIAL
      M4_WRITE_DBG({". enter ReadConfig"})

      file=cat2(pfxconfig,mpi_sfxin)
      

      M4_WRITE_INFO({"opening config (0): ", TRIM(file)})
      open(UNITTMP,FILE=file,STATUS="old", IOSTAT=ios)
      M4_OPEN_ERROR(ios,file)

      err = .false.
      lcount = 1

      do

         call readline(UNITTMP,lcount,eof,line)
         if ( eof ) exit
         
         call getstring(line,string,err)

         M4_SYNTAX_ERROR(err,lcount,{"[STRING]"})

         M4_WRITE_DBG({"got token ",TRIM(string)})
 
         select case ( string )
            
         case( "(GRID" )
            M4_WRITE_INFO({"-> ReadConfigGrid"})
            call ReadConfigGrid(UNITTMP,lcount,string,dim)
            gotgrid = .true.
         case( "(CHECKPOINT" )
            M4_WRITE_INFO({"-> ReadConfigCheckpoint"})
            call ReadConfigCheckpoint(UNITTMP,lcount,string)
         case( "(FDTD" )
            M4_WRITE_INFO({"-> ReadConfigFdtd"})
            call ReadConfigFdtd(UNITTMP,lcount,string)
            gotfdtd = .true.
         case( "(BOUND" )
            M4_WRITE_INFO({"-> ReadConfigBound"})
            call ReadConfigBound(UNITTMP,lcount,string)
         case( "(LUMPED" )
            M4_WRITE_INFO({"-> ReadConfigLumped"})
            call ReadConfigLumped(UNITTMP,lcount)
         case default
            if ( string(1:4) .eq. "(SRC" ) then
               M4_WRITE_INFO({"-> ReadConfigSrc: ",TRIM(string(5:))})
               call ReadConfigSrc(UNITTMP,lcount,string)
               cycle
            endif
               
            if ( string(1:4) .eq. "(MAT" ) then
               M4_WRITE_INFO({"-> ReadConfigMat: ",TRIM(string(5:))})
               call ReadConfigMat(UNITTMP,lcount,string)
               cycle
            endif
               
            if ( string(1:5) .eq. "(DIAG" ) then
               M4_WRITE_INFO({"-> ReadConfigDiag: ",TRIM(string(6:))})
               call ReadConfigDiag(UNITTMP,lcount,string)
               cycle
            endif
            
            M4_BADTOKEN_ERROR(err,lcount,string)

         end select

      end do

      M4_WRITE_INFO({"closing config (0): ", TRIM(file)})

      close(UNITTMP)

      if ( .not. gotgrid ) then
         M4_FATAL_ERROR({"NO GRID SECTION IN ",TRIM(file)})
      endif
      if ( .not. gotfdtd ) then
         M4_FATAL_ERROR({"NO FDTD SECTION IN ",TRIM(file)})
      endif
      
      M4_WRITE_DBG({". exit ReadConfig"})

!OMP END SERIAL

    end subroutine ReadConfig
      

!----------------------------------------------------------------------

  end module config

! Authors:  J.Hamm, R.Crowter 
! Modified: 02/07/2009
!
! =====================================================================






