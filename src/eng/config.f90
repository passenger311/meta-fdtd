!-*- F90 -*------------------------------------------------------------
!
!  module: config / meta3
!
!  config file reader. 
!
!  subs:
!
!    ReadConfig
!
!----------------------------------------------------------------------


! =====================================================================
!


module config

  use strings
  use constant
  use reglist

  use grid
  use fdtd
  use bound
  use mat 
  use diag

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


  subroutine ReadConfig
      
      character(len=STRLNG) :: file, string, line, skiptill
      logical :: gotgrid = .false.
      logical :: gotfdtd = .false.
      integer :: ios
      
      M4_WRITE_DBG({". enter ReadConfig"})

      file=cat2(pfxconfig,mpi_sfxin)
      
      M4_WRITE_DBG({"trying to open ", TRIM(file)})

      open(UNITTMP,FILE=file,STATUS="old", IOSTAT=ios)
      M4_OPEN_ERROR(ios,file)

      skiptill = ""
      do
         read(UNITTMP,*, IOSTAT=ios) line
         string = TRIM(ADJUSTL(line))

         if(ios .ne. 0) exit

         if ( skiptill .ne. "" ) then 
            M4_WRITE_DBG({"skipping line ",TRIM(string)})
            if ( string .eq. skiptill ) skiptill = ""  
            cycle              
         endif
 
         select case ( string )
            
         case( "(GRID" )
            M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigGrid"})
            call ReadConfigGrid(UNITTMP,string)
            gotgrid = .true.
         case( "(FDTD" )
            M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigFdtd"})
            call ReadConfigFdtd(UNITTMP,string)
            gotfdtd = .true.
         case( "(BOUND" )
            M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigPml"})
            call ReadConfigBound(UNITTMP,string)    
         case("")
         case default

            if ( string(1:2) .eq. "(!" ) then
               skiptill = cat2(")",string(3:))
               M4_WRITE_DBG({"got token (! -> skiptill = ", TRIM(skiptill)})  
               cycle
            end if

            if ( string(1:4) .eq. "(MAT" ) then
            M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigMat"})
               call ReadConfigMat(UNITTMP,string)
               cycle
            endif
            
            if ( string(1:5) .eq. "(DIAG" ) then
               M4_WRITE_DBG({"got token ",TRIM(string),"-> invoking ReadConfigDiag"})
               call ReadConfigDiag(UNITTMP,string)
               cycle
            endif

            if ( string(1:1) .ne. "!" ) then
               M4_FATAL_ERROR({"BAD TOKEN ", TRIM(string)})
            endif

         end select
         
      enddo
      close(UNITTMP)

      if ( .not. gotgrid ) then
         M4_FATAL_ERROR({"NO GRID SECTION IN ",TRIM(file)})
      endif
      if ( .not. gotfdtd ) then
         M4_FATAL_ERROR({"NO FDTD SECTION IN ",TRIM(file)})
      endif

      M4_WRITE_DBG({". exit ReadConfig"})

    end subroutine ReadConfig
      

!----------------------------------------------------------------------

  end module config

! =====================================================================






