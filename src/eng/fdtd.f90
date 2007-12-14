!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd / meta3
!
!  fdtd core algorithm.
!
!  subs:
!
!  InitializeFdtd
!    AllocateFields
!    ReadEpsilonField
!  FinalizeFdtd
!  StepE
!  StepH
!
!----------------------------------------------------------------------


!======================================================================
!
! 
!

module fdtd
 
  use constant
  use reglist
  use outlist
  use initlist
  use grid
 
  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), parameter :: modname = 'FDTD'
  logical, private :: modconfigured

  ! --- Public Methods

  public :: ReadConfigFdtd
  public :: InitializeFdtd
  public :: FinalizeFdtd
  public :: StepE
  public :: StepH

  ! --- Public Data

  public :: Ex, Ey, Ez, Hx, Hy, Hz, EPSINV
  public :: pfxepsilon

  ! --- Constants

  character(len=20), parameter :: pfxepsilon = 'epsilon'

  ! --- Data

  M4_FTYPE, allocatable, dimension(:, :, :) :: Ex, Ey, Ez
  M4_FTYPE, allocatable, dimension(:, :, :) :: Hx, Hy, Hz
  real(kind=8), allocatable, dimension(:, :, :) :: EPSINV

contains

!----------------------------------------------------------------------

 subroutine ReadConfigFdtd(funit,string)

   integer :: funit
   character(len=*) :: string 
   character(len=STRLNG) :: line, skiptill
   type (T_OUT) :: out
   type (T_INIT) :: init
   real(kind=8) :: val = 1.0

   M4_WRITE_DBG({". enter ReadConfigFdtd"})
      
   if ( string .ne. "(FDTD" ) then
      M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigFdtd"})
   endif
   
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
      case("(INIT") 
         M4_WRITE_DBG({"got token (INIT -> ReadInitObj"})
         call ReadInitObj(init, fdtdreg, 6, "INIT", funit)
      case("(EPSILON") 
         M4_WRITE_DBG({"got token (EPSILON -> ReadInitObj"})
         call ReadInitObj(init, fdtdreg, 1, "EPSILON", funit)
      case("(OUT") 
	M4_WRITE_DBG({"got token (OUT -> ReadOutObj"})
        call ReadOutObj(out, fdtdreg, modname, funit)
     case default	
        if ( string(1:2) .eq. "(!" ) then
           skiptill = cat2(")",string(3:))
           M4_WRITE_DBG({"got token (! -> skiptill = ", TRIM(skiptill)})  
           cycle
        end if
        M4_WRITE_DBG({"read terminator: ", TRIM(string)})
        if ( string(1:1) .ne. ")" ) then
          M4_FATAL_ERROR({"BAD TERMINATOR: ReadConfigFdtd"})
        end if
        exit
      end select
    enddo	

   modconfigured = .true.

   M4_WRITE_DBG({". exit ReadConfigFdtd"})
   
 end subroutine ReadConfigFdtd
 

!----------------------------------------------------------------------

 subroutine InitializeFdtd
   
   integer :: n
   type(T_INIT) :: init
   M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    
   M4_WRITE_DBG({". enter InitializeFdtd"})

   if ( .not. modconfigured ) then
      M4_FATAL_ERROR({"NOT CONFIGURED: InitializeFdtd"})
   endif
 
   call AllocateFields

   M4_WRITE_DBG({"setting inital values of field components"})   

   Ex = 0.0
   Ey = 0.0
   Ez = 0.0
   Hx = 0.0
   Hy = 0.0
   Hz = 0.0

   do n = 1, numinitobj

      init = initobj(n)
      M4_WRITE_DBG({"found initializer type: ",TRIM(init%type)})
      if ( init%type .eq. "INIT" ) then 
         reg = regobj(init%regidx)
         M4_IFELSE_DBG({call EchoRegObj(regobj(numregobj))})
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,
         Ex(i,j,k) = w*init%val(1)
         Ey(i,j,k) = w*init%val(2)
         Ez(i,j,k) = w*init%val(3)
         Hx(i,j,k) = w*init%val(4)
         Hy(i,j,k) = w*init%val(5)
         Hz(i,j,k) = w*init%val(6)
         )
      endif

   end do

   M4_WRITE_DBG({"setting epsilon = 1.0"})   
   EPSINV = 1.0

   M4_WRITE_DBG({"reading epsilon field"})   
   call ReadEpsilonField

   M4_WRITE_DBG({"overwriting epsilon with data from (EPSILON-)"})   
   ! overwrite epsilon field with values from (EPSILON definition block
   do n = 1, numinitobj

      init = initobj(n)
      if ( init%type .eq. "EPSILON" ) then 
         reg = regobj(init%regidx)
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,
         EPSINV(i,j,k) = 1./(w*init%val(1))
         )
      endif

   end do

   M4_WRITE_DBG({". exit InitializeFdtd"})   

 contains
      
   subroutine AllocateFields

     integer :: err

     M4_WRITE_DBG({". enter InitializeFdtd.AllocateFields"})
     
     allocate(Ex(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")
     
     allocate(Ey(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")
     
     allocate(Ez(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")
     
     allocate(Hx(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")

     allocate(Hy(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")
     
     allocate(Hz(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")
     
     allocate(EPSINV(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields")
     
     M4_WRITE_DBG({". exit InitializeFdtd.AllocateFields"})

   end subroutine AllocateFields
   
   
   subroutine ReadEpsilonField

     integer :: ios, i, j, k
     real(kind=8) :: val
     
     character(len=STRLNG) :: file
     
     M4_WRITE_DBG({". enter InitializeFdtd.ReadEpsilonField"})

     file = cat2(pfxepsilon,mpi_sfxin)
     
     M4_WRITE_DBG({"trying to open ", file})
     
     open(UNITTMP, FILE=file, STATUS='old', IOSTAT=ios)

     if ( ios .eq. 0 ) then

        do k=KBEG,KEND+1
           do j=JBEG,JEND+1
              do i=IBEG, IEND+1
                 read(UNITTMP,*) val
                 EPSINV(i,j,k)=1.0/val
              end do
           end do
        end do
     
     close(UNITTMP)

     else 
        M4_WRITE_WARN({"COULD NOT OPEN ", TRIM(file),"!"})
     endif


     M4_WRITE_DBG({". exit InitializeFdtd.ReadEpsilonField"})
     
   end subroutine ReadEpsilonField
   
 end subroutine InitializeFdtd

!----------------------------------------------------------------------

 subroutine FinalizeFdtd
   
   M4_WRITE_DBG({". enter FinalizeFdtd"})
   
   deallocate(Hz)
   deallocate(Hy)
   deallocate(Hx)
   deallocate(Ez)
   deallocate(Ey)
   deallocate(Ex)
   deallocate(EPSINV)

   M4_WRITE_DBG({". exit FinalizeFdtd"})
   
 end subroutine FinalizeFdtd
 
 !----------------------------------------------------------------------

  subroutine StepH

    implicit none
    
    real(kind=8) :: dtx, dty, dtz
    M4_FTYPE :: Exh,Eyh,Ezh
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz


    ! H in GT+1/2DT

!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh)
    do k=KSIG, KEIG-1
       do j=JSIG, JEIG-1
          do i=ISIG, IEIG-1

             Exh=Ex(i,j,k)
             Eyh=Ey(i,j,k)
             Ezh=Ez(i,j,k)

             Hx(i,j,k) =  Hx(i,j,k) &
                  - dty*( Ez(i,j+1,k) - Ezh ) &
                  + dtz*( Ey(i,j,k+1) - Eyh )
             Hy(i,j,k) = Hy(i,j,k) &
                  - dtz*( Ex(i,j,k+1) - Exh ) &
                  + dtx*( Ez(i+1,j,k) - Ezh )
             Hz(i,j,k) = Hz(i,j,k) &
                  - dtx*( Ey(i+1,j,k) - Eyh ) &
                  + dty*( Ex(i,j+1,k) - Exh )

          enddo
       enddo
    enddo
 !$OMP END PARALLEL DO        
       
  end subroutine StepH

!----------------------------------------------------------------------

  subroutine StepE

    implicit none
    
    real(kind=8) :: dtx, dty, dtz
    M4_FTYPE :: Hxh, Hyh, Hzh
    real(kind=8) :: epsinvx, epsinvy, epsinvz
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz
    
    ! E in GT + DT

 !$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,epsinvx,epsinvy,epsinvz)
   do k=KSIG, KEIG-1
       do j=JSIG, JEIG-1
          do i=ISIG, IEIG-1

             Hxh=Hx(i,j,k)
             Hyh=Hy(i,j,k)
             Hzh=Hz(i,j,k) 

             epsinvx = 0.5*(EPSINV(i,j,k) +  EPSINV(i+1,j,k))
             epsinvy = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j+1,k))
             epsinvz = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j,k+1))

             Ex(i,j,k) =  Ex(i,j,k) +  epsinvx* &
                  ( dty*( Hzh - Hz(i,j-1,k) ) &
                  - dtz*( Hyh - Hy(i,j,k-1) ))
             Ey(i,j,k) =  Ey(i,j,k) +  epsinvy* &
                  ( dtz*( Hxh - Hx(i,j,k-1) ) &
                  - dtx*( Hzh - Hz(i-1,j,k) ))
             Ez(i,j,k) =  Ez(i,j,k) +  epsinvz* &
                  ( dtx*( Hyh - Hy(i-1,j,k) )  &
                  - dty*( Hxh - Hx(i,j-1,k) ))

          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
    
  end subroutine StepE

!----------------------------------------------------------------------

end module fdtd

!
! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
