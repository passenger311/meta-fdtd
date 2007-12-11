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

  type(T_REG) :: fdtdreg

contains

!----------------------------------------------------------------------

 subroutine ReadConfigFdtd(funit,string)

   integer :: funit
   character(len=*) :: string
   type (T_OUT) :: out
   type (T_REG) :: reg
   type (T_INIT) :: init
   real(kind=8) :: val = 1.0

   integer :: ios

   M4_WRITE_DBG({". enter ReadConfigFdtd/fdtd"})
      
   reg = CreateRegObjStart()
   call SetBoxRegObj(reg, IBEG, IEND, 1, JBEG, JEND, 1, KBEG, KEND, 1, val)
   call CreateRegObjEnd(reg,.false.)

   if ( string .ne. "(FDTD" ) then
      M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigFdtd/fdtd"})
   endif
   
   do	
      read(funit,*) string
      select case (string)
      case("(INIT") 
         M4_WRITE_DBG({"got token (INIT -> ReadInitObj"})
         call ReadInitObj(init, reg, 6, "INIT", funit)
      case("(EPSILON") 
         M4_WRITE_DBG({"got token (EPSILON -> ReadInitObj"})
         call ReadInitObj(init, reg, 1, "EPSILON", funit)
      case("(OUT") 
	M4_WRITE_DBG({"got token (OUT -> ReadOutObj"})
        call ReadOutObj(out, reg, modname, funit)
      case default	
        M4_WRITE_DBG({"read terminator: ", TRIM(string)})
        if ( string(1:1) .ne. ")" ) then
          M4_FATAL_ERROR({"BAD TERMINATOR: ReadConfigFdtd/fdtd"})
        end if
        exit
      end select
    enddo	

   modconfigured = .true.

   M4_WRITE_DBG({". exit ReadConfigFdtd/fdtd"})
   
 end subroutine ReadConfigFdtd
 

!----------------------------------------------------------------------

 subroutine InitializeFdtd
   
   integer :: n
   type(T_INIT) :: init
   M4_REGLOOP_DECL(reg,p,i,j,k,w)  
    
   M4_WRITE_DBG({". enter InitializeFdtd/fdtd"})

   if ( .not. modconfigured ) then
      M4_FATAL_ERROR({"NOT CONFIGURED: InitializeFdtd/fdtd"})
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
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,
         Ex(i,j,k) = init%val(1)
         Ey(i,j,k) = init%val(2)
         Ez(i,j,k) = init%val(3)
         Hx(i,j,k) = init%val(4)
         Hy(i,j,k) = init%val(5)
         Hz(i,j,k) = init%val(6)
         )
      endif

   end do

   EPSINV = 1.0

   call ReadEpsilonField

   ! overwrite epsilon field with values from (EPSILON definition block
   do n = 1, numinitobj

      init = initobj(n)
      if ( init%type .eq. "EPSILON" ) then 
         reg = regobj(init%regidx)
         M4_REGLOOP_EXPR(reg,p,i,j,k,w,
         EPSINV(i,j,k) = 1./init%val(1)
         )
      endif

   end do

   M4_WRITE_DBG({". exit InitializeFdtd/fdtd"})   

 contains
      
   subroutine AllocateFields

     integer :: err

     M4_WRITE_DBG({". enter InitializeFdtd/fdtd | AllocateFields"})
     
     allocate(Ex(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Ey(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Ez(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Hx(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")

     allocate(Hy(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(Hz(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     allocate(EPSINV(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), STAT=err)
     M4_ALLOC_ERROR(err, "AllocateFields/fdtd")
     
     M4_WRITE_DBG({". exit InitializeFdtd/fdtd | AllocateFields"})

   end subroutine AllocateFields
   
   
   subroutine ReadEpsilonField

     integer :: ios, i, j, k
     real(kind=8) :: val
     
     character(len=STRLNG) :: file
     
     M4_WRITE_DBG({". enter InitializeFdtd/fdtd | ReadEpsilonField"})

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
        write(6,*) "!WARN COULD NOT OPEN ", TRIM(file),"!"
     endif


     M4_WRITE_DBG({". exit InitializeFdtd/fdtd | ReadEpsilonField"})
     
   end subroutine ReadEpsilonField
   
 end subroutine InitializeFdtd

!----------------------------------------------------------------------

 subroutine FinalizeFdtd
   
   M4_WRITE_DBG({". enter FinalizeFdtd/fdtd"})
   
   deallocate(Hz)
   deallocate(Hy)
   deallocate(Hx)
   deallocate(Ez)
   deallocate(Ey)
   deallocate(Ex)
   deallocate(EPSINV)

   M4_WRITE_DBG({". exit FinalizeFdtd/fdtd"})
   
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
    
  end subroutine StepE

!----------------------------------------------------------------------

end module fdtd

!
! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
