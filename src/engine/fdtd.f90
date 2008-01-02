!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd / meta
!
!  fdtd core algorithm.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------


!======================================================================
!
! Note: MUINV and EPSINV are cell centered with respect to the H and E
! subgrids. This means we need to average in the respective direction
! of the component to retrieve the proper face-centered values.
!

module fdtd
 
  use constant
  use reglist
  use outlist
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
M4_IFELSE_WMU({
  public :: MUINV
})
  public :: pfxepsilon

  ! --- Constants

  character(len=20), parameter :: pfxepsilon = 'epsilon'

  ! --- Data

  M4_FTYPE, allocatable, dimension(:, :, :) :: Ex, Ey, Ez
  M4_FTYPE, allocatable, dimension(:, :, :) :: Hx, Hy, Hz
  real(kind=8), allocatable, dimension(:, :, :) :: EPSINV
  
  M4_IFELSE_WMU({
  real(kind=8), allocatable, dimension(:, :, :) :: MUINV
  })

  integer :: reginitidx = -1, regepsidx = -1, regmuidx = -1
  
contains

!----------------------------------------------------------------------

 subroutine ReadConfigFdtd(funit,string)

   integer :: funit
   character(len=*) :: string 
   character(len=STRLNG) :: line, skiptill
   type (T_OUT) :: out
   type (T_REG) :: reg
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
      end if
      M4_WRITE_DBG({"got token ",TRIM(string)})
 
      select case (string)
      case("(EHFIELDS") 
         call ReadRegObj(reg, fdtdreg, funit, 6)
		 reginitidx = reg%idx
      case("(EPSILON") 
         call ReadRegObj(reg, fdtdreg, funit, 1)
		 regepsidx = reg%idx
      M4_IFELSE_WMU({           
      case("(MU") 
         call ReadRegObj(reg, fdtdreg, funit, 1)
		 regmuidx = reg%idx
      })
      case("(OUT") 
         call ReadOutObj(out, fdtdreg, modname, funit)
      case default	
         if ( string(1:2) .eq. "(!" ) then
            skiptill = cat2(")",string(3:))
            M4_WRITE_DBG({"skiptill = ", TRIM(skiptill)})  
            cycle
         end if
         M4_WRITE_DBG({"read terminator: ", TRIM(string)})
         if ( string(1:1) .ne. ")" ) then
            M4_FATAL_ERROR({"BAD TERMINATOR: ",TRIM(string)," in ReadConfigFdtd"})
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
   M4_REGLOOP_DECL(reg,p,i,j,k,w(6))  
   real(kind=8) :: v(1)
       
   M4_WRITE_DBG({". enter InitializeFdtd"})

   if ( .not. modconfigured ) then
      M4_FATAL_ERROR({"NOT CONFIGURED: InitializeFdtd"})
   endif

   call InitializeEpsilon
M4_IFELSE_WMU({ call InitializeMu })
   call InitializeEHFields

   M4_WRITE_DBG({". exit InitializeFdtd"})   

 end subroutine InitializeFdtd

!----------------------------------------------------------------------

 subroutine InitializeEpsilon

   M4_REGLOOP_DECL(reg,p,i,j,k,v(1))  
   integer :: err

   M4_WRITE_DBG({"allocate epsilon"}) 

   allocate(EPSINV(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeEpsilon")
     
   M4_WRITE_DBG({"setting epsilon = 1.0"})   
   EPSINV = 1.0

   if ( regepsidx .ge. 1 ) then
     ! overwrite epsilon field with values from (EPSILON definition block
	  M4_WRITE_DBG({"initializing epsilon"})   
      reg = regobj(regepsidx)

      M4_REGLOOP_EXPR(reg,p,i,j,k,v, {

         EPSINV(i,j,k) = 1./v(1)
         
      } )
      ! get rid of initialization fields!
      if ( regobj(regepsidx)%numnodes .gt. 0 ) call DestroyRegObj(regobj(regepsidx))
      
   else
      M4_WRITE_WARN({"EPSILON NOT INITIALIZED, USING 1.0!"})   
   end if
   
   call ExtendRealField(EPSINV)
 
 end subroutine InitializeEpsilon
 
!----------------------------------------------------------------------

M4_IFELSE_WMU({

 subroutine InitializeMu

   M4_REGLOOP_DECL(reg,p,i,j,k,v(1))  
   integer :: err

   M4_WRITE_DBG({"allocate mu"})   

   allocate(MUINV(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateMuField")
     
   
   M4_WRITE_DBG({"setting mu = 1.0"})   
   MUINV = 1.0

   if ( regmuidx .ge. 1 ) then
     ! overwrite epsilon field with values from (EPSILON definition block
	  M4_WRITE_DBG({"initializing mu"})   
      reg = regobj(regmuidx)

      M4_REGLOOP_EXPR(reg,p,i,j,k,v, {

         MUINV(i,j,k) = 1./v(1)
         
      } )
      ! get rid of initialization fields!
      if ( regobj(regmuidx)%numnodes .gt. 0 ) call DestroyRegObj(regobj(regmuidx))
      
   else
      M4_WRITE_WARN({"MU NOT INITIALIZED, USING 1.0!"})   
   end if
   
   call ExtendRealField(MUINV)

 end subroutine InitializeMu

})

!----------------------------------------------------------------------

 subroutine InitializeEHFields

   M4_REGLOOP_DECL(reg,p,i,j,k,w(6))  
   integer :: err

   M4_WRITE_DBG({"allocate EH fields"})   

   allocate(Ex(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateEHFields")
   
   allocate(Ey(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateEHFields")
   
   allocate(Ez(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateEHFields")
   
   allocate(Hx(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateEHFields")
   
   allocate(Hy(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateEHFields")
   
   allocate(Hz(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "AllocateEHFields")

   M4_WRITE_DBG({"setting EH fields = 0.0"})   
   Ex = 0.0
   Ey = 0.0
   Ez = 0.0
   Hx = 0.0
   Hy = 0.0
   Hz = 0.0

   if ( reginitidx .ge. 1  ) then
      reg = regobj(reginitidx)
      M4_WRITE_DBG({"initializing EH fields"})
      M4_IFELSE_DBG({call EchoRegObj(reg)})
      M4_REGLOOP_EXPR(reg,p,i,j,k,w, {
         
         Ex(i,j,k) = w(1)
         Ey(i,j,k) = w(2)
         Ez(i,j,k) = w(3)
         Hx(i,j,k) = w(4)
         Hy(i,j,k) = w(5)
         Hz(i,j,k) = w(6)
         
         } )
      ! get rid of initialization fields!
      if ( regobj(reginitidx)%numnodes .gt. 0 ) call DestroyRegObj(regobj(reginitidx))
   end if

 end subroutine InitializeEHFields

!----------------------------------------------------------------------

 subroutine ExtendRealField(field)

   real(kind=8), dimension(:,:,:) :: field

   field(IEND+1,:,:) = field(IEND,:,:) 
   M4_IFELSE_2D({
   field(:,JEND+1,:) = field(:,JEND,:) 
   field(IEND+1,JEND+1,:) = field(IEND,JEND,:)
   })
   M4_IFELSE_3D({
   field(:,:,KEND+1) = field(:,:,KEND) 
   field(IEND+1,:,KEND+1) = field(IEND,:,KEND)
   field(:,JEND+1,KEND+1) = field(:,JEND,KEND)
   field(IEND+1,JEND+1,KEND+1) = field(IEND,JEND,KEND)
   })

 end subroutine ExtendRealField
 
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
M4_IFELSE_WMU({ deallocate(MUINV) })

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
    ! M4_IFELSE_1D({1D},{NOT 1D})
    ! M4_IFELSE_2D({2D},{NOT 2D})
    ! M4_IFELSE_3D({3D},{NOT 3D})

M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh)})
    do k=KSIG, KEIG-1

M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh)})
       do j=JSIG, JEIG-1

M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh)})
          do i=ISIG, IEIG-1

             Exh=Ex(i,j,k)
             Eyh=Ey(i,j,k)
             Ezh=Ez(i,j,k)

             Hx(i,j,k) =  Hx(i,j,k) + M4_MUINVX(i,j,k) * ( &
M4_IFELSE_3D({       + dtz*( Ey(i,j,k+1) - Eyh ) &   },{})
M4_IFELSE_1D({0.&},{ - dty*( Ez(i,j+1,k) - Ezh ) &      })
                     )
             Hy(i,j,k) = Hy(i,j,k) + M4_MUINVY(i,j,k) * ( &
                    + dtx*( Ez(i+1,j,k) - Ezh ) &
M4_IFELSE_3D({      - dtz*( Ex(i,j,k+1) - Exh ) &      })
                  )
             Hz(i,j,k) = Hz(i,j,k) +  M4_MUINVZ(i,j,k) * ( &
                    - dtx*( Ey(i+1,j,k) - Eyh ) &
M4_IFELSE_1D({},{   + dty*( Ex(i,j+1,k) - Exh ) &      })
				    )

          enddo
M4_IFELSE_1D({ !$OMP END PARALLEL DO })     

       enddo
M4_IFELSE_2D({ !$OMP END PARALLEL DO })     

    enddo
M4_IFELSE_3D({ !$OMP END PARALLEL DO })     
       
  end subroutine StepH

!----------------------------------------------------------------------

  subroutine StepE

    implicit none
    
    real(kind=8) :: dtx, dty, dtz
    M4_FTYPE :: Hxh, Hyh, Hzh
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz
    
    ! E in GT + DT

M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh)})
   do k=KSIG, KEIG-1
   
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh)})
       do j=JSIG, JEIG-1

M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh)})
          do i=ISIG, IEIG-1

             Hxh=Hx(i,j,k)
             Hyh=Hy(i,j,k)
             Hzh=Hz(i,j,k) 

             Ex(i,j,k) =  Ex(i,j,k) + M4_EPSINVX(i,j,k) *( &
M4_IFELSE_1D({0.&},{ + dty*( Hzh - Hz(i,j-1,k) ) &             })
M4_IFELSE_3D({       - dtz*( Hyh - Hy(i,j,k-1) ) &          },{})
                     )
             Ey(i,j,k) =  Ey(i,j,k) + M4_EPSINVY(i,j,k) * ( &
                     - dtx*( Hzh - Hz(i-1,j,k) ) &
M4_IFELSE_3D({       + dtz*( Hxh - Hx(i,j,k-1) ) &             })
                    )
             Ez(i,j,k) =  Ez(i,j,k) + M4_EPSINVZ(i,j,k) * ( &
                     + dtx*( Hyh - Hy(i-1,j,k) )  &
M4_IFELSE_1D({},{    - dty*( Hxh - Hx(i,j-1,k) )  &            })
                     )

          enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO })     

       enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO })     

    enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO })     
    
  end subroutine StepE

!----------------------------------------------------------------------

end module fdtd

! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 4/12/2007
!
!======================================================================
