!-*- F90 -*------------------------------------------------------------
!
!  module: fdtd / meta
!
!  fdtd core algorithm.
!
!----------------------------------------------------------------------


!======================================================================
!
! 1. note: (outdated) -> 2.
!
! MUINV and EPSINV are cell centered with respect to the H and E
! subgrids. This means we need to average in the respective direction
! of the component to retrieve the proper face-centered values.
!
! 2. update to 1.: 
!
! Instead of 1., read in 3 fields, epsinvx epsinvy, epsinvz which are
! each centered at the same node than Ex, Ey, Ez. Same for muinv. This
! implies that no average needs to be used. The effort of calculating
! the epsilons on the various faces has to be done when configuring the
! engine for a simulation.
!
! 3. note:
!
! StepH : update eq. H(n+1/2) = H(n-1/2) + ... E(n) ...
! StepE : update eq. E(n+1) = E(n) + ... H(n+1/2) ...
!
! StepE, StepH only include the static (epsilon/mu) medium response.
! The contributions of the dynamic parts (P and Q or J and K) are added
! just after StepE and StepH respectively.
!
!
!
!
module fdtd
 
  use constant
  use checkpoint
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

  public :: Ex, Ey, Ez, Hx, Hy, Hz
  public :: epsinvx, epsinvy, epsinvz
M4_IFELSE_WMU({
  public :: muinvx, muinvy, muinvz
})
  public :: pfxepsilon

  ! --- Constants

  character(len=20), parameter :: pfxepsilon = 'epsilon'

  ! --- Data

  M4_FTYPE, allocatable, dimension(:, :, :) :: Ex, Ey, Ez
  M4_FTYPE, allocatable, dimension(:, :, :) :: Hx, Hy, Hz
  real(kind=8), allocatable, dimension(:, :, :) :: epsinvx
  real(kind=8), allocatable, dimension(:, :, :) :: epsinvy
  real(kind=8), allocatable, dimension(:, :, :) :: epsinvz
  
  M4_IFELSE_WMU({
  real(kind=8), allocatable, dimension(:, :, :) :: muinvx
  real(kind=8), allocatable, dimension(:, :, :) :: muinvy
  real(kind=8), allocatable, dimension(:, :, :) :: muinvz
  })

  integer :: reginitidx = -1, regepsidx = -1, regmuidx = -1
  
contains

!----------------------------------------------------------------------

 subroutine ReadConfigFdtd(funit,lcount,string)

   integer :: funit, lcount
   character(len=*) :: string 
   character(len=LINELNG) :: line
   type (T_OUT) :: out
   type (T_REG) :: reg
   logical :: err, eof

   M4_WRITE_DBG({". enter ReadConfigFdtd"})
      
   M4_SYNTAX_ERROR({string .ne. "(FDTD"},lcount,{"(FDTD"})
   
   err = .false.
   
   do

      call readline(funit,lcount,eof,line)
      call getstring(line,string,err)

      M4_SYNTAX_ERROR(err,lcount,{"[STRING]"})
      M4_WRITE_DBG({"got token ",TRIM(string)})
 
      select case (string)
      case("(EHFIELDS") 
         call readline(funit,lcount,eof,line)
         call getstring(line,string,err)
         M4_SYNTAX_ERROR({string .ne. "(REG"},lcount,{"(REG"})
         call ReadRegObj(reg, fdtdreg, funit, lcount, 6, .false.)
         reginitidx = reg%idx
      case("(EPSILON")
         call readline(funit,lcount,eof,line)
         call getstring(line,string,err)
         M4_SYNTAX_ERROR({string .ne. "(REG"},lcount,{"(REG"})
         call ReadRegObj(reg, fdtdreg, funit, lcount, 3, .false.)
         regepsidx = reg%idx
      M4_IFELSE_WMU({           
      case("(MU") 
         call readline(funit,lcount,eof,line)
         call getstring(line,string,err)
         M4_SYNTAX_ERROR({string .ne. "(REG"},lcount,{"(REG"})
         call ReadRegObj(reg, fdtdreg, funit, lcount, 3, .false.)
         regmuidx = reg%idx
      })
      case("(OUT") 
         call ReadOutObj(out, fdtdreg, funit, lcount, modname, .false.)
      case(")FDTD")
         exit
      case default	
         M4_BADTOKEN_ERROR(err,lcount,string)
      end select

   enddo

   modconfigured = .true.

   M4_WRITE_DBG({". exit ReadConfigFdtd"})
   
 end subroutine ReadConfigFdtd
 

!----------------------------------------------------------------------

 subroutine InitializeFdtd
   
   integer :: i,j,k

   M4_WRITE_DBG({". enter InitializeFdtd"})

   if ( .not. modconfigured ) then
      M4_FATAL_ERROR({"NOT CONFIGURED: InitializeFdtd"})
   endif

   call InitializeEpsilon
M4_IFELSE_WMU({ call InitializeMu })
   call InitializeEHFields

! load from checkpoint file

   if ( load_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

      read(UNITCHK) Ex, Ey, Ez
      read(UNITCHK) Hx, Hy, Hz

   end if


   M4_WRITE_DBG({". exit InitializeFdtd"})   

 end subroutine InitializeFdtd

!----------------------------------------------------------------------

 subroutine InitializeEpsilon

   M4_REGLOOP_DECL(reg,p,i,j,k,v(3))  
   integer :: err

   M4_WRITE_DBG({"allocate epsilon"}) 

   allocate(epsinvx(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeEpsilon")

   allocate(epsinvy(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeEpsilon")

   allocate(epsinvz(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeEpsilon")
     
   M4_WRITE_DBG({"pre-initializing all epsilon = 1.0"})   
   epsinvx = 1.0
   epsinvy = 1.0
   epsinvz = 1.0

   if ( regepsidx .ge. 1 ) then

      M4_WRITE_DBG({"initializing epsilon from input"})   
      reg = regobj(regepsidx)

      M4_REGLOOP_EXPR(reg,p,i,j,k,v, {

         epsinvx(i,j,k) = 1./v(1)
         epsinvy(i,j,k) = 1./v(2)
         epsinvz(i,j,k) = 1./v(3)
         
      } )

      ! free the memory of the reg object
      if ( regobj(regepsidx)%numnodes .gt. 0 ) call DestroyRegObj(regobj(regepsidx))
      
   else
      M4_WRITE_WARN({"EPSILON NOT INITIALIZED, USING 1.0!"})   
   end if

   M4_WRITE_DBG({"extend epsilon towards boundary"})   

   call ExtendField(epsinvx)
   call ExtendField(epsinvy)
   call ExtendField(epsinvz)
    
 end subroutine InitializeEpsilon
 
!----------------------------------------------------------------------

M4_IFELSE_WMU({

 subroutine InitializeMu

   M4_REGLOOP_DECL(reg,p,i,j,k,v(3))  
   integer :: err

   M4_WRITE_DBG({"allocate mu"})   

   allocate(muinvx(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeMu")

   allocate(muinvy(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeMu")

   allocate(muinvz(M4_RANGE(IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)), STAT=err)
   M4_ALLOC_ERROR(err, "InitializeMu")
     
   M4_WRITE_DBG({"pre-initializing all epsilon = 1.0"})   
   muinvx = 1.0
   muinvy = 1.0
   muinvz = 1.0

   if ( regmuidx .ge. 1 ) then

      M4_WRITE_DBG({"initializing mu from input"})   
      reg = regobj(regmuidx)

      M4_REGLOOP_EXPR(reg,p,i,j,k,v, {

         muinvx(i,j,k) = 1./v(1)
         muinvy(i,j,k) = 1./v(2)
         muinvz(i,j,k) = 1./v(3)
         
      } )

      ! free memory allocated by reg object
      if ( regobj(regmuidx)%numnodes .gt. 0 ) call DestroyRegObj(regobj(regmuidx))
      
   else
      M4_WRITE_WARN({"MU NOT INITIALIZED, USING 1.0!"})   
   end if
   
   M4_WRITE_DBG({"extend mu towards boundary"})   

   call ExtendField(muinvx)
   call ExtendField(muinvy)
   call ExtendField(muinvz)

 end subroutine InitializeMu

})

!----------------------------------------------------------------------

subroutine ExtendField(field)

  real(kind=8) :: field(M4_RANGE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX))


  field(IMIN,JBEG:JEND,KBEG:KEND) = field(IBEG,JBEG:JEND,KBEG:KEND)
  field(IMAX,JBEG:JEND,KBEG:KEND) = field(IEND,JBEG:JEND,KBEG:KEND)

  M4_IFELSE_1D({},{

  field(IBEG:IEND,JMIN,KBEG:KEND) = field(IBEG:IEND,JBEG,KBEG:KEND)
  field(IBEG:IEND,JMAX,KBEG:KEND) = field(IBEG:IEND,JEND,KBEG:KEND)

  })

  M4_IFELSE_3D({

  field(IBEG:IEND,JBEG:JEND,KMIN) = field(IBEG:IEND,JBEG:JEND,KBEG)
  field(IBEG:IEND,JBEG:JEND,KMAX) = field(IBEG:IEND,JBEG:JEND,KEND)

  })

end subroutine ExtendField


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

 subroutine FinalizeFdtd
   
   integer :: i,j,k

   M4_WRITE_DBG({". enter FinalizeFdtd"})
   

! save to checkpoint file

   if ( save_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

      write(UNITCHK) Ex, Ey, Ez
      write(UNITCHK) Hx, Hy, Hz

   end if

   deallocate(Hz)
   deallocate(Hy)
   deallocate(Hx)
   deallocate(Ez)
   deallocate(Ey)
   deallocate(Ex)
   deallocate(epsinvx)
   deallocate(epsinvy)
   deallocate(epsinvz)
M4_IFELSE_WMU({ 
   deallocate(muinvx)
   deallocate(muinvy)
   deallocate(muinvz)
})

   M4_WRITE_DBG({". exit FinalizeFdtd"})
   
 end subroutine FinalizeFdtd
 
 !----------------------------------------------------------------------

  subroutine StepH

    M4_FTYPE :: Exh,Eyh,Ezh
    integer :: i, j, k

! !$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) PRIVATE(Exh,Eyh,Ezh) COLLAPSE(3)
! M4_DIM)
M4_IFELSE_3D({!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) &
!$OMP& PRIVATE(Exh,Eyh,Ezh,i,j,k) COLLAPSE(3)})
    do k=KBIG, KEIG

M4_IFELSE_2D({!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) &
!$OMP& PRIVATE(Exh,Eyh,Ezh,i,j) COLLAPSE(2)})
       do j=JBIG, JEIG

M4_IFELSE_1D({!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) PRIVATE(Exh,Eyh,Ezh)})
          do i=IBIG, IEIG

M4_IFELSE_TM({
             Exh=Ex(i,j,k)
             Eyh=Ey(i,j,k)
})
M4_IFELSE_TE({
             Ezh=Ez(i,j,k)
})

M4_IFELSE_TE({
             Hx(i,j,k) =  Hx(i,j,k) + M4_MUINVX(i,j,k) * ( &
M4_IFELSE_3D({       + DT/M4_HSZ(i,j,k)*( Ey(i,j,k+1) - Eyh ) &   },{})
M4_IFELSE_1D({0.&},{ - DT/M4_HSY(i,j,k)*( Ez(i,j+1,k) - Ezh ) &      })
                     )
             Hy(i,j,k) = Hy(i,j,k) + M4_MUINVY(i,j,k) * ( &
                    + DT/M4_HSX(i,j,k)*( Ez(i+1,j,k) - Ezh ) &
M4_IFELSE_3D({      - DT/M4_HSZ(i,j,k)*( Ex(i,j,k+1) - Exh ) &      })
                  )
})

M4_IFELSE_TM({
             Hz(i,j,k) = Hz(i,j,k) +  M4_MUINVZ(i,j,k) * ( &
                    - DT/M4_HSX(i,j,k)*( Ey(i+1,j,k) - Eyh ) &
M4_IFELSE_1D({},{   + DT/M4_HSY(i,j,k)*( Ex(i,j+1,k) - Exh ) &      })
				    )
})


          enddo
M4_IFELSE_1D({ !$OMP END PARALLEL DO })     

       enddo
M4_IFELSE_2D({ !$OMP END PARALLEL DO })     

    enddo
M4_IFELSE_3D({ !$OMP END PARALLEL DO })     
!!$OMP END PARALLEL DO

  end subroutine StepH

!----------------------------------------------------------------------

  subroutine StepE

    implicit none
    
    M4_FTYPE :: Hxh, Hyh, Hzh
    integer :: i, j, k

!!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) PRIVATE(Hxh,Hyh,Hzh) COLLAPSE(3)
! M4_DIM)
M4_IFELSE_3D({!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) &
!$OMP& PRIVATE(Hxh,Hyh,Hzh,i,j,k) COLLAPSE(3)})
   do k=KBIG, KEIG
   
M4_IFELSE_2D({!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) &
!$OMP& PRIVATE(Hxh,Hyh,Hzh,i,j) COLLAPSE(2)})
       do j=JBIG, JEIG

M4_IFELSE_1D({!$OMP PARALLEL DO SHARED(Hx,Hy,Hz,Ex,Ey,Ez,epsinvx,epsinvy,epsinvz) PRIVATE(Hxh,Hyh,Hzh)})
          do i=IBIG, IEIG

M4_IFELSE_TE({
             Hxh=Hx(i,j,k)
             Hyh=Hy(i,j,k)
})
M4_IFELSE_TM({
             Hzh=Hz(i,j,k) 
})

M4_IFELSE_TM({

             Ex(i,j,k) =  Ex(i,j,k) + epsinvx(i,j,k) *( &
M4_IFELSE_1D({0.&},{ + DT/M4_SY(i,j,k)*( Hzh - Hz(i,j-1,k) ) &             })
M4_IFELSE_3D({       - DT/M4_SZ(i,j,k)*( Hyh - Hy(i,j,k-1) ) &          },{})
                     )
             Ey(i,j,k) =  Ey(i,j,k) + epsinvy(i,j,k) * ( &
                     - DT/M4_SX(i,j,k)*( Hzh - Hz(i-1,j,k) ) &
M4_IFELSE_3D({       + DT/M4_SZ(i,j,k)*( Hxh - Hx(i,j,k-1) ) &             })
                    )
})

M4_IFELSE_TE({

             Ez(i,j,k) =  Ez(i,j,k) + epsinvz(i,j,k) * ( &
                     + DT/M4_SX(i,j,k)*( Hyh - Hy(i-1,j,k) )  &
M4_IFELSE_1D({},{    - DT/M4_SY(i,j,k)*( Hxh - Hx(i,j-1,k) )  &            })
                     )
})

          enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO })     

       enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO })     

    enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO })     
!!$OMP END PARALLEL DO

  end subroutine StepE

!----------------------------------------------------------------------

end module fdtd

! Authors:  S.Scholz, A.Klaedtke, J.Hamm 
! Modified: 1/5/2008
!
!======================================================================
