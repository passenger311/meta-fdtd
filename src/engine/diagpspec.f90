!-*- F90 -*------------------------------------------------------------
!
!  module: diagpspec / meta
!
!  Dummy diagnostics module.
!
!----------------------------------------------------------------------


! =====================================================================
!
! The DiagPSpec module calculates the power spectrum of the poynting
! vector over selected spatial points.


module diagpspec

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MODHEAD_DECL({DIAGPSPEC},100,{

  integer :: ns, ne, dn  ! time stepping 

  real(kind=8) :: theta, phi   ! angles of projection

  real(kind=8) :: kinc(3)                   ! normed k vector of plane wave
     
  character(len=80) :: fmt, filename

  integer :: numfreq, stride

  logical :: done
  
  integer :: npointer 


  M4_FTYPE, dimension(:,:) :: field

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagPSpecObj(funit,lcount)

    M4_MODREAD_DECL({DIAGPSPEC}, funit,lcount,diag,reg,out)
    integer :: v(3)
    real(kind=8) :: w(2)
    logical :: err

    M4_WRITE_DBG(". enter ReadMatPSpecObj")
    
    M4_MODREAD_EXPR({DIAGPSPEC}, funit,diag,reg,0,out, {

    err = .false.

    call readline(funit,lcount,eof,line)
    call getstring(line,diag%fmt,err)
    call getstring(line,diag%filename,err)
    M4_SYNTAX_ERROR({eof .or. err .or. line .ne. ""},lcount,"FMT FILENAME")
    M4_WRITE_DBG({"fmt filename: ",TRIM(diag%fmt)," ", TRIM(diag%filename)})


    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)

    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. dn .lt. 1 ) then
       M4_PARSE_ERROR({error in time window specification})
    end if

    call readfloats(funit,lcount,v,2) 
    diag%phi = v(1)
    diag%theta = v(2)
    

    })

    M4_WRITE_DBG(". exit ReadMatPSpecObj")

  end subroutine ReadDiagPSpecObj

!----------------------------------------------------------------------

  subroutine InitializeDiagPSpec

    M4_MODLOOP_DECL({DIAGPSPEC},diag)
    type (REG_T) :: reg
    integer :: err

    M4_WRITE_DBG(". enter InitializeMatPSpec")
    M4_MODLOOP_EXPR({DIAGPSPEC},diag,{
    
    M4_DIAGOBJ_GETREG(diag,reg)

    diag%kinc(1) = sin(DEG*diag%theta)*cos(DEG*diag%phi)
    diag%kinc(2) = sin(DEG*diag%theta)*sin(DEG*diag%phi)
    diag%kinc(3) = cos(DEG*diag%theta)
    
    diag%numfreq = (diag%ne-diag%ns)/diag%dn + 1

    diag%stride = reg%numnodes * 6

    allocate(field(1:diag%numfreq,1:reg%numnodes,6))
    M4_ALLOC_ERROR(err,"InitializeDiagPSpec")

    diag%npointer = 0
    
    diag%done = .false.

    M4_IFELSE_DBG({call EchoDiagPSpecObj(diag)}, stat = err)

    })
    M4_WRITE_DBG(". exit InitializeDiagPSpec")

  end subroutine InitializeDiagPSpec

!----------------------------------------------------------------------

  subroutine FinalizeDiagPSpec

    M4_MODLOOP_DECL({DIAGPSPEC},diag)
    M4_WRITE_DBG(". enter FinalizeDiagPSpec")
    M4_MODLOOP_EXPR({DIAGPSPEC},diag,{

       deallocate(field)

    })
    M4_WRITE_DBG(". exit FinalizeDiagPSpec")

  end subroutine FinalizeDiagPSpec

!----------------------------------------------------------------------

  subroutine StepHDiagPSpec(ncyc)

    integer :: ncyc

    ! nop
  
  end subroutine StepHDiagPSpec


!----------------------------------------------------------------------


  subroutine StepEDiagPSpec(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({DIAGPSPEC},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({DIAGPSPEC},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       if ( ncyc .ge. diag%ns .and. ncyc .le. diag%ne .and. &
            mod(ncyc-diag%ns,diag%dn) .eq. 0) then
          
! record e and h field components
          
          diag%npointer = diag%npointer + 1

          diag%fields(diag%npointer,p,1) = Ex(i,j,k)
          diag%fields(diag%npointer,p,2) = Ey(i,j,k)
          diag%fields(diag%npointer,p,3) = Ez(i,j,k)

          diag%fields(diag%npointer,p,4) = Hx(i,j,k)
          diag%fields(diag%npointer,p,5) = Hy(i,j,k)
          diag%fields(diag%npointer,p,6) = Hz(i,j,k)
          
       end if



       })      
    })

  end subroutine StepEDiagPSpec

!----------------------------------------------------------------------

   subroutine EchoDiagPSpecObj(diag)

    type(T_DIAGPSPEC) :: diag
 
    M4_WRITE_INFO({"--- diagpspec # ",&
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"something = ",diag%something })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))
    

  end subroutine EchoDiagPSpecObj
  
!----------------------------------------------------------------------

end module diagpspec

! =====================================================================


