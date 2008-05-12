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

  integer :: numsteps, numfield, lot

  logical :: done, nohfield
  
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

    ! time window
    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. dn .lt. 1 ) then
       M4_PARSE_ERROR({error in time window specification})
    end if

    ! face normal vector for S.n projection
    call readfloats(funit,lcount,v,2) 
    diag%phi = v(1)
    diag%theta = v(2)
    
    ! assume |H| = |E| / n?
    call readlogical(funit,lcount,diag%nohfield)

    
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
    
    diag%numsteps = (diag%ne-diag%ns)/diag%dn + 1

    if ( mat%nohfield ) then

       diag%numfield = 3

    else

       diag%numfield = 6

    end if

    diag%lot = reg%numnodes * diag%numfield

    allocate(diag%field(1:diag%numsteps,1:reg%numnodes,diag%numfield)) ! test allocation (do we have the memory?)
    M4_ALLOC_ERROR(err,"InitializeDiagPSpec")

    deallocate(field)

    diag%npointer = 0
    
    diag%field = 0

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


    })
    M4_WRITE_DBG(". exit FinalizeDiagPSpec")

  end subroutine FinalizeDiagPSpec

!----------------------------------------------------------------------

  subroutine StepHDiagPSpec(ncyc)

    integer :: ncyc

    ! nop -> do everything in StepE!
  
  end subroutine StepHDiagPSpec


!----------------------------------------------------------------------


  subroutine StepEDiagPSpec(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({DIAGPSPEC},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({DIAGPSPEC},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)

       if ( diag%done ) cycle

       if ( ncyc .eq. diag%ns ) then

          allocate(field(0:diag%numsteps-1,1:reg%numnodes,diag%numfield))
          M4_ALLOC_ERROR(err,"StepEDiagPSpec")

          call InitializeFields

       end if

       if ( ncyc .ge. diag%ns .and. ncyc .le. diag%ne .and. &
            mod(ncyc-diag%ns,diag%dn) .eq. 0) then
          
          
          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
! record e and h field components averaged to yee-cell center


          diag%field(diag%npointer,p,1) = 0.5 * ( Ex(M4_COORD(i-1,j,k)) + Ex(M4_COORD(i,j,k)) )
          diag%field(diag%npointer,p,2) = 0.5 * ( Ey(M4_COORD(i,j-1,k)) + Ey(M4_COORD(i,j,k)) )
          diag%field(diag%npointer,p,3) = 0.5 * ( Ez(M4_COORD(i,j,k-1)) + Ez(M4_COORD(i,j,k)) )

          if ( diag%nohfield ) then

             diag%field(diag%npointer,p,4) = 0.25 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) + Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k-1)) )
             diag%field(diag%npointer,p,5) = 0.25 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i,j,k-1)) + Hy(M4_COORD(i-1,j,k-1)) )
             diag%field(diag%npointer,p,6) = 0.25 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) + Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) )
          
          end if

          diag%npointer = diag%npointer + 1
          
       })      

       end if

       if ( ncyc .ge. diag%ne .and. .not. diag%done ) then

          call FourierTransformFields

          call WriteSpectrum

          deallocate(field)

          diag%done = .true. ! we are done here

       end if


    })

  end subroutine StepEDiagPSpec

!----------------------------------------------------------------------

  subroutine InitializeFields

    integer :: ier = 0

    M4_IFELSE_CF({

    call CFFTMI(diag%numsteps, diag%wsave, diag%lensav, ier)
    M4_FATAL_ERROR({ier .ne. 0},{"CFFTMI failed!"}

    },{

    call RFFTMI(diag%numsteps, diag%wsave, diag%lensav, ier)
    M4_FATAL_ERROR({ier .ne. 0},{"RFFTMI failed!"}

    })

  end subroutine InitializeFields

!----------------------------------------------------------------------

  subroutine FourierTransformFields

    integer :: ier
    
    M4_IFELSE_CF({

! ---- perform complex fft
! SUBROUTINE CFFTMF (LOT, JUMP, N, INC, C, LENC, WSAVE, LENSAV, WORK, LENWRK, IER)
! INTEGER LOT, JUMP, N, INC, LENC, LENSAV, LENWRK, IER 
! COMPLEX C(LENC) 
! REAL WSAVE(LENSAV), WORK(LENWRK)

    call CFFTMB(diag%lot, diag%numsteps, diag%numsteps, 1, fields, diag%lot*diag%numsteps, diag%wsave, diag%lensav, diag%work, diag%lenwrk, ier)
    
    },{
! perform real fft    

    call RFFTMB(diag%lot, diag%numsteps, diag%numsteps, 1, fields, diag%lot*diag%numsteps, diag%wsave, diag%lensav, diag%work, diag%lenwrk, ier)
    
    })

  end subroutine FourierTransformFields


!----------------------------------------------------------------------

  subroutine WriteSpectrum

    



  end subroutine WriteSpectrum

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


