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
  use parse
  use fdtd

  implicit none
  private
  save


  character(len=20), parameter :: sfx = ".pspec"


  M4_MODHEAD_DECL({DIAGPSPEC},100,{

  integer :: ns, ne, dn  ! time stepping 

  real(kind=8) :: theta, phi, psi   ! angles for face normal projection and polarizer

  real(kind=8) :: kinc(3), finc(6,2) ! normal plane vector and field components
     
  character(len=80) :: fmt, filename

  integer :: numsteps, numfield, lot

  logical :: done, nohfield
  
  integer :: npointer 

  real(kind=8), pointer, dimension(:,:,:) :: field

  ! for fourier transform

  real(kind=8), pointer, dimension(:) :: wsave
  integer :: lensav

  real(kind=8), pointer, dimension(:) :: work
  integer :: lenwrk

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagPSpecObj(funit,lcount)

    M4_MODREAD_DECL({DIAGPSPEC}, funit,lcount,diag,reg,out)
    integer :: v(3)
    real(kind=8) :: w(3)
    logical :: err,eof
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadMatPSpecObj")
    
    M4_MODREAD_EXPR({DIAGPSPEC}, funit,lcount,diag,reg,0,out, {

    err = .false.

    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR(eof,lcount)
    call getstring(line,diag%filename,err)
    M4_SYNTAX_ERROR({line .ne. ""},lcount,"FILENAME")

    M4_WRITE_DBG({"filename: ", TRIM(diag%filename)})

    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)

    ! time window
    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. diag%dn .lt. 1 ) then
       M4_FATAL_ERROR({"BAD TIME WINDOW!"})
    end if

    ! face normal vector for S.n projection and polarizer (psi)
    call readfloats(funit,lcount,w,3) 
    diag%phi = w(1)
    diag%theta = w(2)
    diag%psi = w(3)
    
    ! assume |H| = |E| / n?
    call readlogical(funit,lcount,diag%nohfield)

    
    })

    M4_WRITE_DBG(". exit ReadMatPSpecObj")

  end subroutine ReadDiagPSpecObj

!----------------------------------------------------------------------

  subroutine InitializeDiagPSpec

    M4_MODLOOP_DECL({DIAGPSPEC},diag)
    type (T_REG) :: reg
    integer :: ier
    real(kind=8) :: psi

    M4_WRITE_DBG(". enter InitializeMatPSpec")
    M4_MODLOOP_EXPR({DIAGPSPEC},diag,{
    
    M4_MODOBJ_GETREG(diag,reg)

    diag%kinc(1) = sin(DEG*diag%theta)*cos(DEG*diag%phi)
    diag%kinc(2) = sin(DEG*diag%theta)*sin(DEG*diag%phi)
    diag%kinc(3) = cos(DEG*diag%theta)
    
    ! polarisation 1
    psi = diag%psi
    diag%finc(1,1) = cos(DEG*psi)*sin(DEG*diag%phi) - sin(DEG*psi)*cos(DEG*diag%theta)*cos(DEG*diag%phi)
    diag%finc(2,1) = -cos(DEG*psi)*cos(DEG*diag%phi) - sin(DEG*psi)*cos(DEG*diag%theta)*sin(DEG*diag%phi)
    diag%finc(3,1) = sin(DEG*psi)*sin(DEG*diag%theta)
    diag%finc(4,1) = sin(DEG*psi)*sin(DEG*diag%phi) + cos(DEG*psi)*cos(DEG*diag%theta)*cos(DEG*diag%phi)
    diag%finc(5,1) = -sin(DEG*psi)*cos(DEG*diag%phi) + cos(DEG*psi)*cos(DEG*diag%theta)*sin(DEG*diag%phi)
    diag%finc(6,1) = -cos(DEG*psi)*sin(DEG*diag%theta)

    ! polarisation 2, orthogonal to polarisation 1
    psi = diag%psi + 90.
    diag%finc(1,2) = cos(DEG*psi)*sin(DEG*diag%phi) - sin(DEG*psi)*cos(DEG*diag%theta)*cos(DEG*diag%phi)
    diag%finc(2,2) = -cos(DEG*psi)*cos(DEG*diag%phi) - sin(DEG*psi)*cos(DEG*diag%theta)*sin(DEG*diag%phi)
    diag%finc(3,2) = sin(DEG*psi)*sin(DEG*diag%theta)
    diag%finc(4,2) = sin(DEG*psi)*sin(DEG*diag%phi) + cos(DEG*psi)*cos(DEG*diag%theta)*cos(DEG*diag%phi)
    diag%finc(5,2) = -sin(DEG*psi)*cos(DEG*diag%phi) + cos(DEG*psi)*cos(DEG*diag%theta)*sin(DEG*diag%phi)
    diag%finc(6,2) = -cos(DEG*psi)*sin(DEG*diag%theta)

    diag%numsteps = (diag%ne-diag%ns)/diag%dn + 1

    if ( diag%nohfield ) then

       diag%numfield = 2

    else

       diag%numfield = 4

    end if

    diag%lot = reg%numnodes * diag%numfield

    M4_WRITE_DBG({"kinc = ",diag%kinc(1),diag%kinc(2),diag%kinc(3)})

    M4_WRITE_INFO({"points x steps = ",TRIM(i2str(reg%numnodes))," x ",TRIM(i2str(diag%numsteps)) })

    allocate(diag%field(0:diag%numsteps-1,1:reg%numnodes,diag%numfield),stat=ier) ! test allocation (do we have the memory?)
    M4_ALLOC_ERROR({ier},"InitializeDiagPSpec")

    deallocate(diag%field)


    diag%npointer = 0
    
    diag%done = .false.

    M4_IFELSE_DBG({call EchoDiagPSpecObj(diag)})

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
    integer :: ier
    M4_MODLOOP_DECL({DIAGPSPEC},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))
    
    real(kind=8) :: Exh, Eyh, Ezh, Hxh, Hyh, Hzh, Ep1, Ep2, Hp1, Hp2

    M4_MODLOOP_EXPR({DIAGPSPEC},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)

       if ( diag%done ) cycle

       if ( ncyc .eq. diag%ns ) then

          M4_WRITE_INFO({"initializing fft #",TRIM(i2str(diag%idx))})

          allocate(diag%field(0:diag%numsteps-1,1:reg%numnodes,diag%numfield),stat=ier)
          M4_ALLOC_ERROR(ier,"StepEDiagPSpec")

          diag%field = 0.
          diag%npointer = 0

          call InitializeFields(diag)

       end if

       if ( ncyc .ge. diag%ns .and. ncyc .le. diag%ne .and. &
            mod(ncyc-diag%ns,diag%dn) .eq. 0) then
          
          
          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
          ! store E field projections
          
          Exh = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
          Eyh = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )
          Ezh = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )

          Ep1 = diag%finc(1,1)*Exh + diag%finc(2,1)*Eyh + diag%finc(3,1)*Ezh
          Ep2 = diag%finc(1,2)*Exh + diag%finc(2,2)*Eyh + diag%finc(3,2)*Ezh
          
          diag%field(diag%npointer,p,1) = Ep1
          diag%field(diag%npointer,p,2) = Ep2

          if ( .not. diag%nohfield ) then ! store H field projections

             Hxh = real( 0.25 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) + Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k-1)) ) )
             Hyh = real( 0.25 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i,j,k-1)) + Hy(M4_COORD(i-1,j,k-1)) ) )
             Hzh = real( 0.25 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) + Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )

             Hp1 =  diag%finc(4,1)*Hxh + diag%finc(5,1)*Hyh + diag%finc(6,1)*Hzh
             Hp2 =  diag%finc(4,2)*Hxh + diag%finc(5,2)*Hyh + diag%finc(6,2)*Hzh

             diag%field(diag%npointer,p,3) = Hp1
             diag%field(diag%npointer,p,4) = Hp2
          
          end if
          
          })      

          diag%npointer = diag%npointer + 1

       end if

       if ( ncyc .ge. diag%ne .and. .not. diag%done ) then

          call FourierTransformFields(diag)

          call WriteSpectrum(diag)

          deallocate(diag%field)

          diag%done = .true. ! we are done here

       end if


    })

  end subroutine StepEDiagPSpec

!----------------------------------------------------------------------

  subroutine InitializeFields(diag)

    type(T_DIAGPSPEC) :: diag

    integer :: ier = 0

    diag%lensav = 2 * diag%numsteps + int(log(real(diag%numsteps))) + 4 

    ! allocate and initialize wsave array for fft

    allocate(diag%wsave(1:diag%lensav), stat = ier)
    M4_ALLOC_ERROR({ier},{"InitializeFields"})

    ! call fftpack5 real initialization

    call RFFTMI(diag%numsteps, diag%wsave, diag%lensav, ier)
    if ( ier .ne. 0 ) then 
       M4_FATAL_ERROR({"RFFTMI FAILED!"})
    end if

    ! allocate work array for fft

    diag%lenwrk = diag%numsteps * diag%lot

    allocate(diag%work(1:diag%lenwrk), stat = ier)
    M4_ALLOC_ERROR({ier},{"InitializeFields"})


  end subroutine InitializeFields

!----------------------------------------------------------------------

  subroutine FourierTransformFields(diag)
    
    type(T_DIAGPSPEC) :: diag

    integer :: ier = 0

    M4_WRITE_INFO({"performing fft #",TRIM(i2str(diag%idx))})
    

    call RFFTMF(diag%lot, diag%numsteps, diag%numsteps, 1, diag%field, &
         diag%lot*diag%numsteps, diag%wsave, diag%lensav, diag%work, diag%lenwrk, ier)
    if ( ier .ne. 0 ) then
       M4_FATAL_ERROR({"RFFTMF failed!"})
    end if

    deallocate(diag%wsave)
    deallocate(diag%work)

  end subroutine FourierTransformFields


!----------------------------------------------------------------------

  subroutine WriteSpectrum(diag)

    type(T_DIAGPSPEC) :: diag

    type(T_REG) :: reg
    character(len=STRLNG) :: fn
    integer :: w, ios,i,j,k,p

    real(kind=8) :: nrefr
    real(kind=8) :: Ep1c, Ep2c, Hp1c, Hp2c,  Ep1s, Ep2s, Hp1s, Hp2s
    real(kind=8) :: SumUp1, SumUp2
    integer :: nh

    fn = cat2(diag%filename,sfx)
    
    M4_WRITE_INFO({"writing spectrum #",TRIM(i2str(diag%idx)), " -> ", TRIM(fn)})

    M4_MODOBJ_GETREG(diag,reg)

    open(UNITTMP,FILE=fn,STATUS="unknown", IOSTAT=ios)
    M4_OPEN_ERROR(ios,fn)

    write(UNITTMP,*) "! PSPEC"
    write(UNITTMP,*) "! N: ",TRIM(i2str(diag%ns))," ",TRIM(i2str(diag%ne))," ",TRIM(i2str(diag%dn))

    if ( mod(diag%numsteps,2) .eq. 0 ) then 
       nh = diag%numsteps/2 - 1
    else
       nh = (diag%numsteps-1)/2
    end if

    write(UNITTMP,*) "! F: 1 ",TRIM(i2str(nh))

    do w = 1, nh 

       SumUp1 = 0.
       SumUp2 = 0.

! integrate power flux over spatial area       
       do p = 1, reg%numnodes

          Ep1c = diag%field(2*w,p,1)
          Ep2c = diag%field(2*w,p,2)

          Ep1s = diag%field(2*w+1,p,1)
          Ep2s = diag%field(2*w+1,p,2)

          if ( .not. diag%nohfield ) then
            
             Hp1c = diag%field(2*w,p,3)
             Hp2c = diag%field(2*w,p,4)

             Hp1s = diag%field(2*w+1,p,3)
             Hp2s = diag%field(2*w+1,p,4)

          else

             ! that's not quite correct, but does it should not matter too much!

             i = reg%i(p)
             j = reg%j(p)
             k = reg%k(p)

             nrefr = 1./sqrt(epsinvx(M4_COORD(i,j,k))*M4_MUINVX(i,j,k))

             Hp1c = Ep1c/nrefr
             Hp2c = Ep2c/nrefr

             Hp1s = Ep1s/nrefr
             Hp2s = Ep2s/nrefr
             
          end if

          SumUp1 = SumUp1 + Ep1c * Hp1c + Ep1s * Hp1s
          SumUp2 = SumUp2 + Ep2c * Hp2c + Ep1s * Hp1s

       end do

       write(UNITTMP,*) w, SumUp1, SumUp2

    end do

    close(UNITTMP)


  end subroutine WriteSpectrum

!----------------------------------------------------------------------

   subroutine EchoDiagPSpecObj(diag)

    type(T_DIAGPSPEC) :: diag
 
    M4_WRITE_INFO({"--- diagpspec # ",&
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))
    

  end subroutine EchoDiagPSpecObj
  
!----------------------------------------------------------------------

end module diagpspec

! =====================================================================


