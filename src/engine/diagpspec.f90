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

  integer :: ns, ne, dn, nh  ! time stepping 

  real(kind=8) :: theta, phi, psi   ! angles for face normal projection and polarizer

  real(kind=8) :: kinc(3), finc(6,2) ! normal plane vector and field components
     
  character(len=80) :: mode, filename, refname

  integer :: phasefw, phasebw

  integer :: numsteps, numfield, lot, numcomp

  logical :: done
  
  integer :: npointer 

  real(kind=4), pointer, dimension(:,:,:) :: field

  ! for fourier transform

  real(kind=4), pointer, dimension(:) :: wsave
  integer :: lensav

  real(kind=4), pointer, dimension(:) :: work
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

    call readstring(funit,lcount,diag%filename)
    M4_WRITE_DBG({"filename: ", TRIM(diag%filename)})

    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR(eof,lcount)
    call getstring(line,diag%mode,err)

    select case ( diag%mode ) 
    case ( "S" )
       diag%numcomp = 2
    case ( "Ecs" ) 
       diag%numcomp = 4
    case ( "Hcs" ) 
       diag%numcomp = 4
    case ( "Eap" ) 
       diag%numcomp = 4
    case ( "Hap" ) 
       diag%numcomp = 4
    case default
        M4_SYNTAX_ERROR(.true.,lcount,"S | Ecs | Eap | Hcs | Hap")
    end select
         
    M4_WRITE_DBG({"mode: ", TRIM(diag%mode)})

    call getstring(line,diag%refname,err)
    if ( err ) diag%refname = "not_defined"

    M4_SYNTAX_ERROR({line .ne. ""},lcount,"MODE [REF.FILENAME]")

    M4_WRITE_DBG({"ref. filename: ", TRIM(diag%refname)})

    call readints(funit,lcount,v,2)
    diag%phasefw = v(1)
    diag%phasebw = v(2)
    if ( diag%phasefw  .ne. 0 ) diag%phasefw = 1
    if ( diag%phasebw  .ne. 0 ) diag%phasebw = 1

    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)

    diag%ns = max(diag%ns,0)
    diag%ne = min(diag%ne,NCYCMAX)

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
    ! call readlogical(funit,lcount,diag%nohfield)

    
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


    if ( mod(diag%numsteps,2) .eq. 0 ) then 
       diag%nh = diag%numsteps/2 - 1
    else
       diag%nh = (diag%numsteps-1)/2
    end if

    diag%numfield = 4

    diag%lot = reg%numnodes * diag%numfield

    M4_WRITE_DBG({"kinc = ",diag%kinc(1),diag%kinc(2),diag%kinc(3)})

    M4_WRITE_INFO({"comp x points x steps = ",TRIM(i2str(diag%numfield))," x ",TRIM(i2str(reg%numnodes))," x ",TRIM(i2str(diag%numsteps)) })

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

          M4_WRITE_INFO({"recording fft #",TRIM(i2str(diag%idx))})

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

          Hxh = real( 0.25 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) + Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k-1)) ) )
          Hyh = real( 0.25 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i,j,k-1)) + Hy(M4_COORD(i-1,j,k-1)) ) )
          Hzh = real( 0.25 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) + Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )

          Hp1 =  diag%finc(4,1)*Hxh + diag%finc(5,1)*Hyh + diag%finc(6,1)*Hzh
          Hp2 =  diag%finc(4,2)*Hxh + diag%finc(5,2)*Hyh + diag%finc(6,2)*Hzh

          diag%field(diag%npointer,p,3) = Hp1
          diag%field(diag%npointer,p,4) = Hp2
          
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
       M4_FATAL_ERROR({"DFFTMI FAILED!"})
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

    character(len=STRLNG) :: fn
    integer :: l, m, ios
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    real(kind=8) :: nrefr, df
    real(kind=8) :: Ep1c, Ep2c, Hp1c, Hp2c,  Ep1s, Ep2s, Hp1s, Hp2s
    real(kind=8) :: SumS1, SumS2, SumE1c, SumE1s, SumE2c, SumE2s, SumH1c, SumH1s, SumH2c, SumH2s   
    real(kind=8) :: SumEa1,SumEa2,SumHa1,SumHa2,SumEph1,SumEph2,SumHph1,SumHph2
    real(kind=8) :: SumEph1_old,SumEph2_old,SumHph1_old,SumHph2_old
    real(kind=8) :: norm, freq, j1, j2, ph1, ph2, ph1o, ph2o
    real(kind=8) :: rfreq, rv(6), d1, d2, phaseoffs
    logical :: hasref

    integer :: nh

    ! increase phase unwrap sensitivity by factor 2 if either phasefw or phasebw is set (but not both)
    if ( diag%phasefw + diag%phasebw .eq. 1 ) then 
       phaseoffs = 0.
    else
       phaseoffs = PI
    end if

    fn = cat2(diag%filename,sfx)
    
    M4_WRITE_INFO({"writing spectrum #",TRIM(i2str(diag%idx)), " -> ", TRIM(fn)})

    M4_MODOBJ_GETREG(diag,reg)

    df = 1./((diag%ne-diag%ns+1)*DT)

    open(UNITTMP,FILE=fn,STATUS="unknown", IOSTAT=ios)
    M4_OPEN_ERROR(ios,fn)

    write(UNITTMP,*) "# GPL: PSPEC"
    write(UNITTMP,*) "# ",TRIM(diag%mode)," ! mode"
    write(UNITTMP,*) "# ",TRIM(i2str(diag%phasefw))," ",TRIM(i2str(diag%phasefw))," ! phase unwrap fw/bw"
    write(UNITTMP,*) "# ",TRIM(i2str(diag%ns))," ",TRIM(i2str(diag%ne))," ",TRIM(i2str(diag%dn)), " ! tframe"
    write(UNITTMP,*) "# ", DT, " ! dt"
    write(UNITTMP,*) "# 1 ",TRIM(i2str(diag%nh))," ! fframe"
    write(UNITTMP,*) "# ", df, " ! df"

    call OpenRefFile(UNITTMP+1, hasref)

    j1 = 0.
    j2 = 0.

    SumHph1_old = 0.
    SumHph2_old = 0.
    SumEph1_old = 0.
    SumEph2_old = 0.

    select case ( diag%mode ) 
    case( "Ecs" )
       rv = 1.
    case( "Hcs" )
       rv = 1.
    case( "Eap" )
       rv = 1.
       rv(2) = 0.
       rv(4) = 0.
    case( "Hap" ) 
       rv = 1.
       rv(2) = 0.
       rv(4) = 0.
    case ("S")
       rv = 1.
    end select
    
    do l = 1, diag%nh

       SumS1 = 0.
       SumS2 = 0.

       SumE1c = 0.
       SumE1s = 0.
       SumE2c = 0.
       SumE2s = 0.

       SumH1c = 0.
       SumH1s = 0.
       SumH2c = 0.
       SumH2s = 0.

       SumEa1 = 0.
       SumEa2 = 0.
       SumHa1 = 0.
       SumHa2 = 0.

       SumEph1 = 0.
       SumEph2 = 0.
       SumHph1 = 0.
       SumHph2 = 0.

       freq = l*df

! integrate power flux over spatial area
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
          Ep1c = diag%field(2*l-1,p,1) * cos(PI*freq*DT) - diag%field(2*l,p,1) * sin(PI*freq*DT)
          Ep2c = diag%field(2*l-1,p,2) * cos(PI*freq*DT) - diag%field(2*l,p,2) * sin(PI*freq*DT)

          Ep1s = diag%field(2*l-1,p,1) * sin(PI*freq*DT) + diag%field(2*l,p,1) * cos(PI*freq*DT)
          Ep2s = diag%field(2*l-1,p,2) * sin(PI*freq*DT) + diag%field(2*l,p,2) * cos(PI*freq*DT)

          Hp1c = diag%field(2*l-1,p,3)
          Hp2c = diag%field(2*l-1,p,4)

          Hp1s = diag%field(2*l,p,3)
          Hp2s = diag%field(2*l,p,4)

          SumS1 = SumS1 + Ep1c * Hp1c + Ep1s * Hp1s
          SumS2 = SumS2 + Ep2c * Hp2c + Ep2s * Hp2s

          SumE1c = SumE1c + Ep1c 
          SumE1s = SumE1s + Ep1s
          SumE2c = SumE2c + Ep2c 
          SumE2s = SumE2s + Ep2s

          SumH1c = SumH1c + Hp1c 
          SumH1s = SumH1s + Hp1s
          SumH2c = SumH2c + Hp2c 
          SumH2s = SumH2s + Hp2s

          SumEa1 = SumEa1 + Ep1c * Ep1c + Ep1s * Ep1s
          SumEa2 = SumEa2 + Ep2c * Ep2c + Ep2s * Ep2s

          SumHa1 = SumHa1 + Hp1c * Hp1c + Hp1s * Hp1s
          SumHa2 = SumHa2 + Hp2c * Hp2c + Hp2s * Hp2s

          SumEph1 = SumEph1 + atan2(Ep1s,Ep1c)
          SumEph2 = SumEph2 + atan2(Ep2s,Ep2c)

          SumHph1 = SumHph1 + atan2(Hp1s,Hp1c)
          SumHph2 = SumHph2 + atan2(Hp2s,Hp2c)

       })
       norm = reg%numnodes

       if ( hasref ) then
          read(UNITTMP+1,*,iostat=ios) rfreq, (rv(m), m=1, diag%numcomp, 1 )
          if ( ios .ne. 0 ) then 
             M4_WRITE_WARN({"eof in reference file!"})
             close (UNITTMP+1)
             hasref = .false.
          end if
       end if

       select case ( diag%mode ) 
       case( "Ecs" )
          write(UNITTMP,*) freq, SumE1c/norm/rv(1), SumE1s/norm/rv(2), SumE2c/norm/rv(3), SumE2s/norm/rv(4)
       case( "Hcs" ) 
          write(UNITTMP,*) freq, SumH1c/norm/rv(1), SumH1s/norm/rv(2), SumH2c/norm/rv(3), SumH2s/norm/rv(4)
       case( "Eap" )
          SumEph1 = SumEph1/norm
          SumEph2 = SumEph2/norm
          d1 = SumEph1 - SumEph1_old
          d2 = SumEph2 - SumEph2_old
          if ( diag%phasefw .ne. 0 .and. d1 .lt. -phaseoffs ) j1 = j1 + 2.*PI
          if ( diag%phasebw .ne. 0 .and. d1 .gt. phaseoffs ) j1 = j1 - 2.*PI
          if ( diag%phasefw .ne. 0 .and. d2 .lt. -phaseoffs ) j2 = j2 + 2.*PI
          if ( diag%phasebw .ne. 0 .and. d2 .gt. phaseoffs ) j2 = j2 - 2.*PI
          SumEph1_old = SumEph1
          SumEph2_old = SumEph2
          SumEph1 = SumEph1 + PI + j1
          SumEph2 = SumEph2 + PI + j2
          if ( hasref ) then
             SumEph1 = SumEph1 - rv(2)
             SumEph2 = SumEph1 - rv(4)
          end if
          write(UNITTMP,*) freq, sqrt(SumEa1/norm)/rv(1), SumEph1-PI,  &
               sqrt(SumEa2/norm)/rv(3),SumEph2-PI
       case( "Hap" ) 
          SumHph1 = SumHph1/norm
          SumHph2 = SumHph2/norm
          d1 = SumHph1 - SumHph1_old
          d2 = SumHph2 - SumHph2_old
          if ( diag%phasefw .ne. 0 .and. d1 .lt. -phaseoffs ) j1 = j1 + 2.*PI
          if ( diag%phasebw .ne. 0 .and. d1 .gt. phaseoffs ) j1 = j1 - 2.*PI
          if ( diag%phasefw .ne. 0 .and. d2 .lt. -phaseoffs ) j2 = j2 + 2.*PI
          if ( diag%phasebw .ne. 0 .and. d2 .gt. phaseoffs ) j2 = j2 - 2.*PI
          SumHph1_old = SumHph1
          SumHph2_old = SumHph2
          SumHph1 = SumHph1 + PI + j1
          SumHph2 = SumHph2 + PI + j2
          if ( hasref ) then
             SumHph1 = SumHph1 - rv(2)
             SumHph2 = SumHph1 - rv(4)
          end if
          write(UNITTMP,*) freq, sqrt(SumHa1)/norm/rv(1), SumHph1-PI, &
               sqrt(SumHa2)/norm/rv(3), SumHph2-PI
       case ("S")
          write(UNITTMP,*) freq, SumS1/(norm*rv(1)), SumS2/(norm*rv(2))
       end select
          

    end do

    close(UNITTMP)

    if ( hasref ) close(UNITTMP+1)

  contains
    
    subroutine OpenRefFile(unit, hasref)

      integer :: unit
      logical :: hasref
      character(len=STRLNG) :: fn
      integer :: ios
      character(LEN=1) :: mark
      character(LEN=10) :: str
      integer :: ns, ne, dn, fs, fe, phfw, phbw
      real(kind=8) :: dt1, df

      hasref = .false.

      if ( diag%refname .eq. diag%filename ) return

      fn = cat2(diag%refname,sfx)

      open(unit,FILE=fn,STATUS="old", IOSTAT=ios)

      if ( ios .ne. 0 ) then 
         M4_WRITE_WARN({"could not open reference spectrum: ", TRIM(fn),"!"})
         return
      end if
         
      hasref = .true.

      ! compare whether header of reference file matches
      read(unit,*) mark, str
      read(unit,*) mark, str
      if ( str .ne. diag%mode ) hasref = .false. 
      read(unit,*) mark, phfw, phbw 
      read(unit,*) mark, ns, ne, dn 
      if ( (ne-ns) .ne. (diag%ne-diag%ns) .or. dn .ne. diag%dn ) hasref = .false. 
      read(unit,*) mark, dt1
      if ( abs(dt1-dt)/dt .gt. 1.e-5 ) hasref = .false.  
      read(unit,*) mark, fs, fe 
      read(unit,*) mark, df

      if ( .not. hasref ) then 
         M4_WRITE_WARN({"incompatible reference spectrum: ", TRIM(fn),"!"})
         close(unit)
         return
      end if

      M4_WRITE_INFO({"using reference spectrum: ", TRIM(fn)})
     
    end subroutine OpenRefFile

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


