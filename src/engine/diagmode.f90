!-*- F90 -*------------------------------------------------------------
!
!  module: diagmode / meta
!
!  Dummy diagnostics module.
!
!----------------------------------------------------------------------


! =====================================================================
!
! The DiagMode module calculates the power spectrum of the poynting
! vector over selected spatial points.


module diagmode

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

  M4_MODHEAD_DECL({DIAGMODE},100,{

  integer :: ns, ne, dn  ! time stepping 

  integer :: numsteps

  character(len=80) :: filename, outfile, mode

  logical :: done

  ! for DFT
  
  real(kind=8), pointer, dimension(:) :: freqs
  integer :: numfreqs

  real(kind=8), pointer, dimension(:,:,:) :: Fcos
  real(kind=8), pointer, dimension(:,:,:) :: Fsin

  ! for tangential DFT
  integer :: numcomp, face

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagModeObj(funit,lcount)

    M4_MODREAD_DECL({DIAGMODE}, funit,lcount,diag,reg,out)
    integer :: v(3)
    logical :: err,eof
    character(len=LINELNG) :: line

    M4_WRITE_DBG(". enter ReadDiagModeObj")
    
    M4_MODREAD_EXPR({DIAGMODE}, funit,lcount,diag,reg,0,out, {

    err = .false.

    call readstring(funit,lcount,diag%filename)
    M4_WRITE_DBG({"filename: ", TRIM(diag%filename)})

    call readstring(funit,lcount,diag%outfile)
    M4_WRITE_DBG({"outfile: ", TRIM(diag%outfile)})

    call readline(funit,lcount,eof,line)
    M4_EOF_ERROR(eof,lcount)
    call getstring(line,diag%mode,err)
         
    M4_WRITE_DBG({"mode: ", TRIM(diag%mode)})

    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)

    diag%ns = max(diag%ns,0)
    diag%ne = min(diag%ne,NCYCMAX)

    ! time window
    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. diag%dn .lt. 0 ) then
       M4_FATAL_ERROR({"BAD TIME WINDOW!"})
    end if
    
    })

    M4_WRITE_DBG(". exit ReadDiagModeObj")

  end subroutine ReadDiagModeObj

!----------------------------------------------------------------------

  subroutine InitializeDiagMode

    M4_MODLOOP_DECL({DIAGMODE},diag)
    type (T_REG) :: reg
    integer :: ier
    integer :: ios, i
    real(kind=8) :: freqs(1000),maxfreq

    M4_WRITE_DBG(". enter InitializeDiagMode")
    M4_MODLOOP_EXPR({DIAGMODE},diag,{
    
    M4_MODOBJ_GETREG(diag,reg)

    open(UNITTMP,FILE=diag%filename,STATUS="old", IOSTAT=ios)

    if ( ios .ne. 0 ) then 
       M4_WRITE_WARN({"could not open frequency file: ",  TRIM(diag%filename),"!"})
       return
    end if

    diag%numfreqs = 0
    ios = 0

    do while ( ios .eq. 0 )
       diag%numfreqs = diag%numfreqs + 1
       read(UNITTMP,*,iostat=ios) freqs(diag%numfreqs)
    end do

    close(UNITTMP)

    diag%numfreqs = diag%numfreqs - 1

    M4_WRITE_INFO({"recording dft for ",TRIM(i2str(diag%numfreqs))," frequencies"})

    allocate(diag%freqs(diag%numfreqs))

    maxfreq = 0.
    do i = 1, diag%numfreqs
       diag%freqs(i) = freqs(i)
       if ( freqs(i) .gt. maxfreq ) then
          maxfreq = freqs(i)
       end if
    end do

    diag%numsteps = (diag%ne-diag%ns)/diag%dn + 1

    diag%done = .false.

    !check if box is a plane and if so which one

    if ( reg%is .eq. reg%ie ) then
       diag%face = 1
    else if ( reg%js .eq. reg%je ) then
       diag%face = 2
    else if ( reg%ks .eq. reg%ke ) then
       diag%face = 3
    else if ( diag%mode .eq. "EHTG" .or. diag%mode .eq. "EHN" .or. diag%mode .eq. "EHT" ) then
       diag%mode = "F"
       M4_WRITE_WARN({"Mode changed to: ",  TRIM(diag%mode),"!"})
    end if

    !set the number of components needed to evaluate DFT in volume or on plane only

    if ( diag%mode .eq. "EHN" ) then
       diag%numcomp = 2
    elseif ( diag%mode .eq. "EHT" ) then
       diag%numcomp = 4
    else
       diag%numcomp = 6
    end if

    M4_IFELSE_DBG({call EchoDiagModeObj(diag)})

    })
    M4_WRITE_DBG(". exit InitializeDiagMode")

  end subroutine InitializeDiagMode

!----------------------------------------------------------------------

  subroutine FinalizeDiagMode

    M4_MODLOOP_DECL({DIAGMODE},diag)
    M4_WRITE_DBG(". enter FinalizeDiagMode")
    M4_MODLOOP_EXPR({DIAGMODE},diag,{

       deallocate(diag%freqs)

    })

    M4_WRITE_DBG(". exit FinalizeDiagMode")

  end subroutine FinalizeDiagMode

!----------------------------------------------------------------------

  subroutine StepHDiagMode(ncyc)

    integer :: ncyc

    ! nop -> do everything in StepE!
  
  end subroutine StepHDiagMode


!----------------------------------------------------------------------


  subroutine StepEDiagMode(ncyc)

    integer :: ncyc
    integer :: ier, c, l
    M4_MODLOOP_DECL({DIAGMODE},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    real(kind=8) :: Exh, Eyh, Ezh, Hxh, Hyh, Hzh, phi, cosfac, sinfac
    ! tangential field projections
    real(kind=8) :: Ep1, Ep2, Hp1, Hp2, Hp12, Hp22

    M4_MODLOOP_EXPR({DIAGMODE},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)

       if ( diag%done ) cycle

       if ( ncyc .eq. diag%ns ) then

          M4_WRITE_INFO({"recording dft #",TRIM(i2str(diag%idx))})

          allocate(diag%Fcos(diag%numfreqs,reg%numnodes,diag%numcomp),stat=ier) 
          M4_ALLOC_ERROR({ier},"StepEDiagMode")

          allocate(diag%Fsin(diag%numfreqs,reg%numnodes,diag%numcomp),stat=ier) 
          M4_ALLOC_ERROR({ier},"StepEDiagMode")

          diag%Fcos = 0.
          diag%Fsin = 0.

       end if

       if ( ncyc .ge. diag%ns .and. ncyc .le. diag%ne .and. &
            mod(ncyc-diag%ns,diag%dn) .eq. 0) then


          if ( diag%mode .eq. "EHN" ) then

             ! store E and H field projections

             select case ( diag%face )

             case( 1 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
                Hp1 = real( 0.25 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k)) + Hx(M4_COORD(i,j-1,k-1)) ) )

                do c = 1, diag%numfreqs 

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Hp1 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Hp1 * sinfac

                end do

                })

             case( 2 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )
                Hp1 = real( 0.25 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i,j,k-1)) + Hy(M4_COORD(i-1,j,k-1)) ) )

                do c = 1, diag%numfreqs 

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Hp1 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Hp1 * sinfac

                end do

                })

             case( 3 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )
                Hp1 = real( 0.25 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) + Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )

                do c = 1, diag%numfreqs 

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Hp1 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Hp1 * sinfac

                end do

                })

             end select

          elseif ( diag%mode .eq. "EHT" ) then

             ! store E and H field projections

             select case ( diag%face )

             case( 1 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )
                Ep2 = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )

                Hp1 = real( 0.5 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i,j,k-1)) ) )
                Hp2 = real( 0.5 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i,j-1,k)) ) )

                do c = 1, diag%numfreqs

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Ep2 * cosfac
                   diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Hp1 * cosfac
                   diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hp2 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Ep2 * sinfac
                   diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Hp1 * sinfac
                   diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hp2 * sinfac

                end do

                })

             case( 2 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
                Ep2 = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )

                Hp1 = real( 0.5 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j,k-1)) ) )
                Hp2 = real( 0.5 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) ) )

                do c = 1, diag%numfreqs

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Ep2 * cosfac
                   diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Hp1 * cosfac
                   diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hp2 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Ep2 * sinfac
                   diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Hp1 * sinfac
                   diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hp2 * sinfac

                end do

                })

             case( 3 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
                Ep2 = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )

                Hp1 = real( 0.5 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) ) )
                Hp2 = real( 0.5 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) ) )

                do c = 1, diag%numfreqs

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Ep2 * cosfac
                   diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Hp1 * cosfac
                   diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hp2 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Ep2 * sinfac
                   diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Hp1 * sinfac
                   diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hp2 * sinfac

                end do

                })

             end select

          elseif ( diag%mode .eq. "EHTG" ) then

             ! store E and H field projections

             select case ( diag%face )

             case( 1 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )
                Ep2 = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )

                Hp1 = real( 0.5 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i,j,k-1)) ) )
                Hp12 = real( 0.5 * ( Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i-1,j,k-1)) ) )
                Hp2 = real( 0.5 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i,j-1,k)) ) )
                Hp22 = real( 0.5 * ( Hz(M4_COORD(i-1,j,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )

                do c = 1, diag%numfreqs 

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Ep2 * cosfac
                   diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Hp1 * cosfac
                   diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hp2 * cosfac
                   diag%Fcos(c,p,5) = diag%Fcos(c,p,5) + Hp12 * cosfac
                   diag%Fcos(c,p,6) = diag%Fcos(c,p,6) + Hp22 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Ep2 * sinfac
                   diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Hp1 * sinfac
                   diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hp2 * sinfac
                   diag%Fsin(c,p,5) = diag%Fsin(c,p,5) + Hp12 * sinfac
                   diag%Fsin(c,p,6) = diag%Fsin(c,p,6) + Hp22 * sinfac

                end do

                })

             case( 2 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
                Ep2 = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )

                Hp1 = real( 0.5 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j,k-1)) ) )
                Hp12 = real( 0.5 * ( Hx(M4_COORD(i,j-1,k)) + Hx(M4_COORD(i,j-1,k-1)) ) )
                Hp2 = real( 0.5 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) ) )
                Hp22 = real( 0.5 * ( Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )

                do c = 1, diag%numfreqs 

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Ep2 * cosfac
                   diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Hp1 * cosfac
                   diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hp2 * cosfac
                   diag%Fcos(c,p,5) = diag%Fcos(c,p,5) + Hp12 * cosfac
                   diag%Fcos(c,p,6) = diag%Fcos(c,p,6) + Hp22 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Ep2 * sinfac
                   diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Hp1 * sinfac
                   diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hp2 * sinfac
                   diag%Fsin(c,p,5) = diag%Fsin(c,p,5) + Hp12 * sinfac
                   diag%Fsin(c,p,6) = diag%Fsin(c,p,6) + Hp22 * sinfac

                end do

                })

             case( 3 )

                M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

                Ep1 = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
                Ep2 = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )

                Hp1 = real( 0.5 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) ) )
                Hp12 = real( 0.5 * ( Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k-1)) ) )
                Hp2 = real( 0.5 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) ) )
                Hp22 = real( 0.5 * ( Hy(M4_COORD(i,j,k-1)) + Hy(M4_COORD(i-1,j,k-1)) ) )

                do c = 1, diag%numfreqs 

                   phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                   cosfac = cos( phi )
                   sinfac = sin( phi )

                   diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Ep1 * cosfac
                   diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Ep2 * cosfac
                   diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Hp1 * cosfac
                   diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hp2 * cosfac
                   diag%Fcos(c,p,5) = diag%Fcos(c,p,5) + Hp12 * cosfac
                   diag%Fcos(c,p,6) = diag%Fcos(c,p,6) + Hp22 * cosfac

                   diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Ep1 * sinfac
                   diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Ep2 * sinfac
                   diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Hp1 * sinfac
                   diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hp2 * sinfac
                   diag%Fsin(c,p,5) = diag%Fsin(c,p,5) + Hp12 * sinfac
                   diag%Fsin(c,p,6) = diag%Fsin(c,p,6) + Hp22 * sinfac

                end do

                })

             end select

          else

             M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

             ! store E and H field projections

             Exh = real( 0.5 * ( Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) )
             Eyh = real( 0.5 * ( Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) )
             Ezh = real( 0.5 * ( Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) )

             Hxh = real( 0.25 * ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) + Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k-1)) ) )
             Hyh = real( 0.25 * ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i,j,k-1)) + Hy(M4_COORD(i-1,j,k-1)) ) )
             Hzh = real( 0.25 * ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) + Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )

             do c = 1, diag%numfreqs 

                phi = 2. * PI * diag%freqs(c) * ( ncyc - diag%ns ) * DT

                cosfac = cos( phi )
                sinfac = sin( phi )

                diag%Fcos(c,p,1) = diag%Fcos(c,p,1) + Exh * cosfac
                diag%Fcos(c,p,2) = diag%Fcos(c,p,2) + Eyh * cosfac
                diag%Fcos(c,p,3) = diag%Fcos(c,p,3) + Ezh * cosfac
                diag%Fcos(c,p,4) = diag%Fcos(c,p,4) + Hxh * cosfac
                diag%Fcos(c,p,5) = diag%Fcos(c,p,5) + Hyh * cosfac
                diag%Fcos(c,p,6) = diag%Fcos(c,p,6) + Hzh * cosfac

                diag%Fsin(c,p,1) = diag%Fsin(c,p,1) + Exh * sinfac
                diag%Fsin(c,p,2) = diag%Fsin(c,p,2) + Eyh * sinfac
                diag%Fsin(c,p,3) = diag%Fsin(c,p,3) + Ezh * sinfac
                diag%Fsin(c,p,4) = diag%Fsin(c,p,4) + Hxh * sinfac
                diag%Fsin(c,p,5) = diag%Fsin(c,p,5) + Hyh * sinfac
                diag%Fsin(c,p,6) = diag%Fsin(c,p,6) + Hzh * sinfac

             end do

             })

          end if

       end if

       if ( ncyc .ge. diag%ne .and. .not. diag%done ) then

          do c = 1, diag%numfreqs

             M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

             do l = 1, diag%numcomp

                diag%Fcos(c,p,l) = 2. / diag%numsteps * diag%Fcos(c,p,l)
                diag%Fsin(c,p,l) = 2. / diag%numsteps * diag%Fsin(c,p,l)

             end do

             })

             call WriteMode(diag, c)

          end do

          deallocate(diag%Fcos)
          deallocate(diag%Fsin)

          diag%done = .true. ! we are done here

       end if


    })

  end subroutine StepEDiagMode

!----------------------------------------------------------------------

  subroutine WriteMode(diag, c)

    type(T_DIAGMODE) :: diag
    integer :: c

    character(len=STRLNG) :: fn, sfx
    integer :: l, m, ios, mask 
    real(kind=8) :: En, F(12), F_tmp

    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    sfx = ".set"

    fn = cat4(diag%outfile,"_",TRIM(i2str(c)),sfx)

    M4_WRITE_INFO({"writing mode #",TRIM(i2str(diag%idx)), " -> ", TRIM(fn)})

    M4_MODOBJ_GETREG(diag,reg)

    open(UNITTMP,FILE=fn,STATUS="unknown", IOSTAT=ios)
    M4_OPEN_ERROR(ios,fn)

    write(UNITTMP,"(A)") "(SET"
    write(UNITTMP,*) reg%is,reg%ie,reg%di,reg%js,reg%je,reg%dj,reg%ks,reg%ke,reg%dk

    if ( diag%mode .eq. "En" ) then

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

          En = 0.
          En = En + 1./epsinvx(i,j,k) * ( diag%Fcos(c,p,1)**2 + diag%Fsin(c,p,1)**2 )
          En = En + 1./epsinvy(i,j,k) * ( diag%Fcos(c,p,2)**2 + diag%Fsin(c,p,2)**2 )
          En = En + 1./epsinvz(i,j,k) * ( diag%Fcos(c,p,3)**2 + diag%Fsin(c,p,3)**2 )
          En = En + 1./M4_MUINVX(i,j,k)* ( diag%Fcos(c,p,4)**2 + diag%Fsin(c,p,4)**2 )
          En = En + 1./M4_MUINVY(i,j,k)* ( diag%Fcos(c,p,5)**2 + diag%Fsin(c,p,5)**2 )
          En = En + 1./M4_MUINVZ(i,j,k)* ( diag%Fcos(c,p,6)**2 + diag%Fsin(c,p,6)**2 )

          write(UNITTMP,*) real(En,8)
          !write(UNITTMP,"(E15.6E3)") real(En,4)

       })

    elseif ( diag%mode .eq. "EHTG" ) then

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       do l = 1, diag%numcomp
          F(2*l-1) =  diag%Fcos(c,p,l)
          F(2*l) = diag%Fsin(c,p,l)
       end do

       call GeoAverage(F(5:6),F(9:10))
       call GeoAverage(F(7:8),F(11:12))

       write(UNITTMP,*) ( real(F(l),8), l = 1, 2*(diag%numcomp-2), 1 )

       })

    else 

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       do l = 1, diag%numcomp
          F(2*l-1) =  diag%Fcos(c,p,l)
          F(2*l) = diag%Fsin(c,p,l)
       end do

       write(UNITTMP,*) ( real(F(l),8), l = 1, 2*diag%numcomp, 1 )
       !write(UNITTMP,"(12E15.6E3)") ( real(F(l),4), l = 1, 12, 1 )

       })

    endif

    write(UNITTMP,"(A)") ")SET"

    close(UNITTMP)


  end subroutine WriteMode

!----------------------------------------------------------------------

  subroutine GeoAverage(F1,F2)
    real(kind=8) :: F1(2),F2(2),F_tmp
    integer :: mask

    if ( ( F1(1) >= 0. ) .and. ( F2(1) >= 0. ) .and. ( F1(2)*F2(2) <= 0. ) ) then
       mask = -1
    else
       mask = 1
    endif

    call CompSquareRoot(F1)
    call CompSquareRoot(F2)

    F_tmp = F1(1)
    F1(1) = mask * ( F_tmp*F2(1) - F1(2)*F2(2) )
    F1(2) = mask * ( F_tmp*F2(2) + F1(2)*F2(1) )

  end subroutine GeoAverage

!----------------------------------------------------------------------

  subroutine CompSquareRoot(F)
    real(kind=8) :: F(2), F_tmp
    integer :: mask

    if ( F(2) < 0. ) then
       mask = -1
    else
       mask = 1
    endif
    F_tmp = sqrt ( F(1)**2 + F(2)**2 )
    F(2) = sqrt ( 0.5 * ( F_tmp - F(1) ) )
    F(1) = mask * sqrt ( 0.5 * ( F_tmp + F(1) ) )

  end subroutine CompSquareRoot

!----------------------------------------------------------------------

   subroutine EchoDiagModeObj(diag)

    type(T_DIAGMODE) :: diag

    M4_WRITE_INFO({"--- diagmode # ",&
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))


  end subroutine EchoDiagModeObj

!----------------------------------------------------------------------

end module diagmode

! =====================================================================


