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
    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. diag%dn .lt. 1 ) then
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
    integer :: unit, ios, i
    real(kind=8) :: freqs(1000)

    M4_WRITE_DBG(". enter InitializeDiagMode")
    M4_MODLOOP_EXPR({DIAGMODE},diag,{
    
    M4_MODOBJ_GETREG(diag,reg)

    diag%numsteps = (diag%ne-diag%ns)/diag%dn + 1

    open(unit,FILE=diag%filename,STATUS="old", IOSTAT=ios)

    if ( ios .ne. 0 ) then 
       M4_WRITE_WARN({"could not open frequency file: ",  TRIM(diag%filename),"!"})
       return
    end if

    diag%numfreqs = 0
    ios = 0

    do while ( ios .eq. 0 )
       diag%numfreqs = diag%numfreqs + 1
       read(unit,*,iostat=ios) freqs(diag%numfreqs)
    end do

    close(unit)

!    diag%numfreqs = diag%numfreqs - 1
  
    allocate(diag%freqs(diag%numfreqs))

    do i = 1, diag%numfreqs
       diag%freqs(i) = freqs(i)
    end do
    
    diag%done = .false.

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

    M4_MODLOOP_EXPR({DIAGMODE},diag,{

       ! this loops over all diag structures, setting diag

       M4_MODOBJ_GETREG(diag,reg)

       if ( diag%done ) cycle

       if ( ncyc .eq. diag%ns ) then

          M4_WRITE_INFO({"recording dft #",TRIM(i2str(diag%idx))})

          allocate(diag%Fcos(diag%numfreqs,reg%numnodes,6),stat=ier) 
          M4_ALLOC_ERROR({ier},"StepEDiagMode")
          
          allocate(diag%Fsin(diag%numfreqs,reg%numnodes,6),stat=ier) 
          M4_ALLOC_ERROR({ier},"StepEDiagMode")
   
          diag%Fcos = 0.
          diag%Fsin = 0.

       end if

       if ( ncyc .ge. diag%ns .and. ncyc .le. diag%ne .and. &
            mod(ncyc-diag%ns,diag%dn) .eq. 0) then
          
          
          M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
          ! store E field projections
          
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

       if ( ncyc .ge. diag%ne .and. .not. diag%done ) then
          
          do c = 1, diag%numfreqs 
             
             do l = 1, 6
             
                diag%Fcos(c,p,l) = 2./diag%numsteps * diag%Fcos(c,p,l)
                diag%Fsin(c,p,l) = 2./diag%numsteps * diag%Fsin(c,p,l)

             end do

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
    integer :: l, m, ios
    real(kind=4) :: En, F(12)

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

          write(UNITTMP,"(E15.6E3)") real(En,4)

       })

    else 

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
        do l = 1,6
           F(2*l-1) =  diag%Fcos(c,p,l)
           F(2*l) = diag%Fsin(c,p,l)
        end do
            
        write(UNITTMP,"(12E15.6E3)") ( real(F(l),4), l = 1, 12, 1 )
       
       })
     
    endif

    write(UNITTMP,"(A)") ")SET"

    close(UNITTMP)


  end subroutine WriteMode

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


