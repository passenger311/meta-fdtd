!-*- F90 -*------------------------------------------------------------
!
!  module: diagebal / meta
!
!  Energy balance diagnostics module.
!
!----------------------------------------------------------------------


! =====================================================================
!
! The DiagEBal module measures the energy flux in and out of a given
! volume and the amount of energy generated inside this volume.
!
! du/dt + div S = - J E - K H
!
! This equation is centered at timestep n+1/2. du/dt and J dot E must be 
! spatially centered to the center of the Yee cell, while the components 
! of S are sitting on the respective faces of the Yee cell. 
!
! dudt + divs + jekh = def
!
! where def = 0 should hold at each point at each timestep.  

module diagebal

  use constant
  use mpiworld
  use reglist
  use outlist
  use grid
  use fdtd
  use mat

  implicit none
  private
  save

  M4_MODHEAD_DECL({DIAGEBAL},100,{

  integer :: ns, ne, dn  ! time stepping 

  ! spatially integrated energy contributions
  real(kind=8) :: dudt, divs, jekh, res
  ! spatially and time integrated energy contributions
  real(kind=8) :: sumdudt, sumdivs, sumjekh, sumres
  
  ! h field at previous timestep
  real(kind=8), pointer, dimension(:) :: hxo, hyo, hzo
  logical, pointer, dimension(:,:,:) :: mask
  
  ! partial contributions to energy terms
  real(kind=8) :: hb1(3), hb2(3), ed1(3),jekh1, jekh2
  real(kind=8), dimension(0:1) ::  divsx1, divsy1, divsz1
  real(kind=8), dimension(0:1) ::  divsx2, divsy2, divsz2

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagEBalObj(funit)

    M4_MODREAD_DECL({DIAGEBAL}, funit,diag,reg,out)

    M4_WRITE_DBG(". enter ReadMatEBalObj")
    
    M4_MODREAD_EXPR({DIAGEBAL}, funit,diag,reg,0,out, {

    ! read diag parameters here, as defined in diag data structure
    read(funit,*) diag%ns, diag%ne, diag%dn
 
    })

    M4_WRITE_DBG(". exit ReadMatEBalObj")

  end subroutine ReadDiagEBalObj

!----------------------------------------------------------------------

  subroutine InitializeDiagEBal

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({DIAGEBAL},diag)

    M4_WRITE_DBG(". enter InitializeMatEBal")
    M4_MODLOOP_EXPR({DIAGEBAL},diag,{
    
    diag%dudt = 0.
    diag%divs = 0.
    diag%jekh = 0.
    diag%res = 0.

    diag%sumdudt = 0.
    diag%sumdivs = 0.
    diag%sumjekh = 0.
    diag%sumres = 0.

    diag%ed1 = 0
    diag%hb1 = 0
    diag%hb2 = 0

    reg = regobj(diag%regidx)

    M4_WRITE_DBG("allocate h field")
    allocate(diag%hxo(reg%numnodes),diag%hyo(reg%numnodes),diag%hzo(reg%numnodes), stat = err)
    M4_ALLOC_ERROR(err,"InitializeDiagEBal")

    diag%hxo = 0.
    diag%hyo = 0.
    diag%hzo = 0.

    M4_WRITE_DBG("allocate mask")
    allocate(diag%mask(IBEG:IEND,JBEG:JEND,KBEG:KEND), stat = err)
    M4_ALLOC_ERROR(err,"InitializeDiagEBal")
 
    M4_WRITE_DBG("initialize mask")
    call SetMaskRegObj(reg,diag%mask,IBEG,IEND,JBEG,JEND,KBEG,KEND)

    M4_IFELSE_DBG({call EchoDiagEBalObj(diag)})

    })
    M4_WRITE_DBG(". exit InitializeMatEBal")

  end subroutine InitializeDiagEBal

!----------------------------------------------------------------------

  subroutine FinalizeDiagEBal

    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_WRITE_DBG(". enter FinalizeMatEBal")
    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

    deallocate(diag%hxo, diag%hyo, diag%hzo, diag%mask )

    })
    M4_WRITE_DBG(". exit FinalizeMatEBal")

  end subroutine FinalizeDiagEBal

!----------------------------------------------------------------------

  subroutine StepHDiagEBal(ncyc)

    integer :: ncyc, m, mo, moo, s
    real(kind=8) :: hxc, hyc, hzc, bxc, byc, bzc, dsx(0:1), dsy(0:1), dsz(0:1)
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m = mod(ncyc+3,3) + 1
    mo = mod(ncyc+3-1,3) + 1
    moo = mod(ncyc+3-2,3) + 1
    
    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag

       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       diag%hb1(m) = 0.
       diag%hb2(m) = 0.
       diag%jekh2 = 0.

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       ! this loops over all points of the region 

       hxc = 0.25 * real( ( Hx(M4_COORD(i,j,k)) + Hx(M4_COORD(i,j-1,k)) ) &
            + ( Hx(M4_COORD(i,j,k-1)) + Hx(M4_COORD(i,j-1,k-1)) ) )
       hyc = 0.25 * real( ( Hy(M4_COORD(i,j,k)) + Hy(M4_COORD(i,j,k-1)) ) &
            + ( Hy(M4_COORD(i-1,j,k)) + Hy(M4_COORD(i-1,j,k-1)) ) )
       hzc = 0.25 * real( ( Hz(M4_COORD(i,j,k)) + Hz(M4_COORD(i-1,j,k)) ) &
            + ( Hz(M4_COORD(i,j-1,k)) + Hz(M4_COORD(i-1,j-1,k)) ) )


       M4_IFELSE_WMU({

       bxc = 0.25 * real( ( 1./M4_MUINVX(i,j,k)*Hx(M4_COORD(i,j,k)) +  1./M4_MUINVX(i,j-1,k)*Hx(M4_COORD(i,j-1,k)) ) &
            + ( 1./M4_MUINVX(i,j,k-1)*Hx(M4_COORD(i,j,k-1)) + 1./M4_MUINVX(i,j-1,k-1)*Hx(M4_COORD(i,j-1,k-1)) ) )
       byc = 0.25 * real( ( 1./M4_MUINVY(i,j,k)*Hy(M4_COORD(i,j,k)) + 1./M4_MUINVY(i,j,k-1)*Hy(M4_COORD(i,j,k-1)) ) &
            + ( 1./M4_MUINVY(i-1,j,k)*Hy(M4_COORD(i-1,j,k)) + 1./M4_MUINVY(i-1,j,k-1)*Hy(M4_COORD(i-1,j,k-1)) ) )
       bzc = 0.25 * real( ( 1./M4_MUINVZ(i,j,k)*Hz(M4_COORD(i,j,k)) + 1./M4_MUINVZ(i-1,j,k)*Hz(M4_COORD(i-1,j,k)) ) &
            + ( 1./M4_MUINVZ(i,j-1,k)*Hz(M4_COORD(i,j-1,k)) + 1./M4_MUINVZ(i-1,j-1,k)*Hz(M4_COORD(i-1,j-1,k)) ) )

       },{

       bxc = hxc
       byc = hyc
       bzc = hzc

       })


       diag%hb1(m) = diag%hb1(m) + 0.5/DT * ( bxc*hxc + byc*hyc + bzc*hzc )
       diag%hb2(m) = diag%hb2(m) + 0.5/DT * ( bxc*diag%hxo(p) + byc*diag%hyo(p) + bzc*diag%hzo(p) )

       diag%hxo(p) = hxc
       diag%hyo(p) = hyc
       diag%hzo(p) = hzc
       
       })

       ! add up contributions at time step n+3/2 ------------------------------------------------------------------
       ! 

       !diag%dudt = diag%hb2(m) - diag%hb2(mo) + diag%ed1(mo) - diag%ed1(moo)
       diag%dudt = 0.25*diag%hb1(m) + .5 * diag%hb2(m) - .5* diag%hb2(mo) - 0.25 * diag%hb1(moo) + diag%ed1(mo) - diag%ed1(moo)

       diag%divs = diag%divsx1(0) + diag%divsy1(0) + diag%divsz1(0) + diag%divsx2(0) + diag%divsy2(0) + diag%divsz2(0) - ( &
            diag%divsx1(1) + diag%divsy1(1) + diag%divsz1(1) + diag%divsx2(1) + diag%divsy2(1) + diag%divsz2(1) )
       diag%jekh = diag%jekh1 +  diag%jekh2
       if ( diag%dudt .ne. 0.0 ) then
          diag%res = abs(diag%dudt + diag%divs + diag%jekh)/( abs(diag%dudt) + abs(diag%divs) + abs(diag%jekh) )
       end if

       write(6,*) diag%dudt, diag%divs, diag%jekh, diag%dudt + diag%divs + diag%jekh

       ! time integration

       diag%sumdudt = diag%sumdudt + diag%dudt * DT
       diag%sumdivs = diag%sumdivs + diag%divs * DT
       diag%sumjekh = diag%sumjekh + diag%jekh * DT
       diag%sumres =  diag%sumres + diag%res

       write(6,*) diag%sumdudt, diag%sumdivs, diag%sumjekh, diag%sumdudt + diag%sumdivs + diag%sumjekh

       !  --------------------------------------------------------------------------------------------------------

       ! ------ calculate part of div S 

       diag%divsx2 = 0.
       diag%divsy2 = 0.
       diag%divsz2 = 0.

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       dsx = 0.
       dsy = 0.
       dsz = 0.

       ! loop over front / back face
       do s = 0, 1

          dsx(s) = 0.125 * ( &
               real(Ey(M4_COORD(i-s,j,k)) + Ey(M4_COORD(i-s+1,j,k)) + Ey(M4_COORD(i-s,j-1,k)) + Ey(M4_COORD(i-s+1,j-1,k)))* &
               real(Hz(M4_COORD(i-s,j,k)) + Hz(M4_COORD(i-s,j-1,k)) ) -  &
               real(Ez(M4_COORD(i-s,j,k)) + Ez(M4_COORD(i-s+1,j,k)) + Ez(M4_COORD(i-s,j,k-1)) + Ez(M4_COORD(i-s+1,j,k-1)))* &
               real(Hy(M4_COORD(i-s,j,k)) + Hy(M4_COORD(i-s,j,k-1)) ) &
               )
          dsy(s) = 0.125 * ( &
               real(Ez(M4_COORD(i,j-s,k)) + Ez(M4_COORD(i,j-s+1,k)) + Ez(M4_COORD(i,j-s,k-1)) + Ez(M4_COORD(i,j-s+1,k-1)))* &
               real(Hx(M4_COORD(i,j-s,k)) + Hx(M4_COORD(i,j-s,k-1)) ) - &
               real(Ex(M4_COORD(i,j-s,k)) + Ex(M4_COORD(i,j-s+1,k)) + Ex(M4_COORD(i-1,j-s,k)) + Ex(M4_COORD(i-1,j-s+1,k)))* &
               real(Hz(M4_COORD(i,j-s,k)) + Hz(M4_COORD(i-1,j-s,k)) ) &
               )
          dsz(s) = 0.125 * ( &
               real(Ex(M4_COORD(i,j,k-s)) + Ex(M4_COORD(i,j,k-s+1)) + Ex(M4_COORD(i-1,j,k-s)) + Ex(M4_COORD(i-1,j,k-s+1)))* &
               real(Hy(M4_COORD(i,j,k-s)) + Hy(M4_COORD(i-1,j,k-s)) ) - &
               real(Ey(M4_COORD(i,j,k-s)) + Ey(M4_COORD(i,j,k-s+1)) + Ey(M4_COORD(i,j-1,k-s)) + Ey(M4_COORD(i,j-1,k-s+1)))* &
               real(Hx(M4_COORD(i,j,k-s)) + Hx(M4_COORD(i,j-1,k-s)) ) &
               )

          diag%divsx2(s) = diag%divsx2(s) + 0.5/SX * dsx(s)
          diag%divsy2(s) = diag%divsy2(s) + 0.5/SY * dsy(s)
          diag%divsz2(s) = diag%divsz2(s) + 0.5/SZ * dsz(s)
          
       end do


       })

       diag%jekh2 = 0.5 * SumJEKHMat(diag%mask,ncyc)

    })
  
  end subroutine StepHDiagEBal


!----------------------------------------------------------------------


  subroutine StepEDiagEBal(ncyc)

    integer :: ncyc, m, mo, moo, s
    real(kind=8) :: exc, eyc, ezc, dxc, dyc, dzc, dsx(0:1), dsy(0:1), dsz(0:1)
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m = mod(ncyc+3,3) + 1
    mo = mod(ncyc+3-1,3) + 1
    moo = mod(ncyc+3-2,3) + 1

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag
       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       diag%ed1(m) = 0.
       diag%divsx1 = 0.
       diag%divsy1 = 0.
       diag%divsz1 = 0.
       diag%jekh1 = 0.

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       

       exc = 0.5 * real (  Ex(M4_COORD(i,j,k)) + Ex(M4_COORD(i-1,j,k)) ) 
       eyc = 0.5 * real (  Ey(M4_COORD(i,j,k)) + Ey(M4_COORD(i,j-1,k)) ) 
       ezc = 0.5 * real (  Ez(M4_COORD(i,j,k)) + Ez(M4_COORD(i,j,k-1)) ) 

       dxc = 0.5 * real (  1./epsinvx(M4_COORD(i,j,k))*Ex(M4_COORD(i,j,k)) + 1./epsinvx(M4_COORD(i-1,j,k))*Ex(M4_COORD(i-1,j,k)) ) 
       dyc = 0.5 * real (  1./epsinvy(M4_COORD(i,j,k))*Ey(M4_COORD(i,j,k)) + 1./epsinvy(M4_COORD(i,j-1,k))*Ey(M4_COORD(i,j-1,k)) ) 
       dzc = 0.5 * real (  1./epsinvz(M4_COORD(i,j,k))*Ez(M4_COORD(i,j,k)) + 1./epsinvz(M4_COORD(i,j,k-1))*Ez(M4_COORD(i,j,k-1)) ) 

       diag%ed1(m) = diag%ed1(m) + 0.5/DT * ( exc*dxc + eyc*dyc + ezc*dzc )
       
       dsx = 0.
       dsy = 0.
       dsz = 0.
       
       ! loop over front / back face
       do s = 0, 1
          
          dsx(s) = 0.125 * ( &
               real(Ey(M4_COORD(i-s,j,k)) + Ey(M4_COORD(i-s+1,j,k)) + Ey(M4_COORD(i-s,j-1,k)) + Ey(M4_COORD(i-s+1,j-1,k)))* &
               real(Hz(M4_COORD(i-s,j,k)) + Hz(M4_COORD(i-s,j-1,k)) ) -  &
               real(Ez(M4_COORD(i-s,j,k)) + Ez(M4_COORD(i-s+1,j,k)) + Ez(M4_COORD(i-s,j,k-1)) + Ez(M4_COORD(i-s+1,j,k-1)))* &
               real(Hy(M4_COORD(i-s,j,k)) + Hy(M4_COORD(i-s,j,k-1)) ) &
               )
          dsy(s) = 0.125 * ( &
               real(Ez(M4_COORD(i,j-s,k)) + Ez(M4_COORD(i,j-s+1,k)) + Ez(M4_COORD(i,j-s,k-1)) + Ez(M4_COORD(i,j-s+1,k-1)))* &
               real(Hx(M4_COORD(i,j-s,k)) + Hx(M4_COORD(i,j-s,k-1)) ) - &
               real(Ex(M4_COORD(i,j-s,k)) + Ex(M4_COORD(i,j-s+1,k)) + Ex(M4_COORD(i-1,j-s,k)) + Ex(M4_COORD(i-1,j-s+1,k)))* &
               real(Hz(M4_COORD(i,j-s,k)) + Hz(M4_COORD(i-1,j-s,k)) ) &
               )
          dsz(s) = 0.125 * ( &
               real(Ex(M4_COORD(i,j,k-s)) + Ex(M4_COORD(i,j,k-s+1)) + Ex(M4_COORD(i-1,j,k-s)) + Ex(M4_COORD(i-1,j,k-s+1)))* &
               real(Hy(M4_COORD(i,j,k-s)) + Hy(M4_COORD(i-1,j,k-s)) ) - &
               real(Ey(M4_COORD(i,j,k-s)) + Ey(M4_COORD(i,j,k-s+1)) + Ey(M4_COORD(i,j-1,k-s)) + Ey(M4_COORD(i,j-1,k-s+1)))* &
               real(Hx(M4_COORD(i,j,k-s)) + Hx(M4_COORD(i,j-1,k-s)) ) &
               )

          diag%divsx1(s) = diag%divsx1(s) + 0.5/SX * dsx(s)
          diag%divsy1(s) = diag%divsy1(s) + 0.5/SY * dsy(s)
          diag%divsz1(s) = diag%divsz1(s) + 0.5/SZ * dsz(s)
           
       end do
       
      
       })

       diag%jekh1 = 0.5 * SumJEKHMat(diag%mask,ncyc)
      
    })

  end subroutine StepEDiagEBal

!----------------------------------------------------------------------

   subroutine EchoDiagEBalObj(diag)

    type(T_DIAGEBAL) :: diag
 
    M4_WRITE_INFO({"--- diagebal # ",&
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"ns ne dn = ",diag%ns, diag%ne, diag%dn })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))
    

  end subroutine EchoDiagEBalObj
  
!----------------------------------------------------------------------

end module diagebal

!
! Authors:  J.Hamm, E.Kirby
! Modified: 23/01/2007
!
! =====================================================================


