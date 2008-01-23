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
  real(kind=8) :: hb1(3), hb2(3), ed1(3), divs1, divs2, jekh1, jekh2

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
    real(kind=8) :: hxc, hyc, hzc, bxc, byc, bzc, dsx, dsy, dsz
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    write(6,*) "H: ", ncyc

    m = mod(ncyc+3,3) + 1
    mo = mod(ncyc+3-1,3) + 1
    moo = mod(ncyc+3-2,3) + 1
    
    diag%hb1(m) = 0.
    diag%hb2(m) = 0.
    diag%divs2 = 0.
    diag%jekh2 = 0.

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag

       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

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

       bxc = 0.25 * real( ( ( 1./M4_MUINVX(i,j,k)*Hx(M4_COORD(i,j,k)) +  1./M4_MUINVX(i,j-1,k)*Hx(M4_COORD(i,j-1,k)) ) &
            + ( 1./M4_MUINVX(i,j,k-1)*Hx(M4_COORD(i,j,k-1)) + 1./M4_MUINVX(i,j-1,k-1)*Hx(M4_COORD(i,j-1,k-1)) ) )
       byc = 0.25 * real( ( ( 1./M4_MUINVY(i,j,k)*Hy(M4_COORD(i,j,k)) + 1./M4_MUINVY(i,j,k-1)*Hy(M4_COORD(i,j,k-1)) ) &
            + ( 1./M4_MUINVY(i-1,j,k)*Hy(M4_COORD(i-1,j,k)) + 1./M4_MUINVY(i-1,j,k-1)*Hy(M4_COORD(i-1,j,k-1)) ) )
       bzc = 0.25 * real( ( ( 1./M4_MUINVZ(i,j,k)*Hz(M4_COORD(i,j,k)) + 1./M4_MUINVZ(i-1,j,k)*Hz(M4_COORD(i-1,j,k)) ) &
            + ( 1./M4_MUINVZ(i,j-1,k)*Hz(M4_COORD(i,j-1,k)) + 1./M4_MUINVZ(i-1,j-1,k)*Hz(M4_COORD(i-1,j-1,k)) ) )

       },{

       bxc = hxc
       byc = hyc
       bzc = hzc

       })


       diag%hb1(m) = diag%hb1(m) + 0.125/DT * ( bxc*hxc + byc*hyc + bzc*hzc )
       diag%hb2(m) = diag%hb2(m) + 0.125/DT * 2. * ( bxc*diag%hxo(p) + byc*diag%hyo(p) + bzc*diag%hzo(p) )
       
       })

       ! add up contributions at time step n+3/2 ------------------------------------------------------------------
       ! 
       
       diag%dudt = diag%hb1(m) + diag%hb2(m) - diag%hb2(mo) - diag%hb1(moo) + diag%ed1(mo) - diag%ed1(moo)
       diag%divs = diag%divs1 +  diag%divs2
       diag%jekh = diag%jekh1 +  diag%jekh2
       diag%res = diag%dudt + diag%divs + diag%jekh

       write(6,*) diag%dudt, diag%divs, diag%jekh, diag%res

       ! time integration

       diag%sumdudt = diag%sumdudt + diag%dudt * DT
       diag%sumdivs = diag%sumdivs + diag%divs * DT
       diag%sumjekh = diag%sumjekh + diag%jekh * DT
       diag%sumres = diag%sumres + diag%res * DT

       write(6,*) diag%sumdudt, diag%sumdivs, diag%sumjekh, diag%sumres

       !  --------------------------------------------------------------------------------------------------------

       ! ------ calculate part of div S 

       dsx = 0.
       dsy = 0.
       dsz = 0.

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       ! loop over front / back face
       do s = 0, 1

          dsx = dsx + (-1)**s * ( &
               0.25*real(Ey(M4_COORD(i-s,j,k)) + Ey(M4_COORD(i-s+1,j,k)) + Ey(M4_COORD(i-s,j-1,k)) + Ey(M4_COORD(i-s+1,j-1,k)))* &
               0.5*real(Hz(M4_COORD(i-s,j,k)) + Hz(M4_COORD(i-s,j-1,k)) ) -  &
               0.25*real(Ez(M4_COORD(i-s,j,k)) + Ez(M4_COORD(i-s+1,j,k)) + Ez(M4_COORD(i-s,j,k-1)) + Ez(M4_COORD(i-s+1,j,k-1)))* &
               0.5*real(Hy(M4_COORD(i-s,j,k)) + Hy(M4_COORD(i-s,j,k-1)) ) &
               )
          dsy = dsy + (-1)**s * ( &
               0.25*real(Ez(M4_COORD(i,j-s,k)) + Ez(M4_COORD(i,j-s+1,k)) + Ez(M4_COORD(i,j-s,k-1)) + Ez(M4_COORD(i,j-s+1,k-1)))* &
               0.5*real(Hx(M4_COORD(i,j-s,k)) + Hx(M4_COORD(i,j-s,k-1)) ) - &
               0.25*real(Ex(M4_COORD(i,j-s,k)) + Ex(M4_COORD(i,j-s+1,k)) + Ex(M4_COORD(i-1,j-s,k)) + Ex(M4_COORD(i-1,j-s+1,k)))* &
               0.5*real(Hz(M4_COORD(i,j-s,k)) + Hz(M4_COORD(i-1,j-s,k)) ) &
               )
          dsz = dsz +  (-1)**s * ( &
               0.25*real(Ex(M4_COORD(i,j,k-s)) + Ex(M4_COORD(i,j,k-s+1)) + Ex(M4_COORD(i-1,j,k-s)) + Ex(M4_COORD(i-1,j,k-s+1)))* &
               0.5*real(Hy(M4_COORD(i,j,k-s)) + Hy(M4_COORD(i-1,j,k-s)) ) - &
               0.25*real(Ey(M4_COORD(i,j,k-s)) + Ey(M4_COORD(i,j,k-s+1)) + Ey(M4_COORD(i,j-1,k-s)) + Ey(M4_COORD(i,j-1,k-s+1)))* &
               0.5*real(Hx(M4_COORD(i,j,k-s)) + Hx(M4_COORD(i,j-1,k-s)) ) &
               )
       
       end do

       diag%divs2 = diag%divs2 + 0.5/SX * dsx + 0.5/SY * dsy + 0.5/SZ * dsz 

       diag%hxo(p) = hxc
       diag%hyo(p) = hyc
       diag%hzo(p) = hzc

       })

       diag%jekh2 = 0.5 * SumJEKHMat(diag%mask,ncyc)

    })
  
  end subroutine StepHDiagEBal


!----------------------------------------------------------------------


  subroutine StepEDiagEBal(ncyc)

    integer :: ncyc, m, mo, moo, s
    real(kind=8) :: exc, eyc, ezc, dxc, dyc, dzc, dsx, dsy, dsz
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    write(6,*) "E: ", ncyc

    m = mod(ncyc+3,3) + 1
    mo = mod(ncyc+3-1,3) + 1
    moo = mod(ncyc+3-2,3) + 1

    diag%ed1(m) = 0.
    diag%divs1 = 0.
    diag%jekh1 = 0.

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag
       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

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
          
          dsx = dsx + (-1)**s * ( &
               0.25*real(Ey(M4_COORD(i-s,j,k)) + Ey(M4_COORD(i-s+1,j,k)) + Ey(M4_COORD(i-s,j-1,k)) + Ey(M4_COORD(i-s+1,j-1,k)))* &
               0.5*real(Hz(M4_COORD(i-s,j,k)) + Hz(M4_COORD(i-s,j-1,k)) ) -  &
               0.25*real(Ez(M4_COORD(i-s,j,k)) + Ez(M4_COORD(i-s+1,j,k)) + Ez(M4_COORD(i-s,j,k-1)) + Ez(M4_COORD(i-s+1,j,k-1)))* &
               0.5*real(Hy(M4_COORD(i-s,j,k)) + Hy(M4_COORD(i-s,j,k-1)) ) &
               )
          dsy = dsy + (-1)**s * ( &
               0.25*real(Ez(M4_COORD(i,j-s,k)) + Ez(M4_COORD(i,j-s+1,k)) + Ez(M4_COORD(i,j-s,k-1)) + Ez(M4_COORD(i,j-s+1,k-1)))* &
               0.5*real(Hx(M4_COORD(i,j-s,k)) + Hx(M4_COORD(i,j-s,k-1)) ) - &
               0.25*real(Ex(M4_COORD(i,j-s,k)) + Ex(M4_COORD(i,j-s+1,k)) + Ex(M4_COORD(i-1,j-s,k)) + Ex(M4_COORD(i-1,j-s+1,k)))* &
               0.5*real(Hz(M4_COORD(i,j-s,k)) + Hz(M4_COORD(i-1,j-s,k)) ) &
               )
          dsz = dsz +  (-1)**s * ( &
               0.25*real(Ex(M4_COORD(i,j,k-s)) + Ex(M4_COORD(i,j,k-s+1)) + Ex(M4_COORD(i-1,j,k-s)) + Ex(M4_COORD(i-1,j,k-s+1)))* &
               0.5*real(Hy(M4_COORD(i,j,k-s)) + Hy(M4_COORD(i-1,j,k-s)) ) - &
               0.25*real(Ey(M4_COORD(i,j,k-s)) + Ey(M4_COORD(i,j,k-s+1)) + Ey(M4_COORD(i,j-1,k-s)) + Ey(M4_COORD(i,j-1,k-s+1)))* &
               0.5*real(Hx(M4_COORD(i,j,k-s)) + Hx(M4_COORD(i,j-1,k-s)) ) &
               )
       
       end do


       diag%divs1 = diag%divs1 + 0.5/SX * dsx + 0.5/SY * dsy + 0.5/SZ * dsz 
       
      
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
! Authors:  E.Kirby, J.Hamm
! Modified: 23/01/2007
!
! =====================================================================


