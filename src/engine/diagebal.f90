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
  
  logical, pointer, dimension(:,:,:) :: mask
  
  ! partial contributions to energy terms
  real(kind=8) :: en(2)
  real(kind=8) :: divs1, divs2(2), divs3, divs4
  real(kind=8) :: q1, q2(2), q3, q4

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

    diag%en = 0
    diag%divs2 = 0
    diag%q2 = 0

    reg = regobj(diag%regidx)

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

    deallocate(diag%mask)

    })
    M4_WRITE_DBG(". exit FinalizeMatEBal")

  end subroutine FinalizeDiagEBal

!----------------------------------------------------------------------

  subroutine StepHDiagEBal(ncyc)

    integer :: ncyc
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag

       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       M4_MODOBJ_GETREG(diag,reg)

       diag%divs1 = 0.
       diag%divs4 = 0.

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
       
       diag%divs1 = diag%divs1 + real(Hx(i,j,k)) * real( &
            M4_IFELSE_1D({0.&},{ + 1./SY*( Ez(i,j+1,k) - Ez(i,j,k) ) & })
            M4_IFELSE_3D({- 1./SZ*( Ey(i,j,k+1) - Ey(i,j,k) ) & },{+ 0. &})
       ) 
       diag%divs1 = diag%divs1 + real(Hy(i,j,k)) * real( &
            - 1./SX*( Ez(i+1,j,k) - Ez(i,j,k) ) &
            M4_IFELSE_3D({ + 1./SZ*( Ex(i,j,k+1) - Ex(i,j,k) ) & })
       )
       diag%divs1 = diag%divs1 + real(Hz(i,j,k)) * real( &
            + 1./SX*( Ey(i+1,j,k) - Ey(i,j,k) ) &
            M4_IFELSE_1D({},{   - 1./SY*( Ex(i,j+1,k) - Ex(i,j,k) ) & })
       )
       
       diag%divs4 = diag%divs4 + real(Ex(i,j,k)) * real( &
            M4_IFELSE_1D({0.&},{ + 1./SY*( Hz(i,j,k) - Hz(i,j-1,k) ) & })
            M4_IFELSE_3D({  - 1./SZ*( Hy(i,j,k) - Hy(i,j,k-1) ) & },{+ 0. &})
       )
       
       diag%divs4 = diag%divs4 + real(Ey(i,j,k)) * real( &
            - 1./SX*( Hz(i,j,k) - Hz(i-1,j,k) ) &
            M4_IFELSE_3D({  + 1./SZ*( Hx(i,j,k) - Hx(i,j,k-1) ) & })
       )

       diag%divs4 = diag%divs4 + real(Ez(i,j,k)) * real( &
            + 1./SX*( Hy(i,j,k) - Hy(i-1,j,k) ) &
            M4_IFELSE_1D({},{  - 1./SY*( Hx(i,j,k) - Hx(i,j-1,k) )  & })
       )
      
       })

       diag%q1 = SumKHMat(diag%mask,ncyc)
       diag%q4 = SumJEMat(diag%mask,ncyc)

    })
  
  end subroutine StepHDiagEBal


!----------------------------------------------------------------------


  subroutine StepEDiagEBal(ncyc)

    integer :: ncyc, m, mo
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m  = mod(ncyc+1,2) + 1
    mo = mod(ncyc,2) + 1

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag
       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       diag%divs2(m) = 0.
       diag%divs3 = 0.
       diag%en(m) = 0.

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{


       diag%divs2(m) = diag%divs2(m) + real(Hx(i,j,k)) * real( &
            M4_IFELSE_1D({0.&},{ + 1./SY*( Ez(i,j+1,k) - Ez(i,j,k) ) & })
            M4_IFELSE_3D({ - 1./SZ*( Ey(i,j,k+1) - Ey(i,j,k) ) & },{+ 0. &})
       ) 
       diag%divs2(m) = diag%divs2(m) + real(Hy(i,j,k)) * real( &
            - 1./SX*( Ez(i+1,j,k) - Ez(i,j,k) ) &
            M4_IFELSE_3D({ + 1./SZ*( Ex(i,j,k+1) - Ex(i,j,k) ) & })
       )
       diag%divs2(m) = diag%divs2(m) + real(Hz(i,j,k)) * real( &
            + 1./SX*( Ey(i+1,j,k) - Ey(i,j,k) ) &
            M4_IFELSE_1D({},{   - 1./SY*( Ex(i,j+1,k) - Ex(i,j,k) ) & })
       )
       
       diag%divs3 = diag%divs3 + real(Ex(i,j,k)) * real( &
            M4_IFELSE_1D({0.&},{ + 1./SY*( Hz(i,j,k) - Hz(i,j-1,k) ) & })
            M4_IFELSE_3D({  - 1./SZ*( Hy(i,j,k) - Hy(i,j,k-1) ) & },{+ 0. &})
       )
       diag%divs3 = diag%divs3 + real(Ey(i,j,k)) * real( &
            - 1./SX*( Hz(i,j,k) - Hz(i-1,j,k) ) &
            M4_IFELSE_3D({  + 1./SZ*( Hx(i,j,k) - Hx(i,j,k-1) ) & })
       )
       diag%divs3 = diag%divs3 + real(Ez(i,j,k)) * real( &
            + 1./SX*( Hy(i,j,k) - Hy(i-1,j,k) ) &
            M4_IFELSE_1D({},{  - 1./SY*( Hx(i,j,k) - Hx(i,j-1,k) )  & })
       )

       diag%en(m) =  diag%en(m) + &
            epsinvx(i,j,k) * real(Ex(i,j,k))*real(Ex(i,j,k)) + &
            epsinvy(i,j,k) * real(Ey(i,j,k))*real(Ey(i,j,k)) + &
            epsinvz(i,j,k) * real(Ez(i,j,k))*real(Ez(i,j,k)) + &
            M4_MUINVX(i,j,k) * real(Hx(i,j,k))*real(Hx(i,j,k)) + &
            M4_MUINVY(i,j,k) * real(Hy(i,j,k))*real(Hy(i,j,k)) + &
            M4_MUINVZ(i,j,k) * real(Hz(i,j,k))*real(Hz(i,j,k))
             
       })

       diag%q2(m) = SumKHMat(diag%mask,ncyc)
       diag%q3 = SumJEMat(diag%mask,ncyc)
       
       ! ---------------------------------------------------


       diag%dudt = .5/DT * ( diag%en(m) - diag%en(mo) )
       diag%divs = .5 * ( diag%divs1 + diag%divs2(mo) - diag%divs3 - diag%divs4 ) 
       diag%jekh = .5 * ( diag%q1 + diag%q2(mo) + diag%q3 + diag%q4 ) 
       diag%res = diag%dudt + diag%divs + diag%jekh
       

       diag%sumdudt = diag%sumdudt + DT * diag%dudt
       diag%sumdivs = diag%sumdivs + DT * diag%divs
       diag%sumjekh = diag%sumjekh + DT * diag%jekh
       diag%sumres = diag%sumdudt + diag%sumdivs + diag%sumjekh

!       write(6,*) diag%dudt, diag%divs, diag%jekh, diag%res
!       write(6,*) diag%sumdudt, diag%sumdivs, diag%sumjekh, diag%sumres


       ! ---------------------------------------------------
 
      
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


