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
! dudt + ds + je + kh = res
!
! ds = dsxp + dsxm + dsyp + dsym + dszp + dszm
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
  real(kind=8) :: dudt, ds, je, kh, res
  real(kind=8) :: dsxp, dsxm, dsyp, dsym, dszp, dszm
  ! spatially and time integrated energy contributions
  real(kind=8) :: sumdudt, sumds, sumje, sumkh, sumres
  real(kind=8) :: sumdsxp, sumdsxm, sumdsyp, sumdsym, sumdszp, sumdszm
  
  logical, pointer, dimension(:,:,:) :: mask
  
  ! partial contributions to energy terms
  real(kind=8), dimension(3) :: en, skh, sje
  real(kind=8), dimension(3) :: dsx1p, dsx2p, dsx3p, dsx4p
  real(kind=8), dimension(3) :: dsy1p, dsy2p, dsy3p, dsy4p
  real(kind=8), dimension(3) :: dsz1p, dsz2p, dsz3p, dsz4p
  real(kind=8), dimension(3) :: dsx1m, dsx2m, dsx3m, dsx4m
  real(kind=8), dimension(3) :: dsy1m, dsy2m, dsy3m, dsy4m
  real(kind=8), dimension(3) :: dsz1m, dsz2m, dsz3m, dsz4m

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
    diag%ds = 0.
    diag%je = 0.
    diag%kh = 0.
    diag%res = 0.
    diag%dsxp = 0.
    diag%dsxm = 0.
    diag%dsyp = 0.
    diag%dsym = 0.
    diag%dszp = 0.
    diag%dszm = 0.

    diag%sumdudt = 0.
    diag%sumds = 0.
    diag%sumje = 0.
    diag%sumkh = 0.
    diag%sumres = 0.
    diag%sumdsxp = 0.
    diag%sumdsxm = 0.
    diag%sumdsyp = 0.
    diag%sumdsym = 0.
    diag%sumdszp = 0.
    diag%sumdszm = 0.

    reg = regobj(diag%regidx)

    M4_WRITE_DBG("allocate mask")
    allocate(diag%mask(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX), stat = err)
    M4_ALLOC_ERROR(err,"InitializeDiagEBal")
 
    M4_WRITE_DBG("initialize mask")
    call SetMaskRegObj(reg,diag%mask,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX)

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

    integer :: ncyc, m, mo, moo
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m    = mod(ncyc*2+2,3) + 1

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag

       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       M4_MODOBJ_GETREG(diag,reg)

       diag%dsx1p(m) = 0.
       diag%dsx2p(m) = 0.
       diag%dsx3p(m) = 0. 
       diag%dsx4p(m) = 0.

       diag%dsx1m(m) = 0.
       diag%dsx2m(m) = 0.
       diag%dsx3m(m) = 0.
       diag%dsx4m(m) = 0.

       diag%dsy1p(m) = 0.
       diag%dsy2p(m) = 0. 
       diag%dsy3p(m) = 0.
       diag%dsy4p(m) = 0.

       diag%dsy1m(m) = 0.
       diag%dsy2m(m) = 0.
       diag%dsy3m(m) = 0.
       diag%dsy4m(m) = 0.

       diag%dsz1p(m) = 0.
       diag%dsz2p(m) = 0.
       diag%dsz3p(m) = 0. 
       diag%dsz4p(m) = 0.

       diag%dsz1m(m) = 0.
       diag%dsz2m(m) = 0.
       diag%dsz3m(m) = 0.
       diag%dsz4m(m) = 0.
       

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
   
      ! poynting vector
 
       if ( .not. diag%mask(i+1,j,k) ) then 
          diag%dsx1p(m) = diag%dsx1p(m) + 1./SX * real(Ey(i,j,k)) * real(Hz(i,j,k))
          diag%dsx2p(m) = diag%dsx2p(m) + 1./SX * real(Hz(i,j,k)) * real(Ey(i+1,j,k))
          diag%dsx3p(m) = diag%dsx3p(m) - 1./SX * real(Ez(i,j,k)) * real(Hy(i,j,k))
          diag%dsx4p(m) = diag%dsx4p(m) - 1./SX * real(Hy(i,j,k)) * real(Ez(i+1,j,k))
       end if

       if ( .not. diag%mask(i-1,j,k) ) then
          diag%dsx1m(m) = diag%dsx1m(m) - 1./SX * real(Ey(i,j,k)) * real(Hz(i-1,j,k))
          diag%dsx2m(m) = diag%dsx2m(m) - 1./SX * real(Hz(i,j,k)) * real(Ey(i,j,k))        
          diag%dsx3m(m) = diag%dsx3m(m) + 1./SX * real(Ez(i,j,k)) * real(Hy(i-1,j,k))
          diag%dsx4m(m) = diag%dsx4m(m) + 1./SX * real(Hy(i,j,k)) * real(Ez(i,j,k))
       end if


M4_IFELSE_1D({},{


       if ( .not. diag%mask(i,j+1,k) ) then 
          diag%dsy1p(m) = diag%dsy1p(m) + 1./SY * real(Ez(i,j,k)) * real(Hx(i,j,k))
          diag%dsy2p(m) = diag%dsy2p(m) + 1./SY * real(Hx(i,j,k)) * real(Ez(i,j+1,k))
          diag%dsy3p(m) = diag%dsy3p(m) - 1./SY * real(Ex(i,j,k)) * real(Hz(i,j,k))
          diag%dsy4p(m) = diag%dsy4p(m) - 1./SY * real(Hz(i,j,k)) * real(Ex(i,j+1,k))
       end if

       if ( .not. diag%mask(i,j-1,k) ) then 
          diag%dsy1m(m) = diag%dsy1m(m) - 1./SY * real(Ez(i,j,k)) * real(Hx(i,j-1,k))
          diag%dsy2m(m) = diag%dsy2m(m) - 1./SY * real(Hx(i,j,k)) * real(Ez(i,j,k))
          diag%dsy3m(m) = diag%dsy3m(m) + 1./SY * real(Ex(i,j,k)) * real(Hz(i,j-1,k))
          diag%dsy4m(m) = diag%dsy4m(m) + 1./SY * real(Hz(i,j,k)) * real(Ex(i,j,k))
       end if

})

M4_IFELSE_3D({

       if ( .not. diag%mask(i,j,k+1) ) then 
          diag%dsz1p(m) = diag%dsz1p(m) + 1./SZ * real(Ex(i,j,k)) * real(Hy(i,j,k))
          diag%dsz2p(m) = diag%dsz2p(m) + 1./SZ * real(Hy(i,j,k)) * real(Ex(i,j,k+1))
          diag%dsz3p(m) = diag%dsz3p(m) - 1./SZ * real(Ey(i,j,k)) * real(Hx(i,j,k))
          diag%dsz4p(m) = diag%dsz4p(m) - 1./SZ * real(Hx(i,j,k)) * real(Ey(i,j,k+1))
       end if

       if ( .not. diag%mask(i,j,k-1) ) then 
          diag%dsz1m(m) = diag%dsz1m(m) - 1./SZ * real(Ex(i,j,k)) * real(Hy(i,j,k-1))
          diag%dsz2m(m) = diag%dsz2m(m) - 1./SZ * real(Hy(i,j,k)) * real(Ex(i,j,k))
          diag%dsz3m(m) = diag%dsz3m(m) + 1./SZ * real(Ey(i,j,k)) * real(Hx(i,j,k-1))
          diag%dsz4m(m) = diag%dsz4m(m) + 1./SZ * real(Hx(i,j,k)) * real(Ey(i,j,k))
       end if
})

       })

       ! material sources

       diag%skh(m) = SumKHMat(diag%mask,ncyc)
       diag%sje(m) = SumJEMat(diag%mask,ncyc)

    })
  
  end subroutine StepHDiagEBal


!----------------------------------------------------------------------


  subroutine StepEDiagEBal(ncyc)

    integer :: ncyc, m, mo, moo
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m    = mod(ncyc*2+3,3) + 1
    mo   = mod(ncyc*2+2,3) + 1 ! half timestep back
    moo  = mod(ncyc*2+1,3) + 1 ! full timestep back

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag
       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       diag%en(m) = 0.

       diag%dsx1p(m) = 0.
       diag%dsx2p(m) = 0.
       diag%dsx3p(m) = 0. 
       diag%dsx4p(m) = 0.

       diag%dsx1m(m) = 0.
       diag%dsx2m(m) = 0.
       diag%dsx3m(m) = 0.
       diag%dsx4m(m) = 0.

       diag%dsy1p(m) = 0.
       diag%dsy2p(m) = 0. 
       diag%dsy3p(m) = 0.
       diag%dsy4p(m) = 0.

       diag%dsy1m(m) = 0.
       diag%dsy2m(m) = 0.
       diag%dsy3m(m) = 0.
       diag%dsy4m(m) = 0.

       diag%dsz1p(m) = 0.
       diag%dsz2p(m) = 0.
       diag%dsz3p(m) = 0. 
       diag%dsz4p(m) = 0.

       diag%dsz1m(m) = 0.
       diag%dsz2m(m) = 0.
       diag%dsz3m(m) = 0.
       diag%dsz4m(m) = 0.

       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       ! energy density

       diag%en(m) =  diag%en(m) + &
            epsinvx(i,j,k) * real(Ex(i,j,k))*real(Ex(i,j,k)) + &
            epsinvy(i,j,k) * real(Ey(i,j,k))*real(Ey(i,j,k)) + &
            epsinvz(i,j,k) * real(Ez(i,j,k))*real(Ez(i,j,k)) + &
            M4_MUINVX(i,j,k) * real(Hx(i,j,k))*real(Hx(i,j,k)) + &
            M4_MUINVY(i,j,k) * real(Hy(i,j,k))*real(Hy(i,j,k)) + &
            M4_MUINVZ(i,j,k) * real(Hz(i,j,k))*real(Hz(i,j,k))

       ! poynting vector
    
       if ( .not. diag%mask(i+1,j,k) ) then 
          diag%dsx1p(m) = diag%dsx1p(m) + 1./SX * real(Ey(i,j,k)) * real(Hz(i,j,k))
          diag%dsx2p(m) = diag%dsx2p(m) + 1./SX * real(Hz(i,j,k)) * real(Ey(i+1,j,k))
          diag%dsx3p(m) = diag%dsx3p(m) - 1./SX * real(Ez(i,j,k)) * real(Hy(i,j,k))
          diag%dsx4p(m) = diag%dsx4p(m) - 1./SX * real(Hy(i,j,k)) * real(Ez(i+1,j,k))
       end if

       if ( .not. diag%mask(i-1,j,k) ) then
          diag%dsx1m(m) = diag%dsx1m(m) - 1./SX * real(Ey(i,j,k)) * real(Hz(i-1,j,k))
          diag%dsx2m(m) = diag%dsx2m(m) - 1./SX * real(Hz(i,j,k)) * real(Ey(i,j,k))        
          diag%dsx3m(m) = diag%dsx3m(m) + 1./SX * real(Ez(i,j,k)) * real(Hy(i-1,j,k))
          diag%dsx4m(m) = diag%dsx4m(m) + 1./SX * real(Hy(i,j,k)) * real(Ez(i,j,k))
       end if


M4_IFELSE_1D({},{


       if ( .not. diag%mask(i,j+1,k) ) then 
          diag%dsy1p(m) = diag%dsy1p(m) + 1./SY * real(Ez(i,j,k)) * real(Hx(i,j,k))
          diag%dsy2p(m) = diag%dsy2p(m) + 1./SY * real(Hx(i,j,k)) * real(Ez(i,j+1,k))
          diag%dsy3p(m) = diag%dsy3p(m) - 1./SY * real(Ex(i,j,k)) * real(Hz(i,j,k))
          diag%dsy4p(m) = diag%dsy4p(m) - 1./SY * real(Hz(i,j,k)) * real(Ex(i,j+1,k))
       end if

       if ( .not. diag%mask(i,j-1,k) ) then 
          diag%dsy1m(m) = diag%dsy1m(m) - 1./SY * real(Ez(i,j,k)) * real(Hx(i,j-1,k))
          diag%dsy2m(m) = diag%dsy2m(m) - 1./SY * real(Hx(i,j,k)) * real(Ez(i,j,k))
          diag%dsy3m(m) = diag%dsy3m(m) + 1./SY * real(Ex(i,j,k)) * real(Hz(i,j-1,k))
          diag%dsy4m(m) = diag%dsy4m(m) + 1./SY * real(Hz(i,j,k)) * real(Ex(i,j,k))
       end if

})

M4_IFELSE_3D({

       if ( .not. diag%mask(i,j,k+1) ) then 
          diag%dsz1p(m) = diag%dsz1p(m) + 1./SZ * real(Ex(i,j,k)) * real(Hy(i,j,k))
          diag%dsz2p(m) = diag%dsz2p(m) + 1./SZ * real(Hy(i,j,k)) * real(Ex(i,j,k+1))
          diag%dsz3p(m) = diag%dsz3p(m) - 1./SZ * real(Ey(i,j,k)) * real(Hx(i,j,k))
          diag%dsz4p(m) = diag%dsz4p(m) - 1./SZ * real(Hx(i,j,k)) * real(Ey(i,j,k+1))
       end if

       if ( .not. diag%mask(i,j,k-1) ) then 
          diag%dsz1m(m) = diag%dsz1m(m) - 1./SZ * real(Ex(i,j,k)) * real(Hy(i,j,k-1))
          diag%dsz2m(m) = diag%dsz2m(m) - 1./SZ * real(Hy(i,j,k)) * real(Ex(i,j,k))
          diag%dsz3m(m) = diag%dsz3m(m) + 1./SZ * real(Ey(i,j,k)) * real(Hx(i,j,k-1))
          diag%dsz4m(m) = diag%dsz4m(m) + 1./SZ * real(Hx(i,j,k)) * real(Ey(i,j,k))
       end if
})

      
       })

       ! material sources

       diag%skh(m) = SumKHMat(diag%mask,ncyc)
       diag%sje(m) = SumJEMat(diag%mask,ncyc)
       
       ! ---------------------------------------------------
       
        
!       diag%dudt = .5/DT * ( diag%en(m) - diag%en(mo) )
!       diag%divs = .5 * ( diag%divs1 + diag%divs2(mo) - diag%divs3 - diag%divs4 ) 
!       diag%jekh = .5 * ( diag%q1 + diag%q2(mo) + diag%q3 + diag%q4 ) 

       diag%dudt = .5/DT * ( diag%en(m) - diag%en(moo) )

       diag%dsxp = .5 * ( diag%dsx1p(m) + diag%dsx1p(mo) + diag%dsx3p(m) + diag%dsx3p(mo) + &
            diag%dsx2p(mo) + diag%dsx2p(moo) +  diag%dsx4p(mo) + diag%dsx4p(moo) )
       diag%dsxm = .5 * ( diag%dsx1m(m) + diag%dsx1m(mo) + diag%dsx3m(m) + diag%dsx3m(mo) + &
            diag%dsx2m(mo) + diag%dsx2m(moo) +  diag%dsx4m(mo) + diag%dsx4m(moo) )
       diag%dsyp = .5 * ( diag%dsy1p(m) + diag%dsy1p(mo) + diag%dsy3p(m) + diag%dsy3p(mo) + &
            diag%dsy2p(mo) + diag%dsy2p(moo) +  diag%dsy4p(mo) + diag%dsy4p(moo) )
       diag%dsym = .5 * ( diag%dsy1m(m) + diag%dsy1m(mo) + diag%dsy3m(m) + diag%dsy3m(mo) + &
            diag%dsy2m(mo) + diag%dsy2m(moo) +  diag%dsy4m(mo) + diag%dsy4m(moo) )
       diag%dszp = .5 * ( diag%dsz1p(m) + diag%dsz1p(mo) + diag%dsz3p(m) + diag%dsz3p(mo) + &
            diag%dsz2p(mo) + diag%dsz2p(moo) +  diag%dsz4p(mo) + diag%dsz4p(moo) )
       diag%dszm = .5 * ( diag%dsz1m(m) + diag%dsz1m(mo) + diag%dsz3m(m) + diag%dsz3m(mo) + &
            diag%dsz2m(mo) + diag%dsz2m(moo) +  diag%dsz4m(mo) + diag%dsz4m(moo) )

       diag%je = 0.5 * ( diag%sje(m) + diag%sje(mo) )
       diag%kh = 0.5 * ( diag%skh(mo) + diag%skh(moo) )

       diag%ds   = diag%dsxp + diag%dsxm + diag%dsyp + diag%dsym + diag%dszp + diag%dszm
       diag%res = diag%dudt + diag%ds + diag%je + diag%kh
       
       diag%sumdudt = diag%sumdudt + DT * diag%dudt
       diag%sumje = diag%sumje + DT * diag%je
       diag%sumkh = diag%sumkh + DT * diag%kh
       diag%sumdsxp = diag%sumdsxp +  DT * diag%dsxp
       diag%sumdsxm = diag%sumdsxm +  DT * diag%dsxm
       diag%sumdsyp = diag%sumdsyp +  DT * diag%dsyp
       diag%sumdsym = diag%sumdsym +  DT * diag%dsym
       diag%sumdszp = diag%sumdszp +  DT * diag%dszp
       diag%sumdszm = diag%sumdszm +  DT * diag%dszm
       diag%sumds   = diag%sumdsxp + diag%sumdsxm + diag%sumdsyp + diag%sumdsym + diag%sumdszp + diag%sumdszm
       diag%sumres = diag%sumdudt + diag%sumds + diag%sumje + diag%sumkh

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


