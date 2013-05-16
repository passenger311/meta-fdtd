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
  use checkpoint
  use mpiworld
  use reglist
  use outlist
  use grid
  use parse
  use fdtd
  use mat

  implicit none
  private
  save

  M4_MODHEAD_DECL({DIAGEBAL},100,{

  integer :: ns, ne, dn  ! time stepping 

  ! spatially integrated energy contributions
  real(kind=8) :: dudt, ds, res
  real(kind=8) :: je(MAXEBALCH), kh(MAXEBALCH), jedisp, jeabs, khdisp, khabs
  real(kind=8) :: dsx, dsy, dsz
  ! spatially and time integrated energy contributions
  real(kind=8) :: sumdudt, sumds, sumje, sumkh, sumres
  real(kind=8) :: sumdsx, sumdsy, sumdsz
  
  logical, pointer, dimension(:,:,:) :: mask
  
  ! partial contributions to energy terms
  real(kind=8), dimension(MAXEBALCH) :: skh1, skh2, sje1, sje2
  real(kind=8), dimension(3) :: en
  real(kind=8), dimension(3) :: dsx1, dsx2
  real(kind=8), dimension(3) :: dsy1, dsy2
  real(kind=8), dimension(3) :: dsz1, dsz2

  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagEBalObj(funit,lcount)

    M4_MODREAD_DECL({DIAGEBAL}, funit,lcount,diag,reg,out)
    integer :: v(3)

    M4_WRITE_DBG(". enter ReadMatEBalObj")
    
    M4_MODREAD_EXPR({DIAGEBAL}, funit,lcount,diag,reg,0,out, {

    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)
    
    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. diag%dn .lt. 1 ) then
       M4_FATAL_ERROR({"BAD TIME WINDOW!"})
    end if

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
    diag%dsx = 0.
    diag%dsy = 0.
    diag%dsz = 0.
    diag%je = 0.
    diag%kh = 0.
    diag%res = 0.

    diag%sumdudt = 0.
    diag%sumds = 0.
    diag%sumdsx = 0.
    diag%sumdsy = 0.
    diag%sumdsz = 0.
    diag%sumje = 0.
    diag%sumkh = 0.
    diag%sumres = 0.

    reg = regobj(diag%regidx)

    M4_WRITE_DBG("allocate mask")
    allocate(diag%mask(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX), stat = err)
    M4_ALLOC_ERROR(err,"InitializeDiagEBal")
 
    M4_WRITE_DBG("initialize mask")
    call SetMaskRegObj(reg,diag%mask,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX)

! load from checkpoint file

    if ( load_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

       read(UNITCHK) diag%en, diag%skh1, diag%skh2, diag%sje1, diag%sje2
       read(UNITCHK) diag%dsx1, diag%dsx2, diag%dsy1, diag%dsy2, diag%dsz1, diag%dsz2
       read(UNITCHK) diag%sumdudt, diag%sumje, diag%sumkh
       read(UNITCHK) diag%sumds, diag%sumdsx, diag%sumdsy, diag%sumdsz

    end if

    M4_IFELSE_DBG({call EchoDiagEBalObj(diag)})

    })
    M4_WRITE_DBG(". exit InitializeMatEBal")

  end subroutine InitializeDiagEBal

!----------------------------------------------------------------------

  subroutine FinalizeDiagEBal



    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_WRITE_DBG(". enter FinalizeMatEBal")
    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

! save to checkpoint file

    if ( save_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

       write(UNITCHK) diag%en, diag%skh1, diag%skh2, diag%sje1, diag%sje2
       write(UNITCHK) diag%dsx1, diag%dsx2, diag%dsy1, diag%dsy2, diag%dsz1, diag%dsz2
       write(UNITCHK) diag%sumdudt, diag%sumje, diag%sumkh
       write(UNITCHK) diag%sumds, diag%sumdsx, diag%sumdsy, diag%sumdsz

    end if

    deallocate(diag%mask)

    })
    M4_WRITE_DBG(". exit FinalizeMatEBal")

  end subroutine FinalizeDiagEBal

!----------------------------------------------------------------------

  subroutine StepHDiagEBal(ncyc)

    integer :: ncyc, m, mo, moo, idx
    real(kind=8) :: sum
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m    = mod(ncyc*2+2,3) + 1

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag

       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       M4_MODOBJ_GETREG(diag,reg)


       diag%dsx1(m) = 0.
       diag%dsx2(m) = 0.
       diag%dsy1(m) = 0.
       diag%dsy2(m) = 0. 
       diag%dsz1(m) = 0.
       diag%dsz2(m) = 0.

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{
   
      ! poynting vector
 
       diag%dsx1(m) = diag%dsx1(m) &
            + M4_VOLEY(i,j,k)/M4_SX(i,j,k) * ( dble(Ey(i,j,k)) * dble( Hz(i,j,k) - Hz(i-1,j,k) ) ) &
            - M4_VOLEZ(i,j,k)/M4_SX(i,j,k) * ( dble(Ez(i,j,k)) * dble( Hy(i,j,k) - Hy(i-1,j,k) ) )

       diag%dsx2(m) = diag%dsx2(m) &
            + M4_VOLHZ(i,j,k)/M4_HSX(i,j,k) * ( dble(Hz(i,j,k)) * dble( Ey(i+1,j,k) - Ey(i,j,k) ) ) &
            - M4_VOLHY(i,j,k)/M4_HSX(i,j,k) * ( dble(Hy(i,j,k)) * dble( Ez(i+1,j,k) - Ez(i,j,k) ) )

M4_IFELSE_1D({},{


       diag%dsy1(m) = diag%dsy1(m) &
            + M4_VOLEZ(i,j,k)/M4_SY(i,j,k) * ( dble(Ez(i,j,k)) * dble( Hx(i,j,k) - Hx(i,j-1,k ) ) ) &
            - M4_VOLEX(i,j,k)/M4_SY(i,j,k) * ( dble(Ex(i,j,k)) * dble( Hz(i,j,k) - Hz(i,j-1,k) ) )

       diag%dsy2(m) = diag%dsy2(m) &
            + M4_VOLHX(i,j,k)/M4_HSY(i,j,k) * ( dble(Hx(i,j,k)) * dble( Ez(i,j+1,k) - Ez(i,j,k) ) ) &
            - M4_VOLHZ(i,j,k)/M4_HSY(i,j,k) * ( dble(Hz(i,j,k)) * dble( Ex(i,j+1,k) - Ex(i,j,k) ) )

})

M4_IFELSE_3D({

       diag%dsz1(m) = diag%dsz1(m) &
            + M4_VOLEX(i,j,k)/M4_SZ(i,j,k) * ( dble(Ex(i,j,k)) * dble( Hy(i,j,k) - Hy(i,j,k-1) ) ) &
            - M4_VOLEY(i,j,k)/M4_SZ(i,j,k) * ( dble(Ey(i,j,k)) * dble( Hx(i,j,k) - Hx(i,j,k-1) ) )

       diag%dsz2(m) = diag%dsz2(m) &
            + M4_VOLHY(i,j,k)/M4_HSZ(i,j,k) * ( dble(Hy(i,j,k)) * dble( Ex(i,j,k+1) - Ex(i,j,k) ) ) &
            - M4_VOLHX(i,j,k)/M4_HSZ(i,j,k) * ( dble(Hx(i,j,k)) * dble( Ey(i,j,k+1) - Ey(i,j,k) ) )

})

       })

       ! material sources

       diag%sje1 = 0.
       call SumJEMat(diag%mask,ncyc,diag%sje1,idx,.true.)
       diag%skh1 = 0.
       call SumKHMat(diag%mask,ncyc,diag%skh1,idx,.true.)

    })
  
  end subroutine StepHDiagEBal


!----------------------------------------------------------------------


  subroutine StepEDiagEBal(ncyc)

    integer :: ncyc, m, mo, moo, idx
    real(kind=8) :: sum
    M4_MODLOOP_DECL({DIAGEBAL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m    = mod(ncyc*2+3,3) + 1
    mo   = mod(ncyc*2+2,3) + 1 ! half timestep back
    moo  = mod(ncyc*2+1,3) + 1 ! full timestep back

    M4_MODLOOP_EXPR({DIAGEBAL},diag,{

       ! this loops over all diag structures, setting diag
       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       diag%en(m) = 0.
       diag%dsx1(m) = 0.
       diag%dsx2(m) = 0.
       diag%dsy1(m) = 0.
       diag%dsy2(m) = 0. 
       diag%dsz1(m) = 0.
       diag%dsz2(m) = 0.
   
       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       ! energy density

       diag%en(m) =  diag%en(m) + ( &
            M4_VOLEX(i,j,k) / epsinvx(i,j,k) * dble(Ex(i,j,k))*dble(Ex(i,j,k)) + &
            M4_VOLEY(i,j,k) / epsinvy(i,j,k) * dble(Ey(i,j,k))*dble(Ey(i,j,k)) + &
            M4_VOLEZ(i,j,k) / epsinvz(i,j,k) * dble(Ez(i,j,k))*dble(Ez(i,j,k)) + &
            M4_VOLHX(i,j,k) / M4_MUINVX(i,j,k) * dble(Hx(i,j,k))*dble(Hx(i,j,k)) + &
            M4_VOLHY(i,j,k) / M4_MUINVY(i,j,k) * dble(Hy(i,j,k))*dble(Hy(i,j,k)) + &
            M4_VOLHZ(i,j,k) / M4_MUINVZ(i,j,k) * dble(Hz(i,j,k))*dble(Hz(i,j,k)) &
            )
       
       ! poynting vector

    
       diag%dsx1(m) = diag%dsx1(m) &
            + M4_VOLEY(i,j,k)/M4_SX(i,j,k) * ( dble(Ey(i,j,k)) * dble( Hz(i,j,k) - Hz(i-1,j,k) ) ) &
            - M4_VOLEZ(i,j,k)/M4_SX(i,j,k) * ( dble(Ez(i,j,k)) * dble( Hy(i,j,k) - Hy(i-1,j,k) ) )

       diag%dsx2(m) = diag%dsx2(m) &
            + M4_VOLHZ(i,j,k)/M4_HSX(i,j,k) * ( dble(Hz(i,j,k)) * dble( Ey(i+1,j,k) - Ey(i,j,k) ) ) &
            - M4_VOLHY(i,j,k)/M4_HSX(i,j,k) * ( dble(Hy(i,j,k)) * dble( Ez(i+1,j,k) - Ez(i,j,k) ) )

M4_IFELSE_1D({},{


       diag%dsy1(m) = diag%dsy1(m) &
            + M4_VOLEZ(i,j,k)/M4_SY(i,j,k) * ( dble(Ez(i,j,k)) * dble( Hx(i,j,k) - Hx(i,j-1,k ) ) ) &
            - M4_VOLEX(i,j,k)/M4_SY(i,j,k) * ( dble(Ex(i,j,k)) * dble( Hz(i,j,k) - Hz(i,j-1,k) ) )

       diag%dsy2(m) = diag%dsy2(m) &
            + M4_VOLHX(i,j,k)/M4_HSY(i,j,k) * ( dble(Hx(i,j,k)) * dble( Ez(i,j+1,k) - Ez(i,j,k) ) ) &
            - M4_VOLHZ(i,j,k)/M4_HSY(i,j,k) * ( dble(Hz(i,j,k)) * dble( Ex(i,j+1,k) - Ex(i,j,k) ) )

})

M4_IFELSE_3D({

       diag%dsz1(m) = diag%dsz1(m) &
            + M4_VOLEX(i,j,k)/M4_SZ(i,j,k) * ( dble(Ex(i,j,k)) * dble( Hy(i,j,k) - Hy(i,j,k-1) ) ) &
            - M4_VOLEY(i,j,k)/M4_SZ(i,j,k) * ( dble(Ey(i,j,k)) * dble( Hx(i,j,k) - Hx(i,j,k-1) ) )

       diag%dsz2(m) = diag%dsz2(m) &
            + M4_VOLHY(i,j,k)/M4_HSZ(i,j,k) * ( dble(Hy(i,j,k)) * dble( Ex(i,j,k+1) - Ex(i,j,k) ) ) &
            - M4_VOLHX(i,j,k)/M4_HSZ(i,j,k) * ( dble(Hx(i,j,k)) * dble( Ey(i,j,k+1) - Ey(i,j,k) ) )

})
      
       })

       ! material sources
       diag%sje2 = 0.
       call SumJEMat(diag%mask,ncyc,diag%sje2,idx,.false.)
       diag%skh2 = 0.
       call SumKHMat(diag%mask,ncyc,diag%skh2,idx,.false.)
       
       ! ---------------------------------------------------
       
        
!       diag%dudt = .5/DT * ( diag%en(m) - diag%en(mo) )
!       diag%divs = .5 * ( diag%divs1 + diag%divs2(mo) - diag%divs3 - diag%divs4 ) 
!       diag%jekh = .5 * ( diag%q1 + diag%q2(mo) + diag%q3 + diag%q4 ) 

       diag%dudt = .5/DT * ( diag%en(m) - diag%en(moo) )

       diag%dsx = .5 * ( diag%dsx1(m) + diag%dsx1(mo) + diag%dsx2(mo) + diag%dsx2(moo) )
       diag%dsy = .5 * ( diag%dsy1(m) + diag%dsy1(mo) + diag%dsy2(mo) + diag%dsy2(moo) )
       diag%dsz = .5 * ( diag%dsz1(m) + diag%dsz1(mo) + diag%dsz2(mo) + diag%dsz2(moo) )

       diag%je = .5 * ( diag%sje1 + diag%sje2 )
       diag%kh = .5 * ( diag%skh1 + diag%skh2 )

       diag%jedisp = 0.
       diag%jeabs = 0.
       diag%khdisp = 0.
       diag%khabs = 0.

       do i = 1,idx-1,NUMEBALCH
          diag%jedisp = diag%jedisp + diag%je(i) 
          diag%jeabs = diag%jeabs + diag%je(i+1) + diag%je(i+2)
          diag%khdisp = diag%khdisp + diag%kh(i) 
          diag%khabs = diag%khabs + diag%kh(i+1) + diag%kh(i+2)
       enddo

       diag%dudt = diag%dudt + diag%jedisp + diag%khdisp

       diag%ds   = diag%dsx + diag%dsy + diag%dsz 
       diag%res = - ( diag%dudt + diag%ds + diag%jeabs + diag%khabs )
       
       diag%sumdudt = diag%sumdudt + DT * diag%dudt
       diag%sumje = diag%sumje + DT * diag%jeabs
       diag%sumkh = diag%sumkh + DT * diag%khabs
       diag%sumds = diag%sumds + DT * diag%ds
       diag%sumdsx = diag%sumdsx + DT * diag%dsx
       diag%sumdsy = diag%sumdsy + DT * diag%dsy
       diag%sumdsz = diag%sumdsz + DT * diag%dsz
            
       diag%sumres = - ( diag%sumdudt + diag%sumds + diag%sumje + diag%sumkh )

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
! Changed : 7/07/2011 S.Wuestner
!
! =====================================================================


