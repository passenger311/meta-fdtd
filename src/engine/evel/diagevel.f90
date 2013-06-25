!-*- F90 -*------------------------------------------------------------
!
!  module: diagevel / meta
!
!  Energy balance diagnostics module.
!
!----------------------------------------------------------------------


! =====================================================================
!
! The DiagEvel module the center of energy associated with the
! given volume.
!
!

module diagevel

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

  M4_MODHEAD_DECL({DIAGEVEL},100,{

  integer :: ns, ne, dn  ! time stepping 

  ! spatially and time integrated energy contributions
  real(kind=8) :: u, ux, uy, uz, ux2, uy2, uz2
  
  logical, pointer, dimension(:,:,:) :: mask
  
  })

contains

!----------------------------------------------------------------------

  subroutine ReadDiagEvelObj(funit,lcount)

    M4_MODREAD_DECL({DIAGEVEL}, funit,lcount,diag,reg,out)
    integer :: v(3)

    M4_WRITE_DBG(". enter ReadMatEvelObj")
    
    M4_MODREAD_EXPR({DIAGEVEL}, funit,lcount,diag,reg,0,out, {

    call readints(funit,lcount,v,3) 
    diag%ns = v(1)
    diag%ne = v(2)
    diag%dn = v(3)
    
    if ( diag%ns .ge. diag%ne .or. diag%ns .lt. 0 .or. diag%dn .lt. 1 ) then
       M4_FATAL_ERROR({"BAD TIME WINDOW!"})
    end if

    })

    M4_WRITE_DBG(". exit ReadMatEvelObj")

  end subroutine ReadDiagEvelObj

!----------------------------------------------------------------------

  subroutine InitializeDiagEvel

    type (T_REG) :: reg
    integer :: err
    M4_MODLOOP_DECL({DIAGEVEL},diag)

    M4_WRITE_DBG(". enter InitializeMatEvel")
    M4_MODLOOP_EXPR({DIAGEVEL},diag,{
    
    diag%u = 0.
    diag%ux = 0.
    diag%uy = 0.
    diag%uz = 0.
    diag%ux2 = 0.
    diag%uy2 = 0.
    diag%uz2 = 0.

    reg = regobj(diag%regidx)

    M4_WRITE_DBG("allocate mask")
    allocate(diag%mask(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX), stat = err)
    M4_ALLOC_ERROR(err,"InitializeDiagEvel")
 
    M4_WRITE_DBG("initialize mask")
    call SetMaskRegObj(reg,diag%mask,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX)

! load from checkpoint file

    if ( load_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

       read(UNITCHK) diag%u, diag%ux, diag%uy, diag%uz, diag%ux2, diag%uy2, diag%uz2    

    end if

    M4_IFELSE_DBG({call EchoDiagEvelObj(diag)})

    })
    M4_WRITE_DBG(". exit InitializeMatEvel")

  end subroutine InitializeDiagEvel

!----------------------------------------------------------------------

  subroutine FinalizeDiagEvel



    M4_MODLOOP_DECL({DIAGEVEL},diag)
    M4_WRITE_DBG(". enter FinalizeMatEvel")
    M4_MODLOOP_EXPR({DIAGEVEL},diag,{

! save to checkpoint file

    if ( save_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

      write(UNITCHK) diag%u, diag%ux, diag%uy, diag%uz, diag%ux2, diag%uy2, diag%uz2    

    end if

    deallocate(diag%mask)

    })
    M4_WRITE_DBG(". exit FinalizeMatEvel")

  end subroutine FinalizeDiagEvel

!----------------------------------------------------------------------

  subroutine StepHDiagEvel(ncyc)

    integer :: ncyc

    return
  
  end subroutine StepHDiagEvel


!----------------------------------------------------------------------


  subroutine StepEDiagEvel(ncyc)

    integer :: ncyc, m, idx
    real(kind=8) :: u
    M4_MODLOOP_DECL({DIAGEVEL},diag)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(0))

    m    = mod(ncyc*2+3,3) + 1

    M4_MODLOOP_EXPR({DIAGEVEL},diag,{

       ! this loops over all diag structures, setting diag
       if (  ncyc .lt. diag%ns .or. ncyc .gt. diag%ne ) cycle

       diag%u = 1D-300
       diag%ux = 0
       diag%uy = 0
       diag%uz = 0
       diag%ux2 = 0
       diag%uy2 = 0
       diag%uz2 = 0
   
       M4_MODOBJ_GETREG(diag,reg)
       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{

       ! energy density

       u = ( &
            M4_VOLEX(i,j,k) / epsinvx(i,j,k) * dble(Ex(i,j,k))*dble(Ex(i,j,k)) + &
            M4_VOLEY(i,j,k) / epsinvy(i,j,k) * dble(Ey(i,j,k))*dble(Ey(i,j,k)) + &
            M4_VOLEZ(i,j,k) / epsinvz(i,j,k) * dble(Ez(i,j,k))*dble(Ez(i,j,k)) + &
            M4_VOLHX(i,j,k) / M4_MUINVX(i,j,k) * dble(Hx(i,j,k))*dble(Hx(i,j,k)) + &
            M4_VOLHY(i,j,k) / M4_MUINVY(i,j,k) * dble(Hy(i,j,k))*dble(Hy(i,j,k)) + &
            M4_VOLHZ(i,j,k) / M4_MUINVZ(i,j,k) * dble(Hz(i,j,k))*dble(Hz(i,j,k)) &
            )

       diag%u =  diag%u + u
       diag%ux =  diag%ux + u*i 
       diag%uy =  diag%uy + u*j
       diag%uz =  diag%uz + u*k 
       diag%ux2 =  diag%ux2 + u*i*i 
       diag%uy2 =  diag%uy2 + u*j*j 
       diag%uz2 =  diag%uz2 + u*k*k 
      
       })

       diag%ux =  diag%ux/diag%u
       diag%uy =  diag%uy/diag%u
       diag%uz =  diag%uz/diag%u

       diag%ux2 =  (((diag%ux2/diag%u) - (diag%ux**2)))**0.5
       diag%uy2 =  (((diag%uy2/diag%u) - (diag%uy**2)))**0.5
       diag%uz2 =  (((diag%uz2/diag%u) - (diag%uz**2)))**0.5

       ! ---------------------------------------------------
       
       
            
       ! ---------------------------------------------------
 
      
    })

  end subroutine StepEDiagEvel

!----------------------------------------------------------------------

   subroutine EchoDiagEvelObj(diag)

    type(T_DIAGEVEL) :: diag
 
    M4_WRITE_INFO({"--- diagevel # ",& 
         TRIM(i2str(diag%idx))," ", TRIM(diag%type)})

    ! -- write parameters to console 
    M4_WRITE_INFO({"ns ne dn = ",diag%ns, diag%ne, diag%dn })

    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(diag%regidx))
    

  end subroutine EchoDiagEvelObj
  
!----------------------------------------------------------------------

end module diagevel

!
! Authors:  J.Hamm, T. Pickering
! Modified: 23/11/2011
!
! =====================================================================


