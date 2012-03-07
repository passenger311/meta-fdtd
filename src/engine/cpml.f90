!-*- F90 -*------------------------------------------------------------
!
!  module: cpml / meta
!
!  boundary conditions using uniform perfectly matched layers.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------

!======================================================================
!

module cpml

  use constant
  use strings
  use checkpoint
  use grid  
  use fdtd
  use pec

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'CPML'
  logical, private :: modconfigured = .false.

  ! --- Public Methods

  public :: ReadConfigCPml
  public :: InitializeCPml
  public :: FinalizeCPml
  public :: StepEBoundCPml
  public :: StepHBoundCPml

  ! --- Public Data

  ! --- Constants

  ! --- Data

  integer :: planecpml(6)       ! 1 if cpml applies 0 if not
  integer :: cpmlcount          ! cpml count
  integer :: cpmlpart           ! flag (inner/outer partition)
  integer :: CPMLMAX            ! number of CPml layers
  real(kind=8) :: potcpml       ! exponent for sigma and kappa
  real(kind=8) :: sigmamax      ! absorption coefficient 
  real(kind=8) :: kappamax      ! stretching parameter of grid 
  real(kind=8) :: alphamax      ! absorption of evanescent waves
  real(kind=8) :: alphapot      ! exponent for alpha 


  ! numeric CPml coefficients

  real(kind=8), allocatable, dimension(:,:) :: cexcpml,ceycpml,cezcpml
  real(kind=8), allocatable, dimension(:,:) :: chxcpml,chycpml,chzcpml

  ! auxilliary Psi fields on the 6 planes PsiE1 .. PsiE6
  
  M4_FTYPE, allocatable, dimension(:,:,:,:) :: PsiE1,PsiE2,PsiE3,PsiE4,PsiE5,PsiE6
  M4_FTYPE, allocatable, dimension(:,:,:,:) :: PsiH1,PsiH2,PsiH3,PsiH4,PsiH5,PsiH6

contains

!----------------------------------------------------------------------


  subroutine ReadConfigCPml(funit,lcount,string)

    integer :: funit, lcount
    character(len=*) :: string

    character(len=LINELNG) :: line

    M4_WRITE_DBG({". enter ReadConfigCPml"})

    if ( string .ne. "(CPML" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigCPml"})
    endif

    call readint(funit, lcount, CPMLMAX)
    M4_WRITE_DBG({"read CPMLMAX: ", CPMLMAX})

    call readfloat(funit, lcount, potcpml)
    M4_WRITE_DBG({"read potcpml: ", potcpml})

    call readfloat(funit, lcount, sigmamax)
    M4_WRITE_DBG({"read sigmamax: ", sigmamax})

    call readfloat(funit, lcount, kappamax)
    M4_WRITE_DBG({"read kappamax: ", kappamax})

    call readfloat(funit, lcount, alphamax)
    M4_WRITE_DBG({"read alphamax: ", alphamax})

    call readfloat(funit, lcount, alphapot)
    M4_WRITE_DBG({"read alphapot: ", alphapot})

    call readtoken(funit, lcount, ")CPML")

    ! TODO: add some checks on numerical values

    modconfigured = .true.

    M4_WRITE_DBG({". exit ReadConfigCPml"})

  end subroutine ReadConfigCPml

!----------------------------------------------------------------------

  subroutine InitializeCPml(planebound, num)

    integer :: planebound(6), num
    integer :: i

    M4_WRITE_DBG({". enter InitializeCPml"})
    
    M4_WRITE_DBG({"got planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    if ( .not. modconfigured ) then

       M4_WRITE_WARN({"CPMLS NOT CONFIGURED -> USING DEFAULT PARAMETERS!"})
       ! sigmamax (= SigmaOpt, see Tavlove 2, pp 286)
       ! It is Sig = Sig[SI] / (!c*eps0)
       CPMLMAX = 8
       potcpml = 3.2
       sigmamax = (real(POTCPML)+1.0)*0.8/(3.0*DT)
       kappamax = 1.1
       alphamax = 0.0
       alphapot = 1.0
    end if

    planecpml = 0
    cpmlcount = 0
    do i = 1, 6
       if ( planebound(i) .eq. num ) then 
          planecpml(i) = 1
          cpmlcount = cpmlcount + 1
       endif
    end do

    if ( cpmlcount .eq. 0 ) return

    M4_WRITE_DBG({"planecpml : ",(planecpml(i), i = 1,6)})

    ! modify fdtd core ranges for CPML sheets [IBIG,IEIG]
    
    ! the k planes are set 0 if in 2D mode

    if(planecpml(1) .eq. 1) IBIG=IBEG+CPMLMAX+1
    if(planecpml(2) .eq. 1) IEIG=IMAX-CPMLMAX-1
    if(planecpml(3) .eq. 1) JBIG=JBEG+CPMLMAX+1
    if(planecpml(4) .eq. 1) JEIG=JMAX-CPMLMAX-1
    if(planecpml(5) .eq. 1) KBIG=KBEG+CPMLMAX+1
    if(planecpml(6) .eq. 1) KEIG=KMAX-CPMLMAX-1
    
    M4_WRITE_DBG({"set IBIG/IEIG: ", IBIG, IEIG})
    M4_WRITE_DBG({"set JBIG/JEIG: ", JBIG, JEIG})
    M4_WRITE_DBG({"set KBIG/KEIG: ", KBIG, KEIG})

    if((IEIG+1 .le. IBIG) & 
M4_IFELSE_1D({},{    .or. (JEIG+1 .le. JBIG) &    }) 
M4_IFELSE_3D({       .or. (KEIG+1 .le. KBIG) &    })
    ) then
       write(STDERR,*) "OPPOSITE CPML LAYERS OVERLAP!"
       stop
    endif    

    call AllocateFields
    call CalcCoefficients(IBEG, IEND, IBIG-1, IEIG+1, cexcpml, chxcpml)
    call CalcCoefficients(JBEG, JEND, JBIG-1, JEIG+1, ceycpml, chycpml)
    call CalcCoefficients(KBEG, KEND, KBIG-1, KEIG+1, cezcpml, chzcpml)

! load from checkpoint file

   if ( load_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

      read(UNITCHK) PsiE1, PsiH1, PsiE2, PsiH2
M4_IFELSE_1D({},{  
      read(UNITCHK) PsiE3, PsiH3, PsiE4, PsiH4
})
M4_IFELSE_3D({     
      read(UNITCHK) PsiE5, PsiH5, PsiE6, PsiH6
})

   end if

    M4_WRITE_DBG({". exit InitializeCPml"})

  contains
    
    subroutine AllocateFields

      integer :: err, i, j, k, l

      M4_WRITE_DBG({". enter InitializeCPml.AllocateFields"})
      

      M4_WRITE_DBG({". allocating coefficient fields ",IBEG,IEND,JBEG,JEND,KBEG,KEND})
     
      ! numeric coefficient-fields
      

      allocate(cexcpml(1:3,IBEG:IEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      allocate(chxcpml(1:3,IBEG:IEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      cexcpml=1.0
      cexcpml(3,:)=0.0
      chxcpml=1.0
      chxcpml(3,:)=0.0

      allocate(ceycpml(1:3,JBEG:JEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      allocate(chycpml(1:3,JBEG:JEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      ceycpml=1.0
      ceycpml(3,:)=0.0
      chycpml=1.0
      chycpml(3,:)=0.0

      allocate(cezcpml(1:3,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      allocate(chzcpml(1:3,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
          
      cezcpml=1.0
      cezcpml(3,:)=0.0
      chzcpml=1.0
      chzcpml(3,:)=0.0


      M4_WRITE_DBG({". allocating PsiE/PsiH fields"})
     
      !  B and D fields for each of the 6 layers 
      
      allocate(PsiH1(1:6,IBEG:IBIG-1,JBEG:JEND,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiE1(1:6,IBEG:IBIG-1,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiH2(1:6,IEIG+1:IEND,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiE2(1:6,IEIG+1:IEND,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      PsiH1 = 0.0
      PsiH2 = 0.0
      PsiE1 = 0.0
      PsiE2 = 0.0
      
M4_IFELSE_1D({},{      
      allocate(PsiH3(1:6,IBIG:IEIG,JBEG:JBIG-1,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(PsiE3(1:6,IBIG:IEIG,JBEG:JBIG-1,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiH4(1:6,IBIG:IEIG,JEIG+1:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiE4(1:6,IBIG:IEIG,JEIG+1:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      PsiH3 = 0.0
      PsiH4 = 0.0
      PsiE3 = 0.0
      PsiE4 = 0.0
})
      
M4_IFELSE_3D({      
      allocate(PsiH5(1:6,IBIG:IEIG,JBIG:JEIG,KBEG:KBIG-1), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiE5(1:6,IBIG:IEIG,JBIG:JEIG,KBEG:KBIG-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(PsiH6(1:6,IBIG:IEIG,JBIG:JEIG,KEIG+1:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(PsiE6(1:6,IBIG:IEIG,JBIG:JEIG,KEIG+1:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      PsiH5 = 0.0
      PsiH6 = 0.0
      PsiE5 = 0.0
      PsiE6 = 0.0
})      
      

      M4_WRITE_DBG({". exit InitializeCPml.AllocateFields"})
      
    end subroutine AllocateFields
    
    subroutine CalcCoefficients(lbeg, lend, ls, le, ce, ch)
      
      integer, intent(in) :: lbeg, lend, ls, le
      real(kind=8), dimension(1:3,lbeg:lend) :: ce, ch
      real(kind=8), dimension(-1:CPMLMAX) :: sigma, kappa, alpha, sigmap, kappap, alphap!val1,val1p,val2,val2p
      real(kind=8) :: x
      integer :: l

      M4_WRITE_DBG({". enter InitializeCPml.CalcCoefficients"})
      
      ! help values

      sigma = 0.0
      sigmap = 0.0
      kappa = 1.0
      kappap = 1.0
      alpha = alphamax
      alphap = alphamax

      do l=0, CPMLMAX-1
         
         ! points on the grid corners
         
         x = real(l) / real(CPMLMAX-1)
         sigma(l) = sigmamax*x**potcpml
         kappa(l) = 1.0+(kappamax-1.0)*x**potcpml
         alpha(l) = alphamax*(1.0-x)**alphapot

         ! points in the grid center
         
         x = (real(l)+0.5) / real(CPMLMAX-1)
         sigmap(l) = sigmamax*x**potcpml
         kappap(l) = 1.0+(kappamax-1.0)*x**potcpml
         if (x<1.0) then
            alphap(l) = alphamax*(1.0-x)**alphapot
         else
            alphap(l) = 0.0
         endif

      enddo

! coefficients for top cpml sheet
      do l=0, lend - le
         ce(1,le+l) = 1.0/kappa(l)
         ce(2,le+l) = dexp(-(sigma(l)+kappa(l)*alpha(l))*DT/kappa(l))
         if (sigma(l).eq.0.0) then
!write(*,*) "top sheet e", l, "sigma", sigma(l), "alpha", alpha(l)
            ce(3,le+l) = 0.0
         else
            ce(3,le+l) = sigma(l)/(sigma(l)+kappa(l)*alpha(l))/kappa(l)*(ce(2,le+l)-1)
         endif
         
         ch(1,le+l) = 1.0/kappap(l)
         ch(2,le+l) = dexp(-(sigmap(l)+kappap(l)*alphap(l))*DT/kappap(l))
         ch(3,le+l) = sigmap(l)/(sigmap(l)+kappap(l)*alphap(l))/kappap(l)*(ch(2,le+l)-1)
         
      enddo

! coefficients for bottom cpml sheet
      do l=0, ls - lbeg
         ce(1,ls-l) = 1.0/kappa(l)
         ce(2,ls-l) = dexp(-(sigma(l)+kappa(l)*alpha(l))*DT/kappa(l))
         if (sigma(l).eq.0.0) then
!write(*,*) "bottom sheet e", l,"sigma", sigma(l), "alpha", alpha(l)
            ce(3,ls-l) = 0.0
         else
            ce(3,ls-l) = sigma(l)/(sigma(l)+kappa(l)*alpha(l))/kappa(l)*(ce(2,ls-l)-1)
         endif

         ch(1,ls-l) = 1.0/kappap(l-1)
         ch(2,ls-l) = dexp(-(sigmap(l-1)+kappap(l-1)*alphap(l-1))*DT/kappap(l-1))
         if (sigmap(l-1).eq.0) then
!write(*,*) "bottom sheet h", l-1, "sigma", sigmap(l-1), "alpha", alphap(l-1)
            ch(3,ls-l) = 0.0
         else
            ch(3,ls-l) = sigmap(l-1)/(sigmap(l-1)+kappap(l-1)*alphap(l-1))/kappap(l-1)*(ch(2,ls-l)-1)
         endif

      enddo
!write(40,"(E15.7)") sigma
!write(41,"(E15.7)") sigmap
!write(42,"(E15.7)") kappa
!write(43,"(E15.7)") kappap
!write(44,"(E15.7)") alpha
!write(45,"(E15.7)") alphap
!write(46,"(3E15.7)") ce
!write(47,"(3E15.7)") ch


      M4_WRITE_DBG({". exit InitializeCPml.CalcCoefficients"})
      
    end subroutine CalcCoefficients
      
  end subroutine InitializeCPml

!----------------------------------------------------------------------

  subroutine FinalizeCPml

    M4_WRITE_DBG({". FinalizeCPml"})

    if ( cpmlcount .eq. 0 ) return

! save checkpoint file

 if ( save_state .and. ( detail_level .eq. 1 .or. detail_level .eq. 3 ) ) then

      write(UNITCHK) PsiE1, PsiH1, PsiE2, PsiH2
M4_IFELSE_1D({},{  
      write(UNITCHK) PsiE3, PsiH3, PsiE4, PsiH4
})
M4_IFELSE_3D({     
      write(UNITCHK) PsiE5, PsiH5, PsiE6, PsiH6
})

   end if

M4_IFELSE_3D({
    deallocate(PsiE6)
    deallocate(PsiH6)
    deallocate(PsiE5)
    deallocate(PsiH5)
})

M4_IFELSE_1D({},{
    deallocate(PsiE4)
    deallocate(PsiH4)
    deallocate(PsiE3)
    deallocate(PsiH3)
})

    deallocate(PsiE2)
    deallocate(PsiH2)
    deallocate(PsiE1)
    deallocate(PsiH1)

    deallocate(chzcpml)
    deallocate(chycpml)
    deallocate(chxcpml)
    deallocate(cezcpml)
    deallocate(ceycpml)
    deallocate(cexcpml)

  end subroutine FinalizeCPml

!----------------------------------------------------------------------

  ! update the H-fields of all CPml layers
  
  subroutine StepHBoundCPml(i)
    
    integer :: i

    if ( i .gt. 2 .and. M4_IS1D ) return
    if ( i .gt. 4 .and. M4_IS2D ) return

    select case ( i )
    case ( 1 ) 
       call DoStepHCPml(IBEG,IBIG-1,JBEG,JEND,KBEG,KEND,PsiH1)
    case ( 2 )
       call DoStepHCPml(IEIG+1,IEND,JBEG,JEND,KBEG,KEND,PsiH2)
    case ( 3 )
       call DoStepHCPml(IBIG,IEIG,JBEG,JBIG-1,KBEG,KEND,PsiH3)
    case ( 4  )
       call DoStepHCPml(IBIG,IEIG,JEIG+1,JEND,KBEG,KEND,PsiH4)
    case ( 5  )
       call DoStepHCPml(IBIG,IEIG,JBIG,JEIG,KBEG,KBIG-1,PsiH5)
    case ( 6  )
       call DoStepHCPml(IBIG,IEIG,JBIG,JEIG,KEIG+1,KEND,PsiH6)
    end select
    
  contains
    
    subroutine DoStepHCPml(is,ie,js,je,ks,ke,PsiH)
      
      integer is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:6,is:ie,js:je,ks:ke) :: PsiH
      integer i, j, k
  
M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
      do k=ks, ke     
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
            do i=is, ie


               ! Update PsiH
M4_IFELSE_TE({

M4_IFELSE_1D({},{
               PsiH(1,i,j,k) = chycpml(2,j)*PsiH(1,i,j,k) + chycpml(3,j) * &
                     ( ( Ez(i,j+1,k) - Ez(i,j,k) )/M4_SY(i,j,k) )
})
M4_IFELSE_3D({
               PsiH(2,i,j,k) = chzcpml(2,k)*PsiH(2,i,j,k) + chzcpml(3,k) * &
                     ( ( Ey(i,j,k+1) - Ey(i,j,k) )/M4_SZ(i,j,k) )
               PsiH(3,i,j,k) = chzcpml(2,k)*PsiH(3,i,j,k) + chzcpml(3,k) * &
                     ( ( Ex(i,j,k+1) - Ex(i,j,k) )/M4_SZ(i,j,k) )
})
               PsiH(4,i,j,k) = chxcpml(2,i)*PsiH(4,i,j,k) + chxcpml(3,i) * &
                     ( ( Ez(i+1,j,k) - Ez(i,j,k) )/M4_SX(i,j,k) )

})
M4_IFELSE_TM({

               PsiH(5,i,j,k) = chxcpml(2,i)*PsiH(5,i,j,k) + chxcpml(3,i) * &
                     ( ( Ey(i+1,j,k) - Ey(i,j,k) )/M4_SX(i,j,k) )
M4_IFELSE_1D({},{
               PsiH(6,i,j,k) = chycpml(2,j)*PsiH(6,i,j,k) + chycpml(3,j) * &
                     ( ( Ex(i,j+1,k) - Ex(i,j,k) )/M4_SY(i,j,k) )
})

})

               ! Calc H
M4_IFELSE_TE({

               Hx(i,j,k) = Hx(i,j,k) - DT*M4_MUINVX(i,j,k) * ( &
M4_IFELSE_1D({0.&},{ chycpml(1,j)/M4_SY(i,j,k)*( Ez(i,j+1,k) - Ez(i,j,k) ) &
                     + PsiH(1,i,j,k)                                     &     })
M4_IFELSE_3D({       - chzcpml(1,k)/M4_SZ(i,j,k)*( Ey(i,j,k+1) - Ey(i,j,k) ) &
                     - PsiH(2,i,j,k)                                     &     })
                     )
               Hy(i,j,k) = Hy(i,j,k) - DT*M4_MUINVY(i,j,k) * ( &
M4_IFELSE_3D({      chzcpml(1,k)/M4_SZ(i,j,k)*( Ex(i,j,k+1) - Ex(i,j,k) ) &
                    + PsiH(3,i,j,k)                                      &     })
                    - chxcpml(1,i)/M4_SX(i,j,k)*( Ez(i+1,j,k) - Ez(i,j,k) ) &
                    - PsiH(4,i,j,k) &
                    )

})
M4_IFELSE_TM({

               Hz(i,j,k) = Hz(i,j,k) - DT*M4_MUINVZ(i,j,k) * ( &
                    chxcpml(1,i)/M4_SX(i,j,k)*( Ey(i+1,j,k) - Ey(i,j,k) ) &
                    + PsiH(5,i,j,k) &
M4_IFELSE_1D({},{   - chycpml(1,j)/M4_SY(i,j,k)*( Ex(i,j+1,k) - Ex(i,j,k) ) &
                    - PsiH(6,i,j,k)                                      &     })
                    )

})
              
            enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO})
         enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO})
      enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO})

    end subroutine DoStepHCPml
    
  end subroutine StepHBoundCPml

!----------------------------------------------------------------------

  ! update the E-fields of all cpml layers
  
  subroutine StepEBoundCPml(i)
    
    integer :: i
    
    if ( i .gt. 2 .and. M4_IS1D ) return
    if ( i .gt. 4 .and. M4_IS2D ) return

    select case ( i )
    case ( 1 ) 
       call DoStepECPml(IBEG,IBIG-1,JBEG,JEND,KBEG,KEND,PsiE1)
    case ( 2 ) 
       call DoStepECPml(IEIG+1,IEND,JBEG,JEND,KBEG,KEND,PsiE2)
    case ( 3 ) 
       call DoStepECPml(IBIG,IEIG,JBEG,JBIG-1,KBEG,KEND,PsiE3)
    case ( 4 )
       call DoStepECPml(IBIG,IEIG,JEIG+1,JEND,KBEG,KEND,PsiE4)
    case ( 5 ) 
       call DoStepECPml(IBIG,IEIG,JBIG,JEIG,KBEG,KBIG-1,PsiE5)
    case ( 6 ) 
       call DoStepECPml(IBIG,IEIG,JBIG,JEIG,KEIG+1,KEND,PsiE6)
    end select

    call StepEBoundPec(i) ! need to set electric conductor bcs

  contains
    
    subroutine DoStepECPml(is,ie,js,je,ks,ke,PsiE)
      
      integer :: is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:6,is:ie,js:je,ks:ke) :: PsiE
      
      integer :: i, j, k

M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
      do k=ks, ke
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
            do i=is, ie


               ! Update PsiE
M4_IFELSE_TM({

M4_IFELSE_1D({},{
               PsiE(1,i,j,k) = ceycpml(2,j)*PsiE(1,i,j,k) + ceycpml(3,j) * &
                      ( ( Hz(i,j,k) - Hz(i,j-1,k) )/M4_SY(i,j,k) )
})
M4_IFELSE_3D({
               PsiE(2,i,j,k) = cezcpml(2,k)*PsiE(2,i,j,k) + cezcpml(3,k) * &
                      ( ( Hy(i,j,k) - Hy(i,j,k-1) )/M4_SZ(i,j,k) )
               PsiE(3,i,j,k) = cezcpml(2,k)*PsiE(3,i,j,k) + cezcpml(3,k) * &
                      ( ( Hx(i,j,k) - Hx(i,j,k-1) )/M4_SZ(i,j,k) )
})
               PsiE(4,i,j,k) = cexcpml(2,i)*PsiE(4,i,j,k) + cexcpml(3,i) * &
                      ( ( Hz(i,j,k) - Hz(i-1,j,k) )/M4_SX(i,j,k) )

})
M4_IFELSE_TE({

               PsiE(5,i,j,k) = cexcpml(2,i)*PsiE(5,i,j,k) + cexcpml(3,i) * &
                     ( ( Hy(i,j,k) - Hy(i-1,j,k) )/M4_SX(i,j,k) )
M4_IFELSE_1D({},{
               PsiE(6,i,j,k) = ceycpml(2,j)*PsiE(6,i,j,k) + ceycpml(3,j) * &
                     ( ( Hx(i,j,k) - Hx(i,j-1,k) )/M4_SY(i,j,k) )
})

})

               ! Calc E

M4_IFELSE_TM({

               Ex(i,j,k) = Ex(i,j,k) + DT*epsinvx(i,j,k) * ( &
M4_IFELSE_1D({0.&},{ ceycpml(1,j)/M4_SY(i,j,k)*( Hz(i,j,k) - Hz(i,j-1,k) ) &
                     + PsiE(1,i,j,k)                                     &     })
M4_IFELSE_3D({       - cezcpml(1,k)/M4_SZ(i,j,k)*( Hy(i,j,k) - Hy(i,j,k-1) ) &
                     - PsiE(2,i,j,k)                                     &     })
                    )
               Ey(i,j,k) = Ey(i,j,k) + DT*epsinvy(i,j,k) * ( &
M4_IFELSE_3D({      cezcpml(1,k)/M4_SZ(i,j,k)*( Hx(i,j,k) - Hx(i,j,k-1) ) &
                    + PsiE(3,i,j,k)                                &     },{0.&})
                    - cexcpml(1,i)/M4_SX(i,j,k)*( Hz(i,j,k) - Hz(i-1,j,k) ) &
                    - PsiE(4,i,j,k) &
                    )

})
M4_IFELSE_TE({

               Ez(i,j,k) = Ez(i,j,k) + DT*epsinvz(i,j,k) * ( &
                    cexcpml(1,i)/M4_SX(i,j,k)*( Hy(i,j,k) - Hy(i-1,j,k) ) &
                    + PsiE(5,i,j,k) &
M4_IFELSE_1D({},{   - ceycpml(1,j)/M4_SY(i,j,k)*( Hx(i,j,k) - Hx(i,j-1,k) ) &
                    - PsiE(6,i,j,k)                                     &      })
                    )

})


            enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO})
         enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO})
      enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO})
 
    end subroutine DoStepECPml
    
  end subroutine StepEBoundCPml


!----------------------------------------------------------------------

end module cpml

!
! Authors:  S. Wuestner, J.Hamm
! Modified: 17/2/2012
!
!======================================================================
