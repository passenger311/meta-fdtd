!-*- F90 -*------------------------------------------------------------
!
!  module: pml / meta
!
!  boundary conditions using uniform perfectly matched layers.
!
!  CF,1D,2D,3D
!
!----------------------------------------------------------------------

!======================================================================
!

module pml

  use constant
  use strings
  use grid  
  use fdtd
  use pec

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'PML'
  logical, private :: modconfigured = .false.

  ! --- Public Methods

  public :: ReadConfigPml
  public :: InitializePml
  public :: FinalizePml
  public :: StepEBoundPml
  public :: StepHBoundPml

  ! --- Public Data

  ! --- Constants

  ! --- Data

  integer :: planepml(6)       ! 1 if pml applies 0 if not
  integer :: pmlcount          ! pml count
  integer :: pmlpart           ! flag (inner/outer partition)
  integer :: PMLMAX            ! number of Pml layers
  real(kind=8) :: potpml       ! exponent for sigma and kappa
  real(kind=8) :: sigmamax     ! absorption coefficient 
  real(kind=8) :: kappamax     ! absorption of evanescent waves 

  ! numeric Pml coefficients

  real(kind=8), allocatable, dimension(:,:) :: cexpml,ceypml,cezpml
  real(kind=8), allocatable, dimension(:,:) :: cmxpml,cmypml,cmzpml

  ! auxilliary B and D fields on the 6 planes E1 .. E6

  M4_FTYPE, allocatable, dimension(:,:,:,:) :: BE1,BE2,BE3,BE4,BE5,BE6 
  M4_FTYPE, allocatable, dimension(:,:,:,:) :: DE1,DE2,DE3,DE4,DE5,DE6 

contains

!----------------------------------------------------------------------


  subroutine ReadConfigPml(funit,lcount,string)

    integer :: funit, lcount
    character(len=*) :: string

    character(len=LINELNG) :: line

    M4_WRITE_DBG({". enter ReadConfigPml"})

    if ( string .ne. "(PML" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigPml"})
    endif

    call readint(funit, lcount, PMLMAX)
    M4_WRITE_DBG({"read PMLMAX: ", PMLMAX})

    call readfloat(funit, lcount, potpml)
    M4_WRITE_DBG({"read potpml: ", potpml})

    call readfloat(funit, lcount, sigmamax)
    M4_WRITE_DBG({"read sigmamax: ", sigmamax})

    call readfloat(funit, lcount, kappamax)
    M4_WRITE_DBG({"read kappamax: ", kappamax})

    call readtoken(funit, lcount, ")PML")

    ! TODO: add some checks on numerical values

    modconfigured = .true.

    M4_WRITE_DBG({". exit ReadConfigPml"})

  end subroutine ReadConfigPml

!----------------------------------------------------------------------

  subroutine InitializePml(planebound, num)

    integer :: planebound(6), num
    integer :: i

    M4_WRITE_DBG({". enter InitializePml"})

    M4_WRITE_DBG({"got planebound(i): ",  (planebound(i),i=1, M4_SDIM*2)})

    if ( .not. modconfigured ) then

       M4_WRITE_WARN({"PMLS NOT CONFIGURED -> USING DEFAULT PARAMETERS!"})
       ! sigmamax (= SigmaOpt, see Tavlove 2, pp 286)
       ! It is Sig = Sig[SI] / (!c*eps0)
       PMLMAX = 8
       potpml = 3.2
       sigmamax = (real(POTPML)+1.0)*0.8/(3.0*DT)
       kappamax = 1.1

    end if

    planepml = 0
    pmlcount = 0
    do i = 1, 6
       if ( planebound(i) .eq. num ) then 
          planepml(i) = 1
          pmlcount = pmlcount + 1
       endif
    end do

    if ( pmlcount .eq. 0 ) return

    M4_WRITE_DBG({"planepml : ",(planepml(i), i = 1,6)})

    ! modify fdtd core ranges for PML sheets [IBIG,IEIG]
    
    ! the k planes are set 0 if in 2D mode

    if(planepml(1) .eq. 1) IBIG=IBEG+PMLMAX+1
    if(planepml(2) .eq. 1) IEIG=IMAX-PMLMAX-1
    if(planepml(3) .eq. 1) JBIG=JBEG+PMLMAX+1
    if(planepml(4) .eq. 1) JEIG=JMAX-PMLMAX-1
    if(planepml(5) .eq. 1) KBIG=KBEG+PMLMAX+1
    if(planepml(6) .eq. 1) KEIG=KMAX-PMLMAX-1
    
    M4_WRITE_DBG({"set IBIG/IEIG: ", IBIG, IEIG})
    M4_WRITE_DBG({"set JBIG/JEIG: ", JBIG, JEIG})
    M4_WRITE_DBG({"set KBIG/KEIG: ", KBIG, KEIG})

    if((IEIG+1 .le. IBIG) & 
M4_IFELSE_1D({},{    .or. (JEIG+1 .le. JBIG) &    }) 
M4_IFELSE_3D({       .or. (KEIG+1 .le. KBIG) &    })
    ) then
       write(STDERR,*) "OPPOSITE PML LAYERS OVERLAP!"
       stop
    endif    

    call AllocateFields
    call CalcCoefficients(IBEG, IEND, IBIG-1, IEIG+1, cexpml, cmxpml)
    call CalcCoefficients(JBEG, JEND, JBIG-1, JEIG+1, ceypml, cmypml)
    call CalcCoefficients(KBEG, KEND, KBIG-1, KEIG+1, cezpml, cmzpml)

    M4_WRITE_DBG({". exit InitializePml"})

  contains
    
    subroutine AllocateFields

      integer :: err, i, j, k, l

      M4_WRITE_DBG({". enter InitializePml.AllocateFields"})
      

      M4_WRITE_DBG({". allocating coefficient fields ",IBEG,IEND,JBEG,JEND,KBEG,KEND})
     
      ! numeric coefficient-fields
      

      allocate(cexpml(1:4,IBEG:IEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      allocate(cmxpml(1:4,IBEG:IEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      cexpml=1.0
      cmxpml=1.0

      allocate(ceypml(1:4,JBEG:JEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      allocate(cmypml(1:4,JBEG:JEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      ceypml=1.0
      cmypml=1.0

      allocate(cezpml(1:4,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      allocate(cmzpml(1:4,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
          
      cezpml=1.0
      cmzpml=1.0


      M4_WRITE_DBG({". allocating B/D fields"})
     
      !  B and D fields for each of the 6 layers 
      
      allocate(BE1(1:3,IBEG:IBIG-1,JBEG:JEND,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE1(1:3,IBEG:IBIG-1,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE2(1:3,IEIG+1:IEND,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE2(1:3,IEIG+1:IEND,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      BE1 = 0.0
      BE2 = 0.0
      DE1 = 0.0
      DE2 = 0.0
      
M4_IFELSE_1D({},{      
      allocate(BE3(1:3,IBIG:IEIG,JBEG:JBIG-1,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(DE3(1:3,IBIG:IEIG,JBEG:JBIG-1,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE4(1:3,IBIG:IEIG,JEIG+1:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE4(1:3,IBIG:IEIG,JEIG+1:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      BE3 = 0.0
      BE4 = 0.0
      DE3 = 0.0
      DE4 = 0.0
})
      
M4_IFELSE_3D({      
      allocate(BE5(1:3,IBIG:IEIG,JBIG:JEIG,KBEG:KBIG-1), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE5(1:3,IBIG:IEIG,JBIG:JEIG,KBEG:KBIG-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE6(1:3,IBIG:IEIG,JBIG:JEIG,KEIG+1:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(DE6(1:3,IBIG:IEIG,JBIG:JEIG,KEIG+1:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      BE5 = 0.0
      BE6 = 0.0
      DE5 = 0.0
      DE6 = 0.0
})      
      

      M4_WRITE_DBG({". exit InitializePml.AllocateFields"})
      
    end subroutine AllocateFields
    
    subroutine CalcCoefficients(lbeg, lend, ls, le, ce, cm)
      
      integer, intent(in) :: lbeg, lend, ls, le
      real(kind=8), dimension(1:4,lbeg:lend) :: ce, cm
      real(kind=8), dimension(-1:PMLMAX) :: val1,val1p,val2,val2p
      real(kind=8) :: x, sigma, kappa
      integer :: l

      M4_WRITE_DBG({". enter InitializePml.CalcCoefficients"})
      
      ! help values

      val1 = 1.0
      val2 = 1.0
      val1p = 1.0
      val2p = 1.0

      do l=0, PMLMAX-1
         
         ! points on the grid corners
         
         x = real(l) / real(PMLMAX-1)
         sigma = sigmamax*x**potpml
         kappa = 1.0+(kappamax-1.0)*x**potpml
         val1(l) = kappa+0.5*sigma*DT
         val2(l) = kappa-0.5*sigma*DT
         
         ! points in the grid center
         
         x = (real(l)+0.5) / real(PMLMAX-1)
         sigma = sigmamax*x**potpml
         kappa = 1.0+(kappamax-1.0)*x**potpml
         val1p(l) = kappa+0.5*sigma*DT
         val2p(l) = kappa-0.5*sigma*DT
  
      enddo
      


! coefficients for top pml sheet
      do l=0, lend - le
         ce(1,le+l) = val1p(l)
         ce(3,le+l) = val2p(l)
         ce(2,le+l) = 1.0 / val1(l)
         ce(4,le+l) = val2(l) / val1(l)
         
         cm(1,le+l) = val1(l)
         cm(3,le+l) = val2(l)
         cm(2,le+l) = 1.0 / val1p(l)
         cm(4,le+l) = val2p(l) / val1p(l)
      enddo


! coefficients for bottom pml sheet
      do l=0, ls - lbeg
         ce(1,ls-l) = val1p(l-1)
         ce(3,ls-l) = val2p(l-1)
         ce(2,ls-l) = 1.0 / val1(l)
         ce(4,ls-l) = val2(l) / val1(l)
        
         cm(1,ls-l) = val1(l)
         cm(3,ls-l) = val2(l)
         cm(2,ls-l) = 1.0 / val1p(l-1)
         cm(4,ls-l) = val2p(l-1) / val1p(l-1)
      enddo

! note that these are:
!         ce(1,ls) = 1.0
!         ce(3,ls) = 1.0
!         cm(2,ls) = 1.0
!         cm(4,ls) = 1.0

      M4_WRITE_DBG({". exit InitializePml.CalcCoefficients"})
      
    end subroutine CalcCoefficients
      
  end subroutine InitializePml

!----------------------------------------------------------------------

  subroutine FinalizePml

    M4_WRITE_DBG({". FinalizePml"})

    if ( pmlcount .eq. 0 ) return

M4_IFELSE_3D({
    deallocate(DE6)
    deallocate(BE6)
    deallocate(DE5)
    deallocate(BE5)
})

M4_IFELSE_1D({},{
    deallocate(DE4)
    deallocate(BE4)
    deallocate(DE3)
    deallocate(BE3)
})

    deallocate(DE2)
    deallocate(BE2)
    deallocate(DE1)
    deallocate(BE1)

    deallocate(cmzpml)
    deallocate(cmypml)
    deallocate(cmxpml)
    deallocate(cezpml)
    deallocate(ceypml)
    deallocate(cexpml)

  end subroutine FinalizePml

!----------------------------------------------------------------------

  ! update the H-fields of all Pml layers
  
  subroutine StepHBoundPml(i)
    
    integer :: i

    if ( i .gt. 2 .and. M4_IS1D ) return
    if ( i .gt. 4 .and. M4_IS2D ) return

    select case ( i )
    case ( 1 ) 
       call DoStepHPml(IBEG,IBIG-1,JBEG,JEND,KBEG,KEND,BE1)
    case ( 2 )
       call DoStepHPml(IEIG+1,IEND,JBEG,JEND,KBEG,KEND,BE2)
    case ( 3 )
       call DoStepHPml(IBIG,IEIG,JBEG,JBIG-1,KBEG,KEND,BE3)
    case ( 4  )
       call DoStepHPml(IBIG,IEIG,JEIG+1,JEND,KBEG,KEND,BE4)
    case ( 5  )
       call DoStepHPml(IBIG,IEIG,JBIG,JEIG,KBEG,KBIG-1,BE5)
    case ( 6  )
       call DoStepHPml(IBIG,IEIG,JBIG,JEIG,KEIG+1,KEND,BE6)
    end select
    
  contains
    
    subroutine DoStepHPml(is,ie,js,je,ks,ke,B)
      
      integer is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:3,is:ie,js:je,ks:ke) :: B
      integer i, j, k
      M4_FTYPE :: Bxo, Byo, Bzo, Exh, Eyh, Ezh
      
  
M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
      do k=ks, ke     
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
            do i=is, ie
  
M4_IFELSE_TE({
               Exh = Ex(i,j,k)
               Eyh = Ey(i,j,k)
})
M4_IFELSE_TM({             
               Ezh = Ez(i,j,k)
})              
 
M4_IFELSE_TM({
               Bxo = B(1,i,j,k)
               Byo = B(2,i,j,k) 
})
M4_IFELSE_TE({
               Bzo = B(3,i,j,k)
})               
               ! Update B
               
M4_IFELSE_TM({
               B(1,i,j,k) = cmypml(4,j)*B(1,i,j,k) + cmypml(2,j) * ( &
M4_IFELSE_1D({0.&},{ - DT/M4_HSY(i,j,k)*( Ez(i,j+1,k) - Ezh )           &     })
M4_IFELSE_3D({       + DT/M4_HSZ(i,j,k)*( Ey(i,j,k+1) - Eyh )           &     })
                     )
               B(2,i,j,k) = cmzpml(4,k)*B(2,i,j,k) + cmzpml(2,k) * ( &
M4_IFELSE_3D({       - DT/M4_HSZ(i,j,k)*( Ex(i,j,k+1) - Exh )           &     })
                     + DT/M4_HSX(i,j,k)*( Ez(i+1,j,k) - Ezh )           &
                     )
})

M4_IFELSE_TE({
               B(3,i,j,k) = cmxpml(4,i)*B(3,i,j,k) + cmxpml(2,i) * ( &
                     - DT/M4_HSX(i,j,k)*( Ey(i+1,j,k) - Eyh )           &
M4_IFELSE_1D({},{    + DT/M4_HSY(i,j,k)*( Ex(i,j+1,k) - Exh )           &     })
                     )        
             
})  
               ! Calc H
               
M4_IFELSE_TM({
               Hx(i,j,k) = cmzpml(4,k)*Hx(i,j,k) + cmzpml(2,k) * &
                    M4_MUINVX(i,j,k)*(cmxpml(1,i)*B(1,i,j,k) - cmxpml(3,i)*Bxo)
               Hy(i,j,k) = cmxpml(4,i)*Hy(i,j,k) + cmxpml(2,i) * &
                    M4_MUINVY(i,j,k)*(cmypml(1,j)*B(2,i,j,k) - cmypml(3,j)*Byo)
})
M4_IFELSE_TE({
               Hz(i,j,k) = cmypml(4,j)*Hz(i,j,k) + cmypml(2,j) * &
                    M4_MUINVZ(i,j,k)*(cmzpml(1,k)*B(3,i,j,k) - cmzpml(3,k)*Bzo)           
})               
            enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO})
         enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO})
      enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO})

    end subroutine DoStepHPml
    
  end subroutine StepHBoundPml

!----------------------------------------------------------------------

  ! update the E-fields of all pml layers
  
  subroutine StepEBoundPml(i)
    
    integer :: i
    
    if ( i .gt. 2 .and. M4_IS1D ) return
    if ( i .gt. 4 .and. M4_IS2D ) return

    select case ( i )
    case ( 1 ) 
       call DoStepEPml(IBEG,IBIG-1,JBEG,JEND,KBEG,KEND,DE1)
    case ( 2 ) 
       call DoStepEPml(IEIG+1,IEND,JBEG,JEND,KBEG,KEND,DE2)
    case ( 3 ) 
       call DoStepEPml(IBIG,IEIG,JBEG,JBIG-1,KBEG,KEND,DE3)
    case ( 4 )
       call DoStepEPml(IBIG,IEIG,JEIG+1,JEND,KBEG,KEND,DE4)
    case ( 5 ) 
       call DoStepEPml(IBIG,IEIG,JBIG,JEIG,KBEG,KBIG-1,DE5)
    case ( 6 ) 
       call DoStepEPml(IBIG,IEIG,JBIG,JEIG,KEIG+1,KEND,DE6)
    end select

    call StepEBoundPec(i) ! need to set electric conductor bcs

  contains
    
    subroutine DoStepEPml(is,ie,js,je,ks,ke,D)
      
      integer :: is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:3,is:ie,js:je,ks:ke) :: D
      
      integer :: i, j, k
      M4_FTYPE :: Dxo, Dyo, Dzo, Hxh, Hyh, Hzh    
  

M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
      do k=ks, ke
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
            do i=is, ie

M4_IFELSE_TM({
               Hxh = Hx(i,j,k)
               Hyh = Hy(i,j,k)
})
M4_IFELSE_TE({
               Hzh = Hz(i,j,k)
})          
 
M4_IFELSE_TE({
               Dxo = D(1,i,j,k)
               Dyo = D(2,i,j,k)
})
M4_IFELSE_TM({
               Dzo = D(3,i,j,k)
})           
               ! Update D

 
M4_IFELSE_TE({

               D(1,i,j,k) =  ceypml(4,j)*D(1,i,j,k) + ceypml(2,j) * ( &
M4_IFELSE_1D({0.&},{ + DT/M4_SY(i,j,k) * ( Hzh - Hz(i,j-1,k) )           &    })
M4_IFELSE_3D({       - DT/M4_SZ(i,j,k) * ( Hyh - Hy(i,j,k-1) )           &    })
                     )
               D(2,i,j,k) =  cezpml(4,k)*D(2,i,j,k) + cezpml(2,k) * ( &
M4_IFELSE_3D({       + DT/M4_SZ(i,j,k) * ( Hxh - Hx(i,j,k-1) )           &    })
                     - DT/M4_SX(i,j,k) * ( Hzh - Hz(i-1,j,k) )           &
                    )
})

M4_IFELSE_TM({
               D(3,i,j,k) =  cexpml(4,i)*D(3,i,j,k) + cexpml(2,i) * ( &
                    + DT/M4_SX(i,j,k) * ( Hyh - Hy(i-1,j,k) )           &
M4_IFELSE_1D({},{   - DT/M4_SY(i,j,k) * ( Hxh - Hx(i,j-1,k) )           &    })
                    )
})
               ! Calc E

M4_IFELSE_TE({
               Ex(i,j,k) = cezpml(4,k)*Ex(i,j,k) + cezpml(2,k) * &
                    epsinvx(i,j,k)*(cexpml(1,i)*D(1,i,j,k) - cexpml(3,i)*Dxo)
               Ey(i,j,k) = cexpml(4,i)*Ey(i,j,k) + cexpml(2,i) * &
                    epsinvy(i,j,k)*(ceypml(1,j)*D(2,i,j,k) - ceypml(3,j)*Dyo)
})
M4_IFELSE_TM({
               Ez(i,j,k) = ceypml(4,j)*Ez(i,j,k) + ceypml(2,j) * &
                    epsinvz(i,j,k)*(cezpml(1,k)*D(3,i,j,k) - cezpml(3,k)*Dzo)
})

            enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO})
         enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO})
      enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO})
 
    end subroutine DoStepEPml
    
  end subroutine StepEBoundPml


!----------------------------------------------------------------------

end module pml

!
! Authors:  A.Klaedtke, S.Scholz, J.Hamm
! Modified: 1/5/2008
!
!======================================================================
