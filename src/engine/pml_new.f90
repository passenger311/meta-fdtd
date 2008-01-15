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
  integer :: pmlpart           ! flag (inner/outer partition)
  integer :: countpml          ! number of active pml sheets
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

  M4_FTYPE, allocatable, dimension(:,:,:,:) :: Fx,Fy,Fz

contains

!----------------------------------------------------------------------


  subroutine ReadConfigPml(funit,string)

    integer :: funit
    character(len=*) :: string
      
    integer :: ios, i

    M4_WRITE_DBG({". enter ReadConfigPml"})

    M4_WRITE_DBG({"received token: ", TRIM(string)})
    if ( string .ne. "(PML" ) then
       M4_FATAL_ERROR({"BAD SECTION IDENTIFIER: ReadConfigPml"})
    endif

    read(UNITTMP,*) PMLMAX
    M4_WRITE_DBG({"read PMLMAX: ",  PMLMAX})
    read(UNITTMP,*) potpml
    M4_WRITE_DBG({"read potpml: ",  potpml})
    read(UNITTMP,*) sigmamax
    M4_WRITE_DBG({"read sigmamax: ",  sigmamax})
    read(UNITTMP,*) kappamax 
    M4_WRITE_DBG({"read kappamax: ",  kappamax})
    read(UNITTMP,*,iostat=ios) string 
    M4_WRITE_DBG({"read terminator: ", TRIM(string)})

    ! TODO: add some checks on numerical values

    if ( string(1:1) .ne. ")" ) then
       M4_FATAL_ERROR({"BAD SECTION TERMINATOR: ReadConfigPml"})
    endif

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
    countpml = 0
    do i = 1, 6
       if ( planebound(i) .eq. num ) then 
          planepml(i) = 1
          countpml = countpml + 1
       end if
    end do

    M4_WRITE_DBG({"planepml : ",(planepml(i), i = 1,6)})

    ! modify fdtd core ranges for PML sheets [IBIG,IEIG]
    
    ! the k planes are set 0 if in 2D mode

    if(planepml(1) .eq. 1) IBIG=IBEG+PMLMAX
    if(planepml(2) .eq. 1) IEIG=IEND-PMLMAX
    if(planepml(3) .eq. 1) JBIG=JBEG+PMLMAX
    if(planepml(4) .eq. 1) JEIG=JEND-PMLMAX
    if(planepml(5) .eq. 1) KBIG=KBEG+PMLMAX
    if(planepml(6) .eq. 1) KEIG=KEND-PMLMAX
    
    M4_WRITE_DBG({"set IBIG/IEIG: ", IBIG, IEIG})
    M4_WRITE_DBG({"set KBIG/KEIG: ", JBIG, JEIG})
    M4_WRITE_DBG({"set KBIG/KEIG: ", KBIG, KEIG})

    if((IEIG+1.le.IBIG) & 
M4_IFELSE_1D({},{    .or. (JEIG+1.le.JBIG) &    }) 
M4_IFELSE_3D({       .or. (KEIG+1.le.KBIG) &    })
    ) then
       write(STDERR,*) "OPPOSITE PML LAYERS OVERLAP!"
       stop
    endif    

    call AllocateFields
    call CalcCoefficients(IBEG,IEND, IBIG,IEIG+1, cexpml, cmxpml)
    call CalcCoefficients(JBEG,JEND, JBIG,JEIG+1, ceypml, cmypml)
    call CalcCoefficients(KBEG,KEND, KBIG,KEIG+1, cezpml, cmzpml)

    if ( countpml .ne. 0 ) call DisplayPml

    M4_WRITE_DBG({". exit InitializePml"})

  contains
    
    subroutine AllocateFields

      integer :: err, i, j, k, l

      M4_WRITE_DBG({". enter InitializePml.AllocateFields"})
      
      ! numeric coefficient-fields
      
      allocate(cexpml(1:4,IBEG:IEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(ceypml(1:4,JBEG:JEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(cezpml(1:4,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(cmxpml(1:4,IBEG:IEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(cmypml(1:4,JBEG:JEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(cmzpml(1:4,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      cexpml=1.0
      ceypml=1.0
      cezpml=1.0
      cmxpml=1.0
      cmypml=1.0
      cmzpml=1.0

      !  B and D fields for each of the 6 layers 
      
      allocate(BE1(1:3,IBEG:IBIG,JBEG:JEND,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE1(1:3,IBEG:IBIG,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE2(1:3,IEIG+1:IEND,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE2(1:3,IEIG+1:IEND,JBEG:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(Fx(IEIG+1:IEIG+1,JBEG:JEND,KBEG:KEND,3), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      BE1 = 0.0
      BE2 = 0.0
      DE1 = 0.0
      DE2 = 0.0
      Fx = 0.0
      
M4_IFELSE_1D({},{      
      allocate(BE3(1:3,IBIG:IEIG,JBEG:JBIG,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(DE3(1:3,IBIG:IEIG,JBEG:JBIG,KBEG:KEND), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE4(1:3,IBIG:IEIG,JEIG+1:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE4(1:3,IBIG:IEIG,JEIG+1:JEND,KBEG:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(Fy(IBIG:IEIG,JEIG+1:JEIG+1,KBEG:KEND,3), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      BE3 = 0.0
      BE4 = 0.0
      DE3 = 0.0
      DE4 = 0.0
      Fy = 0.0
})
      
M4_IFELSE_3D({      
      allocate(BE5(1:3,IBIG:IEIG,JBIG:JEIG,KBEG:KBIG), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE5(1:3,IBIG:IEIG,JBIG:JEIG,KBEG:KBIG), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE6(1:3,IBIG:IEIG,JBIG:JEIG,KEIG+1:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(DE6(1:3,IBIG:IEIG,JBIG:JEIG,KEIG+1:KEND), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(Fz(IBIG:IEIG,JBIG:JEIG,KEIG+1:KEIG+1,3), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      BE5 = 0.0
      BE6 = 0.0
      DE5 = 0.0
      DE6 = 0.0
      Fz = 0.0
})      
      


      M4_WRITE_DBG({". exit InitializePml.AllocateFields"})
      
    end subroutine AllocateFields
    
    subroutine CalcCoefficients(lbeg, lend, ls, le, ce, cm)
      
      integer, intent(in) :: ls, le, lbeg, lend
      real(kind=8), dimension(1:4,lbeg:lend), intent(inout) :: ce, cm
      real(kind=8), dimension(-1:PMLMAX-1) :: val1,val1p,val2,val2p
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
      
      ! calculate coefficients
      
      do l=0, PMLMAX-1
         ce(1,l+le) = val1p(l)
         ce(3,l+le) = val2p(l)
         ce(2,l+le) = 1.0 / val1(l)
         ce(4,l+le) = val2(l) / val1(l)
         
         cm(1,l+le) = val1(l)
         cm(3,l+le) = val2(l)
         cm(2,l+le) = 1.0 / val1p(l)
         cm(4,l+le) = val2p(l) / val1p(l)

         ce(1,ls-l) = ce(1,l+le)
         ce(3,ls-l) = ce(3,l+le)
         ce(2,ls-l) = ce(2,l+le) 
         ce(4,ls-l) = ce(4,l+le)
         
         cm(1,ls-l-1) = cm(1,l+le)
         cm(3,ls-l-1) = cm(3,l+le)
         cm(2,ls-l-1) = cm(2,l+le)
         cm(4,ls-l-1) = cm(4,l+le)

      enddo

      M4_WRITE_DBG({". exit InitializePml.CalcCoefficients"})
      
    end subroutine CalcCoefficients
      
  end subroutine InitializePml

!----------------------------------------------------------------------

  subroutine FinalizePml

    M4_WRITE_DBG({". FinalizePml"})

M4_IFELSE_3D({
    deallocate(DE6)
    deallocate(BE6)
    deallocate(DE5)
    deallocate(BE5)
    deallocate(Fx)
})
M4_IFELSE_1D({},{
    deallocate(DE4)
    deallocate(BE4)
    deallocate(DE3)
    deallocate(BE3)
    deallocate(Fy)
})
    deallocate(DE2)
    deallocate(BE2)
    deallocate(DE1)
    deallocate(BE1)
    deallocate(Fz)

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
       call SaveRestorePml(IBIG,IBIG,  JBEG,JEND,KBEG,KEND,  Hy,Hz,Ex, Fx, Ey,Ez,Hx )
       call DoStepHPml(IBEG,IBIG,JBEG,JEND,KBEG,KEND,BE1)
       call SaveRestorePml(IBIG,IBIG,  JBEG,JEND,KBEG,KEND,  Ey,Ez,Hx, Fx, Hy,Hz,Ex )
    case ( 2 )
       call DoStepHPml(IEIG+1,IEND,JBEG,JEND,KBEG,KEND,BE2)
   case ( 3 )
       call SaveRestorePml(IBIG,IEIG,  JBIG,JBIG  ,KBEG,KEND, Hx,Hz,Ey, Fy, Ex,Ez,Hy)
       call DoStepHPml(IBIG,IEIG,JBEG,JBIG,KBEG,KEND,BE3)
       call SaveRestorePml(IBIG,IEIG,  JBIG,JBIG  ,KBEG,KEND, Ex,Ez,Hy, Fy, Hx,Hz,Ey )
    case ( 4  )
       call DoStepHPml(IBIG,IEIG,JEIG+1,JEND,KBEG,KEND,BE4)
    case ( 5  )
       call SaveRestorePml(IBIG,IEIG,JBIG,JEIG,  KBIG,KBIG, Hx,Hy,Ez, Fz, Ex,Ey,Hz )
       call DoStepHPml(IBIG,IEIG,JBIG,JEIG,KBEG,KBIG,BE5)
       call SaveRestorePml(IBIG,IEIG,JBIG,JEIG,  KBIG,KBIG, Ex,Ey,Hz, Fz, Hx,Hy,Ez )
    case ( 6  )
       call DoStepHPml(IBIG,IEIG,JBIG,JEIG,KEIG+1,KEND,BE6)
    end select
    
  contains
    
    subroutine DoStepHPml(is,ie,js,je,ks,ke,B)
      
      integer is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:3,is:ie,js:je,ks:ke) :: B
      integer i, j, k
      M4_FTYPE :: Bxo, Byo, Bzo, Exh, Eyh, Ezh
      real(kind=8) dtx, dty, dtz
      
      dtx = DT/SX
      dty = DT/SY
      dtz = DT/SZ   

!      M4_WRITE_DBG({"STEP PML H: ", is, ie, js, je, ks, ke})
        
M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
      do k=ks, ke     
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Exh,Eyh,Ezh,Bxo,Byo,Bzo)})
            do i=is, ie
               
               Exh = Ex(i,j,k)
               Eyh = Ey(i,j,k)
               Ezh = Ez(i,j,k)
               
               Bxo = B(1,i,j,k)
               Byo = B(2,i,j,k) 
               Bzo = B(3,i,j,k)
               
               ! Update B
               
               B(1,i,j,k) = cmypml(4,j)*B(1,i,j,k) + cmypml(2,j) * ( &
M4_IFELSE_1D({0.&},{ - dty*( Ez(i,j+1,k) - Ezh )           &     })
M4_IFELSE_3D({       + dtz*( Ey(i,j,k+1) - Eyh )           &     })
                     )
               B(2,i,j,k) = cmzpml(4,k)*B(2,i,j,k) + cmzpml(2,k) * ( &
M4_IFELSE_3D({       - dtz*( Ex(i,j,k+1) - Exh )           &     })
                     + dtx*( Ez(i+1,j,k) - Ezh )           &
                     )
               B(3,i,j,k) = cmxpml(4,i)*B(3,i,j,k) + cmxpml(2,i) * ( &
                     - dtx*( Ey(i+1,j,k) - Eyh )           &
M4_IFELSE_1D({},{    + dty*( Ex(i,j+1,k) - Exh )           &     })
                     )        
               
               ! Calc H
               
               Hx(i,j,k) = cmzpml(4,k)*Hx(i,j,k) + cmzpml(2,k) * &
                    M4_MUINVX(i,j,k)*(cmxpml(1,i)*B(1,i,j,k) - cmxpml(3,i)*Bxo)
               Hy(i,j,k) = cmxpml(4,i)*Hy(i,j,k) + cmxpml(2,i) * &
                    M4_MUINVY(i,j,k)*(cmypml(1,j)*B(2,i,j,k) - cmypml(3,j)*Byo)
               Hz(i,j,k) = cmypml(4,j)*Hz(i,j,k) + cmypml(2,j) * &
                    M4_MUINVZ(i,j,k)*(cmzpml(1,k)*B(3,i,j,k) - cmzpml(3,k)*Bzo)           
               
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
       call SaveRestorePml(IBIG,IBIG,  JBEG,JEND,KBEG,KEND,  Hy,Hz,Ex, Fx, Ey,Ez,Hx )
       call DoStepEPml(IBEG,IBIG,  JBEG,JEND,KBEG,KEND,DE1)
       call SaveRestorePml(IBIG,IBIG,  JBEG,JEND,KBEG,KEND,  Ey,Ez,Hx, Fx, Hy,Hz,Ex )
    case ( 2 ) 
       call DoStepEPml(IEIG+1,IEND,  JBEG,JEND,KBEG,KEND,DE2)
    case ( 3 ) 
       call SaveRestorePml(IBIG,IEIG,  JBIG,JBIG  ,KBEG,KEND, Hx,Hz,Ey, Fy, Ex,Ez,Hy)
       call DoStepEPml(IBIG,IEIG,  JBEG,JBIG  ,KBEG,KEND,DE3)
       call SaveRestorePml(IBIG,IEIG,  JBIG,JBIG  ,KBEG,KEND, Ex,Ez,Hy, Fy, Hx,Hz,Ey )
    case ( 4 )
       call DoStepEPml(IBIG,IEIG,  JEIG+1,JEND  ,KBEG,KEND,DE4)
    case ( 5 ) 
       call SaveRestorePml(IBIG,IEIG,JBIG,JEIG,  KBIG,KBIG, Hx,Hy,Ez, Fz, Ex,Ey,Hz )
       call DoStepEPml(IBIG,IEIG,JBIG,JEIG,  KBEG,KBIG,  DE5)
       call SaveRestorePml(IBIG,IEIG,JBIG,JEIG,  KBIG,KBIG, Ex,Ey,Hz, Fz, Hx,Hy,Ez )
    case ( 6 ) 
       call DoStepEPml(IBIG,IEIG,JBIG,JEIG,  KEIG+1,KEND,  DE6)
    end select

    call StepEBoundPec(i) ! need to set electric conductor bcs

  contains
    
    subroutine DoStepEPml(is,ie,js,je,ks,ke,D)
      
      integer :: is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:3,is:ie,js:je,ks:ke) :: D
      
      integer :: i, j, k
      real(kind=8) :: dtx, dty, dtz
      M4_FTYPE :: Dxo, Dyo, Dzo, Hxh, Hyh, Hzh    
  
      dtx = DT/SX
      dty = DT/SY
      dtz = DT/SZ


!      M4_WRITE_DBG({"STEP PML E: ", is, ie, js, je, ks, ke})
  
M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
      do k=ks, ke
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
            do i=is, ie

               Hxh = Hx(i,j,k)
               Hyh = Hy(i,j,k)
               Hzh = Hz(i,j,k)
           
               Dxo = D(1,i,j,k)
               Dyo = D(2,i,j,k)
               Dzo = D(3,i,j,k)
           
               ! Update D

               D(1,i,j,k) =  ceypml(4,j)*D(1,i,j,k) + ceypml(2,j) * ( &
M4_IFELSE_1D({0.&},{ + dty * ( Hzh - Hz(i,j-1,k) )           &    })
M4_IFELSE_3D({       - dtz * ( Hyh - Hy(i,j,k-1) )           &    })
                     )
               D(2,i,j,k) =  cezpml(4,k)*D(2,i,j,k) + cezpml(2,k) * ( &
M4_IFELSE_3D({       + dtz * ( Hxh - Hx(i,j,k-1) )           &    })
                     - dtx * ( Hzh - Hz(i-1,j,k) )           &
                    )
               D(3,i,j,k) =  cexpml(4,i)*D(3,i,j,k) + cexpml(2,i) * ( &
                    + dtx * ( Hyh - Hy(i-1,j,k) )           &
M4_IFELSE_1D({},{   - dty * ( Hxh - Hx(i,j-1,k) )           &    })
                    )

               ! Calc E

               Ex(i,j,k) = cezpml(4,k)*Ex(i,j,k) + cezpml(2,k) * &
                    epsinvx(i,j,k)*(cexpml(1,i)*D(1,i,j,k) - cexpml(3,i)*Dxo)
               Ey(i,j,k) = cexpml(4,i)*Ey(i,j,k) + cexpml(2,i) * &
                    epsinvy(i,j,k)*(ceypml(1,j)*D(2,i,j,k) - ceypml(3,j)*Dyo)
               Ez(i,j,k) = ceypml(4,j)*Ez(i,j,k) + ceypml(2,j) * &
                    epsinvz(i,j,k)*(cezpml(1,k)*D(3,i,j,k) - cezpml(3,k)*Dzo)
           
            enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO})
         enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO})
      enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO})

!     M4_WRITE_DBG({"done"})
 
    end subroutine DoStepEPml


  end subroutine StepEBoundPml


  !----------------------------------------------------------------------
 
  ! restore fields f1r,f2r,f3r from buffer f and store f1s,f2s,f3s inside instead

    subroutine SaveRestorePml(is,ie,js,je,ks,ke, f1s, f2s, f3s, F, f1r, f2r,f3r )

      integer :: is, ie, js, je, ks, ke
      M4_FTYPE, dimension(M4_RANGE(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)) :: f1s,f2s,f3s,f1r,f2r,f3r
      M4_FTYPE, dimension(is:ie,js:je,ks:ke,3) :: f
      integer :: i,j,k

      M4_IFELSE_3D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
      do k=ks, ke
M4_IFELSE_2D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
         do j=js, je
M4_IFELSE_1D({!$OMP PARALLEL DO PRIVATE(Hxh,Hyh,Hzh,Dxo,Dyo,Dzo,epsinvx,epsinvy,epsinvz)}) 
            do i=is, ie

               f1r(i,j,k) = f(i,j,k,1)
               f2r(i,j,k) = f(i,j,k,2)
               f3r(i,j,k) = f(i,j,k,3)

               f(i,j,k,1) = f1s(i,j,k)
               f(i,j,k,2) = f2s(i,j,k)
               f(i,j,k,3) = f3s(i,j,k)
 
            enddo
M4_IFELSE_1D({!$OMP END PARALLEL DO})
         enddo
M4_IFELSE_2D({!$OMP END PARALLEL DO})
      enddo
M4_IFELSE_3D({!$OMP END PARALLEL DO})

  end subroutine SaveRestorePml


  !----------------------------------------------------------------------

  subroutine DisplayPml

    M4_WRITE_INFO({"[",TRIM(i2str(IBIG)),":",TRIM(i2str(IEIG)),"]",&
         "[",TRIM(i2str(JBIG)),":",TRIM(i2str(JEIG)),"]", &
         "[",TRIM(i2str(KBIG)),":",TRIM(i2str(KEIG)),"]"
    })
    	
  end subroutine DisplayPml

  !----------------------------------------------------------------------

end module pml

!
! Authors:  A.Klaedtke, S.Scholz, J.Hamm
! Modified: 4/12/2007
!
!======================================================================
