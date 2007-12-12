!-*- F90 -*------------------------------------------------------------
!
!  module: pml / meta3
!
!  boundary conditions using uniform perfectly matched layers.
!
!  subs:
!
!    InitializePml
!      ReadConfig
!      Initialize
!      AllocateFields
!      CalcCoefficients
!    FinalizePml 
!    StepHPml 
!      DoStepHPml
!    StepEPml
!      DoStepEPml
!
!----------------------------------------------------------------------

!======================================================================
!
! m4 macro-preprocessor runs over this file and replaces
! M4 FTYPE -> real(kind=8) or complex(kind=8)
!

module pml

  use constant
  use strings
  use grid  
  use fdtd

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'PML'
  logical, private :: modconfigured = .false.
  logical, private :: modinitialized = .false.

  ! --- Public Methods

  public :: ReadConfigPml
  public :: InitializePml
  public :: FinalizePml
  public :: StepEPml
  public :: StepHPml

  ! --- Public Data

  ! --- Constants

  ! --- Data

  integer :: pmlpart           ! flag (inner/outer partition)
  integer :: PMLMAX            ! number of Pml layers
  real(kind=8) :: PotPml       ! exponent for sigma and kappa
  real(kind=8) :: SigmaMax     ! absorption coefficient 
  real(kind=8) :: KappaMax     ! absorption of evanescent waves 

  ! do planes 1-6 have a Pml (yes=1, no=0)

  integer, dimension(6) :: planepml=(/ 1,1,1,1,1,1 /)

  ! numeric Pml coefficients

  real(kind=8), allocatable, dimension(:,:) :: cexpml,ceypml,cezpml
  real(kind=8), allocatable, dimension(:,:) :: cmxpml,cmypml,cmzpml

  ! auxilliary B and D fields on the 6 planes E1 .. E6

  M4_FTYPE, allocatable, dimension(:,:,:,:) :: BE1,BE2,BE3,BE4,BE5,BE6 
  M4_FTYPE, allocatable, dimension(:,:,:,:) :: DE1,DE2,DE3,DE4,DE5,DE6 

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

    read(UNITTMP,*) (planepml(i),i=1, 6)
    M4_WRITE_DBG({"read planepml(i): ",  (planepml(i),i=1, 6)})
    read(UNITTMP,*) PMLMAX
    M4_WRITE_DBG({"read PMLMAX: ",  PMLMAX})
    read(UNITTMP,*) PotPml
    M4_WRITE_DBG({"read PotPml: ",  PotPml})
    read(UNITTMP,*) SigmaMax
    M4_WRITE_DBG({"read SigmaMax: ",  SigmaMax})
    read(UNITTMP,*) KappaMax 
    M4_WRITE_DBG({"read KappaMax: ",  KappaMax})
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

  subroutine InitializePml

    M4_WRITE_DBG({". enter InitializePml"})

    if ( .not. modconfigured ) then
       M4_WRITE_WARN({"NO PMLS CONFIGURED!"})
       return
    end if

    if(planepml(1) .eq. 1) ISIG=IBEG+PMLMAX
    if(planepml(2) .eq. 1) IEIG=IMAX-PMLMAX
    if(planepml(3) .eq. 1) JSIG=JBEG+PMLMAX
    if(planepml(4) .eq. 1) JEIG=JMAX-PMLMAX
    if(planepml(5) .eq. 1) KSIG=KBEG+PMLMAX
    if(planepml(6) .eq. 1) KEIG=KMAX-PMLMAX
    
    M4_WRITE_DBG({"modified ISIG/IEIG: ", ISIG, IEIG})
    M4_WRITE_DBG({"modified KSIG/KEIG: ", JSIG, JEIG})
    M4_WRITE_DBG({"modified KSIG/KEIG: ", KSIG, KEIG})

    if((IEIG.le.ISIG) .or. (JEIG.le.JSIG) .or. (KEIG.le.KSIG)) then
       write(STDERR,*) "OPPOSITE PML LAYERS OVERLAP!"
       stop
    endif
    

    ! SigmaMax (= SigmaOpt, see Tavlove 2, pp 286)
    ! It is Sig = Sig[SI] / (!c*eps0)

    ! PotPML = 3.2
    ! SigmaMax = (real(POTPml)+1.0)*0.8/(3.0*DT)
    ! KappaMax = 1.1
    

    call AllocateFields
    call CalcCoefficients(IBEG, IMAX, ISIG, IEIG, cexpml, cmxpml)
    call CalcCoefficients(JBEG, JMAX, JSIG, JEIG, ceypml, cmypml)
    call CalcCoefficients(KBEG, KMAX, KSIG, KEIG, cezpml, cmzpml)

    modinitialized = .true.

    M4_WRITE_DBG({". exit InitializePml"})

  contains
    
    subroutine AllocateFields

      integer :: err, i, j, k, l

      M4_WRITE_DBG({". enter InitializePml.AllocateFields"})
      
      ! numeric coefficient-fields
      
      allocate(cexpml(1:4,IBEG:IMAX), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(ceypml(1:4,JBEG:JMAX), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(cezpml(1:4,KBEG:KMAX), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(cmxpml(1:4,IBEG:IMAX), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(cmypml(1:4,JBEG:JMAX), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(cmzpml(1:4,KBEG:KMAX), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      !  B and D fields for each of the 6 layers 
      
      allocate(BE1(1:3,IBEG:ISIG-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE1(1:3,IBEG:ISIG-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE2(1:3,IEIG:IMAX-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE2(1:3,IEIG:IMAX-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE3(1:3,ISIG:IEIG-1,JBEG:JSIG-1,KBEG:KMAX-1), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(DE3(1:3,ISIG:IEIG-1,JBEG:JSIG-1,KBEG:KMAX-1), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE4(1:3,ISIG:IEIG-1,JEIG:JMAX-1,KBEG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE4(1:3,ISIG:IEIG-1,JEIG:JMAX-1,KBEG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE5(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KBEG:KSIG-1), STAT=err) 
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(DE5(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KBEG:KSIG-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      allocate(BE6(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KEIG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")

      allocate(DE6(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KEIG:KMAX-1), STAT=err)
      M4_ALLOC_ERROR(err,"AllocFields")
      
      ! initialize fields
      
      cexpml=1.0; ceypml=1.0; cezpml=1.0
      cmxpml=1.0; cmypml=1.0; cmzpml=1.0

      do l = 1, 3
         do j = JBEG, JMAX-1
            do k = KBEG, KMAX-1
               do i = IBEG, ISIG-1
                  BE1(l,i,j,k) = 0.0; DE1(l,i,j,k) = 0.0;
               enddo
               
               do i = IEIG, IMAX-1
                  BE2(l,i,j,k) = 0.0; DE2(l,i,j,k) = 0.0;
               enddo
            enddo
         enddo
         
         do i = ISIG, IEIG-1
            do k = KBEG, KMAX-1
               do j = JBEG, JSIG-1
                  BE3(l,i,j,k) = 0.0; DE3(l,i,j,k) = 0.0;
               enddo
                 
               do j = JEIG, JMAX-1
                  BE4(l,i,j,k) = 0.0; DE4(l,i,j,k) = 0.0; 
               enddo
            enddo
            
            do j = JSIG, JEIG-1
               do k = KBEG, KSIG-1
                  BE5(l,i,j,k) = 0.0; DE5(l,i,j,k) = 0.0;
               enddo
               
               do k = KEIG, KMAX-1
                  BE6(l,i,j,k) = 0.0; DE6(l,i,j,k) = 0.0; 
               enddo
            enddo
         enddo
      enddo

      M4_WRITE_DBG({". exit InitializePml.AllocateFields"})
      
    end subroutine AllocateFields
    
    subroutine CalcCoefficients(lbeg, lmax, ls, le, ce, cm)
      
      integer, intent(in) :: lbeg, lmax, ls, le
      real(kind=8), dimension(1:4,lbeg:lmax), intent(inout) :: ce, cm
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
         sigma = SigmaMax*x**POTPml
         kappa = 1.0+(KappaMax-1.0)*x**POTPml
         val1(l) = kappa+0.5*sigma*DT
         val2(l) = kappa-0.5*sigma*DT
         
         ! points in the grid center
         
         x = (real(l)+0.5) / real(PMLMAX-1)
         sigma = SigmaMax*x**POTPml
         kappa = 1.0+(KappaMax-1.0)*x**POTPml
         val1p(l) = kappa+0.5*sigma*DT
         val2p(l) = kappa-0.5*sigma*DT
  
      enddo
      
      ! calculate coefficients
      
      do l=lbeg, ls-1
         ce(1,l) = val1p(ls-1-l-lbeg-1)
         ce(3,l) = val2p(ls-1-l-lbeg-1)
         ce(2,l) = 1.0 / val1(ls-1-l-lbeg)
         ce(4,l) = val2(ls-1-l-lbeg) / val1(ls-1-l-lbeg)
         
         cm(1,l) = val1(ls-1-l-lbeg)
         cm(3,l) = val2(ls-1-l-lbeg)
         cm(2,l) = 1.0 / val1p(ls-1-l-lbeg-1)
         cm(4,l) = val2p(ls-1-l-lbeg-1) / val1p(ls-1-l-lbeg-1)
      enddo

      do l=le, lmax-1
         ce(1,l) = val1p(l-le)
         ce(3,l) = val2p(l-le)
         ce(2,l) = 1.0 / val1(l-le)
         ce(4,l) = val2(l-le) / val1(l-le)
         
         cm(1,l) = val1(l-le)
         cm(3,l) = val2(l-le)
         cm(2,l) = 1.0 / val1p(l-le)
         cm(4,l) = val2p(l-le) / val1p(l-le)
      enddo

      M4_WRITE_DBG({". exit InitializePml.CalcCoefficients"})
      
    end subroutine CalcCoefficients
      
  end subroutine InitializePml

!----------------------------------------------------------------------

  subroutine FinalizePml

    M4_WRITE_DBG({". FinalizePml"})

    deallocate(DE6)
    deallocate(BE6)
    deallocate(DE5)
    deallocate(BE5)
    deallocate(DE4)
    deallocate(BE4)
    deallocate(DE3)
    deallocate(BE3)
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
  
  subroutine StepHPml

    if ( .not. modinitialized ) return

    call DoStepHPml(IBEG,ISIG-1,JBEG,JMAX-1,KBEG,KMAX-1,BE1) 
    call DoStepHPml(IEIG,IMAX-1,JBEG,JMAX-1,KBEG,KMAX-1,BE2)
    call DoStepHPml(ISIG,IEIG-1,JBEG,JSIG-1,KBEG,KMAX-1,BE3)
    call DoStepHPml(ISIG,IEIG-1,JEIG,JMAX-1,KBEG,KMAX-1,BE4)
    call DoStepHPml(ISIG,IEIG-1,JSIG,JEIG-1,KBEG,KSIG-1,BE5)
    call DoStepHPml(ISIG,IEIG-1,JSIG,JEIG-1,KEIG,KMAX-1,BE6)

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
        
      do k=ks, ke     
         !SOPTION unroll(1)  
         do j=js, je
            !SOPTION unroll(1)  
            do i=is, ie
               
               Exh = Ex(i,j,k)
               Eyh = Ey(i,j,k)
               Ezh = Ez(i,j,k)
               
               Bxo = B(1,i,j,k)
               Byo = B(2,i,j,k) 
               Bzo = B(3,i,j,k)
               
               ! Update B
               
               B(1,i,j,k) = cmypml(4,j)*B(1,i,j,k) + cmypml(2,j) * &
                    (- dty*( Ez(i,j+1,k) - Ezh )           &
                    +  dtz*( Ey(i,j,k+1) - Eyh ))
               B(2,i,j,k) = cmzpml(4,k)*B(2,i,j,k) + cmzpml(2,k) * &
                    (- dtz*( Ex(i,j,k+1) - Exh )           &
                    +  dtx*( Ez(i+1,j,k) - Ezh ))
               B(3,i,j,k) = cmxpml(4,i)*B(3,i,j,k) + cmxpml(2,i) * &
                    (- dtx*( Ey(i+1,j,k) - Eyh )           &
                    +  dty*( Ex(i,j+1,k) - Exh ))        
               
               ! Calc H
               
               Hx(i,j,k) = cmzpml(4,k)*Hx(i,j,k) + cmzpml(2,k) * &
                    (cmxpml(1,i)*B(1,i,j,k) - cmxpml(3,i)*Bxo)
               Hy(i,j,k) = cmxpml(4,i)*Hy(i,j,k) + cmxpml(2,i) * &
                    (cmypml(1,j)*B(2,i,j,k) - cmypml(3,j)*Byo)
               Hz(i,j,k) = cmypml(4,j)*Hz(i,j,k) + cmypml(2,j) * &
                    (cmzpml(1,k)*B(3,i,j,k) - cmzpml(3,k)*Bzo)           
               
            enddo
         enddo
      enddo
      
    end subroutine DoStepHPml
    
  end subroutine StepHPml

!----------------------------------------------------------------------

  ! update the E-fields of all pml layers
  
  subroutine StepEPml

    if ( .not. modinitialized ) return

    call DoStepEPml(IBEG,ISIG-1,JBEG,JMAX-1,KBEG,KMAX-1,DE1)
    call DoStepEPml(IEIG,IMAX-1,JBEG,JMAX-1,KBEG,KMAX-1,DE2)
    call DoStepEPml(ISIG,IEIG-1,JBEG,JSIG-1,KBEG,KMAX-1,DE3)
    call DoStepEPml(ISIG,IEIG-1,JEIG,JMAX-1,KBEG,KMAX-1,DE4)
    call DoStepEPml(ISIG,IEIG-1,JSIG,JEIG-1,KBEG,KSIG-1,DE5)
    call DoStepEPml(ISIG,IEIG-1,JSIG,JEIG-1,KEIG,KMAX-1,DE6)
    
  contains
    
    subroutine DoStepEPml(is,ie,js,je,ks,ke,D)
      
      integer :: is, ie, js, je, ks, ke
      M4_FTYPE, dimension(1:3,is:ie,js:je,ks:ke) :: D
      
      integer :: i, j, k
      real(kind=8) :: dtx, dty, dtz
      M4_FTYPE :: Dxo, Dyo, Dzo, Hxh, Hyh, Hzh    
      real(kind=8) :: epsinvx, epsinvy, epsinvz
  
      dtx = DT/SX
      dty = DT/SY
      dtz = DT/SZ
  
      do k=ks, ke
!SOPTION unroll(2)  
         do j=js, je
!SOPTION unroll(2)  
            do i=is, ie
           
               Hxh = Hx(i,j,k)
               Hyh = Hy(i,j,k)
               Hzh = Hz(i,j,k)
           
               Dxo = D(1,i,j,k)
               Dyo = D(2,i,j,k)
               Dzo = D(3,i,j,k)
           
               ! Update D

               D(1,i,j,k) =  ceypml(4,j)*D(1,i,j,k) + ceypml(2,j) * &
                    ( dty * ( Hzh - Hz(i,j-1,k) )           &
                    - dtz * ( Hyh - Hy(i,j,k-1) ))
               D(2,i,j,k) =  cezpml(4,k)*D(2,i,j,k) + cezpml(2,k) * &
                    ( dtz * ( Hxh - Hx(i,j,k-1) )           &
                    - dtx * ( Hzh - Hz(i-1,j,k) ))
               D(3,i,j,k) =  cexpml(4,i)*D(3,i,j,k) + cexpml(2,i) * &
                    ( dtx * ( Hyh - Hy(i-1,j,k) )           &
                    - dty * ( Hxh - Hx(i,j-1,k) ))

               ! Calc E

               epsinvx = 0.5*(EPSINV(i,j,k) +  EPSINV(i+1,j,k))
               epsinvy = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j+1,k))
               epsinvz = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j,k+1))

               Ex(i,j,k) = cezpml(4,k)*Ex(i,j,k) + cezpml(2,k) * &
                    epsinvx*(cexpml(1,i)*D(1,i,j,k) - cexpml(3,i)*Dxo)
               Ey(i,j,k) = cexpml(4,i)*Ey(i,j,k) + cexpml(2,i) * &
                    epsinvy*(ceypml(1,j)*D(2,i,j,k) - ceypml(3,j)*Dyo)
               Ez(i,j,k) = ceypml(4,j)*Ez(i,j,k) + ceypml(2,j) * &
                    epsinvz*(cezpml(1,k)*D(3,i,j,k) - cezpml(3,k)*Dzo)
           
            enddo
         enddo
      enddo

    end subroutine DoStepEPml
    
  end subroutine StepEPml
  
!----------------------------------------------------------------------

end module pml

!
! Authors:  A.Klaedtke, S.Scholz, J.Hamm
! Modified: 4/12/2007
!
!======================================================================
