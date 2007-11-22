!----------------------------------------------------------------------
!
!  module: upml(-r) / max3d
!
!  boundary conditions using uniform perfectly matched layers.
!
!----------------------------------------------------------------------

module upml

  use constant
  use grid  
  implicit none
  save

  integer :: pmlpart               ! flag (inner/outer partition)
  integer :: PMLMAX                ! Number of PML layers
  real(8) :: PotPml                   ! Exponent for sigma and kappa
  real(8) :: SigmaMax                 ! Absorption coefficient 
  real(8) :: KappaMax                 ! Absorption of evanescent waves 

  ! Planes 1-6 with a PML (Yes=1, No otherwise)

  integer, dimension(6) :: planepml=(/ 1,1,1,1,1,1 /)

  ! Numeric PML Coefficients

  real(8),allocatable,dimension(:,:) :: cexpml,ceypml,cezpml
  real(8),allocatable,dimension(:,:) :: cmxpml,cmypml,cmzpml

  ! Auxilliary B and D fields on the 6 planes E1 .. E6

  real(8), allocatable, dimension(:,:,:,:) :: BE1,BE2,BE3,BE4,BE5,BE6 
  real(8), allocatable, dimension(:,:,:,:) :: DE1,DE2,DE3,DE4,DE5,DE6 

contains

  subroutine InitPML ()

    implicit none

    ! SigmaMax (= SigmaOpt, see Tavlove 2, pp 286)
    ! It is Sig = Sig[SI] / (!c*eps0)

    ! SigmaMax = (real(POTPML)+1.0)*0.8/(3.0*DT)
    ! KappaMax = 1.1


    call InitMem()
    call CalcKoeff(IBEG, IMAX, ISIG, IEIG, cexpml, cmxpml)
    call CalcKoeff(JBEG, JMAX, JSIG, JEIG, ceypml, cmypml)
    call CalcKoeff(KBEG, KMAX, KSIG, KEIG, cezpml, cmzpml)



    contains

      subroutine InitMem()
        ! Allocates fields of module PML
        implicit none

        !Variables
        integer :: err, i, j, k, l

        ! Numerische Koeffizienten-Felder 
        allocate(cexpml(1:4,IBEG:IMAX), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 1 in AllocFields/PML'
           stop
        endif
        allocate(ceypml(1:4,JBEG:JMAX), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 1 in AllocFields/PML'
           stop
        endif
        allocate(cezpml(1:4,KBEG:KMAX), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 1 in AllocFields/PML'
           stop
        endif
        allocate(cmxpml(1:4,IBEG:IMAX), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 1 in AllocFields/PML'
           stop
        endif
        allocate(cmypml(1:4,JBEG:JMAX), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 1 in AllocFields/PML'
           stop
        endif
        allocate(cmzpml(1:4,KBEG:KMAX), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 1 in AllocFields/PML'
           stop
        endif
        
        ! Felder B und D in den 6 Randschichten 
        allocate(BE1(1:3,IBEG:ISIG-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err) 
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(DE1(1:3,IBEG:ISIG-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(BE2(1:3,IEIG:IMAX-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(DE2(1:3,IEIG:IMAX-1,JBEG:JMAX-1,KBEG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(BE3(1:3,ISIG:IEIG-1,JBEG:JSIG-1,KBEG:KMAX-1), STAT=err) 
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(DE3(1:3,ISIG:IEIG-1,JBEG:JSIG-1,KBEG:KMAX-1), STAT=err) 
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(BE4(1:3,ISIG:IEIG-1,JEIG:JMAX-1,KBEG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(DE4(1:3,ISIG:IEIG-1,JEIG:JMAX-1,KBEG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(BE5(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KBEG:KSIG-1), STAT=err) 
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(DE5(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KBEG:KSIG-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(BE6(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KEIG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        allocate(DE6(1:3,ISIG:IEIG-1,JSIG:JEIG-1,KEIG:KMAX-1), STAT=err)
        if(err .ne. 0) then
           write(6,*) 'Fehler 2 in AllocFields/PML'
           stop
        endif
        
        ! Initialize fields
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

      end subroutine InitMem


      subroutine CalcKoeff(lbeg, lmax, ls, le, ce, cm)

        implicit none

        integer, intent(in) :: lbeg, lmax, ls, le
        real(8), dimension(1:4,lbeg:lmax), intent(inout) :: ce, cm
        real(8), dimension(-1:PMLMAX-1) :: val1,val1p,val2,val2p
        real(8) :: x, sigma, kappa
        integer :: l
        
        ! Hilfsgroessen

        val1 = 1.0
        val2 = 1.0
        val1p = 1.0
        val2p = 1.0

        do l=0, PMLMAX-1
           ! Punkte auf den Gitterecken
           x = real(l) / real(PMLMAX-1)
           sigma = SigmaMax*x**POTPML
           kappa = 1.0+(KappaMax-1.0)*x**POTPML
           val1(l) = kappa+0.5*sigma*DT
           val2(l) = kappa-0.5*sigma*DT
           
           ! Punkte in der Gittermitte
           x = (real(l)+0.5) / real(PMLMAX-1)
           sigma = SigmaMax*x**POTPML
           kappa = 1.0+(KappaMax-1.0)*x**POTPML
           val1p(l) = kappa+0.5*sigma*DT
           val2p(l) = kappa-0.5*sigma*DT
        enddo
        
        ! Koeffizienten berechnen
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

      end subroutine CalcKoeff

  end subroutine InitPML


  subroutine FinalisePML ()

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

  end subroutine FinalisePML


  ! Update der H-Felder in allen PML Schichten
  
  subroutine PMLH

    implicit none

    call StepPmlH(IBEG,ISIG-1,JBEG,JMAX-1,KBEG,KMAX-1,BE1) 
    call StepPmlH(IEIG,IMAX-1,JBEG,JMAX-1,KBEG,KMAX-1,BE2)
    call StepPmlH(ISIG,IEIG-1,JBEG,JSIG-1,KBEG,KMAX-1,BE3)
    call StepPmlH(ISIG,IEIG-1,JEIG,JMAX-1,KBEG,KMAX-1,BE4)
    call StepPmlH(ISIG,IEIG-1,JSIG,JEIG-1,KBEG,KSIG-1,BE5)
    call StepPmlH(ISIG,IEIG-1,JSIG,JEIG-1,KEIG,KMAX-1,BE6)

  contains
    
    subroutine StepPmlH(is,ie,js,je,ks,ke,B)

      implicit none

      integer is, ie, js, je, ks, ke
      real(8), dimension(1:3,is:ie,js:je,ks:ke) :: B
      integer i, j, k
      real(8) Bxo, Byo, Bzo, Exh, Eyh, Ezh
      real(8) dtx, dty, dtz

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

    end subroutine StepPmlH

  end subroutine PMLH


  ! Update der H-Felder in allen PML Schichten


  subroutine PmlE

    implicit none

    call StepPmlE(IBEG,ISIG-1,JBEG,JMAX-1,KBEG,KMAX-1,DE1)
    call StepPmlE(IEIG,IMAX-1,JBEG,JMAX-1,KBEG,KMAX-1,DE2)
    call StepPmlE(ISIG,IEIG-1,JBEG,JSIG-1,KBEG,KMAX-1,DE3)
    call StepPmlE(ISIG,IEIG-1,JEIG,JMAX-1,KBEG,KMAX-1,DE4)
    call StepPmlE(ISIG,IEIG-1,JSIG,JEIG-1,KBEG,KSIG-1,DE5)
    call StepPmlE(ISIG,IEIG-1,JSIG,JEIG-1,KEIG,KMAX-1,DE6)
    
  contains

    subroutine StepPmlE(is,ie,js,je,ks,ke,D)

      implicit none

      integer is, ie, js, je, ks, ke
      real(8), dimension(1:3,is:ie,js:je,ks:ke) :: D

      integer :: i, j, k
      real(8) :: dtx, dty, dtz
      real(8) :: Dxo, Dyo, Dzo, Hxh, Hyh, Hzh    
      real(8) :: epsinvx, epsinvy, epsinvz
  
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

    end subroutine StepPmlE

  end subroutine PMLE

end module upml
