!----------------------------------------------------------------------
!
!  program: pre-obj-1 
!
!  preprocessing program reading <obj-epsilon.in>, <grid.in> files and 
!  writing epsilon 
!
!----------------------------------------------------------------------
 
program preobj1

  use grid
  use fdtd

  implicit none
  save

  integer, parameter :: subgrid=5  ! size of subgrid: (2*subgrid+1)*(2*subgrid+1)

  integer, parameter :: unit=10
  integer :: ios, p


  write(6,*) "ReadGrid"
  call InitializeGrid(sfxin)
!  call EchoGrid
  write(6,*) "InitEpsilon"
  call InitEpsilon
  write(6,*) "SetEpsilon"
  call SetEpsilon
  write(6,*) "WriteEpsilon"
  call WriteEpsilon(-1)
  do p=0, PARTITIONS-1
     write(6,*) "WriteEpsilon(",p,")"
     call WriteEpsilon(p)
  end do


contains

  subroutine InitEpsilon

    implicit none

    integer:: err

    allocate(EPSINV(0:IMAX,0:JMAX,0:KMAX),stat=err)
    
    if ( err .ne. 0 ) then
       write(6,*) "InitEpsilon failed to allocate memory"
       stop
    end if

    EPSINV=1.0

  end subroutine InitEpsilon


  subroutine SetEpsilon

    implicit none
    save
    
    character(len=255),parameter :: file='obj-epsilon.in'

    integer :: i, j, k, n, ii, jj, kk, cntr, ios
    integer :: layerdim, ils, ile, jls, jle, kle, kls
    integer :: zylpuc, zyl
    integer :: ucdimi, ucdimj, ucdimk
    integer :: mi, mj, mk, mmax
    integer :: flen
    real(8) :: eps1, epsb, epsval, val, val1, val2
    real(8) :: epscyl, ni, nj, nk, norm 
    real(8) :: cylh, cylmi, cylmj, cylmk, px, py, pz, h, dsqr, ax, ay, az, radq
    character(len=STRLNG) :: set, str
    real(8) :: radius
    character(len=80) :: datei
    character(len=80) :: fnh

    integer :: sphpuc, sph, sphposi, sphposj, sphposk
    real(8),allocatable:: sphradius(:)
    integer, allocatable:: sphi(:)
    integer, allocatable:: sphj(:)
    integer, allocatable:: sphk(:)
    integer :: sphucdimi, sphucdimj, sphucdimk
    real(8) :: epssph


    open(unit,FILE=file, STATUS='unknown')

    !I: Background Material

    write(6,*) "SetEpsilon: backround material"

    do
       read(unit, IOSTAT=ios, FMT=*) str
       if(ios .ne. 0) exit
       if ( str(1:4) .eq. '#INF') then
          read(unit,*) epsb
          val=1.0/epsb
          do k=0,KMAX
             do j=0,JMAX
                do i=0,IMAX
                   EPSINV(i,j,k) = val
                enddo
             enddo
          enddo
          read(unit,*) layerdim
          read(unit,*) zylpuc
          read(unit,*) sphpuc

       endif
    enddo


    !II. Layer in z-Richtung
   
    write(6,*) "SetEpsilon: layers, ", layerdim

    if (layerdim > 0) then
       rewind unit
       do
          read(unit, IOSTAT=ios, FMT=*) str
          if(ios .ne. 0) exit
          if (str(1:4) .eq. '#LAY') then
             do n=1, layerdim
                read(unit,*) ils, ile, jls, jle, kls, kle, eps1
                ils = max(ils,IBEG)
                ile = min(ile,IEND+1)
                jls = max(jls,JBEG)
                jle = min(jle,JEND+1)
                kls = max(kls,KBEG)
                kle = min(kle,KEND+1)
                val=1.0/eps1
                do k=kls, kle
                   do j= jls, jle
                      do i= ils, ile
                         EPSINV(i,j,k) = val
                      enddo
                   enddo
                enddo
             enddo
          endif
       enddo
    endif

    


    !IV. Dielectric Spheres

    write(6,*) "SetEpsilon: spheres"

    if (sphpuc > 0) then 
       allocate(sphradius(1:sphpuc))
       allocate(sphi(1:sphpuc))
       allocate(sphj(1:sphpuc))
       allocate(sphk(1:sphpuc))

       do sph=1, sphpuc

          write(6,*) "Calculating Sphere #", sph

          write(fnh,*) sph
          fnh=adjustl(fnh)
          datei='#SPH'
          datei((len_trim(datei)+1):20)=fnh      
          
          rewind unit       
          do
             read(unit, IOSTAT=ios, FMT=*) str

             if(ios .ne. 0) exit 

             flen=int(log10(sph*1.0))+6

             if (str(1:flen) .eq. datei) then

                ! READ GEOM DATA
                read(unit,*) sphucdimi, sphucdimj, sphucdimk
                read(unit,*) sphradius(sph)
                read(unit,*) epssph
                read(unit,*) sphi(sph), sphj(sph), sphk(sph)

                ! READ CRYSTAL PATTERNS AND SET EPS
                
                do k=0, KMAX, sphucdimk
                   do j=0, JMAX, sphucdimj
                      do i=0, IMAX, sphucdimi

                         do kk=0,  sphucdimk-1
                            if (kk+k .gt. KMAX) exit

                            do jj=0,  sphucdimj-1
                               if (jj+j .gt. JMAX) exit

                               do ii=0,  sphucdimi-1
                                  if (ii+i .gt. IMAX) exit

                                  epsval=0.0
                                  mmax=subgrid
                                  do mi=-mmax, mmax
                                     do mj=-mmax, mmax
                                        do mk=-mmax, mmax
                                        
                                           ax=(sphi(sph) - ii +0.5*mi/mmax)**2
                                           ay=(sphj(sph) - jj +0.5*mj/mmax)**2
                                           az=(sphk(sph) - kk +0.5*mk/mmax)**2
                                           
                                           radq=ax+ay+az
                                           
                                           if (radq .le. sphradius(sph)**2) then
                                              epsval=epsval+epssph/((2*mmax+1)**3)
                                           else
                                              epsval=epsval+1.0/EPSINV(i+ii,j+jj,k+kk)/((2*mmax+1)**3)
                                           endif
                                        
                                        enddo
                                     enddo
                                  enddo

                                  if ((i+ii .le. IMAX) .and. (j+jj .le. JMAX) & 
                                       .and. (k+kk .le. KMAX)) then
                                     EPSINV(i+ii,j+jj,k+kk)=1.0/epsval
                                  endif

                               enddo
                            enddo
                         enddo

                      enddo
                   enddo
                enddo


             endif
          enddo

       enddo
    endif

    !III. Dielectric Cylinders

    write(6,*) "SetEpsilon: cylinders"

    if (zylpuc > 0) then 

       do zyl=1, zylpuc

          write(6,*) "Calculating Cylinder #", zyl

          write(fnh,*) zyl
          fnh=adjustl(fnh)
          datei='#CYL'
          datei((len_trim(datei)+1):20)=fnh      
          
          rewind unit       
          do
             read(unit, IOSTAT=ios, FMT=*) str

             if(ios .ne. 0) exit 

             flen=int(log10(zyl*1.0))+6

             if (str(1:flen) .eq. datei) then
                
                ! READ GEOM DATA
                read(unit,*) ucdimi, ucdimj, ucdimk
                read(unit,*) radius, cylh
                read(unit,*) epscyl
                read(unit,*) cylmi, cylmj, cylmk
                read(unit,*) ni, nj, nk 

                norm=sqrt(ni**2+nj**2+nk**2)
                ni=ni/norm
                nj=nj/norm
                nk=nk/norm
                
                ! READ CRYSTAL PATTERNS AND SET EPS

                do k=0, KMAX, ucdimk

                   write(6,*)"k:",k

                   do j=0, JMAX, ucdimj

                      read(unit,*) set 
                      write(6,*) set
                      cntr=0

                      do i=0, IMAX, ucdimi

                         cntr=cntr+1  

                         do kk=0,  ucdimk-1
                            if (kk+k .gt. KMAX) exit

                            do jj=0,  ucdimj-1
                               if (jj+j .gt. JMAX) exit

                               do ii=0,  ucdimi-1
                                  if (ii+i .gt. IMAX) exit

                                  if (set(cntr:cntr) .eq. '1') then

                                     epsval=0.0
                                     mmax=subgrid
                                     do mi=-mmax, mmax
                                        do mj=-mmax, mmax
                                           do mk=-mmax, mmax
                                              
                                              px = ii + 0.5*mi/mmax
                                              py = jj + 0.5*mj/mmax
                                              pz = kk + 0.5*mk/mmax
                                           
                                              h=(px-cylmi)*ni+(py-cylmj)*nj+(pz-cylmk)*nk
                                              dsqr=(cylmi+h*ni-px)**2+(cylmj+h*nj-py)**2+(cylmk+h*nk-pz)**2

                                              if (((h .ge. 0) .and. (h .le. cylh)) &
                                                   .and. (dsqr .le. radius**2)) then
                                                 epsval=epsval+epscyl/((2*mmax+1)**3)
 !                                              write(6,*) h, dsqr
                                              else
                                                 epsval=epsval+1.0/EPSINV(i+ii,j+jj,k+kk)/((2*mmax+1)**3)
                                              endif
                                        
                                           enddo
                                        enddo
                                     enddo

                                  
                                     if ((i+ii .le. IMAX) .and. (j+jj .le. JMAX) & 
                                          .and. (k+kk .le. KMAX)) then
                                        EPSINV(i+ii,j+jj,k+kk)=1.0/epsval
                                     endif

                                  endif

                               enddo
                            enddo
                         enddo

                      enddo
                   enddo
                enddo
         
             endif
          end do

       end do
    end if

  end subroutine SetEpsilon

  subroutine WriteEpsilon(p)

    implicit none
    integer :: p
    integer i, j, k

    character(len=255) :: file
    character(len=255) :: sfx

    if ( p .ne. -1 ) then
       sfx = cat3('.',i2str(p),'.in')
    else
       sfx = '.in'
    end if

    file = cat2(pfxepsilon,sfx)

    call InitializeGrid(sfx)

    open(unit,FILE=file,STATUS='unknown')
    
    write(6,*) "... range: ",KBEG,KEND+1

    do k=KBEG, KEND+1
       do j=JBEG, JEND+1
          do i=IBEG,IEND+1

             write(unit,*) 1.0/EPSINV(i,j,k)

          end do
       end do
    end do
    
    close(unit)

    call InitializeGrid(sfxin)

  end subroutine WriteEpsilon

end program preobj1
