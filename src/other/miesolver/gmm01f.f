C  Fortran code gmm01f.f for calculating radiative scattering by an 
C  external aggregate of homogeneous spheres in a fixed orientation or  
C  at an average over individual orientations 
C  10/20/1994    Yu-lin Xu
C  revised at 8/1996, 6/1998, 8/2000, 12/2000
C  released to public at 1/2001
C  ---------------------------------------------------------------------
C
C  The scattering formulation and numerical techniques used in this 
C  code can be found in the following references:
C  (1) Cruzan, Q. Appl. Math. 20, 33 (1962)
C  (2) Bruning and Lo, IEEE Trans. Anten. Prop. AP-19, 378 (1971)
C      Fuller and Kattawar, Opt. Lett. 13, 90 & 1063 (1988)
C      Mackowski, Proc. R. Soc. Lond. A 433, 599 (1991)
C      Wang and van der Hulst, Appl. Opt. 30, 106 (1991)
C      H.A. van der Vorst, SIAM, J. Sci. Stat. Comput. 13, 631 (1992)
C  (3) Xu, Appl. Opt. 34, 4573 (1995), errata, ibid. 37, 6494 (1998)  
C          Appl. Opt. 36, 9496 (1997) 
C          Phys. Lett. A 249, 30 (1998)
C          J. Comput. Phys. 139, 137 (1998)
C      Xu and Wang, Phys. Rev. E 58, 3931 (1998)
C      Xu, Gustafson, Giovane, Blum, and Tehranian,
C          Phys. Rev. E 60, 2347 (1999)
C
C  The normalization factor in field-expansions in this work is slightly
C  different from the one used in references (3), which is now 
C             E_0 i^n [(2n+1)(n-m)!/n/(n+1)/(n+m)!]^{1/2},
C  instead of the previous 
C             E_0 i^n (2n+1)(n-m)!/(n+m)!.
C  The constant factors in all scattering formulas used in this code are 
C  thus different from those found in the references.   
C  ---------------------------------------------------------------------
C
C  For questions/comments/suggestions/bugs/problems please contact 
C  Yu-lin Xu at shu@astro.ufl.edu
C
C  ---------------------------------------------------------------------  
	PROGRAM gmm01f
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
C	parameter (nLp=100,np=20)
C  ---------------------------------------------------------------------
C  Two parameters in "gmm01f.par": nLp, np
C  nLp --- the maximum number of spheres, must be equal to or greater 
C          than the sphere-number in an aggregate actually calculated  
C  np  --- the maximum scattering order in the incident and scattered
C          field expansions, must be equal to or greater than the 
C          highest scattering order required in actual calculations
C  an example of "gmm01f.par":   parameter (nLp=100,np=20)
C  
C  ****** The computer memory required by this code is at the level of 
C                      nLp**2*np**3 + np**4.
C  This code is written in double precision arithmetic. An individual  
C  sphere can have a size parameter way beyond ~200. There is no limit  
C  to the overall dimension of an aggregate. The largest individual  
C  sphere size and the maximum number of spheres that an aggregate can  
C  have depends on the availability of computer memory.  ****** 
C  ---------------------------------------------------------------------
	parameter (nmp=np*(np+2),nmp0=(np+1)*(np+4)/2)
	parameter (NXMAX=3000,nangmax=181,MOR=181,ncmax=180)
C  --------------------------------------------------------------------- 
C  NXMAX - the maximum dimension of an auxiliary array in calculating 
C          Mie scattering coefficients, must not be less than 
C                  1.1*np*|m|+10 
C          where |m| is the largest refractive index of all spheres) 
C  ---------------------------------------------------------------------
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)	
	parameter (nrc=4*np*(np+1)*(np+2)/3+np)
	parameter (nij=nLp*(nLp-1)/2)
	integer u,v,u0,nmax(nLp),uvmax(nLp),ind(nLp)
	double precision k,lnfacd,r0(6,nLp),x(nLp),dang(nangmax),
     +     r00(3,nLp),rsr0(NXMAX),rsi0(NXMAX),
     +     rsx0(NXMAX),px0(NXMAX),w1(np),w2(np),w3(np),w4(np),
     +     rsr(np,nLp),rsi(np,nLp),rsx(np,nLp),px(np,nLp),betar(MOR),
     +     thetr(MOR),phair(MOR),smue(4,4),mue(4,4,ncmax,nangmax),
     +     besj(0:2*np+1),besy(0:2*np+1),i11(nangmax),
     +     i21(nangmax),i22(nangmax),i12(nangmax),inat(nangmax),
     +     pol(nangmax),cscaxi(nLp),cscayi(nLp),cextxi(nLp),
     +     cextyi(nLp),cabsxi(nLp),cabsyi(nLp),cexti(nLp),cabsi(nLp),
     +     cscai(nLp),assymi(nLp),assymxi(nLp),assymyi(nLp),
     +     cprxi(nLp),cpryi(nLp),cpri(nLp),drot(nrc,nij),
     +     c0i(nLp),c1i(nLp)
	complex*16 A,B,cmz,Aj,Bj,A2,B2,Aj2,Bj2,A0,B0,ephi,ci,cin,
     +     atr0(ni0,nij),btr0(ni0,nij),atr(2,np,nmp),at(nmp),bt(nmp),
     +     atr1(ni0,nij),btr1(ni0,nij),ek(np,nij),ref(nLp),p0(nLp,nmp),
     +     q0(nLp,nmp),an(np),bn(np),aMie(nLp,np),bMie(nLp,np),
     +	   as(nLp,nmp),bs(nLp,nmp),as0(nLp,nmp),bs0(nLp,nmp),
     +     asc(nLp,nmp),bsc(nLp,nmp),as1(nLp,nmp),bs1(nLp,nmp),
     +     ast(nLp,nmp),bst(nLp,nmp),asp(nLp,nmp),bsp(nLp,nmp),
     +     asv(nLp,nmp),bsv(nLp,nmp),
     +     s2x(ncmax,nangmax),s4x(ncmax,nangmax),
     +	   s3y(ncmax,nangmax),s1y(ncmax,nangmax),
     +	   atj(nmp),btj(nmp),py0(NXMAX),py(NXMAX),dpy(NXMAX)
	CHARACTER FLNAME*20,fileout*20,fileout1*19,fileout2*21,
     +     tailn*3,cnr2*2,cnr3*3,fileoutA*20,cnr1*1,flout*22
	COMMON/MIESUB/ twopi,pih
        common/rot/bcof(0:np+2),dc(-np:np,0:nmp)
        common/fnr/fnr(0:2*(np+2))
	common/pitau/pi(nmp0),tau(nmp0)
        common/tran/atr
	common/ig0/iga0(ni0)
	common/g0/ga0(ng0)
	common/cofmnv0/cof0(ni0)
	common/crot/cofsr(nmp)
	pih   = dacos(0.d0)
	twopi  = 4.d0*pih
	pione  = 2.d0*pih
	ci=dcmplx(0.d0,1.d0)
	cin=dcmplx(0.d0,-1.d0)
	gcs=0.d0
	gcv=0.d0
	idpq=0	
	OPEN(UNIT=1,FILE='gmm01f.in',STATUS='OLD')
	READ(1,'(a20)') FLNAME
C  ---------------------------------------------------------------------
C  FLNAME - the input file name for the aggregate configuration 
C           the first line in this file is the incident wavelength, 
C           the second line is the number of spheres, 
C           rest lines provide coordinates of sphere-centers, radii, and 
C           refractive indexes of the spheres; each line contains six 
C           real numbers of x, y, z, r, Re(m), Im(m) for a sphere
C  ---------------------------------------------------------------------
	write(6,'(a12,a20)') 'input file: ',FLNAME
	READ(1,*) nbeta,nthet,nphai
C  ---------------------------------------------------------------------
C  The product of nbeta, nthet, and nphai is the total number of 
C  orientations to be averaged for the input sphere-aggregate.
C  In the Cartesian coordinate system, the direction of propagation of 
C  the incident plane wave defines the positive z-axis, the x-z plane is 
C  the scattering plane. 
C  nbeta - # of rotation around the z-axis,  
C          this rotation angle corresponds to the Euler angle of "alpha" 
C  nthet - # of rotation around the y-axis, 
C          corresponding to the Euler angle of "beta"
C  nphai - # of rotation by the Euler angle of "gamma"
C  The definitions of the three Euler angles follow the convention used 
C  by Edmonds ["Angular Momentum in Quantum Mechanics," pp.6-8 (1957)].
C  examples:
C  When only a single orientation needs to be calculated:
C    nbeta=nthet=nphai=1 (the input line will be 1,1,1 or 1 1 1)
C  When an aggregate rotates around the y-axis (i.e., rotates in the  
C    scattering plane), the number of orientations of the aggregate to be 
C    calculated is 19, then nbeta=1, nthet=19, nphai=1 (1,19,1 or 1 19 1) 
C  ---------------------------------------------------------------------
	write(6,'(a19,3i5)') 'nbeta,nthet,nphai: ',nbeta,nthet,nphai
	if(nbeta.gt.MOR.or.nthet.gt.MOR.or.nphai.gt.MOR) then
	   write(6,*) '***  parameter MOR too small  ***'
           write(6,*) ' MOR must be >nbeta,nthet,nphai given above'
           write(6,*) ' Please change MOR in the parameter line of the' 
           write(6,*) '   main code, recompile, then try again'
           stop
        endif 
	if(nbeta*nphai.gt.1) then 
	   idran=1
	else
	   idran=0
	endif
	nram=nbeta*nthet*nphai
	if(nram.lt.1) then 
	   write(6,*) 'please check (nbeta,nthet,nphai) in gmm01f.in'
	   stop
	endif
	READ(1,*) betami,betamx,thetmi,thetmx,phaimi,phaimx
C  ---------------------------------------------------------------------
C  betami,betamx,thetmi,thetmx,phaimi,phaimx are in degrees 
C  betami,betamx -- the range of the Euler angle of rotation (alpha)  
C                   to be calculated
C  thetmi,thetmx -- the range of the Euler angle of rotation (beta)  
C                   to be calculated
C  phaimi,phaimx -- the range of the Euler angle of rotation (beta)  
C                   to be calculated
C  an example:
C  When only a fixed orientation needs to be calculated, betami=betamx,
C  thetmi=thetmx, and phaimi=phaimx, this line could be, e.g.,  
C  (0. 0. 0. 0. 0. 0.), (0. 0. 90. 90. 0. 0.), (30. 30. 40. 40. 45. 45.) 
C  ----------------------------------------------------------------------
	if(nbeta.eq.1) betamx=betami
	if(nthet.eq.1) thetmx=thetmi
	if(nphai.eq.1) phaimx=phaimi
	write(6,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +                betami,betamx,thetmi,thetmx,phaimi,phaimx
	READ(1,*) idMie
C  ----------------------------------------------------------------------
C  idMie=1: calculating only coherent Mie-scattering, no interaction 
C  ----------------------------------------------------------------------
	write(6,'(a7,i3)') 'idMie: ',idMie
	READ(1,*) idd
C  ----------------------------------------------------------------------
C  When calculating a single orientation, put idd=0
C  When idd=1, the number of orientations to be calculated is doubled. 
C  Each orientation is coupled with the orientation that the aggregate  
C  is rotated by 90 degrees around the z axis. This insures that the 
C  averaged polarizations are zero when the scattering angle is either 0 
C  or 180 degrees, as is suppossed to be for an average over random 
C  orientations.
C  ----------------------------------------------------------------------
	if(nram.eq.1.and.idd.eq.1) then
	   write(6,'(/)') 
	   write(6,*) 'Warning: idd=1 while calculating a single'
	   write(6,*) '  particle-orientation, the results will be an'
	   write(6,*) '  average over two particle orientations'
	   write(6,*) '  please check: idd=0 or idd=1?'
	   write(6,'(/)')
	endif	
	READ(1,*) idc,iseed
C  ----------------------------------------------------------------------
C  idc for choosing in the three schemes for an orientation average 
C  idc=1 to devide the range of [thetmi,thetmx] by the cos(theta) scheme
C        (theta is the scattering angle, corresponding to the Euler angle 
C        of ritation "beta") 
C  idc=0 to devide the range of [thetmi,thetmx] by "degrees"
C  idc=-1 to pick up all the three Euler angles of rotation randomly using 
C         a random number generator 
C  iseed is the seed number for the random number generator, used only 
C        when idc=-1 (i.e., when idc=1 or 0, iseed has no function)
C        iseed can be an arbitrary positive integer
C  when calculating only one single orientation, i.e., 
C       nram=nbeta*nthet*nphai=1, idc and iseed have no function
C  ----------------------------------------------------------------------
	if(idpq.eq.1) then
	   nram=nphai+nthet
	   idc=0
	   idd=0
	   nbeta=1
	   betami=0.d0
	   betamx=betami
	   thetr(nthet+1)=0.d0
	   phair(nphai+1)=0.d0
	endif
	write(6,'(a5,i3)') 'idd: ',idd
	if(idd.eq.1) nram=2*nram
	write(6,'(a36,i5)') '# of orientations to be averaged: ',nram	
	if(idc.lt.0) then 
	   write(6,'(a11,i4,i12)') 'idc,iseed: ',idc,iseed
	else
	   write(6,'(a5,i3)') 'idc: ',idc
	endif	
	READ(1,*) factor1,factor2,MXINT
C  ----------------------------------------------------------------------
C  factor1 and factor2 are numerical factors used for improving the 
C  convergence behavior of the iterative solution process in solving 
C  interacting equations, which are in the range of [0,1]
C  factor1 is for x-polarized incident plane wave and factor2 for
C  y-polarized incident wave 
C  MXINT is the maximum number of iterations allowed in the iterative 
C  solution process
C  This code uses two alternative methods in the iterative solution of 
C  interacting equations: iteration scheme [see Fuller and Kattawar, Opt. 
C  Lett. 13, 90 (1988); Xu, Appl. Opt. 34, 4573 (1995)] and BI-CGSTAB, the 
C  stabilized Bi-Conjugate Gradient [see H.A. van der Vorst, SIAM J. Sci. 
C  Stat. Comput. 13, 631, (1992)].  
C  When factor1=0 or factor2=0, the code directly goes to BI-CGSTAB 
C  without using the other iteration scheme.  
C  When factor1=factor2=1 it is equivalent to  Fuller and Kattawar's 
C  order-of-scattering method. 
C  In many cases, a divergence will occur when factor1=factor2=1. Then,  
C  setting factor1,factor2<1 may help to converge to a numerical solution. 
C  When MXINT is exceeded, it will automatically switch to BI-CGSTAB.
C  ----------------------------------------------------------------------
	write(6,'(a,2f6.2)') 'Numerical factors for convergence:',
     +                        factor1,factor2
	write(6,'(a37,i5)') 'Maximum iterations allowed:',MXINT     
	READ(1,*) NADD
C  ----------------------------------------------------------------------
C  NADD is the number of terms to add to the scattering orders required  
C  by the Wiscombe's criterion, which can be negative or positive in the 
C  range of [-9,99] 
C  Normally, set NADD=0
C  ----------------------------------------------------------------------
	write(6,'(a41,i3)') 'Scat. orders added to Wiscombe criterion:',
     +                      NADD
   	READ(1,*) eps,small
C  ----------------------------------------------------------------------
C  eps: error tolerance for determining single-sphere field-expansion 
C       truncation, default: 1.d-20 
C       (the default value of 1.d-20 allows to use Wiscombi's criterion) 
C  small: error tolerance for the iterative solution process for solving 
C         the interacting scattering coefficients (1.d-6)
C  ----------------------------------------------------------------------
	write(6,'(a35,e10.2)') 'error tolerance for Mie-expansions:',eps
	write(6,'(a22,e10.2)') 'Convergence criterion:',small
	READ(1,*) fint
C  ----------------------------------------------------------------------
C  fint is the interaction index in the range of [0,1] (default: 0.02)
C  In the scattering calculations, a quantity "f" for a pair of component 
C  spheres is defined by f=(r_i+r_j)/d_{ij}, where (r_i,r_j) are the radii 
C  of the two spheres and d_{ij} is the separation distance of the two 
C  sphere-centers. When f<fint, interaction between the pair is considered
C  to be negligible and no interaction will be calculated between the two
C  spheres. When fint=0, no sphere is excluded in interaction calculations.
C  ----------------------------------------------------------------------
	if(fint.lt.0.d0.or.fint.gt.1.d0) then
	   fint=0.02d0
	   write(6,'(a37)') 'Interaction index: using default 0.02'
	else
	   write(6,'(a18,f7.3)') 'Interaction index:',fint
	endif
	READ(1,*) sang,pang
C  ----------------------------------------------------------------------
C  sang -- the scattering angle interval for output  
C  example: when sang=1, the results in output will be written for every 
C  degree of scattering angle from 0 to 180 
C  pang -- the azimuth angle interval for the two-dimensional Mueller 
C          matrix output
C      (1) For each azimuth angle, the number of scattering angles that 
C          will be calculated is the same as that determined by "sang", 
C          i.e., the number of scattering angles that will calculated is 
C          the same as calculated for the scattering plane of the azimuth 
C          angle=0 (180/sang +1).
C      (2) When pang = 0., no additional calculations for the scattering 
C          matrix map, i.e., calculating only the scattering plane of the 
C          azimuth angle = 0. 
C      (3) When pang>0, the number of azimuth angles in the range of 
C          (0,360) degrees that will be calculated in addition to 0 (360)  
C          degrees is npng-1 with npng=360/pang. For example, when 
C          pang=180, npng=2, the azimuth angles 0 (360) and 180 degrees 
C          will be calculated. In the output for the Mueller matrix at 
C          each azimuth angle, there are nang2 (=2*nang-1) sets of the 
C          16 elements. 
C  ----------------------------------------------------------------------
	write(6,'(a32,f7.3)') 'scat.-angle-interval in output: ',sang
	if(sang.le.0.d0) sang=1.d0
	nang=90.d0/sang+1.d0
	nang2=2*nang-1
	if(nang2.gt.nangmax) then
	   write(6,*) 'sang too small'
	   write(6,*) 'please increase sang in the input file gmm01f.in'
	   write(6,*) '  and try again, or'
	   write(6,*) '  increase nangmax in the parameter line of the'
	   write(6,*) '  main code, recompile, then try again'
	   stop
	endif
	write(6,'(a33,a15,f8.3)') 'azimuth-angle-interval in Mueller ', 
     +     'matrix output: ',pang
	if(pang.lt.0.0001d0) then
	   npng=1
	else
	   npng=360.0d0/pang
	endif
	if(npng.gt.ncmax) then
	   write(6,*) 'pang too small'
	   write(6,*) 'please increase pang in the input file gmm01f.in'
	   write(6,*) 'and try again, or increase ncmax in the parameter'
	   write(6,*) ' line of the main code, recompile, then try again'
	   stop
	endif
	close(1)
	write(6,'(/)')
	betami=betami*pih/90.d0
	betamx=betamx*pih/90.d0
	thetmi=thetmi*pih/90.d0
	thetmx=thetmx*pih/90.d0
	phaimi=phaimi*pih/90.d0
	phaimx=phaimx*pih/90.d0
	if(idc.gt.0) then
	   call orientcd(betami,betamx,thetmi,thetmx,phaimi,phaimx,
     +                MOR,MOR,MOR,nbeta,nthet,nphai,betar,thetr,phair)
        else
           call orientud(betami,betamx,thetmi,thetmx,phaimi,phaimx,
     +                MOR,MOR,MOR,nbeta,nthet,nphai,betar,thetr,phair)
        endif
	fileout='gmm01f.out'
	fileoutA='gmm01f.Aout'
	fileout1=fileout
	fileout2='pq'//fileout1
	if(idMie.eq.1) then
	   fileout='M'//fileout1
	   fileout2='Mpq'//fileout1
           write(6,*) '*** Calculating only coherent Mie-scattering ***'
           write(6,*) '*** No interaction included ********************'
	endif
	OPEN(UNIT=2,FILE=FLNAME,STATUS='OLD')
	READ(2,*) w
C  ----------------------------------------------------------------------
C  w -- incident wavelength
C  ----------------------------------------------------------------------
	READ(2,*) nL
C  ----------------------------------------------------------------------
C  nL -- number of spheres in the aggregate
C  ----------------------------------------------------------------------
	if(nL.gt.nLp) then
	  write(6,*) 'Parameter nLp too small, must be >', nL
	write(6,*) 'Change nLp in gmm01f.par, recompile, then try again'
	  stop
	endif
	if(nL.eq.1) idMie=1
C  ----------------------------------------------------------------------
C  input the configuration and particle parameters for the aggregate
C  each line includes 6 numbers:
C  x-, y-, z-coordinates of the sphere-center, the radius of the sphere 
C  in the same unit of the incident wavelength, the real and imaginary 
C  parts of the refractive index 
C  ----------------------------------------------------------------------
	do 1 i=1,nL
	   read(2,*,err=10) (r0(j,i),j=1,6)
	   x0=r0(1,i)
	   y0=r0(2,i)
	   z0=r0(3,i)
	   r00(1,i)=x0
	   r00(2,i)=y0
	   r00(3,i)=z0
	   if(r0(6,i).gt.0.d0) r0(6,i)=-r0(6,i)
	   if(r0(5,i).eq.1.d0.and.r0(6,i).eq.0.d0) goto 1 
	   gcs=gcs+r0(4,i)*r0(4,i)
	   gcv=gcv+r0(4,i)*r0(4,i)*r0(4,i)
 1	continue
        close(2)
	gcsr=dsqrt(gcs)
	gcvr=gcv**(1.d0/3.d0)
	goto 11
 10	write(6,*) 'fatal error in the input file'
	stop
 11	k=twopi/w
        xv=k*gcvr
        xs=k*gcsr
        write(6,'(a,f7.3,a,f7.3)') ' volume-equiv. xv: ',xv,
     +      '  surface-equiv. xs: ',xs
        write(6,'(/)')
	do i=1,nL
	   x(i)=k*r0(4,i)
	   ref(i)=dcmplx(r0(5,i),r0(6,i))
	   temp1=ref(i)
	enddo
	do j=1,np
	   do i=1,nL
	      aMie(i,j)=0.d0
	      bMie(i,j)=0.d0
	   enddo
	enddo
	nmax0=1
	do i=1,nL
	   if(i.eq.1) goto  12
	   if(x(i).eq.x(i-1).and.ref(i).eq.ref(i-1)) then
	      nmax(i)=nmax(i-1)
	      uvmax(i)=uvmax(i-1)
	      goto 15
	   endif
 12	   write(6,'(a,i3,a,f7.2)') 'sphere #',i, 
     +         '   individual size parameter: ',x(i)
c
c  calculating Mie-scattering coefficients for each spheres
c  the ratio method of Wang and van der Hulst is used in calculating  
c  Riccati-Bessel functions [see Wang and van der Hulst, Appl. Opt. 
c  30, 106 (1991), Xu, Gustafson, Giovane, Blum, and Tehranian, 
c  Phys. Rev. E 60, 2347 (1999)]
c 
 	   call abMiexud(x(i),ref(i),np,NXMAX,nmax(i),an,bn,NADD,
     +                   rsr0,rsi0,rsx0,px0,w1,w2,w3,w4,eps)
	   if(nmax(i).gt.np) then
	      write(6,*) ' Parameter np too small, must be >',nmax(i)
	      write(6,*) ' Please change np in gmm01f.par, recompile,' 
	      write(6,*) '   then try again'
	      stop
	   endif
	   uvmax(i)=nmax(i)*(nmax(i)+2)
	   write(6,'(a,1x,i4)') 
     +        ' Actual single-sphere expansion truncation:',nmax(i)
	   do j=1,nmax(i)
	      rsr(j,i)=rsr0(j)
	      rsi(j,i)=rsi0(j)
	      rsx(j,i)=rsx0(j)
	      px(j,i)=px0(j)
	      temp1=an(j)
	      temp2=bn(j)
	      if(j.eq.1.or.j.eq.nmax(i)) 
     +           write(6,'(i10,4e15.7)') j,temp1,
     +                 dimag(an(j)),temp2,dimag(bn(j))
	   enddo
 15	   do j=1,nmax(i)
	      aMie(i,j)=an(j)
	      bMie(i,j)=bn(j)
	      rsr(j,i)=rsr0(j)
	      rsi(j,i)=rsi0(j)
	      rsx(j,i)=rsx0(j)
	      px(j,i)=px0(j)
	   enddo
	   if(nmax(i).gt.nmax0) nmax0=nmax(i)
	enddo
	cextx=0.d0
	cexty=0.d0
	cabsx=0.d0
	cabsy=0.d0
	cscax=0.d0
	cscay=0.d0
	cprx=0.d0
	cpry=0.d0
	cbakx=0.d0
	cbaky=0.d0
	do i=1,nL
	   cextxi(i)=0.d0
	   cextyi(i)=0.d0
	   cabsxi(i)=0.d0
	   cabsyi(i)=0.d0
	   cscaxi(i)=0.d0
	   cscayi(i)=0.d0
	   cprxi(i)=0.d0
	   cpryi(i)=0.d0
	enddo
	do i=1,nang2
	   i11(i)=0.d0
	   i21(i)=0.d0
           i22(i)=0.d0
	   i12(i)=0.d0
	   do jc=1,npng
	      do j=1,4
	         do m=1,4
	            mue(j,m,jc,i)=0.d0
	         enddo
	      enddo
	   enddo
	enddo
	iram=0
	if(idpq.eq.1) then 
	   open(12,file=fileout2,status='unknown')
	   write(12,*) 'input file: ',FLNAME 
	endif
	write(6,'(/)')
	write(6,*) 'original input sphere-positions: '
	i=1
	write(6,'(i5,3f14.5)') i,r0(1,1),r0(2,1),r0(3,i)
	i=nL
	write(6,'(i5,3f14.5)') i,r0(1,i),r0(2,i),r0(3,i)
	if(idpq.eq.1) then
	   nphaic=nphai+1
	   phair(nphaic)=0.d0
	else
	   nphaic=nphai
	endif
	
	do ibeta=1,nbeta
	   do iphai=1,nphaic
	      if(idpq.eq.1.and.iphai.lt.nphaic) then
	         nthetc=1
	      else
	         nthetc=nthet
	      endif
	      do ithet=1,nthetc
	         if(idc.lt.0) then
	            betar(ibeta)=(betamx-betami)*ran1d(iseed)
	            phair(iphai)=(phaimx-phaimi)*ran1d(iseed)
	            thetr(ithet)=(thetmx-thetmi)*ran1d(iseed)
	         endif 
	         do 19 irot=1,2
	            if(idd.ne.1.and.irot.eq.2) goto 19
	            iram=iram+1
	            if(irot.eq.1) then
	               alph=0.d0
	            else
	               alph=pih
	            endif
	            ca=dcos(alph)
	            sa=dsin(alph)
	      	    beta=betar(ibeta)
	            cb=dcos(beta)
	            sb=dsin(beta)	   
	            do i=1,nL
	               x0=r00(1,i)
	               y0=r00(2,i)
	               r0(1,i)=cb*x0-sb*y0
	               r0(2,i)=sb*x0+cb*y0
	            enddo
	         phai=phair(iphai)
	         thet=thetr(ithet)
	         if(idpq.eq.1.and.nthetc.eq.1) thet=0.d0
	         cb=dcos(phai)
	         sb=dsin(phai)
	   	 cz=dcos(thet)
	         sz=dsin(thet)
	         if(iram.eq.1.or.iram/50*50.eq.iram) then 
	            write(6,'(a,2i5)') 'iram & nram: ', iram,nram
	         endif     
	         do i=1,nL
	            x0=r0(1,i)
	            y0=r0(2,i)
	            z0=r00(3,i)
	            r0(1,i)=ca*cz*x0-(ca*sz*sb+sa*cb)*y0
     +                              +(ca*sz*cb-sa*sb)*z0
	            r0(2,i)=sa*cz*x0-(sa*sz*sb-ca*cb)*y0
     +                              +(sa*sz*cb+ca*sb)*z0
	            r0(3,i)=-sz*x0-cz*sb*y0+cz*cb*z0	
	         enddo
	         if(iram.eq.1.or.iram/50*50.eq.iram) then 
	            i=1	         
	            write(6,'(i5,3f14.5)') i,r0(1,i),r0(2,i),r0(3,i)
	            i=nL
	            write(6,'(i5,3f14.5)') i,r0(1,i),r0(2,i),r0(3,i)
	         endif
C  ----------------------------------------------------------------------
C  calculating constants, Gaunt, rotational and translation coefficients
C  ----------------------------------------------------------------------
	n0=nmax0+2
	fnr(0)=0.d0
        do n=1,2*n0
          fnr(n)=dsqrt(dble(n))
        enddo
        bcof(0)=1.d0
        do n=0,n0-1
         bcof(n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n)/fnr(n+1)/fnr(n+1)
        enddo
C
C  the formulation used here for the calculation of Gaunt coefficients 
C  can be found in Bruning and Lo, IEEE Trans. Anten. Prop. Ap-19, 378 
C  (1971) and Xu, J. Comput. Appl. Math. 85, 53 (1997), J. Comput. Phys. 
C  139, 137 (1998)
C
        call cofsrd(nmax0)	
	call cofd0(nmax0)
        call cofnv0(nmax0)
        call gau0(nmax0)
 	do i=1,nL
           do j=i+1,nL
              ij=(j-1)*(j-2)/2+j-i
              x0=r0(1,i)-r0(1,j)
              y0=r0(2,i)-r0(2,j)
              z0=r0(3,i)-r0(3,j)
              call carsphd(x0,y0,z0,d,xt,sphi,cphi)
              temp=(r0(4,i)+r0(4,j))/d
              if(temp.le.fint) goto 16             
              ephi=dcmplx(cphi,sphi)
              nlarge=max(nmax(i),nmax(j))
              do m=1,nlarge
                 ek(m,ij)=ephi**m
              enddo
              xd=k*d
              nbes=2*nlarge+1
	      call besseljd(nbes,xd,besj)
	      call besselyd(nbes,xd,besy)
C
C  subroutine rotcoef.f calculates (-1)^{m+u)d_{mu}^n where d_{mu}^n is
C  the "reduced rotation matrix elements"
C  rotcoef.f is originally written by Mackowski (taken from scsmfo1b.for 
C  developed by Mackowski, Fuller, and Mishchenko)   
C
              call rotcoef(xt,nlarge)
              irc=0
              do n=1,nlarge
                 n1=n*(n+1)
                 do u=-n,n
                    do m=-n,n
                       imn=n1+m
                       irc=irc+1
                       drot(irc,ij)=dc(u,imn)
                    enddo
                 enddo
              enddo
              itrc=0
              nsmall=min(nmax(i),nmax(j))
C
C  the formulation used here for the calculation of vector translation 
C  coefficients are from Cruzan, Q. Appl. Math. 20, 33 (1962) and  
C  Xu, J. Comput. Phys. 139, 137 (1998)
C
              do m=-nsmall,nsmall
                 n1=max(1,iabs(m))
                 do n=n1,nlarge
                    do v=n1,nlarge
                       itrc=itrc+1
                       call cofxuds0(nmax0,m,n,v,besj,besy,
     +                         atr0(itrc,ij),btr0(itrc,ij),
     +                         atr1(itrc,ij),btr1(itrc,ij))
                    enddo
                 enddo
              enddo
 16        continue
           enddo
        enddo
	indpol=0
        factor=factor1
        write(6,'(a8,i4,a32)') 
     +     ' orien.#',iram,'  Solving for x-pol. inci. state'
 18	do imn=1,nmp
	   do i=1,nL
	      p0(i,imn)=0.d0
	      q0(i,imn)=0.d0
	      as(i,imn)=0.d0
	      bs(i,imn)=0.d0
	   enddo
	enddo
C  ------------------------------------------------
C  calculating incident wave expansion coefficients
C  ------------------------------------------------
	do 20 i=1,nL
	   cz=dcos(k*r0(3,i))
	   sz=dsin(k*r0(3,i))
	   cmz=0.5d0*dcmplx(cz,sz)
	   do 21 n=1,nmax(i)
              imn=n*n+n+1
              A=fnr(2*n+1)*cmz
	      p0(i,imn)=aMie(i,n)*A
	      q0(i,imn)=bMie(i,n)*A
	      p0(i,imn-2)=-p0(i,imn)
	      q0(i,imn-2)=q0(i,imn)
	      if(indpol.gt.1) then
	         p0(i,imn)=p0(i,imn)*cin
	         q0(i,imn)=q0(i,imn)*cin
	         p0(i,imn-2)=p0(i,imn)
	         q0(i,imn-2)=-q0(i,imn)
	      endif
	      as(i,imn)=p0(i,imn)
	      bs(i,imn)=q0(i,imn)
	      as(i,imn-2)=p0(i,imn-2)
	      bs(i,imn-2)=q0(i,imn-2)
 21	   continue
 20	continue
c  ------------------------------------------------------------
c  begins iteration process to solve the interaction equations
c  for partial interacting scatteing coefficients
c  factor1,factor2=0: BI-CGSTAB [see H.A. van der Vorst, SIAM,  
c  J. Sci. Stat. Comput. 13, 631 (1992)]
c  factor1=factor2=1: order-of-scattering [Fuller and Kattawar, 
c  Opt. Lett. 13, 90 & 1063 (1988)]
c  0<factor1,factor2<1: an iteration scheme [Xu, Appl. Opt. 34,
c  4573 (1995)] 
c  When the number of iterations exceeds the maximum MXINT, it 
c  switches to use BI-CGSTAB.
c  In the iteration solution of interacting equations, a vector
c  of wave expansion referred to a sphere-center needs to be 
c  translated to other sphere-centers, which demands the storage 
c  and evaluation of a large number of vector translation
c  coefficients. This code uses Mackowski's three-step technique
c  (rotation-translation-rotation method) [see Mackowski, Proc. 
c  R. Soc. Lond. A 433, 599 (1991)] for the translation process, 
c  which significantly reduces both computer memory usage and 
c  computing time requirements by decomposition of the vector 
c  translation coefficients into rotational and axial 
c  translational parts. By using Mackowski's three-step method, 
c  the calculations involve only translation along z-axis.     
c  -------------------------------------------------------------
	if(idMie.eq.1) goto 490
	if(nL.eq.1) goto 490
	
	do i=1,nL
	   ind(i)=0
	   c0i(i)=0.d0
	   do n=1,nmax(i)
	      imn=n*n+n+1
	      c0i(i)=c0i(i)+p0(i,imn)*dconjg(p0(i,imn)) 
	      c0i(i)=c0i(i)+q0(i,imn)*dconjg(q0(i,imn))
	      c0i(i)=c0i(i)+p0(i,imn-2)*dconjg(p0(i,imn-2)) 
	      c0i(i)=c0i(i)+q0(i,imn-2)*dconjg(q0(i,imn-2))
 	   enddo
 	enddo
 	niter=1
	if(factor1.lt.0.001d0.or.factor2.lt.0.001d0) goto 61
	if(iram.eq.1) 
     +     write(6,*) 'Starting iteration solution process'	
 	do i=1,nL
	   do imn=1,uvmax(i)
	      as0(i,imn)=p0(i,imn)
	      bs0(i,imn)=q0(i,imn)
 	   enddo
 	enddo
 60	call trans(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,
     +              drot,as0,bs0,as1,bs1,ind)
 	do i=1,nL
 	   if(ind(i).gt.0) goto 601
 	   c1i(i)=0.d0
 	   do imn=1,uvmax(i)
 	      n=dsqrt(dble(imn))
 	      as0(i,imn)=p0(i,imn)-aMie(i,n)*as1(i,imn)
 	      bs0(i,imn)=q0(i,imn)-bMie(i,n)*bs1(i,imn)
 	      A=as0(i,imn)-as(i,imn)
 	      B=bs0(i,imn)-bs(i,imn)
 	      c1i(i)=c1i(i)+A*dconjg(A)
 	      c1i(i)=c1i(i)+B*dconjg(B)
 	      as0(i,imn)=as(i,imn)+factor*A
 	      bs0(i,imn)=bs(i,imn)+factor*B
 	      as(i,imn)=as0(i,imn)
 	      bs(i,imn)=bs0(i,imn) 
 	   enddo
 601       continue
 	enddo
 	cext0=0.d0
 	cext1=0.d0
 	do i=1,nL
 	   if(ind(i).gt.0) goto 602
 	   cext0=cext0+c0i(i)
 	   cext1=cext1+c1i(i)
 	   temp=c1i(i)/c0i(i)
 	   if(temp.lt.small) ind(i)=1
 602       continue
	enddo
	temp=cext1/cext0
	if(temp.lt.small) goto 490
	if(iram.eq.1.or.iram.eq.nram) then
     	         write(6,'(a11,i4,2x,e15.7)') 'iteration #',
     +                                        niter,temp
	endif
	if(niter.gt.MXINT) then
	   write(6,*) '  *** Maximum iterations exceeded ***'
	   write(6,*) '  *** Switched to Bi-CGSTAB method***'
	   do i=1,nL
	      ind(i)=0
	      do imn=1,uvmax(i)
	         as(i,imn)=p0(i,imn)
	         bs(i,imn)=q0(i,imn)
 	      enddo
 	   enddo
 	   niter=1
 	   goto 61
 	endif
	niter=niter+1
	goto 60
 61	if(iram.eq.1) 
     +     write(6,*) 'Starting Bi-CGSTAB solution process'
  	call trans(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,
     +              drot,as,bs,as1,bs1,ind)
 	do i=1,nL
 	   c1i(i)=0.d0
 	   do imn=1,uvmax(i)
 	      n=dsqrt(dble(imn))
 	      as1(i,imn)=-aMie(i,n)*as1(i,imn)
 	      bs1(i,imn)=-bMie(i,n)*bs1(i,imn)
 	      c1i(i)=c1i(i)+as1(i,imn)*dconjg(as1(i,imn))
 	      c1i(i)=c1i(i)+bs1(i,imn)*dconjg(bs1(i,imn))
 	   enddo
 	enddo
 	temp=0.d0
 	do i=1,nL
 	   cext0=c1i(i)/c0i(i)
 	   if(cext0.lt.small) ind(i)=1
 	   if(cext0.gt.temp) temp=cext0
	enddo
	if(temp.lt.small) goto 490
        A0=0.d0
        do i=1,nL
           if(ind(i).gt.0) goto 611
	   do imn=1,uvmax(i)
	      asp(i,imn)=as1(i,imn)
	      bsp(i,imn)=bs1(i,imn)
	      as0(i,imn)=as1(i,imn)
	      bs0(i,imn)=bs1(i,imn)
	      A0=A0+as1(i,imn)*as1(i,imn)
	      A0=A0+bs1(i,imn)*bs1(i,imn)
 	   enddo
 611       continue
 	enddo
  62	call trans(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,
     +              drot,asp,bsp,asv,bsv,ind)
 	do i=1,nL
 	   if(ind(i).gt.0) goto 621
 	   do imn=1,uvmax(i)
 	      n=dsqrt(dble(imn))
 	      asv(i,imn)=aMie(i,n)*asv(i,imn)+asp(i,imn)
 	      bsv(i,imn)=bMie(i,n)*bsv(i,imn)+bsp(i,imn)
 	   enddo
 621       continue
 	enddo
 	A=0.d0
 	do i=1,nL
 	   if(ind(i).gt.0) goto 622
 	   do imn=1,uvmax(i)
 	      A=A+asv(i,imn)*as1(i,imn)
 	      A=A+bsv(i,imn)*bs1(i,imn)
 	   enddo
 622       continue
 	enddo
	Aj=A0/A
	do i=1,nL
	   if(ind(i).gt.0) goto 623
 	   do imn=1,uvmax(i)
 	      asc(i,imn)=as0(i,imn)-Aj*asv(i,imn)
 	      bsc(i,imn)=bs0(i,imn)-Aj*bsv(i,imn)
 	   enddo
 623       continue
 	enddo
	call trans(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,
     +              drot,asc,bsc,ast,bst,ind)
 	do i=1,nL
 	   if(ind(i).gt.0) goto 624
 	   do imn=1,uvmax(i)
 	      n=dsqrt(dble(imn))
 	      ast(i,imn)=aMie(i,n)*ast(i,imn)+asc(i,imn)
 	      bst(i,imn)=bMie(i,n)*bst(i,imn)+bsc(i,imn)
 	   enddo
 624       continue
 	enddo
 	A2=0.d0
 	B2=0.d0
 	do i=1,nL
 	   if(ind(i).gt.0) goto 625
 	   do imn=1,uvmax(i)
 	      A2=A2+ast(i,imn)*asc(i,imn)
 	      A2=A2+bst(i,imn)*bsc(i,imn)
 	      B2=B2+ast(i,imn)*ast(i,imn)
 	      B2=B2+bst(i,imn)*bst(i,imn)
 	   enddo
 625       continue
 	enddo
	Bj=A2/B2
	do i=1,nL
	   if(ind(i).gt.0) goto 626
	   c1i(i)=0.d0
 	   do imn=1,uvmax(i)
 	      Aj2=Aj*asp(i,imn)+Bj*asc(i,imn)
 	      Bj2=Aj*bsp(i,imn)+Bj*bsc(i,imn)
 	      c1i(i)=c1i(i)+Aj2*dconjg(Aj2)
 	      c1i(i)=c1i(i)+Bj2*dconjg(Bj2)
 	      as(i,imn)=as(i,imn)+Aj2
 	      bs(i,imn)=bs(i,imn)+Bj2
 	   enddo
 626       continue
 	enddo
 	cext0=0.d0
 	cext1=0.d0
 	do i=1,nL
 	   if(ind(i).gt.0) goto 627
 	   cext0=cext0+c0i(i)
 	   cext1=cext1+c1i(i)
 627       continue
	enddo
	temp=cext1/cext0
	if(temp.lt.small) goto 490
	if(niter.gt.MXINT) then
	   write(6,*) 'Caution:'
	   write(6,*) '  *** Maximum iterations exceeded ***'	
	   goto 490
	endif
	if(iram.eq.1.or.iram.eq.nram) then
     	   write(6,'(a11,i4,2x,e15.7)') 'iteration #',niter,temp
	endif
	do i=1,nL
 	   if(ind(i).gt.0) goto 628	   
	   do imn=1,uvmax(i)
 	      as0(i,imn)=asc(i,imn)-Bj*ast(i,imn)
 	      bs0(i,imn)=bsc(i,imn)-Bj*bst(i,imn)
 	   enddo
 628       continue
 	enddo
 	A2=0.d0
 	do i=1,nL
 	   ast(i,1)=0.d0
 	   if(ind(i).gt.0) goto 629
 	   do imn=1,uvmax(i)
 	      ast(i,1)=ast(i,1)+as0(i,imn)*as1(i,imn)
 	      ast(i,1)=ast(i,1)+bs0(i,imn)*bs1(i,imn)
 	   enddo
 	   A2=A2+ast(i,1)
 629       continue
 	enddo
  	B0=A2/A0
 	B0=B0*Aj/Bj
 	do i=1,nL
 	   if(ind(i).gt.0) goto 630
 	   do imn=1,uvmax(i)
 	      asp(i,imn)=as0(i,imn)+B0*(asp(i,imn)-Bj*asv(i,imn))
 	      bsp(i,imn)=bs0(i,imn)+B0*(bsp(i,imn)-Bj*bsv(i,imn))
 	   enddo
 630       continue
 	enddo 
 	A0=0.d0
 	do i=1,nL
 	   if(ind(i).gt.0) goto 631
 	   cext0=c1i(i)/c0i(i)
 	   if(cext0.lt.small) then
 	      ind(i)=1
 	      goto 631
 	   endif
 	   A0=A0+ast(i,1)
 631       continue
	enddo
	niter=niter+1
 	goto 62
c  --------------------------------------------------------------
c  computing total and differential scattering cross sections and 
c  efficiencies for radiation pressure, using the formulas given  
c  by Xu, Physics Letters A 249, 30 (1998)
c  --------------------------------------------------------------
 490    do i=1,nL
	   ind(i)=0
	enddo
c	open(17,file='scacof.dat',status='unknown')
c	do 282 i=1,nL
c	   write(17,'(a2,i3)') 'i=',i
c	   do 281 j=1,uvmax(i)
c	      n=dsqrt(dble(j))
c	      m=j-n*n-n
c	write(17,'(2i5,4e17.8)') m,n,dble(as(i,j)),dimag(as(i,j)),
c     +                               dble(bs(i,j)),dimag(bs(i,j))
c 281      continue
c 282   continue
c	close(17)
 	call trans(nL,r0,nmax,uvmax,fint,atr1,btr1,ek,
     +              drot,as,bs,as1,bs1,ind)
        do i=1,nL
 	   do imn=1,uvmax(i)
	      at(imn)=as(i,imn)+as1(i,imn)
	      bt(imn)=bs(i,imn)+bs1(i,imn)
	   enddo
	   do n=1,nmax(i)
	      n1=n+1
	      n2=2*n
	      rn=1.0d0/dble(n*n1)
	      p=fnr(n)*fnr(n+2)/fnr(n2+1)/fnr(n2+3)/dble(n1)
	      t=fnr(n-1)*fnr(n+1)/fnr(n2-1)/fnr(n2+1)/dble(n)
	      sc=0.d0
	      temp=0.d0
	      do m=-n,n
	         iL=n*(n+1)+m
                 sc=sc+dble(dconjg(as(i,iL))*at(iL))                   
 	         sc=(sc+dble(dconjg(bs(i,iL))*bt(iL)))
 	         rm=dble(m)*rn
 	         A0=rm*bt(iL)
 	         B0=rm*at(iL)              
 	         if(n.eq.nmax(i)) goto 51
 	         u=(n+1)*(n+2)+m
 	         fnp=fnr(n+m+1)*fnr(n-m+1)*p
 	         A0=A0+fnp*at(u)
 	         B0=B0+fnp*bt(u)
 51	         if(n.eq.1.or.iabs(m).gt.n-1) goto 52
 	         u=(n-1)*n+m
 	         fn=fnr(n+m)*fnr(n-m)*t
 	         A0=A0+fn*at(u)
 	         B0=B0+fn*bt(u)
 52	         temp=temp+dble(dconjg(as(i,iL))*A0)
	         temp=temp+dble(dconjg(bs(i,iL))*B0)
              enddo
              if(indpol.lt.1) then
	         cscaxi(i)=cscaxi(i)+sc
	         cscax=cscax+sc
	         cprxi(i)=cprxi(i)+temp
	         cprx=cprx+temp
	      else
	         cscayi(i)=cscayi(i)+sc
	         cscay=cscay+sc
	         cpryi(i)=cpryi(i)+temp
	         cpry=cpry+temp
	      endif
	   enddo
	enddo
c  ---------------------------------------------------------------
c  computing total and differential extinction and absorption 
c  cross-sections, using the formulas given by Mackowski, Proc. R. 
c  Soc. Lond. A 433, 599 (1991) and Xu, Appl. Opt. 36, 9496 (1997) 
c  ---------------------------------------------------------------
 	do j=1,nL
	   cz=dcos(k*r0(3,j))
	   sz=dsin(k*r0(3,j))
	   cmz=dcmplx(cz,-sz)
	   A=0.d0
	   B=0.d0
	   do n=1,nmax(j)
	      rn=fnr(2*n+1)
	      m0=n*n+n+1
	      u0=n*n+n-1
	      A=A+rn*(as(j,m0)+bs(j,m0))
	      B=B+rn*(as(j,u0)-bs(j,u0))
	   enddo
	   if(indpol.lt.1) then
	      cextxi(j)=cextxi(j)+dble((A-B)*cmz)
	      cextx=cextx+dble((A-B)*cmz)
	   else
	      cextyi(j)=cextyi(j)-dimag((A+B)*cmz)
	      cexty=cexty-dimag((A+B)*cmz)
	   endif
	enddo
        do j=1,nL
	   do n=1,nmax(j)
	      A=ref(j)*dcmplx(rsr(n,j),-rsi(n,j))
	      temp1=-dimag(A)
	      A=px(n,j)*(ref(j)*rsx(n,j)
     +          -dcmplx(rsr(n,j),rsi(n,j)))
	      temp=cdabs(A)*cdabs(A)
	      if(temp.eq.0.d0) then 
	         dn=0.d0
              else
	         dn=temp1/temp
              endif
              A=dcmplx(r0(5,j),-r0(6,j))
     +          *dcmplx(rsr(n,j),-rsi(n,j))
	      temp1=-dimag(A)
              A=px(n,j)*(rsx(n,j)
     +          -ref(j)*dcmplx(rsr(n,j),rsi(n,j)))
	      temp=cdabs(A)*cdabs(A)
	      if(temp.eq.0.d0) then 
	         cn=0.d0
	      else
	         cn=temp1/temp
              endif
	      do m=-n,n
	         i=n*n+n+m      
	         temp1=dn*cdabs(as(j,i))*cdabs(as(j,i))
     +	               +cn*cdabs(bs(j,i))*cdabs(bs(j,i))
	         if(indpol.lt.1) then
	            cabsxi(j)=cabsxi(j)+temp1
	            cabsx=cabsx+temp1
	         else
	            cabsyi(j)=cabsyi(j)+temp1
	            cabsy=cabsy+temp1
	         endif
	      enddo
	   enddo
	enddo
c  ------------------------------------------------------------
c  computing two-dimensional angular distribution of the total 
c  scattered field in an array of (npng+1)Xnang2 specified by
c  the input of (sang,pang)
c  the formulas used here for the scattering amplitude matrix 
c  elements are from Xu, Appl. Opt. 34, 4573 (1995); 36, 9496 
c  (1997) [see also Xu and Wang, Phys. Rev. E 58, 3931 (1998)]
c  ------------------------------------------------------------
	do i=1,nang
	   iang=2*nang-i
	   dang(i)=sang*dble(i-1)
	   dang(iang)=180.0d0-dang(i)
	   theta=dang(i)*pione/180.0d0
	   xt=dcos(theta)
	   st=dsin(theta)
	   call tipitaud(nmax0,xt)
	   do jc=1,npng
	      azphi=pang*pih*dble(jc-1)/90.0d0
	      sphi=dsin(azphi)
	      cphi=dcos(azphi)
	      do imn=1,nmp
	         at(imn)=0.d0
	         bt(imn)=0.d0
	         atj(imn)=0.d0
	         btj(imn)=0.d0
	      enddo
  	      do 31 j=1,nL
	         sb=r0(1,j)*cphi+r0(2,j)*sphi
	         sb=sb*st
	         cb=r0(3,j)*xt
	         cz=k*(sb+cb)          ! for theta
	         sz=k*(sb-cb)          ! for pi-theta
	         A=dcmplx(dcos(cz),-dsin(cz))  ! for theta
	         B=dcmplx(dcos(sz),-dsin(sz))  ! for pi-theta
	         do 32 imn=1,uvmax(j)
                    n=dsqrt(dble(imn))
	            if(n.gt.nmax(j)) goto 31
	            if(idMie.gt.0) then 
	               m=imn-n*n-n
	               if(iabs(m).ne.1) goto 32
	            endif
	            at(imn)=at(imn)+A*as(j,imn)
	            bt(imn)=bt(imn)+A*bs(j,imn)
	            atj(imn)=atj(imn)+B*as(j,imn)
	            btj(imn)=btj(imn)+B*bs(j,imn)
 32	         continue
 31	      continue
	      if(indpol.lt.1) then
	         s2x(jc,i)=0.d0
	         s4x(jc,i)=0.d0
	         s2x(jc,iang)=0.d0
	         s4x(jc,iang)=0.d0
	      else
	         s3y(jc,i)=0.d0
	         s1y(jc,i)=0.d0
	         s3y(jc,iang)=0.d0
	         s1y(jc,iang)=0.d0
	      endif
	      A=0.d0
	      B=0.d0
	      Aj=0.d0
	      Bj=0.d0
	      do j=1,nmax0
	         imn=(j-1)*(j+2)/2+1 !sequence # for pi(0,n),tau(0,n)
	         u=j*j+j
	         A=A+at(u)*tau(imn)
	         B=B+bt(u)*tau(imn)
	         if(i.ne.iang) then
	            t=(-1)**(j+1)
	            Aj=Aj+atj(u)*tau(imn)*t
	            Bj=Bj+btj(u)*tau(imn)*t
	         endif
	      enddo
	      if(indpol.lt.1) then
	         s2x(jc,i)=s2x(jc,i)+A*cphi
	         s4x(jc,i)=s4x(jc,i)+B*dcmplx(0.d0,-1.d0)*cphi
	         if(i.ne.iang) then
	            s2x(jc,iang)=s2x(jc,iang)+Aj*cphi
	            s4x(jc,iang)=s4x(jc,iang)+Bj*dcmplx(0.d0,-1.d0)*cphi
	         endif
	      else
	         s3y(jc,i)=s3y(jc,i)-A*cphi
	         s1y(jc,i)=s1y(jc,i)+B*dcmplx(0.d0,1.d0)*cphi
	         if(i.ne.iang) then
	            s3y(jc,iang)=s3y(jc,iang)-Aj*cphi
	            s1y(jc,iang)=s1y(jc,iang)+Bj*dcmplx(0.d0,1.d0)*cphi
	         endif
	      endif
	      rm=1.d0
	      do 302 m=1,nmax0
	         A=0.d0
	         B=0.d0
	         A2=0.d0
	         B2=0.d0
	         Aj=0.d0
	         Bj=0.d0
	         Aj2=0.d0
	         Bj2=0.d0
	         rm=-rm
	         do 303 j=m,nmax0
	            imn=(j-1)*(j+2)/2+m+1 !sequence # for pi(m,n),tau(m,n)
	            u=j*j+j+m             !sequence # for at(m,n)
	            v=u-2*m               !sequence # for at(-m,n)
	            A0=at(u)*tau(imn)+bt(u)*pi(imn)
	            B0=rm*(at(v)*tau(imn)-bt(v)*pi(imn))
	            A=A+A0+B0
	            A2=A2+A0-B0
  	            A0=at(u)*pi(imn)+bt(u)*tau(imn)
	            B0=rm*(at(v)*pi(imn)-bt(v)*tau(imn))
	            B=B+A0-B0
	            B2=B2+A0+B0
	            if(i.ne.iang) then
	               t=(-1)**(j+m+1)
	               p=-t
	               A0=atj(u)*tau(imn)*t+btj(u)*pi(imn)*p
	               B0=rm*(atj(v)*tau(imn)*t-btj(v)*pi(imn)*p)
	               Aj=Aj+A0+B0
	               Aj2=Aj2+A0-B0
	               A0=atj(u)*pi(imn)*p+btj(u)*tau(imn)*t
	               B0=rm*(atj(v)*pi(imn)*p-btj(v)*tau(imn)*t)
	               Bj=Bj+A0-B0
	               Bj2=Bj2+A0+B0
	            endif
 303             continue	         
	         temp=dble(m-1)*azphi
	         sb=dsin(temp)
	         cb=dcos(temp)
	         if(indpol.lt.1) then
	            s2x(jc,i)=s2x(jc,i)+A*cb+A2*dcmplx(0.d0,1.d0)*sb
	            s4x(jc,i)=s4x(jc,i)+B*dcmplx(0.d0,-1.d0)*cb+B2*sb
	            if(i.ne.iang) then
	               s2x(jc,iang)=s2x(jc,iang)+Aj*cb
	               s2x(jc,iang)=s2x(jc,iang)+Aj2*dcmplx(0.d0,1.d0)*sb
	               s4x(jc,iang)=s4x(jc,iang)+Bj*dcmplx(0.d0,-1.d0)*cb
	               s4x(jc,iang)=s4x(jc,iang)+Bj2*sb
	            endif
	         else
	            s3y(jc,i)=s3y(jc,i)-A*cb
	            s3y(jc,i)=s3y(jc,i)+A2*dcmplx(0.d0,-1.d0)*sb
	            s1y(jc,i)=s1y(jc,i)+B*dcmplx(0.d0,1.d0)*cb
	            s1y(jc,i)=s1y(jc,i)-B2*sb
	            if(i.ne.iang) then
	               s3y(jc,iang)=s3y(jc,iang)-Aj*cb
	               s3y(jc,iang)=s3y(jc,iang)+Aj2*dcmplx(0.d0,-1.d0)*sb
	               s1y(jc,iang)=s1y(jc,iang)+Bj*dcmplx(0.d0,1.d0)*cb
	               s1y(jc,iang)=s1y(jc,iang)-Bj*sb
	            endif
	         endif  
 302	      continue
	   enddo
	enddo
	indpol=indpol+2
        factor=factor2
	if(indpol.lt.3) then
           write(6,'(a8,i4,a32)') 
     +        ' orien.#',iram,'  Solving for y-pol. inci. state'
	   goto 18
	endif
	do i=1,nang2
	   i22(i)=i22(i)+cdabs(s2x(1,i))*cdabs(s2x(1,i))
	   i21(i)=i21(i)+cdabs(s4x(1,i))*cdabs(s4x(1,i))
           i11(i)=i11(i)+cdabs(s1y(1,i))*cdabs(s1y(1,i))
	   i12(i)=i12(i)+cdabs(s3y(1,i))*cdabs(s3y(1,i))
	   do jc=1,npng
	      call mueller(s1y(jc,i),s2x(jc,i),s3y(jc,i),s4x(jc,i),smue)
	      do j=1,4
	         do m=1,4
	            mue(j,m,jc,i)=mue(j,m,jc,i)+smue(j,m)
	         enddo
	      enddo
	   enddo
	enddo
	cbakx=cbakx+cdabs(s2x(1,nang2))*cdabs(s2x(1,nang2))
	cbaky=cbaky+cdabs(s1y(1,nang2))*cdabs(s1y(1,nang2))
        cz=4.0d0/(gcs*k*k)
	if(idpq.eq.1) then
	   temp1=s2x(1,1)*cz
	   temp2=s1y(1,1)*cz
	   write(12,'(i6,1x,2f6.1,1x,4e14.6)') iram,
     +           thet/pih*90.d0,phai/pih*90.0d0,temp1,
     +           dimag(s2x(1,1))*cz,temp2,dimag(s1y(1,1))*cz
          if(idc.gt.0) then 
	      write(6,'(i7,3f7.1)') iram,dang(1),dcos(thet),
     +	                          phai*90.d0/pih
	   else 
	      write(6,'(i7,3f7.1)') iram,dang(1),
     +	                     thet*90.d0/pih,phai*90.d0/pih
	   endif
	endif
	if(idMie.eq.1) goto 19
	if(iram.lt.10) then
	   write(cnr1,'(i1)') iram
	   tailn='00'//cnr1
	else
	   if(iram.lt.100) then
	      write(cnr2,'(i2)') iram
	      tailn='0'//cnr2
	   else
	      write(cnr3,'(i3)') iram
	      tailn=cnr3
	   endif
	endif
	if(nram.eq.1) then
	  if(iram.eq.1) then
            write(6,'(a45,a16)') 
     +        'Output file for scattering amplitude matrix: ',fileoutA
          endif
	  open(23,file=fileoutA,status='unknown')
	  write(23,'(a20,a30)')fileoutA,'(Scattering amplitude matrix)'
	write(23,'(a)') 'The results are for the x-z plane of phi=0 only' 
	  write(23,'(a,f7.4,3x,a,a20)')  
     +        'wavelength: ',w,'input filename: ',FLNAME
	write(23,'(a)') 'sphere#,x,y,z,radius,complex refractive index:'
	  do i=1,nL
	     write(23,'(i5,6f11.4)') i,r0(1,i),r0(2,i),r0(3,i),
     +                               r0(4,i),r0(5,i),r0(6,i)
	  enddo
	  write(23,'(a)') 'scattering angle, s2x(complex), s3y(complex)'
	  write(23,'(a)') '                  s4x(complex), s1y(complex)'
	  do i=1,nang2 
	     write(23,'(f8.2,4e15.6)') dang(i),dble(s2x(1,i)),
     +                 dimag(s2x(1,i)),dble(s3y(1,i)),dimag(s3y(1,i))
	     write(23,'(8x,4e15.6)') dble(s4x(1,i)),dimag(s4x(1,i)),
     +                               dble(s1y(1,i)),dimag(s1y(1,i))
	  enddo
	  close(23)
	endif
c
 19           continue
	      enddo
	   enddo
	enddo
	if(iram.ne.nram) then
	   write(6,*) 'Note: iram is not equal to nram!'
	   write(6,'(a,i6,a,i6)') 'iram: ',iram,'nram: ',nram
	endif
C
C  ending doloop for nbeta,nphi,ntheta
C
	if(idpq.eq.1) close(12)
	cz=nram
	do i=1,nang2
	   i11(i)=i11(i)/cz
	   i21(i)=i21(i)/cz
           i22(i)=i22(i)/cz
	   i12(i)=i12(i)/cz
	   inat(i)=i11(i)+i22(i)+i12(i)+i21(i)
	   pol(i)=(i11(i)-i22(i))/inat(i)
	   do jc=1,npng
	      do j=1,4
	         do m=1,4
	            mue(j,m,jc,i)=mue(j,m,jc,i)/cz
	         enddo
	      enddo
	   enddo
	enddo
	cz=cz*k*k
	cscax=2.0d0*twopi*cscax/cz
	cscay=2.0d0*twopi*cscay/cz
	csca=0.5d0*(cscax+cscay)
	cextx=twopi*cextx/cz
	cexty=twopi*cexty/cz
	cext=0.5d0*(cextx+cexty)
	cabsx=2.0d0*twopi*cabsx/cz
	cabsy=2.0d0*twopi*cabsy/cz
	cabs=0.5d0*(cabsx+cabsy)
        cprx=2.0d0*twopi*cprx/cz
        cpry=2.0d0*twopi*cpry/cz
        cpr=0.5d0*(cprx+cpry)
	assym=(cprx+cpry)/(cscax+cscay)
	assym0=0.5d0*(cprx/cscax+cpry/cscay)
	cbakx=2.0d0*twopi*cbakx/cz
	cbaky=2.0d0*twopi*cbaky/cz
	cbak=0.5d0*(cbakx+cbaky)
	write(6,'(5x,2e15.6)') assym,assym0
	write(6,'(6x,a4,9x,a4,9x,a4,9x,a4,9x,a3,8x,a12)')
     +    'Cext','Cabs','Csca','Cbak','Cpr','<cos(theta)>'
	write(6,'(2x,6e13.5)') cext,cabs,csca,cbak,cext-cpr,assym
	cscax=0.d0
	cscay=0.d0
	cextx=0.d0
	cexty=0.d0
	cabsx=0.d0
	cabsy=0.d0
	cprx=0.d0
	cpry=0.d0
	do i=1,nL
	   cscax=cscax+cscaxi(i)
	   cscay=cscay+cscayi(i)
           cextx=cextx+cextxi(i)
           cabsx=cabsx+cabsxi(i)
           cexty=cexty+cextyi(i)
           cabsy=cabsy+cabsyi(i)
           cprx=cprx+cprxi(i)
           cpry=cpry+cpryi(i)
	enddo
	assymx=cprx/cscax
	assymy=cpry/cscay
	assym0=0.5d0*(assymx+assymy)
	cscax=2.0d0*twopi*cscax/cz
	cscay=2.0d0*twopi*cscay/cz
	csca=0.5d0*(cscax+cscay)
	cextx=twopi*cextx/cz
	cexty=twopi*cexty/cz
	cext=0.5d0*(cextx+cexty)
	cabsx=2.0d0*twopi*cabsx/cz
	cabsy=2.0d0*twopi*cabsy/cz
	cabs=0.5d0*(cabsx+cabsy)
        cprx=2.0d0*twopi*cprx/cz
        cpry=2.0d0*twopi*cpry/cz
        cpr=0.5d0*(cprx+cpry)
	assym=cpr/csca	
	write(6,'(2x,6e13.5)') cext,cabs,csca,cbak,cext-cpr,assym
	do i=1,nL
	   cabsxi(i)=4.0d0*pione*cabsxi(i)/cz
	   cabsyi(i)=4.0d0*pione*cabsyi(i)/cz
	   cextxi(i)=2.0d0*pione*cextxi(i)/cz
	   cextyi(i)=2.0d0*pione*cextyi(i)/cz
	   cscaxi(i)=4.0d0*pione*cscaxi(i)/cz
	   cscayi(i)=4.0d0*pione*cscayi(i)/cz
	   cprxi(i)=4.0d0*pione*cprxi(i)/cz
	   cpryi(i)=4.0d0*pione*cpryi(i)/cz
	   cscai(i)=0.5d0*(cscaxi(i)+cscayi(i))
	   cexti(i)=0.5d0*(cextxi(i)+cextyi(i))
	   cabsi(i)=0.5d0*(cabsxi(i)+cabsyi(i))
	   cpri(i)=0.5d0*(cprxi(i)+cpryi(i))
	   cpri(i)=cscai(i)+cabsi(i)-cpri(i)
	   assymxi(i)=cprxi(i)/cscaxi(i)
	   assymyi(i)=cpryi(i)/cscayi(i)
	   assymi(i)=0.5d0*(cprxi(i)+cpryi(i))/csca
	   cprxi(i)=cscaxi(i)+cabsxi(i)-cprxi(i)
	   cpryi(i)=cscayi(i)+cabsyi(i)-cpryi(i)
	   write(6,'(i5,5e15.6)') i,cexti(i),cabsi(i),
     +           cscai(i),cpri(i),assymi(i)
	enddo
	write(6,'(/,a)') 'efficiencies for radiation pressure'
	write(6,'(5x,6e13.5)') assym,assym0
	write(6,'(5x,6e13.5)') assym,cext-cpr,assymx,cextx-cprx,
     +                                        assymy,cexty-cpry
	do i=1,nL
	   write(6,'(i5,6e13.5)') i,assymi(i),cpri(i),
     +         assymxi(i),cprxi(i),assymyi(i),cpryi(i)  	   
	enddo
	betami=betami*90.d0/pih
	betamx=betamx*90.d0/pih
	thetmi=thetmi*90.d0/pih
	thetmx=thetmx*90.d0/pih
	phaimi=phaimi*90.d0/pih
	phaimx=phaimx*90.d0/pih
	flout='cr'//fileout
	open(33,file=flout,status='unknown')
	write(33,'(a20,a47)') 
     +     flout,'(Total and individual-particle cross sections)'
	write(33,'(a32,2x,a22)') 
     +     'input sphere-aggregate filename:',FLNAME
	write(33,'(a19,3i5)') 'nbeta,nthet,nphai: ',nbeta,nthet,nphai
	write(33,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +                betami,betamx,thetmi,thetmx,phaimi,phaimx
	write(33,'(a28,i5)') '# of orientations averaged: ',nram
	write(33,'(12x,a4,11x,a4,11x,a4,11x,a3,8x,a12)')
     +    'Cext','Cabs','Csca','Cpr','<cos(theta)>'
	write(33,'(a5,5e15.6)') 'total',cext,cabs,csca,cext-cpr,
     +     assym
	do i=1,nL
	   write(33,'(i5,5e15.6)') i,cexti(i),cabsi(i),cscai(i),
     +           cpri(i),assymi(i)
	enddo
	close(33)       
	cz=pione*gcvr*gcvr
	assym=cpr/csca
	assymx=cprx/cscax
	assymy=cpry/cscay
	cabsxv=cabsx/cz
	cabsyv=cabsy/cz
	cextxv=cextx/cz
	cextyv=cexty/cz
	cscaxv=cscax/cz
	cscayv=cscay/cz
	cprxv=cprx/cz
	cprxv=cextxv-cprxv
	cpryv=cpry/cz
	cpryv=cextyv-cpryv
	cscav=0.5d0*(cscaxv+cscayv)
	cextv=0.5d0*(cextxv+cextyv)
	cabsv=0.5d0*(cabsxv+cabsyv)
	cprv=0.5d0*(cprxv+cpryv)
	cbakxv=cbakx/cz
	cbakyv=cbaky/cz
	cbakv=0.5d0*(cbakxv+cbakyv)
	temp=gcvr*gcvr/gcs
	cabsxs=cabsxv*temp
	cabsys=cabsyv*temp
	cextxs=cextxv*temp
	cextys=cextyv*temp
	cscaxs=cscaxv*temp
	cscays=cscayv*temp
	cprxs=cprxv*temp
	cprys=cpryv*temp
	cscas=cscav*temp
	cexts=cextv*temp
	cabss=cabsv*temp
	cprs=cprv*temp
	cbakxs=cbakxv*temp
	cbakys=cbakyv*temp
	cbaks=cbakv*temp
 222    format(1x,a1,6e13.5)
 221    format(6x,a5,8x,a5,8x,a5,8x,a5,8x,a4,5x,a12)        
	write(6,221)
     +    'Qextv','Qabsv','Qscav','Qbakv','Qprv','<cos(theta)>'	
        write(6,222) 't',cextv,cabsv,cscav,cbakv,cprv,assym
	write(6,222) 'x',cextxv,cabsxv,cscaxv,cbakxv,cprxv,assymx
        write(6,222) 'y',cextyv,cabsyv,cscayv,cbakyv,cpryv,assymy
	write(6,221)
     +    'Qexts','Qabss','Qscas','Qbaks','Qprs','<cos(theta)>'
        write(6,222) 't',cexts,cabss,cscas,cbaks,cprs,assym
	write(6,222) 'x',cextxs,cabsxs,cscaxs,cbakxs,cprxs,assymx
        write(6,222) 'y',cextys,cabsys,cscays,cbakys,cprys,assymy
        temp=-(cabs+csca-cext)/cext
	write(6,'(/,a37,e14.5)') 
     +        'Accuracy of this numerical solution: ',temp
	open(12,file=fileout,status='unknown')
	write(12,'(a20,a16,a18,a4,f8.3,a5,f8.3)') fileout, 
     +        '--- input file: ',FLNAME,' xv:',xv,'  xs:',xs
	write(12,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +                betami,betamx,thetmi,thetmx,phaimi,phaimx        
	write(12,'(a19,3i5,6x,a28,i5)') 
     +        'nbeta,nthet,nphai: ',nbeta,nthet,nphai,
     +        '# of orientations averaged: ',nram
	write(12,221)
     +    'Cext','Cabs','Csca','Cbak','Cpr','<cos(theta)>'
	write(12,222) 
     +          't',cext,cabs,csca,cbak,cext-cpr,assym
	write(12,222)
     +          'x',cextx,cabsx,cscax,cbakx,cextx-cprx,assymx
	write(12,222)
     +          'y',cexty,cabsy,cscay,cbaky,cexty-cpry,assymy
	write(12,221)
     +    'Qextv','Qabsv','Qscav','Qbakv','Qprv','<cos(theta)>'
        write(12,222) 't',cextv,cabsv,cscav,cbakv,cprv,assym
	write(12,222) 'x',cextxv,cabsxv,cscaxv,cbakxv,cprxv,assymx
        write(12,222) 'y',cextyv,cabsyv,cscayv,cbakyv,cpryv,assymy
	write(12,221)
     +    'Qexts','Qabss','Qscas','Qbaks','Qprs','<cos(theta)>'
        write(12,222) 't',cexts,cabss,cscas,cbaks,cprs,assym
	write(12,222) 'x',cextxs,cabsxs,cscaxs,cbakxs,cprxs,assymx
        write(12,222) 'y',cextys,cabsys,cscays,cbakys,cprys,assymy
        write(12,'(1x,a4,4x,a7,5x,a4,7x,a3,10x,a3,10x,a3,10x,a3)')
     +    's.a.','i11+i22','pol.','i11','i21','i12','i22'               
	do i=1,nang2
	   write(12,'(f7.1,e13.5,f8.4,4e13.5)') 
     +        dang(i),inat(i),pol(i),i11(i),i21(i),i12(i),i22(i) 
	enddo
	close(12)
	open(11,file='mueller.out',status='unknown')
	write(11,'(a11,6x,a16)') 'mueller.out','(Mueller matrix)'
	write(11,'(a16,a20)') 'Input filename: ',FLNAME
	write(11,'(a19,3i5)') 'nbeta,nthet,nphai: ',nbeta,nthet,nphai
	write(11,'(a24,6f7.2)') 'Ranges of Euler angles: ',
     +                betami,betamx,thetmi,thetmx,phaimi,phaimx
	write(11,'(a28,i5)') '# of orientations averaged: ',nram
	do jc=1,npng
	   t=pang*dble(jc-1)
	   write(11,'(a18,f8.3)') 'phi (in degrees): ',t
	   do i=1,nang2
	      write(11,'(f7.1,4e16.7)') 
     +          dang(i),mue(1,1,jc,i),mue(1,2,jc,i),mue(1,3,jc,i),
     +                  mue(1,4,jc,i)
              write(11,'(7x,4e16.7)') 
     +                  mue(2,1,jc,i),mue(2,2,jc,i),mue(2,3,jc,i),
     +                  mue(2,4,jc,i)
              write(11,'(7x,4e16.7)') 
     +                  mue(3,1,jc,i),mue(3,2,jc,i),mue(3,3,jc,i),
     +                  mue(3,4,jc,i)
              write(11,'(7x,4e16.7)') 
     +                  mue(4,1,jc,i),mue(4,2,jc,i),mue(4,3,jc,i),
     +                  mue(4,4,jc,i)
	   enddo
	enddo
	write(11,*) 'phi (in degrees): 360'
	do i=1,nang2
	   write(11,'(f7.1,4e16.7)') 
     +        dang(i),mue(1,1,1,i),mue(1,2,1,i),mue(1,3,1,i),
     +                mue(1,4,1,i)
           write(11,'(7x,4e16.7)') 
     +                mue(2,1,1,i),mue(2,2,1,i),mue(2,3,1,i),
     +                mue(2,4,1,i)
              write(11,'(7x,4e16.7)') 
     +                mue(3,1,1,i),mue(3,2,1,i),mue(3,3,1,i),
     +                mue(3,4,1,i)
              write(11,'(7x,4e16.7)') 
     +                mue(4,1,1,i),mue(4,2,1,i),mue(4,3,1,i),
     +                mue(4,4,1,i)
	enddo
	close(11)
	STOP
	END

      SUBROUTINE orientcd(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,
     &                  MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI,
     &                  BETA,THETA,PHI)
C Arguments:
      INTEGER MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI
      double precision BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX
      double precision BETA(MXBETA),THETA(MXTHET),PHI(MXPHI)
C Local variables:
      INTEGER J
      double precision DELTA
C***********************************************************************
C Given: BETAMI=minimum value of beta (radians)
C        BETAMX=maximum value of beta (radians)
C        THETMI=minimum value of theta (radians)
C        THETMX=maximum value of theta (radians)
C        PHIMIN=minimum value of phi (radians)
C        PHIMAX=maximum value of phi (radians)
C        MXBETA,MXTHET,MXPHI=dimensions of the arrays BETA,THETA,PHI
C        NBETA=number of values of beta
C        NTHETA=number of values of theta
C        NPHI=number of values of PHI
C Returns: BETA(1-NBETA)=beta values (radians)
C          THETA(1-NTHETA)=theta values (radians)
C          PHI(1-NPHI)=phi values (radians)
C Purpose: to generate a sequence of desired target orientations
C Present version assumes:
C        beta to be uniformly distributed between BETAMI and BETAMX
C        cos(theta) to be uniformly distributed between cos(THETMI) and
C                   cos(THETMX)
C        phi to be uniformly distributed between PHIMIN and PHIMAX
C        If NPHI=1, first angle is THETMI, last angle is THETMX
C        If NPHI>1, then angles are midpoints of intervals of equal
C            range in theta subdividing range from THETMI to THETMX
C***********************************************************************
      BETA(1)=BETAMI            
      IF(NBETA.GT.1)THEN
         DELTA=(BETAMX-BETAMI)/DBLE(NBETA-1)
         DO 1000 J=2,NBETA
           BETA(J)=BETA(1)+DELTA*DBLE(J-1)
 1000    CONTINUE
      ENDIF
      IF(NPHI.EQ.1.AND.NTHETA.GT.1)THEN
         DELTA=(DCOS(THETMX)-DCOS(THETMI))/DBLE(NTHETA-1)
         THETA(1)=THETMI
      ELSE
         DELTA=(DCOS(THETMX)-DCOS(THETMI))/DBLE(NTHETA)
         THETA(1)=DACOS(DCOS(THETMI)+0.5d0*DELTA)
      ENDIF
      IF(NTHETA.GT.1)THEN
         DO 2000 J=2,NTHETA
            THETA(J)=DACOS(DCOS(THETA(1))+DELTA*DBLE(J-1))
 2000    CONTINUE
      ENDIF
c      DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI)
c      PHI(1)=PHIMIN+0.5D0*DELTA
      PHI(1)=PHIMIN
      IF(NPHI.GT.1)THEN
         DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI-1)
         DO 3000 J=2,NPHI
            PHI(J)=PHI(1)+DELTA*DBLE(J-1)
 3000    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE orientud(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,
     &                  MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI,
     &                  BETA,THETA,PHI)
C Arguments:
      INTEGER MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI
      double precision BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX
      double precision BETA(MXBETA),THETA(MXTHET),PHI(MXPHI)
C Local variables:
      INTEGER J
      double precision DELTA
C***********************************************************************
C Given: BETAMI=minimum value of beta (radians)
C        BETAMX=maximum value of beta (radians)
C        THETMI=minimum value of theta (radians)
C        THETMX=maximum value of theta (radians)
C        PHIMIN=minimum value of phi (radians)
C        PHIMAX=maximum value of phi (radians)
C        MXBETA,MXTHET,MXPHI=dimensions of the arrays BETA,THETA,PHI
C        NBETA=number of values of beta
C        NTHETA=number of values of theta
C        NPHI=number of values of PHI
C Returns:  BETA(1-NBETA)=beta values (radians)
C           THETA(1-NTHETA)=theta values (radians)
C           PHI(1-NPHI)=phi values (radians)
C Note: it is assumed that target orientation weight function
C       can be factored into WGTA*WGTB -- i.e., that rotations
C       around a1 are decoupled from orientation of a1.
C Purpose: to generate a sequence of desired target orientations
C Present version assumes:
C        beta to be uniformly distributed between BETAMI and BETAMX
C        theta to be uniformly distributed between THETMI and THETMX
C        phi to be uniformly distributed between PHIMIN and PHIMAX
C            first angle is THETMI, last angle is THETMX
C***********************************************************************
      BETA(1)=BETAMI
      IF(NBETA.GT.1)THEN
         DELTA=(BETAMX-BETAMI)/DBLE(NBETA-1)
         DO 1000 J=2,NBETA
           BETA(J)=BETA(1)+DELTA*DBLE(J-1)
 1000    CONTINUE
      ENDIF
      THETA(1)=THETMI
      IF(NTHETA.GT.1)THEN
         DELTA=(THETMX-THETMI)/DBLE(NTHETA-1)
         DO 2000 J=2,NTHETA
            THETA(J)=THETA(1)+DELTA*DBLE(J-1)
 2000    CONTINUE
      ENDIF
      PHI(1)=PHIMIN
      IF(NPHI.GT.1)THEN
         DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI-1)
         DO 3000 J=2,NPHI
            PHI(J)=PHI(1)+DELTA*DBLE(J-1)
 3000    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE abMiexud(X,REFREL,NP,NMAX,NM,AN,BN,NADD,
     +                    RSR,RSI,RSX,PX,AR,AI,BR,BI,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NMAX,NM,NADD,NSTOP,NX,K,N
      COMPLEX*16 REFREL,AN(NP),BN(NP)
      DOUBLE PRECISION AR(NP),AI(NP),BR(NP),BI(NP),
     +   RSR(NMAX),RSI(NMAX),RSX(NMAX),PX(NMAX)
      IF(EPS.GT.1.D0.OR.EPS.LT.0.D0) EPS=1.0D-20
      IF(NADD.NE.0) EPS=0.D0
      CTC=EPS
      XM=DBLE(REFREL)
      YM=DIMAG(REFREL)
      XMX=X*XM
      YMX=X*YM
      RP2=XMX*XMX+YMX*YMX
      NSTOP=X+4.D0*X**.3333D0
      NSTOP=NSTOP+2
      NM=NSTOP+NADD
      XN=DSQRT(XM**2+YM**2)*X
      NX=1.1D0*XN+10.D0
      if(NX-NM.lt.10) NX=NM+10
      write(6,*) 'Wiscombe criterion:',NSTOP 
      write(6,*) 'NADD:',NADD
      write(6,*) 'NX:',NX
      IF(NX.GT.NMAX) THEN
         WRITE(6,*) 'Parameter NXMAX too small'
         WRITE(6,*) '  NXMAX must be greater than', NX
         WRITE(6,*) 'Please correct NXMAX in main code,' 
         WRITE(6,*) '  recompile, then try again'
         STOP
      ENDIF
      IF(NM.GT.NP) THEN
         WRITE(6,*) 'Parameter np too small'
         WRITE(6,*) '  np must be greater than', NM
         WRITE(6,*) 'Please correct np in gmm01f.par,' 
         WRITE(6,*) '  recompile the code, then try again'
         STOP
      ENDIF
C  DOWN RECURSION FOR RATIOS RSR,RSI,RSX,PNR,PNI,PX
      PNX=X/DBLE(2*NX+3)
      PNR=XMX/DBLE(2*NX+3)
      PNI=YMX/DBLE(2*NX+3)
      DO 5 K=1,NX
         N=NX-K+1
         CN=DBLE(N)
         ALN=(2.D0*CN+1.D0)*XMX/RP2-PNR
         BEN=(2.D0*CN+1.D0)*YMX/RP2+PNI
         RSR(N)=-CN*XMX/RP2+ALN
         RSI(N)=CN*YMX/RP2-BEN
         PZD=ALN*ALN+BEN*BEN
         PNR=ALN/PZD
         PNI=BEN/PZD
         RSX(N)=(CN+1.D0)/X-PNX
         IF(N.EQ.1) GO TO 20
         PNX=X/(2.D0*CN+1.D0-PNX*X)
         PX(N)=PNX
    5    CONTINUE
   20 SNM1X=DSIN(X)
      CNM1X=DCOS(X)
      IF(X-0.1D0) 21,22,22
   21 SNX=X**2./3.D0-X**4./30.D0+X**6./840.D0-X**8./45360.D0
      GO TO 23
   22 SNX=SNM1X/X-CNM1X
   23 CNX=CNM1X/X+SNM1X
      DO 10 N=1,NX
         PX(N)=SNX      !preparing for the calculation of Cabs
         C=DBLE(N)
         DCNX=CNM1X-C*CNX/X
         DSNX=RSX(N)*SNX
C  CALCULATION OF EXTERIOR COEFFICIENTS AN AND BN
         ANNR=RSR(N)*SNX-XM*DSNX
         ANNI=RSI(N)*SNX-YM*DSNX
         TA1=RSR(N)*SNX-RSI(N)*CNX
         TA2=RSI(N)*SNX+RSR(N)*CNX
         ANDR=TA1-XM*DSNX+YM*DCNX
         ANDI=TA2-XM*DCNX-YM*DSNX
         AND=ANDR*ANDR+ANDI*ANDI
         BNNR=(XM*RSR(N)-YM*RSI(N))*SNX-DSNX
         BNNI=(XM*RSI(N)+YM*RSR(N))*SNX
         TB1=RSR(N)*SNX-RSI(N)*CNX
         TB2=RSR(N)*CNX+RSI(N)*SNX
         BNDR=XM*TB1-YM*TB2-DSNX
         BNDI=XM*TB2+YM*TB1-DCNX
         BND=BNDR*BNDR+BNDI*BNDI
         AR(N)=(ANNR*ANDR+ANNI*ANDI)/AND
         AI(N)=(ANNI*ANDR-ANNR*ANDI)/AND
         BR(N)=(BNNR*BNDR+BNNI*BNDI)/BND
         BI(N)=(BNNI*BNDR-BNNR*BNDI)/BND
C  MIE SERIES CONVERGENCE TEST IS MADE BY TESTING AN'S AND BN'S
         TI=AR(N)*AR(N)+AI(N)*AI(N)+BR(N)*BR(N)+BI(N)*BI(N)
         TI=TI/(AR(1)*AR(1)+AI(1)*AI(1)+BR(1)*BR(1)+BI(1)*BI(1))
         IF(TI-CTC) 16,18,18
   18    IF(NM-N) 15,15,6
    6    IF(N-NX)7,15,15
    7    M=N+1
         SNX=PX(M)*SNX
         CNM2X=CNM1X
         CNM1X=CNX
         CNX=(2.D0*C+1.D0)*CNM1X/X-CNM2X
   10 CONTINUE
      GO TO 15
   16 WRITE(6,*) '*** NOTE THAT THE FIELD-EXPANSION TRANCATION'
      WRITE(6,*) '*** IS DETERMINED BY eps GIVEN IN THE INPUT'
      WRITE(6,*) '*** FILE gmm01f.in'
      WRITE(6,*) '*** IN CASE YOU NEED A HIGHER ORDER, eps MUST'
      WRITE(6,'(a,e9.1)')
     +          ' *** BE SMALLER THAN THE CURRENT VALUE',EPS
   15 NM=N
      DO I=1,NM
         AN(I)=DCMPLX(AR(I),-AI(I))
         BN(I)=DCMPLX(BR(I),-BI(I))
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION ran1d(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,
     *IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,
     *EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1d=dmin1(AM*iy,RNMX)
      return
      END

	subroutine cofsrd(nmax)
	include 'gmm01f.par'
	parameter (nmp=np*(np+2))	
	double precision cofsr(nmp),lnfacd,c
	common/crot/cofsr
        i=0
	do n=1,nmax
	   do m=-n,n
	      i=i+1       
	      c=lnfacd(dble(n-m))-lnfacd(dble(n+m))
	      cofsr(i)=0.5d0*c
c	      c=0.5d0*c
c	      cofsr(i)=dexp(c)
           enddo
        enddo
	return
	end

	subroutine cofd0(nmax)
	implicit double precision (a-h,o-z)	
	include 'gmm01f.par'
	parameter (nmp=np*(np+2))	
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)	
	integer v
	double precision lnfacd
	common/cofmnv0/cof0(ni0)
	common/crot/cofsr(nmp)	
        common/fnr/fnr(0:2*(np+2))
	i=0
	sm=-0.5d0*dble((-1)**nmax)
	do m=-nmax,nmax
	   ns=max(1,iabs(m))
	   sm=-sm
	   do n=ns,nmax
	      inm=n*(n+1)-m
	      do v=ns,nmax
	         i=i+1
	         ivm=v*(v+1)+m
	         c=cofsr(inm)+cofsr(ivm)	            
	         c=sm*dexp(c)
	         c0=fnr(2*n+1)*fnr(2*v+1)
	         c1=fnr(n)*fnr(v)*fnr(n+1)*fnr(v+1)
	         c0=c0/c1
                 cof0(i)=c*c0
              enddo
           enddo
        enddo
	return
	end

	subroutine cofnv0(nmax)
	include 'gmm01f.par'
	integer n,v
	double precision c1,lnfacd,cnv(np,np)
	common/cfnv/cnv
	do n=1,nmax
	   do v=n,nmax
	      c1=lnfacd(dble(2*n))+lnfacd(dble(2*v))
	      c1=c1-lnfacd(dble(2*n+2*v))
	      c1=c1+2.d0*lnfacd(dble(n+v))
	      c1=c1-lnfacd(dble(n))-lnfacd(dble(v))
	      cnv(n,v)=c1
	   enddo
	enddo
	return
	end

C  subroutine gau0.f generates tabulated values for  
C  Gaunt coefficients up to n=v=n_max
	subroutine gau0(nmax)
	include 'gmm01f.par'
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
	integer v,qmax,uvmax,iga0(ni0)
	double precision ga0(ng0)
	common/g0/ga0
	common/ig0/iga0
        na=0
        uvmax=nmax*(nmax+2)
        i=0
	do m=-nmax,nmax
	   ns=max(1,iabs(m))
	   do n=ns,nmax
	      do v=ns,nmax
	         call gxurcd0(-m,n,v,qmax,na)
	         i=i+1
                 iga0(i)=na
	         na=na+qmax+1
              enddo
	   enddo
        enddo
	return
	end

c  transforms the rectangular coordinates (x,y,z)
c  to spherical coordinates (r,theta,phi)
	subroutine carsphd(x,y,z,r,xt,sphi,cphi)
	double precision x,y,z,r,xt,sphi,cphi
	r=dsqrt(x*x+y*y+z*z)
	if(r.eq.0.d0) then
	   xt=1.d0
	   sphi=0.d0
	   cphi=1.d0
	   return
	endif
	xt=z/r
	if(y.eq.0.d0.and.x.eq.0.d0) then
	   sphi=0.d0
	   cphi=1.d0
	   return
	endif
	sphi=dsqrt(x*x+y*y)
	cphi=x/sphi
	sphi=y/sphi
	return
	end

c  subroutine besseljd.f  (in double precision arithmetic)
c  returns an array of the spherical Bessel function of the  
c  first kind with a real argument: j_0,j_1,j_2,...,j_{NC}
c  uses Ru Wang's ratio method for the downward recursive 
c  calculation of the Riccati-Bessel function Psi_n(z)=z j_n(z) 
c  [see Xu et al., Physical Review E, v.60, 2347-2365 (1999)]
      SUBROUTINE besseljd(NC,X,BESJ)
      INTEGER NC,NX,K,N
      DOUBLE PRECISION X,BESJ(0:NC),PN,CN,X2
      DO K=1,NC
         BESJ(K)=0.D0
      ENDDO
      IF(DABS(X).LT.1.D-12) THEN
         BESJ(0)=1.D0
         BESJ(1)=1.D0/3.D0*X
         RETURN
      ENDIF
c  down-recursively calculating an array of the ratio functions  
c  P_{NC},P_{NC-1},...,P(2) stored in the same array for j_n, 
c  starting with an asymptotic value P_{NX+1}=X/(2 NX+3) at the 
c  highest order NX+1, where NX=NC+1.1X+10 
      NX=1.1D0*X+10.D0
      NX=NC+NX   
      PN=X/DBLE(2*NX+3)
      DO 5 K=1,NX-1
         N=NX-K+1
         CN=DBLE(N)
         PN=X/(DBLE(2*N+1)-PN*X)
         IF(N.GT.NC) GOTO 5
         BESJ(N)=PN
    5 CONTINUE
C  calculating j_0(x) and j_1(x)
      IF(DABS(X)-0.1D0) 10,11,11
   10 X2=X*X
      BESJ(0)=1.D0-X2/72.D0
      BESJ(0)=1.D0-X2/42.D0*BESJ(0)
      BESJ(0)=1.D0-X2/20.D0*BESJ(0)
      BESJ(0)=1.D0-X2/6.D0*BESJ(0)
      BESJ(1)=1.D0/45360.D0-X2/3991680.D0
      BESJ(1)=1.D0/840.D0-X2*BESJ(1)
      BESJ(1)=1.D0/30.0d0-X2*BESJ(1)
      BESJ(1)=X*(1.D0/3.0d0-X2*BESJ(1))
      GOTO 12
   11 BESJ(0)=DSIN(X)/X
      BESJ(1)=(DSIN(X)/X-DCOS(X))/X      
c  calculating j_2,...,j_{NC} 
   12 DO 20 N=2,NC
         BESJ(N)=BESJ(N)*BESJ(N-1)
   20 CONTINUE
      RETURN
      END

c  sub. besselyd.f  (in double precision arithmetic)
c  returns an array of the spherical Bessel function of
c  the second kind with a real argument: y_0,y_1,...,y_n
       subroutine besselyd(n,x,besy)
       integer i,n
       double precision x,besy(0:n),besyn,x2
       if(x.eq.0.d0) then
         write(6,*) 'bad argument in sub. besselyd'
         stop
       endif
       if(dabs(x)-0.1d0)10,11,11
  10   x2=x*x
       besyn=1.d0-x2/72.d0
       besyn=1.d0-x2/42.d0*besyn
       besyn=1.d0-x2/20.d0*besyn
       besyn=1.d0-x2/6.d0*besyn
       besy(0)=1.d0-x2/56.d0
       besy(0)=1.d0-x2/30.d0*besy(0)
       besy(0)=1.d0-x2/12.d0*besy(0)
       besy(0)=1.d0-0.5d0*x2*besy(0)
       besy(0)=-besy(0)/x
       goto 12
  11   besyn=dsin(x)/x
       besy(0)=-dcos(x)/x
  12   besy(1)=besy(0)/x-besyn
       do i=2,n
         besy(i)=dble(2*i-1)/x*besy(i-1)-besy(i-2)
       enddo
       return
       end

c
c  This subroutine is originally written by D.W. Mackowski   
c  (taken from Mackowski's multisphere-scattering code scsmfo1b.for    
c  released to public by the author at 8/1999).     
c  Very slightly modified for fitting into this code.           
c  Yu-lin Xu   12/2000
c
      subroutine rotcoef(cbe,nmax)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2))
      implicit double precision (a-h,o-z)
      double precision dk0(-2*np:2*np),dk01(-2*np:2*np)
      common/rot/bcof(0:np+2),dc(-np:np,0:nmp)
      common/fnr/fnr(0:2*(np+2))
      sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
      cbe2=.5d0*(1.d0+cbe)
      sbe2=.5d0*(1.d0-cbe)
      in=1
      dk0(0)=1.d0
      sben=1.d0
      dc(0,0)=1.d0
      dk01(0)=0.d0
      do n=1,nmax
         nn1=n*(n+1)
         in=-in
         sben=sben*sbe/2.d0
         dk0(n)=dble(in)*sben*bcof(n)
         dk0(-n)=dble(in)*dk0(n)
         dk01(n)=0.d0
         dk01(-n)=0.d0
         dc(0,nn1+n)=dk0(n)
         dc(0,nn1-n)=dk0(-n)
         do k=-n+1,n-1
            kn=nn1+k
            dkt=dk01(k)
            dk01(k)=dk0(k)
            dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)
     1             /(fnr(n+k)*fnr(n-k))
            dc(0,kn)=dk0(k)
         enddo
         im=1
         do m=1,n
            im=-im
            fmn=1.d0/fnr(n-m+1)/fnr(n+m)
            m1=m-1
            dkm0=0.d0
            do k=-n,n
               kn=nn1+k
               dkm1=dkm0
               dkm0=dc(m1,kn)
               if(k.eq.n) then
                  dkn1=0.d0
               else
                  dkn1=dc(m1,kn+1)
               endif
               dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1
     1           -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1
     1              -dble(k)*sbe*dc(m1,kn))*fmn
               dc(-m,nn1-k)=dble((-1)**(k)*im)*dc(m,kn)
            enddo
         enddo
      enddo
      return
      end

c  subroutine cofxuds0.f returns the two classes of vector 
c  (axial) translation coefficients for a given combination of 
c  (m,n,m,v) and a given dimensionless translation distance kd 
cu uses subroutine gid0.f 
	subroutine cofxuds0(nmax,m,n,v,sja,sya,A,B,Aj,Bj)
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
	integer v,p,qmax
	double precision sja(0:n+v+1),sya(0:n+v+1)
	complex*16 A,B,Aj,Bj,signz
	common/ig0/iga0(ni0)
	common/g0/ga0(ng0)
	common/cofmnv0/cof0(ni0)
	fa(m,p)=dble(-2*m*p*(p-1))
	fb(n,v,p)=dble(p*p-(n+v+1)*(n+v+1))
     +            *dble(p*p-(n-v)*(n-v))/dble(4*p*p-1)
	if(iabs(m).gt.n.or.iabs(m).gt.v) then 
	   write(6,*) '|m|>n or v in subroutine cofxuds0.f'
	   stop
	endif
	A=0.d0
	B=0.d0
	Aj=0.d0
	Bj=0.d0
	call gid0(nmax,m,n,v,id)
	c=cof0(id)
	ig=iga0(id)
	nv2=v*(v+1)+n*(n+1)
	signz=dcmplx(0.d0,1.d0)**(n+v)
	p=n+v+2
	qmax=min(n,v)
	do i=1,qmax+1
	   p=p-2
	   cp=dble(nv2-p*(p+1))*ga0(ig+i)
	   sj=sja(p)
	   sy=sya(p)
	   A=A+dcmplx(sj,sy)*signz*cp
	   Aj=Aj+sj*signz*cp
	   signz=-signz
	enddo
	A=A*c
	Aj=Aj*c
	if(m.eq.0) return
 	signz=dcmplx(0.d0,1.d0)**(n+v+1)
 	p=n+v
	do 20 i=1,qmax
	   p=p-2
	   signz=-signz
	   if(i.eq.1) then
              cp=dble(2*p+3)*fa(m,p+3)      
              cp=cp*ga0(ig+1)/dble((p+3)*(p+2))
              goto 21
           endif
           if(i.eq.qmax) then
              if(p.eq.0) goto 22
              nv2=p*(p+1)
              cp=dble(2*p+3)*fa(m,p+1)
              cp=-cp*ga0(ig+i+1)/dble(nv2)
              goto 21
           endif
 22	   c4=fa(m,p+2)
           cp=-dble((p+1)*(p+2))*fb(n,v,p+2)*ga0(ig+i)
	   cp=cp+dble((p+2)*(p+1))*fb(n,v,p+1)*ga0(ig+i+1)
 	   cp=cp*dble(2*p+3)/c4
 21	   sj=sja(p+1)
	   sy=sya(p+1)
 	   B=B+dcmplx(sj,sy)*signz*cp
 	   Bj=Bj+sj*signz*cp
 20	continue
	B=B*c
	Bj=Bj*c
	return
	end

	subroutine trans(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,
     +                   drot,as,bs,as1,bs1,ind)
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
	parameter (nmp=np*(np+2))
	parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)	
	parameter (nrc=4*np*(np+1)*(np+2)/3+np)
	parameter (nij=nLp*(nLp-1)/2)
	integer v,nmax(nLp),uvmax(nLp),ind(nL)
	double precision r0(6,nLp),drot(nrc,nij)
	complex*16 atr(2,np,nmp),atr0(ni0,nij),btr0(ni0,nij),
     +     ek(np,nij),as(nLp,nmp),bs(nLp,nmp),as1(nLp,nmp),
     +     bs1(nLp,nmp),at1(2,nmp),bt1(2,nmp)
        common/tran/atr
 	do i=1,nL
 	   if(ind(i).gt.0) goto 11
           do imn=1,uvmax(i)
              as1(i,imn)=dcmplx(0.d0,0.d0)
	      bs1(i,imn)=dcmplx(0.d0,0.d0)
	   enddo
	   do 10 j=1,nL
	      if(j.eq.i) goto 10
	      x0=r0(1,i)-r0(1,j)
	      y0=r0(2,i)-r0(2,j)
	      z0=r0(3,i)-r0(3,j) 
	      d=dsqrt(x0*x0+y0*y0+z0*z0)
	      temp=(r0(4,i)+r0(4,j))/d
	      if(temp.le.fint) goto 10
	      if(i.lt.j) then
                 ij=(j-1)*(j-2)/2+j-i
              else
                 ij=(i-1)*(i-2)/2+i-j
              endif
              nlarge=max(nmax(i),nmax(j))
              itrc=0
              nsmall=min(nmax(i),nmax(j))
              do m=-nsmall,nsmall
                 n1=max(1,iabs(m))
                 do n=n1,nlarge
                    do v=n1,nlarge
                       itrc=itrc+1
                       iuv=v*(v+1)+m
                       atr(1,n,iuv)=atr0(itrc,ij)
                       atr(2,n,iuv)=btr0(itrc,ij)
                       if(x0.eq.0.d0.and.y0.eq.0.d0) then
                          if(z0.lt.0.d0) goto 20
                       else
                          if(j.lt.i) goto 20
                       endif
                       goto 21
 20                    sic=dble((-1)**(n+v))
                       atr(1,n,iuv)=sic*atr(1,n,iuv)
                       atr(2,n,iuv)=-sic*atr(2,n,iuv)
 21                    continue
                    enddo
                 enddo
              enddo
              do iuv=1,uvmax(j)
                 at1(1,iuv)=as(j,iuv)
                 at1(2,iuv)=bs(j,iuv)
              enddo
              if(x0.eq.0.d0.and.y0.eq.0.d0) then
                 call trv(at1,nmax(j),nmax(i))
              else 
                 call rtr(at1,nmax(j),nmax(i),ek(1,ij),drot(1,ij))
              endif
     	      do imn=1,uvmax(i)
     	         as1(i,imn)=as1(i,imn)+at1(1,imn)
     	         bs1(i,imn)=bs1(i,imn)+at1(2,imn)
     	      enddo
 10	   continue
 11	   continue
	enddo 
	return
	end


c  subroutine tipitaud.f 
c  calculates pi(m,n) & tau(m,n) up to a specified nmax for all 
c  m=0,1,...n at a given x
c  pi(m,n) & tau(m,n) calculated are normalized by  
c         C_mn=[(2n+1)(n-m)!/n/(n+1)/(n+m)!]^(1/2)
c  Yu-lin Xu    12/2000

	subroutine tipitaud(nmax,x)
	include 'gmm01f.par'
	parameter (nmp0=(np+1)*(np+4)/2)
	implicit double precision (a-h,o-z)
	common/fnr/fnr(0:2*(np+2))
	common/pitau/pi(nmp0),tau(nmp0)	

	nt=(nmax+1)*(nmax+4)/2          ! calculates pi up to nmax+1
	if(nt.gt.nmp0.or.dabs(x).gt.1.d0) then
	   write(6,*) 'dimension or argument wrong in sub. tipitaud'
	   write(6,*) 'argument: ',x
	   stop
	endif
	sx=dsqrt(1.d0-x*x)
	pi(1)=0.d0                     ! pi(0,1)  pi(0,n)=0 when m=0
	pi(2)=dsqrt(.75d0)             ! pi(1,1)
	pi(3)=0.d0                     ! pi(0,2)
	t125=dsqrt(1.25d0) 
	pi(4)=t125*x                      ! pi(1,2)  
	pi(5)=t125*sx                     ! pi(2,2)
	imn=5
	do i=3,nmax+1
	   n=i
	   imn=imn+1
	   pi(imn)=0.d0                ! pi(0,n)=0
	   do 11 j=2,n
	      m=j-1
	      imn=imn+1
	      i1=(n-2)*(n+1)/2+m+1
	      if(m.eq.n-1) then
	         pi(imn)=fnr(n-1)*fnr(2*n+1)/fnr(n+1)*x*pi(i1)
	         goto 11
	      endif
	      i2=(n-3)*n/2+m+1
	      t=fnr(n)*fnr(2*n-3)
	      t=fnr(n-2)*fnr(n-m-1)*fnr(n+m-1)/t
	      pi(imn)=fnr(2*n-1)*x*pi(i1)-t*pi(i2)
	      t=fnr(n+1)*fnr(n-m)*fnr(n+m)
	      t=fnr(n-1)*fnr(2*n+1)/t
	      pi(imn)=t*pi(imn)
 11	   continue
	   imn=imn+1
	   i1=(n-2)*(n+1)/2+n
	   t=fnr(n-1)*fnr(n+1)
	   t=dsqrt(.5d0)*fnr(n)*fnr(2*n+1)/t
	   pi(imn)=t*sx*pi(i1)
	enddo
	tx=x*sx
	tau(1)=-dsqrt(1.5d0)*sx          ! tau(0,1)
	tau(2)=pi(2)*x                   ! tau(1,1)
	tau(3)=-dsqrt(7.5d0)*tx          ! tau(0,2)   
	tau(4)=t125*(2.d0*x*x-1.d0)      ! tau(1,2)
	tau(5)=t125*tx                   ! tau(2,2)
	imn=5
	do i=3,nmax
	   n=i
	   do 21 j=1,n+1
	      m=j-1
	      imn=imn+1
	      if(m.eq.0) then
	         i1=(n-2)*(n+1)/2+1
	         i2=(n-3)*n/2+1
	         t=fnr(2*n-3)
	         t=fnr(n-2)*fnr(n)/t
	         tau(imn)=fnr(2*n-1)*x*tau(i1)-t*tau(i2)
	         t=fnr(n-1)*fnr(n+1)
	         t=fnr(2*n+1)/t
	         tau(imn)=t*tau(imn)
	         goto 21
	      endif
	      i1=n*(n+3)/2+m+1
	      t=fnr(n)*fnr(2*n+3)
	      t=fnr(n+2)*fnr(2*n+1)*fnr(n-m+1)*fnr(n+m+1)/t
	      tau(imn)=t*pi(i1)-dble(n+1)*x*pi(imn)
	      tau(imn)=tau(imn)/dble(m)
 21	   continue
	enddo
	return
	end

c  subroutine mueller.f
c  returns the values of 16 elements of the 4x4 
c  Mueller matrix calculated from the known 2x2 
c  amplitude scattering matrix 
c  using Bohren & Huffman's formulas (p.65)
	subroutine mueller(s1,s2,s3,s4,s)
	double precision s(4,4),s1s,s2s,s3s,s4s
	complex*16 s1,s2,s3,s4
	complex*16 s2s3c,s1s4c,s2s4c,s1s3c,s1s2c
	complex*16 s3s4c,s2s1c,s4s3c,s2cs4,s3cs1
	s1s=cdabs(s1)**2
	s2s=cdabs(s2)**2
	s3s=cdabs(s3)**2
	s4s=cdabs(s4)**2
	s2s3c=s2*dconjg(s3)
	s1s4c=s1*dconjg(s4)
	s2s4c=s2*dconjg(s4)
	s1s3c=s1*dconjg(s3)
	s1s2c=s1*dconjg(s2)
	s3s4c=s3*dconjg(s4)
	s2s1c=s2*dconjg(s1)
	s4s3c=s4*dconjg(s3)
	s2cs4=dconjg(s2)*s4
	s3cs1=dconjg(s3)*s1
	s(1,1)=0.5d0*(s1s+s2s+s3s+s4s)
	s(1,2)=0.5d0*(s2s-s1s+s4s-s3s)
	s(1,3)=s2s3c+s1s4c
	s(1,4)=dimag(s2s3c-s1s4c)
	s(2,1)=0.5d0*(s2s-s1s-s4s+s3s)
	s(2,2)=0.5d0*(s2s+s1s-s4s-s3s)
	s(2,3)=s2s3c-s1s4c
	s(2,4)=dimag(s2s3c+s1s4c)
	s(3,1)=s2s4c+s1s3c
	s(3,2)=s2s4c-s1s3c
	s(3,3)=s1s2c+s3s4c
	s(3,4)=dimag(s2s1c+s4s3c)
	s(4,1)=dimag(s2cs4+s3cs1)
	s(4,2)=dimag(s2cs4-s3cs1)
	s(4,3)=dimag(s1s2c-s3s4c)
	s(4,4)=s1s2c-s3s4c
	return
	end

c  lnfacd.f  (double precision function)
c  returns ln(z!)  z>-1.0
c  based on Lanczos' method [see Xu, Journal of Computational 
c  Physics, v.139, 137-165 (1998)]
      double precision function lnfacd(z)
      integer i
      double precision z,a,b,cp,c0(11)
      data c0/0.16427423239836267d5, -0.48589401600331902d5,
     +        0.55557391003815523d5, -0.30964901015912058d5,
     +        0.87287202992571788d4, -0.11714474574532352d4,
     +        0.63103078123601037d2, -0.93060589791758878d0,
     +        0.13919002438227877d-2,-0.45006835613027859d-8,
     +        0.13069587914063262d-9/ 
      a=1.d0
      cp=2.5066282746310005d0
      b=z+10.5d0
      b=(z+0.5d0)*dlog(b)-b
      do i=1,11
        z=z+1.d0
        a=a+c0(i)/z
      enddo
      lnfacd=b+dlog(cp*a)
      return
      end
 
c  gxurcd0.f to compute Gaunt coefficients a(-m,n,m,v,p)
cu uses lnfacd.f to compute ln(z!)
	subroutine gxurcd0(m,n,v,qmax,na)
	implicit double precision (a-h,o-z)
	include 'gmm01f.par'
	parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
	integer v,qmax,p
	double precision lnfacd,cnv(np,np),ga0(ng0)
	common/cfnv/cnv
	common/g0/ga0
	fb(n,v,p)=dble(p-(n+v+1))*dble(p+(n+v+1))
     +           *dble(p-(n-v))*dble(p+(n-v))
     +           /(dble(2*p+1)*dble(2*p-1))	
	if(iabs(m).gt.n.or.iabs(m).gt.v) then
	   write(6,*) 'warning: |m|>n or v in gxurcd0'
	   qmax=-1	   
	   return
	endif
	qmax=min(n,v)
	nq=qmax+1
	if(n.le.v) then 
	   c1=cnv(n,v)
	else
	   c1=cnv(v,n)
	endif
	c1=c1-lnfacd(dble(n-m))-lnfacd(dble(v+m))
	ga0(na+1)=dexp(c1)
	if(qmax.lt.1) return	
 	p=n+v
	do 8 i=2,nq
	   p=p-2
	   if(m.eq.0) then
	      c1=fb(n,v,p+1)
	      c2=fb(n,v,p+2)
	      goto 2
	   endif
	   c1=fb(n,v,p+1)
	   c2=dble(4*m*m)+fb(n,v,p+2)+fb(n,v,p+3)
	   if(i.eq.2) goto 2
	   c3=-fb(n,v,p+4)
	   goto 4
  2	   ga0(na+i)=c2*ga0(na+i-1)/c1
	   goto 8
  4	   ga0(na+i)=(c2*ga0(na+i-1)+c3*ga0(na+i-2))/c1
  8     continue
	return
	end

	subroutine gid0(nmax,m,n,iv,id)
	nt=nmax*(nmax+1)*(2*nmax+1)/3+nmax*nmax
	ns=max(1,iabs(m))
 	nc0=nmax-iabs(m)
 	id=nc0*(nc0+1)*(2*nc0+1)/6
	if(m) 10,11,12
 10	id=id+(n-ns)*(nc0+1)+iv-ns+1  
	return
 11	id=id+(n-ns)*nmax+iv
	return
 12	id=id+(nmax-n)*(nc0+1)+nmax-iv
	id=nt-id
	return
	end

c
c  This subroutine is the implementation of Mackowski's three-step  
c  method for decomposition of translation matrix into rotational   
c  and axial translational parts. 
c  The rotation-translation-rotation formulation is described 
c  in Mackowski, Proc. R. Soc. Lond. A 433, pp.611-612 (1991).
c  Yu-lin Xu   12/2000
c
      subroutine rtr(anpt,nodrj,nodri,ekt,drot)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter(nmp=np*(np+2))
      double precision drot(*)
      complex*16 anpt(2,*),ant(2,2*np),amt(2,-np:np),a,b,
     +           ekt(*),ek(-np:np),atr(2,np,nmp)
      common/tran/atr
c
      ek(0)=1.d0
      nmax=max(nodrj,nodri)
      do m=1,nmax
         ek(m)=ekt(m)
         ek(-m)=dconjg(ek(m))
      enddo
      irc=0
      do n=1,nodrj
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.d0
            amt(2,m)=0.d0
         enddo
         do k=-n,n
            kn=n1+k
            a=ek(k)*anpt(1,kn)
            b=ek(k)*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
         enddo
         do m=-n,n
            imn=n1+m
            anpt(1,imn)=amt(1,m)
            anpt(2,imn)=amt(2,m)
         enddo
      enddo
      mmax=min(nodrj,nodri)
      do m=-mmax,mmax
         n1=max(1,iabs(m))
         do n=n1,nodrj
            imn=n*(n+1)+m
            do ip=1,2
               ant(ip,n)=anpt(ip,imn)
            enddo
         enddo
         do n=n1,nodri
            imn=n*(n+1)+m
            a=0.d0
            b=0.d0
            do l=n1,nodrj
               ml=l*(l+1)+m
               a=a+atr(1,n,ml)*ant(1,l)
     1            +atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)
     1            +atr(2,n,ml)*ant(1,l)
            enddo
            anpt(1,imn) = a
            anpt(2,imn) = b
         enddo
      enddo
      in=1
      irc=0
      do n=1,nodri
         in=-in
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.d0
            amt(2,m)=0.d0
         enddo
         sik=-in
         do k=-n,n
	    sik=-sik
	    kn=n1+k
	    a=sik*anpt(1,kn)
            b=sik*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
         enddo
         sik=-in
         do m=-n,n
            sik=-sik
            imn=n1+m
            anpt(1,imn)=amt(1,m)*ek(-m)*sik
            anpt(2,imn)=amt(2,m)*ek(-m)*sik
         enddo
      enddo
      return
      end

      subroutine trv(anpt,nodrj,nodri)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter(nmp=np*(np+2))
      complex*16 anpt(2,*),ant(2,2*np),a,b,atr(2,np,nmp)
      common/tran/atr
      mmax=min(nodrj,nodri)
      do m=-mmax,mmax
         n1=max(1,iabs(m))
         do n=n1,nodrj
            imn=n*(n+1)+m
            do ip=1,2
               ant(ip,n)=anpt(ip,imn)
            enddo
         enddo
         do n=n1,nodri
            imn=n*(n+1)+m
            a=0.d0
            b=0.d0
            do l=n1,nodrj
               ml=l*(l+1)+m
               a=a+atr(1,n,ml)*ant(1,l)
     1            +atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)
     1            +atr(2,n,ml)*ant(1,l)
            enddo
            anpt(1,imn) = a
            anpt(2,imn) = b
         enddo
      enddo
      return
      end
