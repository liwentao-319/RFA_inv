      subroutine kerf(a,b,d,h,mmax,cg,ags,pr,dt,itim,u0,w0,u1,w1,data)
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of langston (1979)
c  *************************************************************************
c
	include 'apinv.inc'
c      common /rfinp/ iow,blen,agauss,c,tdelay
c	common /iput/  dt,ta,pr,npts,nft,nfpts,fny,delf
c	common /oput/  iowsy,iowrf
c      COMMON /COEF1/ cofoo,cofos,cofss
c	COMMON /coef2/ nl,alfm(mlay),betm(mlay),vpvsm(mlay),qpm(mlay),
c     *      qsm(mlay),rhom(mlay),thikm(mlay),dum1(mlay),dum2(mlay)
c      common /synth/ u0,w0
	real*4 sr(MAXP),sz(MAXP)
	real*4 a(mlay),b(mlay),d(mlay),h(mlay)
	real*4 alfm(mlay),betm(mlay),vpvsm(mlay),rhom(mlay),thikm(mlay)
      dimension d2(MAXP)
      complex data(MAXP,2),zero1,u0(maxp),w0(maxp),temp,u1(maxp),w1(maxp),datt(MAXP)
      double precision gnorm
      logical yes,yesno,rel
      integer blank,ounit
	integer nly,mmax,Lags,nAp,ITags,m,ip0
	integer itim(2)
	real*4 A0,ags,Ar,Az,Tags2
	real*4 t,dt,ta,pr,fny,delf
	real*4 blen,agauss,c,tdelay
	double precision cg(ntim),cg0,cgp(ntim),cgs(ntim)
c **********************************************************************
c      integer year,jday,hour,min,isec,msec
      character*8 sta,cmpnm,evnm
c      common /tjocm/ dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
c     *            cmpnm,az,cinc,evnm,baz,delta,rayp,depth,decon,agauss,
c     *              c,tq,instr,dlen,begin,t0,t1,t2
c
c **************************************************************************
c
c   parameter definitions may be found in sacio comments
c
c      common /win/ data
c      common /innout/ inunit,ounit
c      inunit=5
c      ounit=6
      zero1=cmplx(0.,0.)
      pi=3.141592654

!c    initialize parameters used in calculation of P-RFs
	iow=0
	blen = 5.
	c = 0.001
	tdelay = 20.
	ta = 40.
c	dt = 0.1	
c	pr = 0.06
c	if(ip .lt. 5) goto 1000

!c     parameters used during calculation of P-RFs
	npts=ifix(ta/dt+1.5)
	nft=npowr2(npts)
	nfpts=nft/2+1
	fny=1./(2.*dt)
	delf=2.*fny/float(nft)
	t=dt*nft
!c     parameters to control the selection of amplitude
!	cbeg = -1.
!	cend = 1.
!	nbeg = floor((cbeg+tdelay)/dt+0.5)
!	nend = floor((cend+tdelay)/dt+0.5)
!	nlen = nend-nbeg+1

	nly=mmax
	do 1 k=1,nly
	alfm(k)=a(k)
	betm(k)=b(k)
	rhom(k)=d(k)
	thikm(k)=h(k)
	vpvsm(k)=alfm(k)/betm(k)
1	continue
c	write(*,*)(alfm(k),k=1,nl)

	call caltim(alfm,betm,thikm,nly-1,pr,dtim)

	if(dtim .le. tdelay)then
	    atim = tdelay - dtim
	    call addly(alfm,betm,rhom,thikm,vpvsm,nly,atim,pr)
            call caltim(alfm,betm,thikm,nly-1,pr,dtim)

	endif

	call kesy(nly,alfm,betm,rhom,vpvsm,thikm,
     *dt,pr,ta,dtim,u0,w0,u1,w1)


	nbeg = ifix(0.5*(dtim-tdelay)/dt+0.5)
	if(nbeg .lt. 0)nbeg = 0

c	write(*,*)nbeg,dtim,tdelay
c	write(*,*)u0(1)
c	if(ip .eq. 2) tdelay = ta - tdelay
	nAp = int(tdelay/dt+1)
c	write(*,*)dtim,nbeg,nAp,nl
c	do 2 k=1,Lags
	agauss = ags
c	open(5,file='kenrf.dat',status='old')
c	read(5,*)iow,npts
c        read(5,*) dt0
c      c=ask('Trough filler, c =  ')
c      agauss=ask('Gaussian scale, a = ')
c      tdelay=ask('Enter phase shift: ')
c	read(5,*) blen,c,agauss,tdelay
c	close (5)
        
c      isyntp=2
c      if(rel) isyntp=1
c      nft=npowr2(npts)
c      nfpts=nft/2 + 1
c      fny=1./(2.*dt0)
c      delf=fny/float(nft/2)
c      write(*,102) npts,nft,fny,delf,dt
c  102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
c     *          'delf=',f8.4,1x,'dt=',f6.3)

c      do 1 i=1,3
c         call zero(data(1,i),1,MAXP)
c    1 continue
c        write(*,*) 'npts   ',npts,dt,tdelay,pr
	do i =1, 2
         call zero(data(1,i),1,maxp)
	enddo

	    do 11 j=1,nfpts
	       data(j,1)=u0(j+nbeg)
	       data(j,2)=w0(j+nbeg)
11	    continue


      if(iow.eq.0) go to 4
c         blen=ask('Length of window (in secs): ')
         call window(data,blen,npts,dt,3)
   4  continue
c
c  change delf to normalize deconvolution so dfftr works
c  dt from forward step is cancelled in decon, so no delf
c  on inverse, just 1/nft
c
      cdelf = 1. / float( nft )
   14 do 5 i=1,2
         call dfftr(data(1,i),nft,'forward',dt)
	
         if(i.ne.1) go to 5
         d2max=0.
         do 7 j=1,nfpts
            d2(j)=real(data(j,i)*conjg(data(j,i)))
            if(d2(j).gt.d2max) d2max=d2(j)
    7    continue
    5 continue
      decon=1.
      t0=tdelay
      phi1=c*d2max
      do 8 i=1,1
         do 9 j=1,nfpts
            freq=float(j-1)*delf
            w=2.*pi*freq
            phi=phi1
            if(d2(j).gt.phi) phi=d2(j)
            gauss=-w*w/(4.*agauss*agauss)
            data(j,i+1)=data(j,i+1)*conjg(data(j,1))*
     *                   cmplx(exp(gauss)/phi,0.)
            data(j,i+1)=data(j,i+1)*exp(cmplx(0.,-w*tdelay))
    9    continue
         call dfftr(data(1,i+1),nft,'inverse',cdelf)
    8 continue
c
c     deconvolve the vertical from itself using the
c         specified water-level parameter and gaussian
c
*     also compute the area under the gaussian filter
*
      gnorm = 0.0d0
      do 19 j=1,nfpts
         freq=float(j-1)*delf
         w=2.*pi*freq
         phi=phi1
         gnorm = gnorm + exp( gauss )
         gauss=-w*w/(4.*agauss*agauss)
         if(d2(j).gt.phi) phi=d2(j)
         data(j,1)=data(j,1)*conjg(data(j,1))*
     &   cmplx(exp(gauss)/phi,0.)
         data(j,1)=data(j,1)*exp(cmplx(0.,-w*tdelay))
19    continue

*
*     Finish the are integration
*
      gnorm = 2 * gnorm * delf

c
c*************************************************************
c
c     inverse transform the equalized vertical component
c
      call dfftr(data(1,1),nft,'inverse',cdelf)

c
c     compute the maximum value of the vertical component
c        to be used later in normailzation
c
      call minmax(data(1,1),npts,dmin,dmax,dmean)
      
*     Output the averaging function
*
*     This can be confusing, but bear with me.
*     To normalize to unit amplitude we must divide by the
*     area under the gaussian filter.  Also we must multiply
*     by dt since the result of a deconvolution of a function
*     from itself is a unit AREA spike (max amp = 1/dt).
*
      gnorm = gnorm * dt
      do 110 j = 1,npts
          data(j,1) = data(j,1) / gnorm
 110  continue
*
*
c!      t0=0.0
c!      begin = -tdelay
c
c     normalize the dmax for the transforms and gaussian
c     Not really necessary, horizontals have the same factors.
c
c     dmax = dmax *  float(nfpts) / gnorm
c
c************************************************************* 
c
c     gnorm = dmax * gnorm / nfpts
c
c     note that (dmax * gnorm / nfpts) = the unormalized dmax
c
      gnorm = dmax
      do 111 i = 2,2
      do 111 j = 1,npts
        data(j,i) = data(j,i) / gnorm
 111  continue

      do 81 i=2,2
          call minmax(data(1,i),npts,dmin,dmax,dmean)
   81 continue

      do 181 j=1,nfpts
        sr(2*(j-1)+1)=real(data(j,2))
        sr(2*j)=aimag(data(j,2))
 181  continue


	  do idt = 1, itim(2)
	    idt0 = nAp+itim(1)+idt-1
	    A0=sr(idt0)
	    cg(idt)=A0
	  enddo


      return
      end subroutine kerf
c
c************************************************************* 
c

      subroutine keat(a,b,d,h,mmax,cg,ags,pr,dt,itim,u0,w0,u1,w1,data)
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of langston (1979)
c    Out: cg -- RFHV Ratios(Svenningsen & Jacobsen, 2007) 
c	If ip = 1: P-RF; if ip = 2: S-RF. by Xu Wang at USTC, 1st/OCT,2019
c  *************************************************************************
c
	include 'apinv.inc'
	real*4 sr(MAXP),sz(MAXP)
	real*4 a(mlay),b(mlay),d(mlay),h(mlay)
	real*4 alfm(mlay),betm(mlay),vpvsm(mlay),rhom(mlay),thikm(mlay)
      dimension d2(MAXP)
      complex data(MAXP,2),zero1,u0(maxp),w0(maxp),temp,u1(maxp),w1(maxp)
      double precision gnorm
      logical yes,yesno,rel
      integer blank,ounit
	integer nly,mmax,Lags,nAp,ITags,m,ip0,it
	integer itim(2)
	real*4 A0,ags,Ar,Az,Tags2
	real*4 t,dt,ta,pr,fny,delf
	real*4 blen,agauss,c,tdelay
	double precision cg(ntim),cg0(ntim),cgp(ntim),cgs(ntim)
c **********************************************************************
      character*8 sta,cmpnm,evnm
c **************************************************************************
c
c   parameter definitions may be found in sacio comments
c
      zero1=cmplx(0.,0.)
      pi=3.141592654
!c    initialize parameters used in calculation of P-RFs
	iow=0
	blen = 5.
	c = 0.001
	tdelay = 20.
	ta = 40.
c	dt = 0.1	
c	pr = 0.06
c	if(ip .lt. 5) goto 1000

!c     parameters used during calculation of P-RFs
	npts=ifix(ta/dt+1.5)
	nft=npowr2(npts)
	nfpts=nft/2+1
	fny=1./(2.*dt)
	delf=2.*fny/float(nft)
	t=dt*nft

	nly=mmax
	do 1 k=1,nly
	alfm(k)=a(k)
	betm(k)=b(k)
	rhom(k)=d(k)
	thikm(k)=h(k)
	vpvsm(k)=alfm(k)/betm(k)
1	continue


	call caltim(alfm,betm,thikm,nly-1,pr,dtim)
	if(dtim .le. tdelay)then
	    atim = tdelay - dtim
	    call addly(alfm,betm,rhom,thikm,vpvsm,nly,atim,pr)
            call caltim(alfm,betm,thikm,nly-1,pr,dtim)
	endif
	call kesy(nly,alfm,betm,rhom,vpvsm,thikm,
     *dt,pr,ta,dtim,u0,w0,u1,w1)


	nbeg = ifix(0.5*(dtim-tdelay)/dt+0.5)
	if(nbeg .lt. 0)nbeg = 0

	nAp = int(tdelay/dt+1)
	agauss = ags

	do i =1, 2
         call zero(data(1,i),1,maxp)
	enddo

	    do 11 j=1,nfpts
	       data(j,1)=u0(j+nbeg)
	       data(j,2)=w0(j+nbeg)
11	    continue

c      if(iow.eq.0) go to 4
c         call window(data,blen,npts,dt,3)
c   4  continue
c
c  change delf to normalize deconvolution so dfftr works
c  dt from forward step is cancelled in decon, so no delf
c  on inverse, just 1/nft
c
      cdelf = 1. / float( nft )
   14 do 5 i=1,2
         call dfftr(data(1,i),nft,'forward',dt)
	
         if(i.ne.1) go to 5
         d2max=0.
         do 7 j=1,nfpts
            d2(j)=real(data(j,i)*conjg(data(j,i)))
            if(d2(j).gt.d2max) d2max=d2(j)
    7    continue
    5 continue
      decon=1.
      t0=tdelay
      phi1=c*d2max
      do 8 i=1,1
         do 9 j=1,nfpts
            freq=float(j-1)*delf
            w=2.*pi*freq
            phi=phi1
            if(d2(j).gt.phi) phi=d2(j)
            gauss=-w*w/(4.*agauss*agauss)
            data(j,i+1)=data(j,i+1)*conjg(data(j,1))*
     *                   cmplx(exp(gauss)/phi,0.)
            data(j,i+1)=data(j,i+1)*exp(cmplx(0.,-w*tdelay))
    9    continue
         call dfftr(data(1,i+1),nft,'inverse',cdelf)
    8 continue
c
c     deconvolve the vertical from itself using the
c         specified water-level parameter and gaussian
c
*     also compute the area under the gaussian filter
*
      gnorm = 0.0d0
      do 19 j=1,nfpts
         freq=float(j-1)*delf
         w=2.*pi*freq
         phi=phi1
         gnorm = gnorm + exp( gauss )
         gauss=-w*w/(4.*agauss*agauss)
         if(d2(j).gt.phi) phi=d2(j)
         data(j,1)=data(j,1)*conjg(data(j,1))*
     &   cmplx(exp(gauss)/phi,0.)
         data(j,1)=data(j,1)*exp(cmplx(0.,-w*tdelay))
19    continue

*
*     Finish the are integration
*
      gnorm = 2 * gnorm * delf

c
c*************************************************************
c
c     inverse transform the equalized vertical component
c
      call dfftr(data(1,1),nft,'inverse',cdelf)

c
c     compute the maximum value of the vertical component
c        to be used later in normailzation
c
      call minmax(data(1,1),npts,dmin,dmax,dmean)

*     Output the averaging function
*
*     This can be confusing, but bear with me.
*     To normalize to unit amplitude we must divide by the
*     area under the gaussian filter.  Also we must multiply
*     by dt since the result of a deconvolution of a function
*     from itself is a unit AREA spike (max amp = 1/dt).
*
      gnorm = gnorm * dt
      do 110 j = 1,npts
          data(j,1) = data(j,1) / gnorm
 110  continue
*
*
c      t0=0.0
c      begin = -tdelay
c
c     normalize the dmax for the transforms and gaussian
c     Not really necessary, horizontals have the same factors.
c
c     dmax = dmax *  float(nfpts) / gnorm
c
c************************************************************* 
c
c     gnorm = dmax * gnorm / nfpts
c
c     note that (dmax * gnorm / nfpts) = the unormalized dmax
c
      gnorm = dmax
      do 111 i = 2,2
      do 111 j = 1,npts
        data(j,i) = data(j,i) / gnorm
 111  continue

      do 81 i=2,2
          call minmax(data(1,i),npts,dmin,dmax,dmean)
   81 continue

      do 181 j=1,nfpts
        sr(2*(j-1)+1)=real(data(j,2))
        sr(2*j)=aimag(data(j,2))
        sz(2*(j-1)+1)=real(data(j,1))
        sz(2*j)=aimag(data(j,1))
 181  continue

c cutting out the RF recordings before itim(2)*dt seconds
c modified for total points according to input par.
	ict = itim(2) + itim(1)

	  do idt = 1, ict
	    idt0 = nAp+idt-1  !for P positive
	    idt1 = nAp-idt+1  !for P negative

	    sr0 = sr(idt0)
	    sz0 = sz(idt0)
		if(idt1 .gt. 0)then
		    sr1 = sr(idt1)	
	 	    sz1 = sz(idt1)
		else
		    sr1 = 0.	
	 	    sz1 = 0.
		endif
	    cgp(idt)=sr0+sr1
	    cgs(idt)=sz0+sz1
	  enddo

	    cgp(1)=cgp(1)/2
	    cgs(1)=cgs(1)/2	

c end cutting

c calculating the RF-HV-Ratio before (ict-1)*dt seconds 
c	    pinc = atan(cgp(1)/cgs(1))
c	    cg(1) = sin(pinc*0.5)/pr
	    do idt = 1, ict
		Ar = 0.
		Az = 0.
		JTags = idt
		Tags2 = real(JTags*2)-1
		do ia0 = 1, JTags
c		    ib0 = ia0 - 1
		    cosia = cos(real(ia0-1)/Tags2*pi)
		    cosia = cosia*cosia
		    Ar = Ar + cgp(ia0)*cosia
		    Az = Az + cgs(ia0)*cosia
		enddo
c	        pinc = atan(Ar/Az)
c		cg(idt) = sin(pinc*0.5)/pr
		cg0(idt) = Ar/Az
c		write(*,*)Jtags,cg(idt)
	    enddo
	it = 0
	do idt = itim(1)+1, ict
	    it = it + 1
	    cg(it) = cg0(idt)
	enddo
      return
      end


c
c************************************************************* 
c


      subroutine window(data,b,npts,dt,ndat)
      parameter(MAXP=8000, MAXP2=MAXP*2)
      dimension data(MAXP2,3)
c      common /win/ data
      data pi/3.1415926/
      bb=pi/b
      nend=ifix((b/dt)+.5) +1
      do 1 i=1,nend
      t=float(i-1)*dt
      windo=.5*(1. + cos(bb*t+pi))
         do 2 j=1,ndat
         data(i,j)=data(i,j)*windo
         data(npts+1-i,j)=data(npts+1-i,j)*windo
    2    continue
    1 continue
      return
      end

      subroutine minmax(x,npts,min,max,mean)
      dimension x(1)
      real min,max,mean
      min=9.0e+19
      max=-9.0e+19
      mean=0.
      do 1 i=1,npts
           if(x(i).gt.max) max=x(i)
           if(x(i).lt.min) min=x(i)
           mean=mean + x(i)
    1 continue
      mean=mean/float(npts)
      return
      end

      subroutine minmaxI(ID,npts,min,max)
      PARAMETER (NP=100)
      INTEGER ID(NP)
      integer min,max
      min=900000
      max=-900000
      do 1 i=1,npts
           if(ID(i).gt.max) max=ID(i)
           if(ID(i).lt.min) min=ID(i)
    1 continue
      return
      end

	subroutine caltim(alf,bet,thk,nl,p,dtim)
	real*4 alf(1),bet(1),thk(1),p,p2,dtim,vsl,vel2
	integer nl,i

	    dtim = 0.
	    p2 = p * p

	        do i = 1, nl
	            vel2 = alf(i)*alf(i)
	            vel2 = 1./vel2
	            vsl=sqrt(vel2 - p2)
		    dtim = dtim + thk(i)*vsl
	        enddo


	end

	subroutine addly(alfm,betm,rhom,thikm,vpvsm,nl,atim,p)
	real*4 alfm(1),betm(1),rhom(1),thikm(1),vpvsm(1)
	real*4 vel2,atim,p,p2
	integer nl

	p2 = p * p

		vel2 = alfm(nl)*alfm(nl)
	        vel2 = 1./vel2
	        vsl=sqrt(vel2 - p2)
c		thikm(nl) = 0.
		thikm(nl) = int(atim / vsl+1.5)
c		write(*,*)vsl,atim,thikm(nl)
		alfm(nl+1) = alfm(nl)
		betm(nl+1) = betm(nl)
		rhom(nl+1) = rhom(nl)
		thikm(nl+1) = 0.
		vpvsm(nl+1) = vpvsm(nl)
		nl = nl + 1
	end

      subroutine zero(x,start,end)
      dimension x(1)
      integer start,end
      do 1 i=start,end
    1 x(i)=0.
      return
      end
