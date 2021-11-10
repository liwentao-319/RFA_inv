!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!C    Program to construct the Vp, Vs, and Rho kernel of visual velocity which   C
!C    was derived from the direct P-Wave amplitude in P-RFs                      C
!C                 -- written by Xu Wang, at USTC, Jun/16th,2019                 C
!C    Insert into the program of Apinv.    Aug/6th, 2019                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   

	subroutine Apkernel(iphc,a,b,d,h,Lmax,Lags,Mpd,pr,pr0,dt,itb,sen_Vs)
	use omp_lib
	implicit none
	include 'apinv.inc'
	real,parameter :: pi=3.1415926535898
	real*4 a(mlay),b(mlay),d(mlay),h(mlay)
	real*4 ags,pd
	real*4 sen_Vs(NO,MLAY*2)
	INTEGER IA,IP,IZ,im,in,ipha,iphc,it,it0
	integer Lmax,Lags,mmax
        complex data(MAXP,2),u0(maxp),w0(maxp),u1(maxp),w1(maxp)
	real*4 dlnVp,dlnVs,dlnrho
	real*4 Mpd(NP),prd(NP),pr(NP),pr0(NP),dt(NP)
	real*4 vpm(mlay),vsm(mlay),rhom(mlay),thkm(mlay),dpt(mlay)
	double precision cg0(NTIM),cg1(NTIM),cg2(NTIM),cg
        double precision dlncg_dlnvs(NO,MLAY),dlncg_dlnvp(NO,MLAY),dlncg_dlnrho(NO,MLAY)
	integer i,k,kk,nbeg,nend,npowr2
	integer itb(NP,2),itim(2),ntps(NP+1)
	real coe_a,vpft,coe_rho

!c    set perturbation of Vs, Vp, and rho (e.g., dVs/Vs) for 
!c    the calculation of sensitivity kernel
	dlnVs = 0.001
	dlnVp = 0.001
	dlnrho = 0.001

	ntps(1) = 0
	do ia=1,Lags
	    prd(ia) = pi/Mpd(ia)
	    ntps(ia+1) = ntps(ia) + itb(ia,2)
	enddo

	ipha = iphc
  !$omp parallel num_threads(10)&
  !$omp default(private) &
  !$omp shared(a,b,h,d,mpd,lags,lmax,dlnVs,dlnVp,dlnrho) &
  !$omp shared(dlncg_dlnvs,dlncg_dlnvp,dlncg_dlnrho) &
  !$omp shared(itb,ntps,dt,pr,pr0,ipha)
  !$omp do


	    do ip=1, lmax       !loop for layer
		do ia = 1, lags !loop for agauss
		do iz=1,lmax
              	   vsm(iz) = b(iz)
             	   vpm(iz) = a(iz)
                   thkm(iz) = h(iz)
                   rhom(iz) = d(iz)
		enddo

!ccc ---- calculate sensitivity kernels for each layer

		pd = Mpd(ia)
		itim(1:2) = itb(ia,1:2)
		vsm(ip) = b(ip) - 0.5*dlnVs*b(ip)
	  	call selke(ipha,vpm,vsm,rhom,thkm,lmax,cg1,pd,pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)
		vsm(ip) = b(ip) + 0.5*dlnVs*b(ip)
	  	call selke(ipha,vpm,vsm,rhom,thkm,lmax,cg2,pd,pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)
		vsm(ip)=b(ip)
		do it=1,itim(2)
		it0 = it + ntps(ia)
                dlncg_dlnvs(it0,ip) = (cg2(it)-cg1(it))/vsm(ip)/dlnVs
		enddo
		vpm(ip) = a(ip) - 0.5*dlnVp*a(ip)
	  	call selke(ipha,vpm,vsm,rhom,thkm,lmax,cg1,pd,pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)
		vpm(ip) = a(ip) + 0.5*dlnVp*a(ip)
	  	call selke(ipha,vpm,vsm,rhom,thkm,lmax,cg2,pd,pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)		
		vpm(ip)=a(ip)
		do it=1,itim(2)
		  it0 = it + ntps(ia)
                  dlncg_dlnvp(it0,ip) = (cg2(it)-cg1(it))/vpm(ip)/dlnVp
		enddo
		rhom(ip) = d(ip) - 0.5*dlnrho*d(ip)
	  	call selke(ipha,vpm,vsm,rhom,thkm,lmax,cg1,pd,pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)
		rhom(ip) = d(ip) + 0.5*dlnrho*d(ip)
	  	call selke(ipha,vpm,vsm,rhom,thkm,lmax,cg2,pd,pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)
		rhom(ip)=d(ip)
		do it=1,itim(2)
		it0 = it + ntps(ia)
                dlncg_dlnrho(it0,ip) = (cg2(it)-cg1(it))/rhom(ip)/dlnrho
		enddo
		enddo   !end for agauss
	    enddo	!end for loop layer
  !$omp end do
  !$omp end parallel


	do im=1,ntps(Lags+1)
	    do in=1,Lmax
		coe_a=(2.0947-0.8206*2*b(in)+&
                     0.2683*3*b(in)**2-0.0251*4*b(in)**3)
                vpft=0.9409 + 2.0947*b(in) - 0.8206*b(in)**2+&
                     0.2683*b(in)**3 - 0.0251*b(in)**4
                coe_rho=coe_a*(1.6612-0.4721*2*vpft+&
                     0.0671*3*vpft**2-0.0043*4*vpft**3+&
                     0.000106*5*vpft**4)
        	sen_Vs(im,in)=sngl(dlncg_dlnvs(im,in))+&
			      sngl(dlncg_dlnvp(im,in))*coe_a+&
			      sngl(dlncg_dlnrho(im,in))*coe_rho
	    enddo
	enddo

	end subroutine Apkernel




!c   IdPS,a,b,d,h,lmax,cg,Mpd,Lags,pr,dt,itb,ntps
	subroutine syn(Iphc,a,b,d,h,lmax,cg,Mpd,Lags,pr,pr0,dt,itb)
	include 'apinv.inc'
	complex u0(maxp),w0(maxp),u1(maxp),w1(maxp),data(maxp,2)
	integer Lags,Lmax,iphc,ia,it,it0
	real*4 a(mlay),b(mlay),d(mlay),h(mlay)
	integer itim(2),itb(NP,2),ntps(NP+1)
	real*4 pr(NP),pr0(NP),dt(NP),cg(NO),Mpd(NP),cgt(Ntim)
	double precision cg0(NTIM)

	ntps(1) = 0
	do i=1,Lags
	    ntps(i+1) = ntps(i) + itb(i,2)
	enddo

  !$omp parallel num_threads(10)&
  !$omp default(private) &
  !$omp shared(a,b,h,d,mpd,lags,lmax) &
  !$omp shared(itb,ntps,dt,pr,pr0,Iphc,cg)
  !$omp do
	do ia = 1, Lags
	  itim(1:2) = itb(ia,1:2)
	  call selke(iphc,a,b,d,h,lmax,cg0,Mpd(ia),pr(ia),pr0(ia),dt(ia),itim,u0,w0,u1,w1,data)
	    do it=1,itim(2)
		it0 = it + ntps(ia)
		cg(it0) = sngl(cg0(it))
	    enddo
	enddo
  !$omp end do
  !$omp end parallel

	end subroutine syn


	subroutine selke(IdPS,a,b,d,h,lmax,cg0,Mpd,pr,pr0,dt,itim,u0,w0,u1,w1,data)
	include 'apinv.inc'
	complex u0(maxp),w0(maxp),u1(maxp),w1(maxp),data(maxp,2)
	integer Lmax
	real*4 a(mlay),b(mlay),d(mlay),h(mlay)
	integer itim(2),IdPS
	real*4 pr,pr0,dt,Mpd
	double precision cg0(NTIM)

	if(IdPS .le. 2)then
	  call kerf(a,b,d,h,lmax,cg0,Mpd,pr,dt,itim,u0,w0,u1,w1,data)
	elseif(IdPS .le. 4)then
	  call keat(a,b,d,h,lmax,cg0,Mpd,pr,dt,itim,u0,w0,u1,w1,data)
	endif

	end subroutine selke



