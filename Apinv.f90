!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!C         Program to invert direct P-waves in P-RFs for Vs structures           C
!C 1)        -- Written by Xu Wang, at USTC, Aug/6th, 2019                       C
!C 2)        -- Trans. from matlab to fortran                                    C
!C 3)        -- Modified to parallel calculation Sep/24th, 2019                  C
!C 4)        -- From sigle- to multiple-time-point  Sep/28th, 2019               C
!C 5)        -- Add the RFH/V Ratios(Svenningsen & Jacobsen, 2007) Oct/8th, 2019 C
!C 6)        -- Simply the architecture for saving time. Oct/31th, 2019. IGGCAS  C
!C 7)        -- Objective function is determined through a hybrid sigmal(std)-   C
!C              weighted residuals, i.e., SUM((d_i-f_i)/sigm_i). Feb/6th, 2020   C
!C 8)        -- More info. refers to Wang et al., unpublished                    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   

	program Apinv
        use lsmrModule, only:lsmr
        use  lsmrblasInterface, only : dnrm2
	implicit none
	include 'apinv.inc'

	real*4 a(mlay),b(mlay),d(mlay),h(mlay),dep(mlay)
	real*4 cbst(NOTT),sen_Vs(NOTT,MLAY*2*NKD),cbst0(NOTT)
!c for int-p
	real*4 aa(mlay),bb(mlay),dd(mlay),hh(mlay),depn(mlay),SS1(mlay)
	integer Nmax,ithk

! for sloving Ax = b
	real*4 cg(NO),cgt(NTIM)
	real*4 rw(NL),col(NL)
	integer Lmax,Lags
	integer dall,count0
	integer m,n,status
	integer leniw,lenrw,iw(NL)
	integer ir
	real*4 dv(mlay*2)
	real*4 dpitv
	Integer Idmin,Idmax,Idme,Id
	real acond,anorm,atol,btol,damp,rnorm,xnorm,arnorm,conlim
	integer itn,itnlim,localsize,nar
	integer ibo,ido,il,im,in,index,istop,iter,k,Lags0
	integer MCOL,nb,nout,indx,i,inc

! data P-RF waveform
	real*4 Mppw(NP),prpw(NP),dtpw(NP)
	real*4 Obspw(NO),senpw(NO,MLAY*2),sigpw(No)
	integer itbpw(NP,2),Lagpw,ntpspw(NP+1),dallpw
! data P-RF ratio
	real*4 Mppr(NP),prpr(NP),dtpr(NP)
	real*4 Obspr(NO),senpr(NO,MLAY*2),sigpr(No)
	integer itbpr(NP,2),Lagpr,ntpspr(NP+1),dallpr

! used for inv.
	integer  itermax,isig,iterNo
	real*4 smooth,resid,resid0
	real*4 thick,Thkw(10)
	integer Tly
	character*80 chatmp,modfil
	character*80 datfilPw,datfilPr
	character*80 sigfilPw,sigfilPr
	integer IdPw,IdPr
	real*4 	Wall,WdPw,WdPr
	integer maxalldat
!c Handle NO. for each file
	real*4 temp,residm,smooth0,sumvs
!c output 
	integer it,it0,itim(2),ia,icg

	open(unit=15,file='Apinv.in',status='old')
	    read(15,301) chatmp
	    read(15,301) modfil !modelfile 
	    read(15,301) chatmp
	    read(15,*) itermax,isig
	    read(15,301) chatmp
	    read(15,*) smooth0,damp
	    read(15,301) chatmp
	    read(15,*) thick,Tly,(Thkw(i),i=1,Tly)
	    read(15,301) chatmp
	    read(15,*) IdPw,WdPw
	    read(15,301) chatmp
	    read(15,*) IdPr,WdPr
	close(15)
301  	format(A)
	
	datfilPw = '../data/synPw.dat'
	datfilPr = '../data/synPr.dat'
	sigfilPw = '../data/stdPw.dat'
	sigfilPr = '../data/stdPr.dat'
!c re-set the weight for each dataset, which if not used will be set to zero
	call formId(IdPw,IdPr)
	if(IdPw .ne. 1) WdPw = 0.
	if(IdPr .ne. 1) WdPr = 0.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c below for real data
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	Lagpw = 0
	Lagpr = 0
	dallpw = 0 
	dallpr = 0 

	Wall = WdPw + WdPr
	WdPw = WdPw/Wall
	WdPr = WdPr/Wall

	write(*,*)'Weights of Pw and Pr are'
	write(*,*)WdPw,WdPr

	if(IdPw .eq. 1)then
	    call readat(datfilPw,Mppw,prpw,dtpw,itbpw,Obspw,Lagpw)
	    ntpspw(1) = 0
	    do i = 1, Lagpw
	        dallpw = dallpw + itbpw(i,2)
	        ntpspw(i+1) = ntpspw(i) + itbpw(i,2)
	    enddo
	    sigpw(1:dallpw) = 1.0
	    if(isig .eq. 1)call readsigma(sigfilPw,itbpw,sigpw)
	    write(*,*)'Lagpw = ',Lagpw,'; dallpw = ',dallpw
	endif

	if(IdPr .eq. 1)then
	    call readat(datfilPr,Mppr,prpr,dtpr,itbpr,Obspr,Lagpr)
	    ntpspr(1) = 0
	    do i = 1, Lagpr
	        dallpr = dallpr + itbpr(i,2)
	        ntpspr(i+1) = ntpspr(i) + itbpr(i,2)
	    enddo
	    sigpr(1:dallpr) = 1.0
	    if(isig .eq. 1)call readsigma(sigfilPr,itbpr,sigpr)
	    write(*,*)'Lagpr = ',Lagpr,'; dallpr = ',dallpr
	endif

	maxalldat = dallpw
	if(dallpr .gt. maxalldat)maxalldat=dallpr
	if(dallpw .gt. 0)WdPw = WdPw * maxalldat / dallpw
	if(dallpr .gt. 0)WdPr = WdPr * maxalldat / dallpr

	smooth0 = smooth0 * maxalldat
!	write(*,*)'smooth',smooth0
	write(*,*)'Weights of Pw and Pr based on data No. pre-sigh weights are'
	write(*,*)WdPw,WdPr
! read in reference velocity model
	call get_mod(modfil,lmax,a,b,d,h)
!	if(Lmax .gt. mmax)Lmax=mmax
! end here
! total number of observations
	Lags = Lagpw + Lagpr 
	write(*,'(A,I10)')'Period No. =',Lags
	dall = dallpw + dallpr
	write(*,'(A,I10)')'Data Points No. =',dall

!c Mcol: column of G matrix
!	damp = 0.1
	iterNo = 1
	if(thick .lt. 0) goto 2001
	do 2000 ithk = 1, Tly
!c Interplotation
	dep(1) = h(1)/2
        do i = 2, Lmax-1
	    dep(i) = dep(i-1) + 0.5*(h(i)+h(i-1))
        enddo
	if(Lmax .ge. 2) dep(Lmax) = dep(Lmax-1) + 2*h(Lmax-1)
	dpitv = Thkw(ithk)
	Nmax = ifix(thick/dpitv+1.5)
	depn(1)=dpitv/2
	hh(1) = dpitv
	do i = 2, Nmax
	    depn(i) = depn(i-1) + dpitv
	    hh(i) = dpitv
	enddo
	call spline(Nmax,Lmax,dep,a,depn,aa,SS1)
	call spline(Nmax,Lmax,dep,b,depn,bb,SS1)
	call spline(Nmax,Lmax,dep,d,depn,dd,SS1)
	WRITE(*,*)LMAX,NMAX,DPITV

	Lmax = Nmax
	a(1:Nmax) = aa(1:Nmax)
	b(1:Nmax) = bb(1:Nmax)
	d(1:Nmax) = dd(1:Nmax)
	h(1:Nmax-1) = dpitv
	h(Nmax) = 0.
!c Ending intp...
!	itermax = 0
2001	continue
! Out model par. during each iteration
	open(28,file='../OUT/model_iter.dat',status='unknown')
	open(29,file='../OUT/mod_data_resid.dat',status='unknown')

	sumvs = 0.
	do i = 1, Lmax
	    sumvs = sumvs + b(i)
	enddo
	sumvs = sumvs / Lmax
	Mcol = lmax

	do 1000 iter = 1, itermax
	inc = 0

	if(IdPw .eq. 1)then

	  call syn(1,a,b,d,h,lmax,cg,Mppw,Lagpw,prpw,prpw,dtpw,itbpw)
	  call Apkernel(1,a,b,d,h,lmax,Lagpw,Mppw,prpw,prpw,dtpw,itbpw,senpw)

          do i = 1,dallpw
	      inc = inc + 1
	      sen_Vs(inc,1:Mcol) = senpw(i,1:Mcol) * WdPw / sigpw(i)
              cbst(inc) = (obspw(i) - cg(i)) * WdPw / sigpw(i)
	      cbst0(inc) = obspw(i) - cg(i)
          enddo
	endif

	if(IdPr .eq. 1)then
	  call syn(3,a,b,d,h,lmax,cg,Mppr,Lagpr,prpr,prpr,dtpr,itbpr)
	  call Apkernel(3,a,b,d,h,lmax,Lagpr,Mppr,prpr,prpr,dtpr,itbpr,senpr)
          do i = 1,dallpr
	      inc = inc + 1
	      sen_Vs(inc,1:Mcol) = senpr(i,1:Mcol)  * WdPr / sigpr(i)
              cbst(inc) = (obspr(i) - cg(i))  * WdPr / sigpr(i)
              cbst0(inc) = obspr(i) - cg(i)
          enddo
	endif

!! out residuals of data and initial model
	if(iterNo .eq. 1)then
	    open(13,file='../OUT/residual_first',status='unknown')
            do k = 1,dall
	        write(13,*)cbst(k),cbst0(k)
            enddo	
	    close(13)

	    open(121,file='../OUT/model_ini',status='unknown')
            do k = 1,Lmax
	        write(121,*)a(k),b(k),d(k),h(k),a(k)/b(k)
            enddo	
	    close(121)

	    iterNo = 0
	endif
!! end out

! parameter counts for matrix
	temp = 0.
	do im = 1, dall
	    do in = 1, Mcol
		if(abs(sen_Vs(im,in)).ge.temp)temp=sen_Vs(im,in)
	    enddo
	enddo
	nar=0
	do im = 1, dall
	    do in = 1, Mcol
		if(abs(sen_Vs(im,in)) .gt. abs(0.001*temp))then
			nar=nar+1
			iw(nar+1)=im
			col(nar)=in
			rw(nar)=sen_Vs(im,in)
		endif
	    enddo
	enddo
	write(*,*)'NAR before and after decrease SEN',dall*Mcol,NAR,NL
	write(*,*)'Max sensitivity',temp

	! ADDING REGULARIZATION TERM
	smooth = smooth0 / (Lmax-2) / sumvs
	write(*,*)'smooth after weighting',smooth0
	count0 = 0 ! index for the term
	if(Lmax .gt. 2 .and. smooth .gt. 0.0)then
	  do ir = 1, Lmax-2
	    count0 = count0 + 1
	    iw(1+nar+1) = dall + count0
	    col(nar+1) = ir
	    rw(nar+1) = smooth * (-1.)
	    iw(1+nar+2) = dall + count0
	    col(nar+2) = ir + 1
	    rw(nar+2) = smooth * 2.
	    iw(1+nar+3) = dall + count0
	    col(nar+3) = ir + 2
	    rw(nar+3) = smooth * (-1.)
	    cbst(dall+count0) = 0.
	    nar=nar+3
	  enddo
	endif !ENDING REGULARIZATION

	lenrw = nar
	leniw = 2*nar+1
	iw(1) = nar
	iw(nar+2:nar*2+1)=col(1:nar)
        atol = 1e-4
        btol = 1e-4
        conlim = 100
        itnlim = 400
        istop = 0
        anorm = 0.0
        acond = 0.0
        arnorm = 0.0
        xnorm = 0.0
        localSize = 10
	m = dall + count0
	n = Mcol
        nout=-1 ! control to output LMQR info. or not, if > 0, yes.
	if(nout .gt. 0)open(nout,file='lsmr.txt',status='unknown')
	dv(1:Mcol) = 0.0
        call LSMR(m, n, leniw, lenrw,iw,rw,cbst, damp,&
                atol, btol, conlim, itnlim, localSize, nout,&
                dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)

	call cha_mod(a,b,d,lmax,dv)
	resid = 0.
        do i = 1,dall
	    resid = resid + cbst(i)*cbst(i)
        enddo
	resid = sqrt(resid)
	resid0 = 0.
        do i = 1,dall
	    resid0 = resid0 + cbst0(i)*cbst0(i)
        enddo
	resid0 = sqrt(resid0)

!	close(nout)
	    write(28,*)iter,Lmax,-1,-1,-1
	    residm = 0.
        do k = 1,Lmax
	    write(*,*)a(k),b(k),d(k),h(k),a(k)/b(k)
	    write(28,*)a(k),b(k),d(k),h(k),a(k)/b(k)
	    if(k .gt.1 .and. k .lt. Lmax)then
		residm = residm + (b(k-1)+b(k+1)-2*b(k))**2
	    endif
        enddo
	    residm = sqrt(residm/(Lmax-2))
	write(*,*)'--iter: ',iter,'--DATA Res aft, bef weighted:',resid,resid0
	write(*,*)'--iter: ',iter,'--MOD Res:',residm
	write(29,*)iter,residm,resid,resid0

1000	continue
2000	continue
	close(28)
	close(29)

	open(12,file='../OUT/model_final',status='unknown')
	write(12,*)Lmax
        do k = 1,Lmax
	    write(12,*)a(k),b(k),d(k),h(k),a(k)/b(k)
        enddo	
	close(12)

	if(IdPw .eq. 1)then
	  call syn(1,a,b,d,h,lmax,cg,Mppw,Lagpw,prpw,prpw,dtpw,itbpw)
	  open(12,file='../OUT/synthe.pw',status='unknown')
	    do ia = 1, Lagpw
	      itim(1:2) = itbpw(ia,1:2)
	      do it=1,itim(2)
		it0 = it + ntpspw(ia)
		cgt(it) = cg(it0)
	      enddo
	      write(12,*)Mppw(ia),prpw(ia),dtpw(ia),itbpw(ia,1),itbpw(ia,2),(cgt(icg),icg=1,itbpw(ia,2))
	   enddo
	  close(12)
	endif

	if(IdPr .eq. 1)then
	  call syn(3,a,b,d,h,lmax,cg,Mppr,Lagpr,prpr,prpr,dtpr,itbpr)
	  open(12,file='../OUT/synthe.pr',status='unknown')
	    do ia = 1, Lagpr
	      itim(1:2) = itbpr(ia,1:2)
	      do it=1,itim(2)
		it0 = it + ntpspr(ia)
		cgt(it) = cg(it0)
	      enddo
	      write(12,*)Mppr(ia),prpr(ia),dtpr(ia),itbpr(ia,1),itbpr(ia,2),(cgt(icg),icg=1,itbpr(ia,2))
	   enddo
	  close(12)
	endif

	open(23,file='../OUT/residual',status='unknown')
        do k = 1,dall
	    write(23,*)cbst(k),cbst0(k)
        enddo	
	close(23)

	stop
	end program Apinv

!c    subroutine for getting model parameters from input file
	subroutine get_mod(mod0,mmax0,aa,bb,dd,hh)
	include 'apinv.inc'
	integer mmax0,k0
	real*4 aa(mlay),bb(mlay),dd(mlay),hh(mlay),vr(mlay)
	character*80 mod0

	open(unit=12,file=mod0,status='old')
	read(12,*) mmax0
        do k0=1,mmax0
	  read(12,*) bb(k0),hh(k0)
          aa(k0) = 0.9409 + 2.0947*bb(k0) - 0.8206*bb(k0)**2+ &
                   0.2683*bb(k0)**3 - 0.0251*bb(k0)**4
          dd(k0) = 1.6612*aa(k0) - 0.4721*aa(k0)**2 + &
                   0.0671*aa(k0)**3 - 0.0043*aa(k0)**4 + & 
                   0.000106*aa(k0)**5
	  vr(k0) = aa(k0)/bb(k0)
      enddo

! Not use since now
!c	do 1 k0=1,mmax0
!c	    read(12,*) vr(k0),bb(k0),hh(k0)
!c	    aa(k0) = vr(k0) * bb(k0)
!c	    dd(k0) = aa(k0) * 0.32 + 0.77
!c	write(*,*)  aa(k0),bb(k0),dd(k0),hh(k0),vr(k0)
!c1	continue

	close(12)
	end


