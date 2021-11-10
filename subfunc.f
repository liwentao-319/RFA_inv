      SUBROUTINE SPLINE(M,N,X,Y,T,SS,SS1) 
	real*4 X(N),Y(N),S2(N),H(N),H1(N),B(N),C(N),
     &DELY(N),S3(N),DELSQY(N),T(M),SS(M),SS1(M) 
	N1=N-1 
	DO I=1,N1 
	    H(I)=X(I+1)-X(I) 
	    DELY(I)=(Y(I+1)-Y(I))/H(I) 
	enddo
	DO I=2,N1 
	    H1(I)=H(I-1)+H(I) 
	    B(I)=0.5*H(I-1)/H1(I) 
	    DELSQY(I)=(DELY(I)-DELY(I-1))/H1(I) 
	    S2(I)=2.0*DELSQY(I) 
	    C(I)=3.0*DELSQY(I) 
	enddo 
	S2(1)=0.0 
	S2(N)=0.0 
	DO 30 JJ=1,26 
	DO 30 I=2,N1 
	    S2(I)=(C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*1.0717968+S2(I) 
30      CONTINUE 
	DO 40 I=1,N1 
40          S3(I)=(S2(I+1)-S2(I))/H(I) 
	DO 50 J=1,M 
	    I=1 
	    IF((T(J)-X(I)).LE.0.0) GOTO 17 
	    IF((T(J)-X(N)).LT.0.0) GOTO 57 
	    GOTO 59 
56          IF((T(J)-X(I)).LT.0.0) GOTO 60 
	IF((T(J)-X(I)).EQ.0.0) GOTO 17 
57      I=I+1 
    	GOTO 56 
59      I=N 
60      I=I-1 
17      HT1=T(J)-X(I) 
	HT2=T(J)-X(I+1) 
	PROD=HT1*HT2 
	DELSQS=(2.0*S2(I)+S2(I+1)+HT1*S3(I))/6.0 
	SS(J)=Y(I)+HT1*DELY(I)+PROD*DELSQS 
	SS1(J)=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3(I)/6.0 
50      CONTINUE 
	RETURN 
	END


	subroutine cha_mod(a,b,d,lmax,dv)
	real a(1),b(1),d(1),dv(1)
	integer i,lmax
	real*4 kap,vso,dkapn,kapn,kapa

        do i = 1, Lmax
	    kapa = a(i)/b(i)
	    if(dv(i) .gt. 0.5) dv(i) = 0.5
	    if(dv(i) .lt. -0.5) dv(i) = -0.5
	    b(i) = b(i) + dv(i)
	    if(b(i) .gt. 5.0) b(i) = 5.0
	    if(b(i) .lt. 0.5) b(i) = 0.5
	    a(i) = 0.9409 + 2.0947*b(i) - 0.8206*b(i)**2+ 
     &0.2683*b(i)**3 - 0.0251*b(i)**4
	    d(i) = 1.6612*a(i) - 0.4721*a(i)**2 + 
     &0.0671*a(i)**3 - 0.0043*a(i)**4 +  
     &0.000106*a(i)**5
	    kapa = a(i)/b(i)

! below is not used since now
!	    a(i) = b(i) * kapa
!	    d(i) = 0.32*a(i)+0.77
!	    write(*,*)b(i),a(i),d(i),a(i)/b(i),kapa
	enddo

	end subroutine cha_mod

	subroutine formId(IdPw,IdPr)
	integer IdPw,IdPr
	    if(IdPw .ne. 1) IdPw = 0
	    if(IdPr .ne. 1) IdPr = 0
	end subroutine formId

c read input data file
	subroutine readat(datf,Mpd,pr,dt,itb,Obst,Lags)
	include 'apinv.inc'
	integer i,n,nb,Lags
	integer status
	real*4 Mpd(NP),pr(NP),dt(NP),Obst(NO)
	integer itb(NP,2)
	character*80 datf
	    n=0
	    nb=0
		open(5,file=datf,status='old')
		do while(.true.)
		n=n+1
		read(5,*,iostat=status)Mpd(n),pr(n),dt(n),itb(n,1),itb(n,2),
     &     (Obst(i),i=nb+1,nb+itb(n,2))
		dt(n) = real(int(dt(n)*1000))/1000
		nb = nb + itb(n,2)
		if(status.lt.0)goto 100
	        end do
100	    close(5)
	    Lags = n - 1
	end subroutine readat

c Program to read in the sigmal (standard deviation) for weighting the data residuals
c Written by Xu WANG, Feb/6th,2020
	subroutine readsigma(datf,itb,Obst)
	include 'apinv.inc'
	integer i,n,nb
	integer status
	real*4 Obst(NO)
	integer itb(NP,2)
	character*80 datf
	    n=0
	    nb=0
		open(5,file=datf,status='old')
		do while(.true.)
		n=n+1
		read(5,*,iostat=status)(Obst(i),i=nb+1,nb+itb(n,2))
		nb = nb + itb(n,2)
		if(status.lt.0)goto 100
	        end do
100	    close(5)
	return
	end subroutine readsigma

