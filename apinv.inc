!c
!c	model parameters
!c
!c      MLAY: max layer No. of velocity model
!c      NP: max period No. for each data set
!c	NTIM: max data point No. for one dataset with sigle period
!c	NO: total observation sample points
!c	NL: for indexes of G matrix
!c	MAXP: for forwarding receiver functions

	integer MLAY,NP,MAXP,NTIM,NO,NL,NKD,NOTT
	PARAMETER (MLAY=100)
	PARAMETER (NP=15)
	PARAMETER (NTIM=201)
	PARAMETER (NO=NP*NTIM)
	
	PARAMETER (NKD=2)
	PARAMETER (NOTT=NKD*NO)
	PARAMETER (NL=NOTT*MLAY)

!c 	If you do not know the subroutine well, do not change this unless
!c      warning in the terminal appears, because I found that if it is too
!c      big, the exe will be wrong under parallel condition but well under
!c      single-core condition.
	PARAMETER (MAXP=8000)
