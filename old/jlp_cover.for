C++********************************************************************
C Set of routines for uv-coverage, spectral and bispectral lists:
C
C Contains:
C COVERA, COVER, AFFECTE
C
C Subroutine COVERA
C To compute the elements of the uv coverage, and of the A matrix
C (To solve the equation A*X=Y, that is to invert the bispectral relations)
C
C Output:
C IR: Maximum radius of the uv-coverage
C NBETA: Number of elements of the spectral list (Number of columns of A)
C NGAMMA: Number of elements of the bispectral list (Number of rows of A)
C
C JLP
C Version 23-11-91
C--********************************************************************
	SUBROUTINE COVERA(IR,NBETA,NGAMMA)
 
C  POUR IR=25
	PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C For IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)
 
	INTEGER*4 NBCOUV,IXY,NGT
C NBCOUV( X from -IRMAX to IRMAX,  Y from 0 to IRMAX )
C IXY( 1 for X; 2 for Y,  NB from 0 to NBMAX )
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)
 
C Computing the A matrix as generated by the uv-coverage (only defined by IR)
	CALL COUV(IR,NBETA,NGAMMA)
 
	WRITE(6,28) IR,NBETA,NGAMMA
C	WRITE(2,28) IR,NBETA,NGAMMA
28	FORMAT(' UV_coverage:   IR =',I4,/,
     1	' NBETA (spectral list) =',I6,/,
     1	' NGAMMA (bispect. list)  =',I8)
 
C Computing the number of closure relations for some values of NB:
	DO 8 NB=3,NBETA,100
	   NFERM = NGT(NB) - NGT(NB-1)
	   WRITE(6,29) NB,NFERM
C	   WRITE(2,29) NB,NFERM
29	   FORMAT(' NB =',I6,' Closure relations:',I8)
8	CONTINUE
 
	RETURN
	END
C*******************************************************************
C COUV defines the uv-coverage and the A matrix:
C Input:
C IR: maximum radius of the uv-coverage
C
C Output:
C NBETA, NGAMMA: number of elements of the spectral and bispectral lists
C
C In common blocks (output): the uv-coverage is accessible from 2 sides:
C
C NBCOUV(I,J): uv-coverage (i.e. number of the spectral list for
C              the pixel(I,J) of the spectrum (I=0,J=0 for null frequency)
C
C IXY(1,I) and IXY(2,I) coordinates of the pixel number I in the spectral list
C              (this allows another entry for the uv-coverage)
C
C************************************************************************
	SUBROUTINE COUV(IR,NBETA,NGAMMA)
 
C  POUR IR=25
	PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C For IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)
 
C NBCOUV( X from -IRMAX to IRMAX,  Y from 0 to IRMAX )
C IXY( 1 for X; 2 for Y,  NB from 0 to NBMAX )
	INTEGER*4 NBCOUV,IXY
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
	INTEGER NCARRE(NULL:50),J2,IR2
 
C Computing the squares:
	DO 1 J=0,IR
	  NCARRE(J) = J**2
1	CONTINUE
	IR2MAX=NCARRE(IR)
 
C Easy cases:
	NB=0
	IXY(1,NB)=0
	IXY(2,NB)=0
	NBCOUV(0,0)=NB
 
	NB=1
	IXY(1,NB)=1
	IXY(2,NB)=0
	NBCOUV(1,0)=NB
 
	NB=2
	IXY(1,NB)=0
	IXY(2,NB)=1
	NBCOUV(0,1)=NB
 
C Reseting the total number of elements the  bispectral list
	NG=0
 
C Main loop: work with successive iterations on circles with
C increasing radii
C Squared radius: IR2 = 2, 3, ... ,IR2MAX
 
	DO 2 IR2=2,IR2MAX
 
C Searching for the couples (I,J) such as: I**2 + J**2 = IR2 with I>=J
	  DO 3 J=0,IR
	    J2=NCARRE(J)
	    I2=IR2-J2
	    IF (I2.LT.J2) GO TO 2
	    I = INT(SQRT(FLOAT(I2)))
	    II2=I*I
C Selecting the points defined by each couple (I,J):
	    IF (II2.EQ.I2)THEN
	      IRS=INT(SQRT(FLOAT(IR2)))
	      CALL AFFECTE(I,J,IRS,NB,NG)
C Now use the symmetry relations:
	      IF (I.NE.J) CALL AFFECTE(J,I,IRS,NB,NG)
	      IF (J.NE.0) CALL AFFECTE(-J,I,IRS,NB,NG)
	      IF (I.NE.J.AND.J.NE.0) CALL AFFECTE(-I,J,IRS,NB,NG)
	    ENDIF
3	  CONTINUE
 
2 	CONTINUE
 
C NBETA: Total number of the spectral list (Number of columns of the X matrix)
C NGAMMA: Total number of the bispectral list (Number of rows of the X matrix)
C (Remember, we have to solve    A*X = Y)
	NBETA=NB
	NGAMMA=NG
 
	RETURN
	END
C*******************************************************************
C AFFECTE:
C Gives a structure to the S group of the A matrix
C Input:
C ISX, ISY : coordinates
C IRS: radius
C
C Output:
C NB: index in the spectral list
C NG: index in the bispectral list
C*******************************************************************
 
	SUBROUTINE AFFECTE(ISX,ISY,IRS,NB,NG)
 
C  POUR IR=25
	PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C For IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)
 
	INTEGER*4 NBCOUV,IXY,NGT,KLM
C NBCOUV( X from -IRMAX to IRMAX,  Y from 0 to IRMAX )
C IXY( 1 for X; 2 for Y,  NB from 0 to NBMAX )
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
C KLM(1,.) for K ; KLM(2,.) for L ; KLM(3,.) for M ;
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C5/NGT(NBMAX)
 
C Recording the new point of the uv coverage (spectral list):
	NB=NB+1
	IXY(1,NB)=ISX
	IXY(2,NB)=ISY
	NBCOUV(ISX,ISY)=NB
 
C Searching for the couples associated with the point NBS=NB
C Generating the rows of the A matrix:
	NBS=NB
 
C Loop on all the possible points (Q):
	DO 1 NBQ=1,NBS-1
 
C Coordinates of the vector T = S - Q
	  ITX=ISX-IXY(1,NBQ)
	  ITY=ISY-IXY(2,NBQ)
 
C Work within the circle of radius IRS, so we can a priori reject
C the points outside the window [-IRS,+IRS]:
	IF((ITX.GE.-IRS.AND.ITX.LE.IRS)
     1	.AND.(ITY.GE.-IRS.AND.ITY.LE.IRS))THEN
 
	  IF((ITY.GT.0).OR.
     1	(ITY.EQ.0.AND.ITX.GE.0))THEN
 
C Case number 1 (which could be : k=t, l=q, k+l=s)
	    NBK=NBCOUV(ITX,ITY)
 
C We select this couple (U,V) if the vector NBK is in [0,NBQ]
	    IF (NBK.NE.0.AND.NBK.LE.NBQ)THEN
	      NG=NG+1
	      KLM(1,NG)=NBK
	      KLM(2,NG)=NBQ
	      KLM(3,NG)=NBS
	    ENDIF
 
	ELSE
 
C Case number 2 (which could be : r=-t, s=s, r+s=m=q)
 
	  NBR=NBCOUV(-ITX,-ITY)
C We select this couple (R,S) if the vector NBR is in [0,NBS]
	  IF (NBR.NE.0.AND.NBR.LE.NBS)THEN
	    NG=NG+1
	    KLM(1,NG)=NBR
	    KLM(2,NG)=NBS
	    KLM(3,NG)=NBQ
	  ENDIF
 
C Nota: we can only have L=NBS or M=NBS (never K=NBS)

	ENDIF
 
	ENDIF
1	CONTINUE
 
C NGT(NB) is the number of the last U,V couple of the group NB=NBS=S:
	NGT(NBS)=NG
 
	RETURN
	END
