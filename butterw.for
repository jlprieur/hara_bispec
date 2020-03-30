C****************************************************************
C Program BUTTERW
C Applies a Butterworth filter to a spectrum, to reduce
C oscillations during inverse Fourier Transformation
C
C JLP
C Version: 06-02-90
C****************************************************************
	PROGRAM BUTTERW
	PARAMETER(IDIM=256)
	CHARACTER NAME*40,COMMENTS*80
	REAL IN(IDIM,IDIM),OUT(IDIM,IDIM)
 
C Format of the files:
	PRINT *,' WARNING: all the files should not have odd',
     1	' numbers for NX and NY (size in X and Y)'
	CALL JLP_INQUIFMT
 
C Input of the squared modulus:
	WRITE(6,37)
37	FORMAT(' Applies Butterworth filter,',/,
     1	' Filter(f)= 1/(1+(f/fc)**2K)',/,
     1	' Cut Frequency fc (in pixel units), and',
     1	' order K (3 is not bad) :')
	READ(5,*) FC,KK
 
	WRITE(6,38)
38	FORMAT(' (The Fourier images are supposed to be',
     1	' centered in the frame)')
	NAME=' '
	CALL JLP_READIMAG(IN,NX,NY,IDIM,NAME,COMMENTS)
 
	CALL BUTTERWORTH(IN,OUT,NX,NY,IDIM,FC,KK)
 
	WRITE(COMMENTS,11) FC,KK,NAME(1:15)
11	FORMAT(' Butterworth ',F9.2,I3,' of:',A15)
	NAME=' '
	CALL JLP_WRITEIMAG(OUT,NX,NY,IDIM,NAME,COMMENTS)
 
	CALL JLP_END
	END
 
C*******************************************************************
C BUTTERWORTH
C
C*******************************************************************
	SUBROUTINE BUTTERWORTH(IN,OUT,NX,NY,IDIM,FC,KK)
	INTEGER NX,NY,IDIM,IXC,IYC,KK,KK2
	REAL IN(IDIM,*),OUT(IDIM,*),FC
 
C Center:
	IXC=(NX/2)+1
	IYC=(NY/2)+1
 
C Main loop:
	KK2=2*KK
	DO 91 J=1,NY
	  DJ=IYC-J
	  DO 92 I=1,NX
	    DI=IXC-I
	    FF=SQRT(DJ*DJ+DI*DI)
	    OUT(I,J)=IN(I,J)/(1.+(FF/FC)**KK2)
92	  CONTINUE
91	CONTINUE
	
	RETURN
	END
