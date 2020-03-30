C++***************************************************************
C POSITIV
C Program to project a complex image to a positiv real image.
C
C JLP
C Version of 26-04-90
C--***************************************************************
	PROGRAM POSITIV
	PARAMETER (IDIM=256)
C Index 0: Initial (and final) data:
	REAL RE0(IDIM,IDIM),IM0(IDIM,IDIM)
     	INTEGER IRMAX,KMAX
	CHARACTER NAME*40,COMMENTS1*80,COMMENTS2*80
 
	CALL JLP_BEGIN
	OPEN(2,FILE='positiv.log',STATUS='UNKNOWN',ERR=998)
 
C Inquire the format (input/output) :
	CALL JLP_INQUIFMT
 
C Options :
	WRITE(6,15)
15	FORMAT(' IRMAX and number of iterations ?')
	READ(5,*) IRMAX,KMAX
	WRITE(6,24) IRMAX,KMAX
24	FORMAT(' IRMAX, KMAX :',I5,X,I5)
 
C Input of the real part :
	WRITE(6,*)' Input image (real part):'
	NAME=' '
	CALL JLP_READIMAG(RE0,NX,NY,IDIM,NAME,COMMENTS1)
	WRITE(2,*)' Input image (real part):',NAME
 
C Input of the imaginary part :
	WRITE(6,*)' Input image (imaginary part):'
	NAME=' '
	CALL JLP_READIMAG(IM0,NX1,NY1,IDIM,NAME,COMMENTS2)
	WRITE(2,*)' Input image (imaginary part):',NAME
 
C Check the size:
	  IF((NX1.NE.NX).OR.(NY1.NE.NY))THEN
	    WRITE(2,21)
	    WRITE(6,21)
21	   FORMAT(' FATAL ERROR: Wrong size for the input images',/,
     1	' (should be the same for the real/imaginary parts)')
	    GOTO 999
	  ENDIF
 
	DO K=1,KMAX
	  WRITE(2,*) ' Iteration:',K
	  WRITE(6,*) ' Iteration:',K
	  CALL POSITIV1(RE0,IM0,NX,NY,IRMAX)
	END DO
 
C Output :
	WRITE(6,*)' Real part:'
	NAME=' '
	CALL JLP_WRITEIMAG(RE0,NX,NY,IDIM,NAME,COMMENTS1)
	WRITE(2,*)' Output image (real part):',NAME
 
	WRITE(6,*)' Imaginary part:'
	NAME=' '
	CALL JLP_WRITEIMAG(IM0,NX,NY,IDIM,NAME,COMMENTS2)
	WRITE(2,*)' Output image (imaginary part):',NAME
 
999	CLOSE(2)
	CALL JLP_END
	STOP
998 	PRINT *,' Fatal error opening fft2.log'
	CALL JLP_END
	STOP
	END
C********************************************************************
C POSITIV1
C Input:
C RE0, IM0: real and imaginary part
C
C Output:
C RE0, IM0: real and imaginary part
C*******************************************************************
	SUBROUTINE POSITIV1(RE0,IM0,NX,NY,IRMAX)
	PARAMETER (IDIM=256)
C Index 0: Initial (and final) data:
	REAL RE0(IDIM,*),IM0(IDIM,*),MOD0(IDIM,IDIM)
C After "positive" operation
	REAL RE1(IDIM,IDIM),IM1(IDIM,IDIM),PHA1_R,PHA1_I,MOD1
	INTEGER NX,NY,IRMAX
 
C Copies the input array to RE1, IM1, and makes the real part positive:
	DO J=1,NY
	  DO I=1,NX
	    RE1(I,J)=MAX(0.,RE0(I,J))
	    IM1(I,J)=IM0(I,J)
	  END DO
	END DO
 
C Computes real and imaginary part and modulus of the FFT
C of the initial image:
	CALL FFT2D(RE0,IM0,NX,NY,IDIM,1)
 
	DO J=1,NY
	  DO I=1,NX
	    MOD0(I,J)=SQRT(RE0(I,J)**2+IM0(I,J)**2)
	  END DO
	END DO
 
C Computes phase and modulus of the FFT of the "positiv" image:
	CALL FFT2D(RE1,IM1,NX,NY,IDIM,1)
 
C Computing now the real and imaginary part where needed:
	IC=NX/2+1
	JC=NY/2+1
	RMAX=FLOAT(IRMAX)
 
	DO J=1,NY
	  XJ=(J-JC)**2
	  DO I=1,NX
	    XI=(I-IC)**2
	    XRAD=SQRT(XI+XJ)
	    IF(XRAD.GT.RMAX)THEN
	      MOD1=SQRT(RE1(I,J)**2+IM1(I,J)**2)
	      PHA1_R=RE1(I,J)/MOD1
	      PHA1_I=IM1(I,J)/MOD1
	      RE0(I,J)=MOD0(I,J)*PHA1_R
	      IM0(I,J)=MOD0(I,J)*PHA1_I
	    ENDIF
	  END DO
	END DO
 
C Inverse FFT now:
	CALL FFT2D(RE0,IM0,NX,NY,IDIM,-1)
 
	RETURN
	END
C********************************************************************
C FFT 2D routine
C Input:
C TR, TI:  real part and imaginary part  (input/output) of the image
C
C KOD=1, or 2 direct
C KOD=-1, or -2 inverse
C KOD=-2, or 2: power spectrum, and phase
C********************************************************************
	SUBROUTINE FFT2D(TR,TI,NX,NY,IDIM,KOD)
	PARAMETER (IDIM1=256)
	REAL TR(IDIM,*),TI(IDIM,*)
	DOUBLE PRECISION DR(IDIM1*IDIM1),DI(IDIM1*IDIM1),WORK(2000)
	INTEGER NX,NY,KOD
	INTEGER LWORK,ND(2),NDIM
 
C Copy the arrays:
	II=0
	DO J=1,NY
	  DO I=1,NX
	    II=II+1
	    DR(II)=TR(I,J)
	    DI(II)=TI(I,J)
	  END DO
	END DO
 
C Compute the FFT:
	NDIM=2
	ND(1)=NX
	ND(2)=NY
	LWORK=2000
	NN=NX*NY
	IFAIL=0
	IF(KOD.GT.0)THEN
	   CALL C06FJF(NDIM,ND,NN,DR,DI,WORK,LWORK,IFAIL)
	ELSE
	   CALL C06GCF(DI,NN,IFAIL)
	   CALL C06FJF(NDIM,ND,NN,DR,DI,WORK,LWORK,IFAIL)
	   CALL C06GCF(DI,NN,IFAIL)
	ENDIF
 
C Copy the arrays:
	II=0
	DO J=1,NY
	  DO I=1,NX
	    II=II+1
	    TR(I,J)=DR(II)
	    TI(I,J)=DI(II)
	  END DO
	END DO
 
	RETURN
	END
