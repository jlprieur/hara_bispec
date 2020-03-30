C****************************************************************
C BISP_TO_MOD.FOR
C Program to extract the square modulus from the bispectrum
C
C JLP
C Version: 20-11-92
C****************************************************************
	PROGRAM BISP_TO_MOD
 
C  POUR IR=25
        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)
 
C  POUR IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
 
	REAL IMAGE(IDIM,IDIM),MODSQ(NBMAX),BISP(NGMAX,3)
	INTEGER NGT,NG1,NG2,NB,NBETA,NGAMMA
	CHARACTER BISP_NAME*60,MODSQ_NAME,COMMENTS*80
        INTEGER KLM,NBCOUV,IXY
        INTEGER NX,NY,NX_MIN,NY_MIN
        REAL XR(NULL:NBMAX),XI(NULL:NBMAX)


C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A NBMAX )
        COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)

C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
        COMMON /C2/KLM(3,NGMAX)

	COMMON /C5/NGT(NBMAX)

	CALL JLP_BEGIN
	CALL JLP_INQUIFMT
C
	WRITE(6,22) IDIM,IDIM
22	FORMAT(' Program BISP_TO_MOD Version 04-11-92',/,
     1   ' To extract the square modulus from the bispectrum ',/,
     1   '  (images up to ',I4,'x',I4,')')

C UV coverage:
	WRITE(6,*) ' Radius (IR) of uv-coverage?'
	READ(5,*) IR
	CALL COVERA(IR,NBETA,NGAMMA)

C Bispectrum: 
	WRITE(6,41)
41	FORMAT(' Input bispectrum (phasor and snr)')
	BISP_NAME=' '
	CALL JLP_READIMAG(BISP,NX1,NY1,NGMAX,BISP_NAME,COMMENTS)
	IF((NY1.LT.3).OR.(NX1.NE.NGAMMA))THEN
	   WRITE(6,45)
45	   FORMAT(' FATAL ERROR: Size of bispectrum inconsistent',
     1	' with IRMAX')
	   STOP
        ENDIF

C Output size of square modulus:
        NX_MIN=2*IR+2
        NY_MIN=NX_MIN
        WRITE(6,51) NX_MIN,NY_MIN
51	FORMAT(' Size of the output 2-D file containing the',
     1  ' square modulus: NX,NY ',/,
     1  '(larger than ',I5,'x',I5,')')
        READ(5,*) NX,NY

C First erases image:
        DO J=1,NY
          DO I=1,NX
	   IMAGE(I,J)=0.
          END DO
	END DO

C Computes the square modulus for each value of  the spectral list: 
	NG1=1
	DO 1 NB=3,NBETA
	  NG2=NGT(NB)
	  INDX=0
	  DO NG=NG1,NG2
            K1=KLM(1,NG)
            K2=KLM(2,NG)
            K3=KLM(3,NG)
C Here K1<K2 so:
            IF(K1.EQ.1)PRINT *,' NB=',NB,' NG=',NG,' NG1=',NG1,'NG2 =',NG2
C K1=0 is not allowed!!!!!!!!
C i.e. zero frequency is not allowed!!!!!!!!
          END DO
         NG1=NG2+1
1	CONTINUE

        IF(I.NE.123)STOP
C Copies back the square modulus (spectral list) to 2-D image:
        DO NB=0,NBETA
C IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
C element NB of the spectral list.
C As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
C we do the following transformation:
          IIX=MOD(IXY(1,NB)+NX,NX)+1
          IIY=MOD(IXY(2,NB)+NY,NY)+1
          IMAGE(IIX,IIY)=MODSQ(NB)
C Writes symmetrical also:
          IMAGE(IXC-IIX,IYC-IIY)=MODSQ(NB)
        END DO

C Output of the 2D image:
	WRITE(6,42)
42	FORMAT(' Output of the 2-D file with the square modulus')
	WRITE(COMMENTS,43) BISP_NAME(1:40)
43	FORMAT(' Extracted square modulus from ',A)
	MODSQ_NAME=' '
	CALL JLP_WRITEIMAG(IMAGE,NX,NY,IDIM,MODSQ_NAME,COMMENTS)
 
	CALL JLP_END
	STOP
	END
C**********************************************************************
	INCLUDE 'jlp_cover.for'
