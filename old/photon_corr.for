C****************************************************************
C PROGRAM  PHOTON_CORR
C To correct from photon noise
C
C Other programs: phot_noise.c and phot_noise_mask.c
C (should be used if uv-mask or very big uv-coverage)
C
C SYNTAX:
C     runs photon_corr ir_max input_extension ffield-name output_extension
C Example:
C     runs photon_corr 12 _5 ffield1 _5p
C Put 0 if no flat field, and give the mean number of photons for correction:
C     runs photon_corr 12 _5 0 34.5 _5p
C
C JLP
C Version: 06-10-93
C****************************************************************
	PROGRAM PHOTON_CORR
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)
 
C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
	REAL YCE1(NGMAX,3)
	REAL MODSQ(IDIM,IDIM),FFIELD(IDIM,IDIM),LONG(IDIM,IDIM)
	CHARACTER NAME*40,COMMENTS*80,BISP_NAME*40,MODSQ_NAME*40
	CHARACTER DATE*24,FFIELD_NAME*80,LONG_NAME*80,IN_EXT*30,OUT_EXT*30
        LOGICAL FFIELD_CORRECTION
	INTEGER NBCOUV,IXY,NX,NY
 
C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A  NBMAX)
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
	CALL JLP_BEGIN
	OPEN(2,FILE='photon_corr.log',STATUS='UNKNOWN',ERR=998)
	CALL JLP_INQUIFMT
 
	CALL JLP_DATE_TIME(DATE)
	WRITE(2,52) IDIM,IDIM,DATE
	WRITE(6,52) IDIM,IDIM,DATE
52	FORMAT(' Program  Photon_corr. Version 06-10-93 (IR_UVcoverage=25)',
     1         /,' Max size of input images:',I4,1X,I4,/,A)
 
10      FORMAT(A)

	WRITE(6,92)
92	FORMAT(' Radius (IR) of uv-coverage:')
	READ(5,*) IR
	WRITE(2,59) IR
59	FORMAT(' IR = ',I4)

	WRITE(6,17)
17	FORMAT(' Input file extension :')  
        READ(5,10) IN_EXT

C UV coverage:
	CALL COVERA(IR,NBETA,NGAMMA)

C Input of the biased square modulus:
	WRITE(6,37)
37	FORMAT(' Mean squared modulus of the FFT',
     1	' (centered in the frame) :')
	WRITE(MODSQ_NAME,18) IN_EXT
18      FORMAT('modsq',A)
	CALL JLP_READIMAG(MODSQ,NX,NY,IDIM,MODSQ_NAME,COMMENTS)
	WRITE(2,57) MODSQ_NAME(1:14)
57	FORMAT(' Mean square modulus of the fft: ',A14)
 
C Input of the flat field:
	WRITE(6,38)
38	FORMAT(' Corresponding flat field :', 
     1        ' (Enter 0 if no flat field correction) ')
	READ(5,10) FFIELD_NAME
        IF(FFIELD_NAME(1:1).EQ.'0')THEN
           FFIELD_CORRECTION=.FALSE.
        ELSE
           FFIELD_CORRECTION=.TRUE.
	   CALL JLP_READIMAG(FFIELD,NX,NY,IDIM,FFIELD_NAME,COMMENTS)
	   WRITE(2,39) FFIELD_NAME(1:14)
39	   FORMAT(' Flat field: ',A14)
 	   WRITE(6,68)
68	   FORMAT(' Long integration :') 
	   WRITE(LONG_NAME,81) IN_EXT
81         FORMAT('long',A)
	   CALL JLP_READIMAG(LONG,NX,NY,IDIM,LONG_NAME,COMMENTS)
	   WRITE(2,69) LONG_NAME(1:14)
69	   FORMAT(' Long integration: ',A14)
        ENDIF
 
C Bispectrum:
	WRITE(6,41)
41	FORMAT(' Input mean bispectrum:')
	WRITE(BISP_NAME,19) IN_EXT
19      FORMAT('bisp1',A)
	CALL JLP_READIMAG(YCE1,NX1,NY1,NGMAX,BISP_NAME,COMMENTS)
	WRITE(2,58) BISP_NAME(1:14)
58	FORMAT(' Input mean bispectrum: ',A14)
	IF((NY1.LT.2).OR.(NX1.NE.NGAMMA))THEN
	   WRITE(6,45)
	   WRITE(2,45)
45	   FORMAT(' Photon_corr/FATAL ERROR: Size of bispectrum',
     1    ' inconsistent with IRMAX')
	   GOTO 999
	ENDIF

	CALL CORREC1(MODSQ,NX,NY,IDIM,YCE1,IR,NBETA,NGAMMA,
     1  XPHOTONS,FFIELD,LONG,FFIELD_CORRECTION)

C Output of the corrected frames:
	WRITE(6,87)
87	FORMAT(' Output file extension :')  
        READ(5,10) OUT_EXT

	WRITE(COMMENTS,79)MODSQ_NAME(1:14),BISP_NAME(1:14),XPHOTONS
79	FORMAT('Photon corr of ',A14,' and ',
     1    A14,' Nph=',G12.5)
	WRITE(6,67)
67	FORMAT(' ***** Output the corrected mean square modulus')
	WRITE(NAME,82) OUT_EXT
82      FORMAT('modsq',A)
	CALL JLP_WRITEIMAG(MODSQ,NX,NY,IDIM,NAME,COMMENTS)
	WRITE(6,61)
61	FORMAT(' ***** Output the corrected mean bispectrum')
	WRITE(NAME,83) OUT_EXT
83      FORMAT('bisp1',A)
	CALL JLP_WRITEIMAG(YCE1,NX1,NY1,NGMAX,NAME,COMMENTS)

999     CLOSE(2)
	CALL JLP_END
	STOP
998	PRINT *,' Fatal error opening photon_corr.log'
	CALL JLP_END
	STOP
	END
C*****************************************************************
C Only for non-normalized files!
C*****************************************************************
	SUBROUTINE CORREC1(MODSQ,NX,NY,IDIM,YCE1,IR,NBETA,
     1        NGAMMA,XPHOTONS,FFIELD,LONG,FFIELD_CORRECTION)
 
C  POUR IR=25
C	PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL MODSQ(IDIM,*),FFIELD(IDIM,*),LONG(IDIM,*)
	REAL YCE1(NGMAX,3),RO(NULL:NBMAX)
	REAL XPHOTONS,XPHOT1
        REAL*8 SUM1,SUM2,SUM3
        LOGICAL NORMALIZED,FFIELD_CORRECTION
	INTEGER IR,NBETA,NGAMMA
	INTEGER K1,K2,K3
	INTEGER NBCOUV,IXY,NX,NY
 
C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A NBMAX )
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
	COMMON /C2/KLM(3,NGMAX)
 
C Central position: 
	IXC = (NX/2) + 1
	IYC = (NY/2) + 1

C At the center, should be XPHOT**2 (or 1. if normalized), so:
        NORMALIZED=((MODSQ(IXC,IYC)-1.).LE.0.01)
        IF(NORMALIZED.AND.(.NOT.FFIELD_CORRECTION))THEN
           WRITE(6,45)
           WRITE(2,45)
45	   FORMAT(' CORREC1/Fatal error: input modsq normalized!')
           STOP
        ENDIF
 
C Biased square modulus:
	DO 1 NB = 0,NBETA
	  IIX = IXC + IXY(1,NB)
	  IIY = IYC + IXY(2,NB)
	  RO(NB)=MODSQ(IIX,IIY)
1	CONTINUE
 
C Debug:
C	DO I=0,4
C	   WRITE(6,23) I,RO(I)
C	   WRITE(2,23) I,RO(I)
C23	   FORMAT(' RO(',I2,') = ',G12.5)
C	END DO

C Check the FFT (zero frequency is at (1,1)):
	XPHOT1=SQRT(RO(0))
	PRINT *,' CORREC1/XPHOT1 (from central value)= ',XPHOT1
 
C Now compute the sums implied in photon noise correction:
        IF(FFIELD_CORRECTION)THEN
          CALL SUM_WEIGHTED(FFIELD,LONG,NX,NY,IDIM,SUM1,SUM2,SUM3)
        ELSE
          SUM2=XPHOT1
	  PRINT *,' Enter the number of photons (for correction): '
          READ(5,*) SUM2
        ENDIF

        WRITE(6,19) SUM2
19      FORMAT(' Bispectrum correction with SUM2=',G12.5)

C Phase factor of the bispectrum (with bispectral list):
C and correction from photon noise effects:
	DO NG=1,NGAMMA
	  K1=KLM(1,NG)
	  K2=KLM(2,NG)
	  K3=KLM(3,NG)
C Debug:
	  IF(NG.LT.4)THEN
	   PRINT *,' K1,K2,K3,RO(K1,K2,K3),YCE1(R,I) (before)',
     1     K1,K2,K3,RO(K1),RO(K2),RO(K3),YCE1(NG,1),YCE1(NG,2)
	  ENDIF
	  YCE1(NG,1)=YCE1(NG,1)-(RO(K1)+RO(K2)+RO(K3))+2.*SUM2
	  WORK=RO(K1)+RO(K2)+RO(K3)-2.*XPHOT1
C Debug:
	  IF(NG.LT.4)THEN
	   PRINT *,' CORR,YCE1(R,I) (after)',
     1      WORK,YCE1(NG,1),YCE1(NG,2)
	  ENDIF
	END DO
 
C Phase factor of the bispectrum:
          DO NG=1,NGAMMA
            XMOD=YCE1(NG,1)*YCE1(NG,1)+YCE1(NG,2)*YCE1(NG,2)
            XMOD=SQRT(XMOD)
	    YCE1(NG,1)=YCE1(NG,1)/XMOD
	    YCE1(NG,2)=YCE1(NG,2)/XMOD
          END DO

C Photon noise correction and normalization of the squared modulus 
        WRITE(6,20) SUM2
20      FORMAT(' Modulus correction with SUM2=',G12.5)

	WORK=MODSQ(IXC,IYC)-SUM2
	DO J=1,NY
	  DO I=1,NX
	    MODSQ(I,J)=(MODSQ(I,J)-SUM2)/WORK
	  END DO
	END DO
 
C Returns nphot for the output comments:
        XPHOTONS = SUM2
	RETURN
	END
C**********************************************************************
C Compute sum of square values of the flat field weighted by
C the long exposure
C**********************************************************************
        SUBROUTINE SUM_WEIGHTED(FFIELD,LONG,NX,NY,IDIM,SUM1,SUM2,SUM3)
        REAL*4 FFIELD(IDIM,*),LONG(IDIM,*)
        REAL*8 SUM1,SUM2,SUM3,SUM_WEIGHTS
        REAL*4 MEAN1,MEAN2 
        INTEGER*4 NX,NY

        SUM_WEIGHTS=0.
        SUM1=0.
        SUM2=0.
        SUM3=0.
        DO J=1,NY
          DO I=1,NX
C First determines how many events have occurred:
C number_of_photo_events(i,j) = cst * long(i,j) / ffield(i,j)
C Thus there is one order less than expected for the sums
C Weighted_sum_of_squares = long(i,j)*ffield(i,j)
            SUM1=SUM1+LONG(I,J)
            SUM2=SUM2+LONG(I,J)*FFIELD(I,J)
            SUM3=SUM3+LONG(I,J)*FFIELD(I,J)*FFIELD(I,J)
            SUM_WEIGHTS=SUM_WEIGHTS+LONG(I,J)
          END DO
        END DO

        WRITE(6,20) SUM1,SUM2,SUM3,SUM_WEIGHTS
20      FORMAT(' Sum of FFIELD weighted by LONG: ',G12.5,
     1     /,' Sum of FFIELD**2 weighted by LONG: ',G12.5,
     1     /,' Sum of FFIELD**3 weighted by LONG: ',G12.5,
     1     /,' Sum of weights (LONG): ',G12.5)

        IF(SUM_WEIGHTS.EQ.0) THEN
           WRITE(6,21)
21         FORMAT(' Fatal error: sum of weigths (LONG) is null !')
           STOP
        ENDIF
         
C Compute mean ffield value (weighted by long exposure):
        MEAN1=SUM1/SUM_WEIGHTS
        WRITE(6,22) MEAN1 
22      FORMAT(' Mean ffield value (weighted by long exposure): ',G12.5) 

C Compute mean square value (weighted by long exposure):
        MEAN2=SUM2/SUM_WEIGHTS
        WRITE(6,23) MEAN2 
23      FORMAT(' Mean ffield square value (weighted by long exposure): ',G12.5) 

        RETURN
        END
C**********************************************************************
	INCLUDE 'jlp_cover_30.for'
