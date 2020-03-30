C****************************************************************
C PROGRAM INV_BISPEC.FOR : 2DPA (2D PHASE A)
C To compute the phase of the spectrum from the bispectrum,
C assuming a full pupil.
C
C Older version than INV_BISPEC1
C Look also at INV_BISPEC1 which uses jlp_cover_mask...
C
C Contains CREYCE_SIMU, CREYCE1, CREYCE2, TRANS, RECURSIVE, EPLUS,
C          WEIGHTS_SIGM1, WEIGHTS_SIGM2
C          NOYAU, PROJ_EPLUS, ATRANS, ATRANS_A, CGRADIENT, ERROR, SORTIE
C          BISP_WEIGHT,BISP_WEIGHT2
C
C Syntax:
C
C Example1:
C RUNS INV_BISPEC 0 12,20,0.08,220,0.004,0.8 real_fft ima_fft 
C
C   with 0=simulation option
C   12=radius of uv cover., 20 = weights or phase error, 0.08=test_to_stop,
C   220=nbeta, 0.004=low_modulus, 0.8= max sigma
C   real_fft and ima_fft= FFT of image to be reconstructed
C
C Example2:
C RUNS INV_BISPEC 1 12,20,0.08,220,0.004,0.8 modsq bisp1
C
C   with 1=real data, 
C   12=radius of uv cover., 20 = weights or phase error, 0.08=test_to_stop,
C   220=nbeta, 0.004=low_modulus, 0.8 = max sigma
C   modsq=mean squared modulus of FFT
C   bisp1=mean bispectrum
C
C Example3:
C RUNS INV_BISPEC 2 12,20,0.08,220,0.004,0.8 modsq real_fft ima_fft bisp1
C
C   with 2=simulated data, 
C   12=radius of uv cover., 20 = weights or phase error, 0.08=test_to_stop,
C   220=nbeta, 0.004=low_modulus, 0.8= max sigma
C   modsq=mean squared modulus of FFT (simulated) 
C   real_fft and ima_fft= FFT of image to be reconstructed
C   bisp1=mean bispectrum (simulated)
C
C JLP
C Version: 02-02-94
C From A. Lannes P.FOR (January 1988)
C****************************************************************
	PROGRAM INV_BISPEC
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)
 
C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
 
C XC: Spectrum (list): phase factor
C RO: Spectrum (list): modulus
C YCE: Bispectrum (list)
C XCR: Reference spectrum (list), phase factor  (not known for real observ.)
	COMPLEX XC,YCE,XCR,CC
	INTEGER IXY,KLM,IFERMAX,NBCOUV,NX,NY,NBX,NBY
        REAL SIGM(NBMAX),SIG_MAX,X,BMASK,XXC,XXCR
	REAL RE,IM,RO,QMOY,QMOY0,WEIGHT,LOWER_MODULUS
	CHARACTER FNAME*40,COMMENTS*80,DATE*24,ERRNAME*40
 
        INTEGER*4 MEMSIZE,PNTR_ARRAY,MADRID(1)

C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A  NBMAX)
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C4/X(NULL:NBMAX,NULL:3)
	COMMON /C5/NGT(NBMAX)
	COMMON /C6/XCR(NULL:NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
C Simply to give a limit to the number of closure relations to
C be taken into account (when this number would be really too big...)
	COMMON /C11/IFERMAX
C For virtual memory:
        COMMON /VMR/MADRID
 
	CALL JLP_BEGIN
        CALL JLP_RANDOM_INIT(1)
	OPEN(2,FILE='inv_bispec_30.log',STATUS='UNKNOWN',ERR=998)
C Max. number of closure phase relations allowed for each pixel of
C the uv-coverage:
	IFERMAX = NBMAX
 

C FORMAT FICHIER
	PRINT *,' WARNING: all the files should not have odd',
     1	' numbers for NX and NY (size in X and Y)'
	CALL JLP_INQUIFMT
 
	WRITE(6,22)
22	FORMAT(' Menu : ',/,' 0=simulation (re, im),'/,
     1	' 1= observations (mosq, bisp)',/,
     1	' 2=simulation (modsq, bisp, re, im)',/,
     1	' Enter the number of the option :')
	READ(5,*) IOPT
	CALL JLP_DATE_TIME(DATE)
	WRITE(2,52) DATE,IOPT,IFERMAX
	WRITE(6,52) DATE,IOPT,IFERMAX
52	FORMAT(' Program inv_bispec_30 Version 15-04-92 ',
     1	A,/,' Option = ',I3,/,
     1  ' Maximum number of closure relations allowed here:',I8)
 
C Angular reference for bispectral perturbation (degrees) 
C Cette reference est associee au carre de la norme de la frequence IR
	CTE = 20.
 
C Test for the main iteration (in degrees)
C This means that when the largest angular correction
C is larger than TEST, we exit from the main loop:
	TEST = 0.1
 
	WRITE(6,92)
92	FORMAT(' Radius (IR) of uv-coverage, phase error ',
     1    '(0 if all weights=1, <0 if modulus weighted), test (degrees),'
     1    ' nbeta, lower_modulus, max sigma, sigma_value_when_modulus_is_null',
     1    '  (Ex: 12,20,0.1,220,0.,0.8,0.7)')
	READ(5,*) IR,CTE,TEST,NBETA,LOWER_MODULUS,SIG_MAX,
     1    SIGMA_NULL
	WRITE(2,59) IR,CTE,TEST,LOWER_MODULUS,SIG_MAX,SIGMA_NULL
59	FORMAT(' IR = ',I4,' CTE = ',F8.2,
     1  ' TEST = ',G12.4,' Lower_modulus: ',G12.5,/,
     1  ' Sig_max:',G10.3,' Sig_null:',G10.3)

C UV coverage:
	CALL COVERA(IR,NBETA_MAX,NGAMMA_MAX)
 
        IF(NBETA.GT.NBETA_MAX)THEN
	  WRITE(6,229) NBETA,NBETA_MAX
	  WRITE(2,229) NBETA,NBETA_MAX
229	  FORMAT(' JLP_BISPEC/Error: NBETA (wanted) =',I4,
     1         ' whereas NBETA_MAX =',I5,/,' I correct it to BETA_MAX')
          NBETA = NBETA_MAX
        ENDIF

	NGAMMA=NGT(NBETA)
	WRITE(2,29) NBETA,NGAMMA
29	FORMAT(' NBETA (used here):',I5,' Corresponding NGAMMA:',I8)

C Max number of iterations in the main loop:
	ITTMAX = 15
 
	WRITE(6,53) CTE,TEST
53	FORMAT(' Angular constant for bispectral noise estimation',
     1 F8.3,' (degrees)',/,' TEST (for exit check) :',E12.3)
 
C Creates YCE : bispectrum phasor,
C and initial set of weights 
	IF(IOPT.EQ.0)THEN
 	  CALL CREYCE_SIMU(IR,CTE,NBETA,NGAMMA,FNAME)
	ELSEIF(IOPT.EQ.1) THEN
 	  CALL CREYCE1(IR,CTE,NBETA,NGAMMA,NGAMMA_MAX,FNAME,
     1    LOWER_MODULUS)
	ELSE
 	  CALL CREYCE2(IR,CTE,NBETA,NGAMMA,NGAMMA_MAX,FNAME,
     1    LOWER_MODULUS)
	ENDIF
 
C******************************************************
C To initialize the solution,
C we solve the problem with the recursive method (Weigelt,...) 
C The initial spectral phase is stored in XC(.)
	CALL RECURSIVE(SIGM,NBETA,NGAMMA,SIGMA_NULL)
 
C New version of the weights:
C JLP93:
        IF(CTE.EQ.0)THEN
	  CALL WEIGHTS_SIGM1(SIGM,NBETA,NGAMMA,SIG_MAX,LOWER_MODULUS)
        ELSE
	  CALL WEIGHTS_SIGM2(SIGM,NBETA,NGAMMA,SIG_MAX,LOWER_MODULUS)
        ENDIF

C Output of SNR map:
        CALL OUTPUT_SNR1(SIGM,NBETA,FNAME)

C Output the errors of initial bispectrum:
        IF(IOPT.LT.-1045)THEN
	  ERRNAME='bisp_error1.dat'
          NBX=NBETA
          NBY=NBETA/2
          MEMSIZE=NBX*NBY*4
          CALL JLP_GETVM(PNTR_ARRAY,MEMSIZE)
	  CALL ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,MADRID(PNTR_ARRAY),
     1   NBX,NBY)
	   ERRNAME='bisperr1'
	   COMMENTS='Quadratic errors'
           CALL JLP_WRITEIMAG(MADRID(PNTR_ARRAY),NBX,NBY,NBX,
     1  ERRNAME,COMMENTS)
           CALL JLP_FREEVM(PNTR_ARRAY,MEMSIZE)
        ENDIF

C Initial error estimation:
	IF(IOPT.NE.1)THEN
	  ERRNAME='error1.dat'
	  CALL ERROR_SIMU(NBETA,EI,EPI,ERRNAME)
	  WRITE(6,64) EI,EPI
	  WRITE(2,64) EI,EPI
64	  FORMAT(' Comparison with the model (since simulation)',/,
     1  ' Initial rms error of the spectrum',1PG12.5,/,
     1	' Initial rms error of the phase factor of the spect.',1PG12.5)
	ENDIF
 
C Projection of the initial solution onto E^+
	CALL EPLUS(NBETA)
 
C Other alternative (Null phases)
C	DO 2 NB=0,NBETA
C2	XC(NB)=(1.,0.)
 
C SORTIE DES PARTIES REELLES ET IMAGINAIRES DE XC
C In rei and imi (i for initial)
	CALL SORTIE(NBETA,1,FNAME)

C LA SUITE PRINCIPALE
 
	DEGRAD = 3.141592654/180.
	TEST = TEST*DEGRAD
 
	QMOY0 = 0.
 
	DO 3 ITT=1,ITTMAX
C Computes phase error term X(.,0) for PHI^+
C Iterative solution by conjugate gradients.
 
	  CALL CGRADIENT(NBETA,NGAMMA,IT,QMOY)
	  WRITE(6,77) QMOY,ITT,IT
	  WRITE(2,77) QMOY,ITT,IT
77	  FORMAT(' rms bisp error',1PG12.5,
     1    ' ITT =',I3,' Internal IT =',I4,' done')
 
C SUP : Maximum of the solution PHI^+ X(.,0) (L1 norm).
	  SUP=0.
	  DO 4 NB=1,NBETA
              R = X(NB,0)*BMASK(NB)
	      R=ABS(R)
	      SUP=AMAX1(R,SUP)
4	  CONTINUE
 
C Test of convergence: 
	  IF (SUP.LT.TEST) THEN
	   WRITE(6,*) ' Normal exit: SUP .LT. TEST '
	   WRITE(2,*) ' Normal exit: SUP .LT. TEST '
	   GO TO 6
	  ENDIF
 
C Generating the new value of the spectral term XC(.)
	  DO 5 NB=0,NBETA
            R=X(NB,0)*BMASK(NB)
	    CALL JLP_COS(COSR,R)
	    CALL JLP_SIN(SINR,R)
	    XC(NB)=XC(NB)*CMPLX(COSR,SINR)
5	  CONTINUE
 
C TEST SUR L'EVOLUTION DE QMOY
	DELQ = ABS(QMOY0 - QMOY) / QMOY
	IF (DELQ.LT.2.E-4.AND.ITT.GT.10)THEN
	   WRITE(6,49)
	   WRITE(2,49)
49	   FORMAT(' Warning: Exit because the relative DELQ is too small')
	   GO TO 6
	ELSE
	   QMOY0 = QMOY
	ENDIF
3	CONTINUE
 
C Sortie de la suite principale
6	CONTINUE
 
C Files ref et imf (f for final)
	CALL SORTIE(NBETA,2,FNAME)
 
C   Calage en translation
	CALL TRANS(IR,NBETA)

C ERREUR FINALE RF
 
C ERREURS FINALES GLOBALES ET BILAN
	IF(IOPT.NE.1)THEN
	  ERRNAME='error2.dat'
	  CALL ERROR_SIMU(NBETA,EF,EPF,ERRNAME)
 
	  WRITE(2,65) EI,EF,EPI,EPF
	  WRITE(6,65) EI,EF,EPI,EPF
65	  FORMAT(' Comparison with the model (since simulation)',/,
     1  ' Initial & final rms error of the spectrum:',2(1X,F8.3),/,
     1  ' Initial & final rms error of the phase factor (spect):',
     1   2(1X,F8.3))
 
C ECART DE PHASE ANGULAIRE POINT PAR POINT EXPRIME EN RADIAN
	  DO 10 NB = 0,NBETA
	    CC = XCR(NB)
	    CALL JLP_ATAN2(XXCR,IMAG(CC),REAL(CC))
	    CC = XC(NB)
	    CALL JLP_ATAN2(XXC,IMAG(CC),REAL(CC))
	    X(NB,0) = BMASK(NB)*(XXCR - XXC)
10	  CONTINUE
 
C Projection of this error X(.,0) onto E^+
	  CALL PROJ_EPLUS(NBETA,0,0)
 
C ECART ANGULAIRE CORRESPONDANT AUX U CROISSANTS EN NORME LE LONG
C D'UN RAYON; FINALEMENT EXPRIME EN DEGRE
	  RADDEG = 180./3.141592654
 
	  DO 11 IRS = 1,IR
	    NBS = NBCOUV(IRS,0)
	    RR = BMASK(NBS)*RADDEG*X(NBS,0)
	    WRITE(6,72) IRS,NBS,RR,SIGM(NBS)
	    WRITE(2,72) IRS,NBS,RR,SIGM(NBS)
72	    FORMAT(' #',I2,' Spect. index ',I4,' Phase err (deg) ',
     1      F8.3,' Sigma ',F8.3)
11	  CONTINUE
 
	ENDIF
 
999	WRITE(6,*) ' Log File in "inv_bispec_30.log"'
C
C Output the errors of final bispectrum:
        IF(IOPT.LT.-1045)THEN
	  ERRNAME='bisp_error2.dat'
          NBX=NBETA
          NBY=NBETA/2
          MEMSIZE=NBX*NBY*4
          CALL JLP_GETVM(PNTR_ARRAY,MEMSIZE)
	  CALL ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,MADRID(PNTR_ARRAY),
     1   NBX,NBY)
	  ERRNAME='bisperr2'
	  COMMENTS='Quadratic errors'
          CALL JLP_WRITEIMAG(MADRID(PNTR_ARRAY),NBX,NBY,NBX,
     1  ERRNAME,COMMENTS)
          CALL JLP_FREEVM(PNTR_ARRAY,MEMSIZE)
        ENDIF

C Eigenvalues:
        I=1234
        IF(I.EQ.1234)THEN
	CALL EIGEN_VALUE1(NBETA,NGAMMA,XLAMBDA1)
	CALL EIGEN_VALUE2(NBETA,NGAMMA,XLAMBDA1,XLAMBDA2)
	XLAMBDA2=MAX(XLAMBDA2,1.E-12)
        WRITE(6,48) XLAMBDA1,XLAMBDA2,XLAMBDA1/XLAMBDA2
        WRITE(2,48) XLAMBDA1,XLAMBDA2,XLAMBDA1/XLAMBDA2
48	FORMAT(' Largest eigen value: ',G12.5,
     1  /,' Smallest eigen value: ',G12.5,
     1  /,' Conditionning number: ',G12.5)
        ENDIF

	CLOSE(2)
	CALL JLP_END
	STOP
998	WRITE(6,*) ' Fatal error opening inv_bispec_30.log'
	CALL JLP_END
	STOP
	END
C*******************************************************************
C CREYCE_SIMU : CRE LES TERMES DE PHASE DU SPECTRE ET DU BISPECTRE
C Output in real/imaginary (from artificial bispectrum)
C Centers the Fourier transform
C*******************************************************************
	SUBROUTINE CREYCE_SIMU(IR,CTE,NBETA,NGAMMA,FNAME)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
 
	INTEGER KLM,NBCOUV,IXY
	REAL RE,IM,WEIGHT,RO,BMASK,COS1,SIN1
	COMPLEX XCR,XC,YCE
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)
	COMMON /C6/XCR(NULL:NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
	
C Reading RE and IM :
	WRITE(6,37)
37	FORMAT(' REAL PART OF THE FFT OF THE IMAGE',
     1	' (centered in the frame)')
	FNAME=' '
	CALL JLP_READIMAG(RE,NX,NY,IDIM,FNAME,COMMENTS)
	WRITE(2,54) NAME(1:14),COMMENTS(1:30)
54	FORMAT(' Real part of the fft: ',A14,' comments: ',A30)
 
C Central position: 
	IXC = (NX/2) + 1
	IYC = (NY/2) + 1
 
	WRITE(6,38)
38	FORMAT(' IMAGINARY PART OF THE FFT OF THE IMAGE',
     1	' (centered in the frame)')
	NAME=' '
	CALL JLP_READIMAG(IM,NX,NY,IDIM,NAME,COMMENTS)
	WRITE(2,56) NAME(1:14),COMMENTS(1:30)
56	FORMAT(' Imag. part of the fft: ',A14,' comments: ',A30)
 
C Mise en place du module et de la phase dans RO et XC
	INULL=0
	DO 1 NB = 0,NBETA
	  IIX = IXC + IXY(1,NB)
	  IIY = IYC + IXY(2,NB)
	  XR = RE(IIX,IIY)
	  XI = IM(IIX,IIY)
C Modulus:
	  XM = SQRT(XR*XR + XI*XI)
	  RO(NB) = XM
C Phase factor:
	  IF(XM.EQ.0.)THEN
	    XC(NB) = (1.,0.)
            INULL=INULL+1
C	    WRITE(6,78) NB
78	    FORMAT(' CREYCE_SIMU/Warning: XC(',I5,') is null!') 
	  ELSE
	    XC(NB) = CMPLX(XR/XM,XI/XM)
	  ENDIF
1	CONTINUE
	    IF(INULL.GT.0)THEN
	      WRITE(6,79) INULL
	      WRITE(2,79) INULL 
79	      FORMAT(' CREYCE_SIMU/Warning: Modulus is null for'
     1               ,I5,' values') 
	    ENDIF

	DO I=0,4
	   WRITE(2,*) ' RO(',I,') =',RO(I)
	   WRITE(6,*) ' RO(',I,') =',RO(I)
	   WRITE(2,*) ' XC(',I,') =',XC(I)
	   WRITE(6,*) ' XC(',I,') =',XC(I)
	END DO
 
C Mask to discard some values of the spectral list:
        DO 91 NB=0,NBETA
          BMASK(NB)=1.
C Check that RO(NB) greater than LOWER_MODULUS:
C Remember that LOWER_MODULUS can be negative...
          RO(NB)=RO(NB)/RO(0)
	  IF(RO(NB).LT.LOWER_MODULUS.OR.
     1       RO(NB).LE.0)BMASK(NB)=0.
91      CONTINUE
 
C FACTEUR DE PHASE DU SPECTRE CALE EN TRANSLATION : XC(.)
	CALL TRANS(IR,NBETA)

C Reference (since simulation)
	DO 3 NB=0,NBETA
	  XCR(NB)=XC(NB)
3	CONTINUE
 
C FACTEUR DE PHASE DU BISPECTRE : YCE
	DO 2 NG=1,NGAMMA
	  YCE(NG)=XC(KLM(1,NG))*XC(KLM(2,NG))*CONJG(XC(KLM(3,NG)))
2	CONTINUE
 
	DO I=1,5
	   WRITE(2,*) ' YCE(',I,') =',YCE(I)
	   WRITE(6,*) ' YCE(',I,') =',YCE(I)
	END DO
 
C Computing the weights:
	IF(CTE.GT.0)THEN
	   CALL BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSEIF(CTE.LT.0)THEN
	   CALL BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSE
	   WRITE(6,34)
	   WRITE(2,34)
34         FORMAT(' Weights set to unity, and then normalized')
	   DO 96 NG=1,NGAMMA
	     WEIGHT(NG)=1.
	     DO 95 KK=1,3
	       K = KLM(KK,NG)
	       IF(BMASK(K).EQ.0.)WEIGHT(NG)=0.
95           CONTINUE
96         CONTINUE
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ENDIF
 
C Perturbation of YCE
	DEGRAD = 3.141592654/180.
	CTE1 = CTE * DEGRAD /FLOAT(IR**2)
 
	DO 4 NG=1,NGAMMA
	  IS2=0
	    DO 6 KK=1,3
	       K = KLM(KK,NG)
 	       IS2 = IS2 + IXY(1,K)**2 + IXY(2,K)**2
6	    CONTINUE
 
C Random generation (Gaussian law, (1.,0.)) of DGAMMA = DELTA GAMMA
          CALL JLP_RANDOM_GAUSS(WORK)
	  DGAMMA = WORK*CTE1*FLOAT(IS2)
          CALL JLP_COS(COS1,DGAMMA)
          CALL JLP_SIN(SIN1,DGAMMA)
 	  YCE(NG) = YCE(NG)*CMPLX(COS1,SIN1)
4	CONTINUE
 
	RETURN
	END
C*******************************************************************
C TRANS.FOR  fait le calage en translation du facteur de phase spectral
C XC(.) en entree, XC(.) en sortie
C*******************************************************************
	SUBROUTINE TRANS(IR,NBETA)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	INTEGER NBCOUV,IXY
	REAL RO,BMASK
	COMPLEX XC,YCE
	COMPLEX TX(MIRMAX:IRMAX),TY(NULL:IRMAX),CX,CCX,CY,CCY
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
 
	CX=CONJG(XC(1))
	CY=CONJG(XC(2))
 
C TABLEAU DE LA FORME LINEAIRE CONCERNE DU NOYAU
	CCX=(1.,0.)
	CCY=(1.,0.)
	TX(0)=(1.,0.)
	TY(0)=(1.,0.)
 
	DO 1 I=1,IR
C
	 CCX=CCX*CX
	 TX(I)=CCX
	 TX(-I)=CONJG(CCX)
C
	 CCY=CCY*CY
	 TY(I)=CCY
1	CONTINUE
 
C CALAGE EN TRANSLATION
	DO 2 NB=0,NBETA
	    XC(NB) = XC(NB) * TX(IXY(1,NB)) * TY(IXY(2,NB))
            XC(NB)=XC(NB)*BMASK(NB)
2	CONTINUE
 
	RETURN
	END
C*******************************************************************
C RECURSIVE: Initial solution of the phasor.
C                      YCE(.)  to  XC(.)
C
C exp(i*YCE(NG)) = exp (i * phase_K) * exp (i * phase_L) * exp (i * phase_M)
C with L=NB or M=NB
C*******************************************************************
	SUBROUTINE RECURSIVE(SIGM,NBETA,NGAMMA,SIGMA_NULL)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE,CC1
C Double precision is necessary for large sums !!!
        REAL*8 SUMXR1,SUMXI1,SUMSQR,SUMSQI,SIGR,SIGI,SUMWEIGHTS
        REAL*8 RNORM
        REAL*4 XR1,XI1,XW,SIGMA_NULL
	REAL RO,BMASK,WEIGHT,SIGM(*)
	INTEGER K,L,M,NG1,KLM,INDX,NGT,IFERMAX
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C11/IFERMAX
 
	XC(0)=(1.,0.)
	XC(1)=(1.,0.)
	XC(2)=(1.,0.)
	NG1=1
	INULL=0
 
	DO 1 NB=3,NBETA
          SUMXI1=0.
          SUMXR1=0.
          SUMSQR=0.
          SUMSQI=0.
          SUMWEIGHTS=0.
          NVAL = 0
	  NG2=NGT(NB)
	  INDX=0
C
	  DO 2 NG=NG1,NG2
	    INDX=INDX+1
	     IF (INDX.LE.IFERMAX)THEN
	       K=KLM(1,NG)
	       L=KLM(2,NG)
	       M=KLM(3,NG)
	       XW=SQRT(WEIGHT(NG))
C Weight should be strictly positive to increment NVAL
               IF(XW.GT.0.) THEN
               NVAL = NVAL + 1
               SUMWEIGHTS=SUMWEIGHTS+XW
	        IF (M.EQ.NB) THEN
C CAS 1 (M = NB)
	          CC1=XC(K)*XC(L)*CONJG(YCE(NG))
	        ELSE
C CAS 2 (L = NB)
	          CC1=YCE(NG)*CONJG(XC(K))*XC(M)
	        ENDIF
                XR1=REAL(CC1)
                XI1=IMAG(CC1)
                SUMSQR=SUMSQR+XR1*XR1*XW
                SUMSQI=SUMSQI+XI1*XI1*XW
                SUMXR1=SUMXR1+XR1*XW
                SUMXI1=SUMXI1+XI1*XW
C End of XW > 0:
               ENDIF
	     ENDIF
2	  CONTINUE
 
C NORMALISATION DU TERME DE PHASE
	  RNORM=SQRT(SUMXR1*SUMXR1+SUMXI1*SUMXI1)
	  IF(RNORM.GT.0)THEN
	    XC(NB)=CMPLX(SUMXR1/RNORM,SUMXI1/RNORM)

C Sigma should have been computed with more than one value...
             IF(NVAL.GT.1) THEN
C SIGR**2 = Sum of (weight_i *(X_i - mean)**2) / Sum of (weight_i) 
C or Sum of (weight_i X_i**2)/Sum of weight_i - mean**2
C Here the weighted mean of real(CC) is SUMXR1/SUMWEIGHTS:
               SIGR=SUMSQR/SUMWEIGHTS-(SUMXR1/SUMWEIGHTS)**2
               SIGI=SUMSQI/SUMWEIGHTS-(SUMXI1/SUMWEIGHTS)**2
C Note that double precision variables are needed for large sums !!!
               SIGM(NB)=SQRT(SIGI+SIGR)
             ELSE
               SIGM(NB)=SIGMA_NULL
             ENDIF
	  ELSE
C The recursive process has been interrupted:
	    INULL=INULL+1
	    XC(NB)=(1.,0.)
C Maximum sigma is one, but I set it to SIGMA_NULL since we assume that
C phase is not important if modulus is too small
            SIGM(NB)=SIGMA_NULL
C Allows breaks only if in the list of null modulus:
	    IF(BMASK(NB).NE.0.)THEN
C	      WRITE(6,78) NB
78	      FORMAT(' RECURSIVE/Warning: Unexpected null modulus,'
     1     ' XC(',I5,') undetermined!') 
	      WRITE(2,82) NB,SUMXR1,SUMXI1,SUMWEIGHTS,SUMSQR,SUMSQI
	      WRITE(6,82) NB,SUMXR1,SUMXI1,SUMWEIGHTS,SUMSQR,SUMSQI
82	   FORMAT(' RECURSIVE/Warning: Recursive process interrupts at',I5,
     1      /,' SUMXR1,SUMXI1,SUMWEIGHTS,SUMSQR,SUMSQI: ',5(G12.5,1X),/)
C	      STOP
	    ENDIF
	  ENDIF
	  NG1=NG2+1
1	CONTINUE
 
	  IF(INULL.GT.0)THEN
	    WRITE(2,79) INULL,SIGMA_NULL,1./SIGMA_NULL 
	    WRITE(6,79) INULL,SIGMA_NULL,1./SIGMA_NULL 
79	    FORMAT(' RECURSIVE/Warning: Null modulus (i.e., spectral'
     1     ,'list is reduced) for',I5,' values',/,
     1      ' WARNING: In that case SIGMA was set to ',G12.5,
     1      '  (and SNR to ',G12.5,')') 
	  ENDIF
 
         RETURN
	END
C*******************************************************************
C WEIGHTS_SIGM2: Weights according to Sigma of recursive solution
C Specially designed with weights according to snr_bispectrum
C*******************************************************************
	SUBROUTINE WEIGHTS_SIGM2(SIGM,NBETA,NGAMMA,SIG_MAX,
     1      LOWER_MODULUS)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE,CC,CW
        REAL*8 SUMSQR,SUMSQI
        REAL*4 LOWER_MODULUS
	REAL RO,BMASK,WEIGHT,SIGR,SIGI,SIGM(*),XNUMB,SIG_MAX
	INTEGER K,L,M,NG1,KLM,NGT,N_OUT,IFERMAX
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C11/IFERMAX
 
	NG1=1
        N_OUT=0
        SUM_WEIGHTS=0.
        SIGM(3)=5.E-2
        SIGM(4)=5.E-2
 
        OPEN(3,FILE='sigma.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
	DO 1 NB=3,NBETA
	  NG2=NGT(NB)
          WRITE(3,*) NB,SIGM(NB)
C Caution to avoid SIGM = 0, and problems when computing 1/SIGM in SNR map:
          SIGM(NB)=MAX(5.E-2,SIGM(NB))
C Troncation to reduce the range of weights:
          IF(SIGM(NB).GE.SIG_MAX.OR.BMASK(NB).EQ.0.) THEN
              N_OUT=N_OUT+1
              BMASK(NB)=0.
	      DO 3 NG=NG1,NG2
                WEIGHT(NG)=0.
3             CONTINUE
            ELSE
	      DO 2 NG=NG1,NG2
C Previous weight divided by this sigma: 
C I have tried without dividing: gives worse results 
C I have tried with **1, **0.8, **0.2 instead: gives worse results 
C JLP93
               WEIGHT(NG)=WEIGHT(NG)/SQRT(SIGM(NB))
C Normalisation to the number of bispectral terms:
C I have tried without dividing: gives better results 
CJLP93               WEIGHT(NG)=WEIGHT(NG)/FLOAT(NG2-NG1+1)
2             CONTINUE
          ENDIF
	  NG1=NG2+1
1	CONTINUE
 
        CLOSE(3)

C Normalisation of the weights:
	 CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)

        IF(N_OUT.NE.0)THEN
	 WRITE(2,79) N_OUT,LOWER_MODULUS,SIG_MAX
	 WRITE(6,79) N_OUT,LOWER_MODULUS,SIG_MAX
79	 FORMAT(' SIGM_WEIGHTS/Warning: ',I5,' discarded BETA values',
     1          ' because of modulus less than',G12.5,
     1          ' or sigma greater than:',G10.3) 
        ENDIF

	RETURN
	END
C*******************************************************************
C WEIGHTS_SIGM1: Weights according to Sigma of recursive solution
C specially designed for unity weights
C*******************************************************************
	SUBROUTINE WEIGHTS_SIGM1(SIGM,NBETA,NGAMMA,SIG_MAX,
     1      LOWER_MODULUS)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE,CC,CW
        REAL*8 SUMSQR,SUMSQI
        REAL*4 LOWER_MODULUS
	REAL RO,BMASK,WEIGHT,SIGR,SIGI,SIGM(*),XNUMB,SIG_MAX
	INTEGER K,L,M,NG1,KLM,NGT,N_OUT,IFERMAX
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C11/IFERMAX
 
	NG1=1
        N_OUT=0
        SUM_WEIGHTS=0.
        SIGM(3)=5.E-2
        SIGM(4)=5.E-2
 
        OPEN(3,FILE='sigma.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
	DO 1 NB=3,NBETA
	  NG2=NGT(NB)
          WRITE(3,*) NB,SIGM(NB)
C Caution to avoid SIGM = 0, and problems when computing 1/SIGM in SNR map:
          SIGM(NB)=MAX(5.E-2,SIGM(NB))
C Troncation to reduce the range of weights:
          IF(SIGM(NB).GE.SIG_MAX.OR.BMASK(NB).EQ.0.) THEN
              N_OUT=N_OUT+1
              BMASK(NB)=0.
	      DO 3 NG=NG1,NG2
                WEIGHT(NG)=0.
3             CONTINUE
            ELSE
	      DO 2 NG=NG1,NG2
C JLP93
C (Unity weights: gives worse results when not dividing by SIGM)
               WEIGHT(NG)=WEIGHT(NG)/SIGM(NB)
C Normalisation to the number of bispectral terms:
C (Unity weights: gives worse results when dividing by nber_terms**2 only)
               WEIGHT(NG)=WEIGHT(NG)/FLOAT(NG2-NG1+1)**2
2             CONTINUE
          ENDIF
	  NG1=NG2+1
1	CONTINUE
 
        CLOSE(3)

C Normalisation of the weights:
	 CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)

        IF(N_OUT.NE.0)THEN
	 WRITE(2,79) N_OUT,LOWER_MODULUS,SIG_MAX
	 WRITE(6,79) N_OUT,LOWER_MODULUS,SIG_MAX
79	 FORMAT(' SIGM_WEIGHTS/Warning: ',I5,' discarded BETA values',
     1          ' because of modulus less than',G12.5,
     1          ' or sigma greater than:',G10.3) 
        ENDIF

	RETURN
	END
C*******************************************************************
C EPLUS projects ALPHA_0 onto E^+
C*******************************************************************
	SUBROUTINE EPLUS(NBETA)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)
 
	COMPLEX XC,CC,YCE
	REAL RO,X,BMASK,COS1,SIN1,XXC
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C4/X(NULL:NBMAX,NULL:3)
 
C VECTEURS DE BASE DU NOYAU (in COMMON/C7)
	CALL NOYAU(NBETA)
 
C ALPHA_0 EN X(.,0)
	DO 1 NB=0,NBETA
	    CC=XC(NB)
	    CALL JLP_ATAN2(XXC,IMAG(CC), REAL(CC) )
	    X(NB,0) = BMASK(NB)*XXC
1	CONTINUE
 
C ALPHA_0^+  EN X(.,0)
	CALL PROJ_EPLUS(NBETA,0,0)
 
C FACTEUR DE PHASE INITIAL CORRESPONDANT
	DO 2 NB=0,NBETA
	  ALPHA=X(NB,0)
          CALL JLP_COS(COS1,ALPHA)
          CALL JLP_SIN(SIN1,ALPHA)
          COS1=BMASK(NB)*COS1
          SIN1=BMASK(NB)*SIN1
	  XC(NB) = CMPLX(COS1,SIN1)
2	CONTINUE
 
	RETURN
	END
C*******************************************************************
C NOYAU:
C DEFINIT LES VECTEURS DE BASE DU NOYAU
C*******************************************************************
	SUBROUTINE NOYAU(NBETA)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	INTEGER NBCOUV,IXY
	COMPLEX XC,YCE
	REAL VECXY,XNORM,RO,BMASK
	REAL*8 SUM,SUM2
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C7/VECXY(NULL:NBMAX,2)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
 
C The base vectors are formed with the coordinates of the uv vectors. 
	DO 1 NB=0,NBETA
	    VECXY(NB,1)=BMASK(NB)*IXY(1,NB)
	    VECXY(NB,2)=BMASK(NB)*IXY(2,NB)
1	CONTINUE

C Normalization of the first vector:
	CALL NORMALIZE_L2(VECXY(1,1),NBETA,XNORM)

C Orthogonalization:
	SUM=0.
	DO NB=0,NBETA
        SUM=SUM+VECXY(NB,1)*VECXY(NB,2)
	END DO

	DO 2 NB=0,NBETA
	    VECXY(NB,2)=VECXY(NB,2)-SUM*VECXY(NB,1)
2	CONTINUE

C Normalization of the second vector:
	CALL NORMALIZE_L2(VECXY(1,2),NBETA,XNORM)
	
	RETURN
	END
C*******************************************************************
C PROJ_EPLUS
C Projects X(.,N1) onto E^+ and stores result in X(.,N2)
C*******************************************************************
	SUBROUTINE PROJ_EPLUS(NBETA,N1,N2)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	INTEGER NBETA,N1,N2
	REAL X,VECXY,DELTA
	REAL*8 XX,YY
	COMMON /C4/X(NULL:NBMAX,NULL:3)
	COMMON /C7/VECXY(NULL:NBMAX,2)
 
C Scalar product of X(.,N1) with the base vectors of the kernel
	XX=0.
	YY=0.
 
	DO 1 NB=1,NBETA
	  XX=XX+X(NB,N1)*VECXY(NB,1)
	  YY=YY+X(NB,N1)*VECXY(NB,2)
1	CONTINUE
 
	DELTA = SQRT(XX**2 + YY**2)
C JLP95: I remove the comment
	WRITE(6,*) ' Distance to E+: DELTA =',DELTA
C	WRITE(2,*) ' Distance to E+: DELTA =',DELTA
 
C Projection onto E^+ and stored in X(.,N2)
	DO NB=1,NBETA
	    X(NB,N2) = X(NB,N1) - XX*VECXY(NB,1) - YY*VECXY(NB,2)
        END DO
 
	RETURN
	END
C*******************************************************************
C ATRANS compute the second member of the system to be solved.
C This member is stored in X(.,2) which is the initial
C residual R_0
C
C  R_0 = AT ( PSI - A PHI_0)
C Here the residual is IMAG(exper.bispectrum * CONJ(estimated bispectrum))) :
C*******************************************************************
	SUBROUTINE ATRANS(NBETA,NGAMMA,QMOY)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE,C1,C2,C3
	INTEGER KLM,NGT,NBETA,NGAMMA
	REAL*8 Q2
	REAL RO,QMOY,X,WEIGHT,BMASK
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C4/X(NULL:NBMAX,NULL:3)
	COMMON /C5/NGT(NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C11/IFERMAX
 
C Initialization to zero:
	DO 1 NB=1,NBETA
	  X(NB,2)=0.
1	CONTINUE
 
C Weighted quadratic measure:
	Q2=0.
 
	NG1=1
	DO 3 NB=3,NBETA
	  NG2=NGT(NB)
	  INDEX=0
	    DO 2 NG=NG1,NG2
	      INDEX=INDEX+1
	        IF(INDEX.LE.IFERMAX.AND.BMASK(NB).NE.0.)THEN
	          K=KLM(1,NG)
	          L=KLM(2,NG)
	          M=KLM(3,NG)
C Computed bispectrum (from the spectrum of the previous iteration):
	          C1=XC(K)*XC(L)*CONJG(XC(M))
C Experimental bispectrum:
	          C2=YCE(NG)
C Corresponding error (between experimental bispectrum and computed
C bispectrum from previous iteration):
	          C3=C2 - C1
C JLP93: Put WEIGHT=g**2 for AT A and for A (Use also g**2 for the
C estimation of convergence...)
	          Q2=Q2+(REAL(C3)**2+IMAG(C3)**2)*WEIGHT(NG)
C Residual is IMAG(YCE(NG)*CONJ(C1)) :
C JLP93: Put WEIGHT=g**2 for AT A and for A: 
	          YY=IMAG(C2*CONJG(C1))*WEIGHT(NG)
C Computing AT [YY] by adding the contribution to X(.,2):
	          X(K,2)=X(K,2)+YY
	          X(L,2)=X(L,2)+YY
	          X(M,2)=X(M,2)-YY
	        ENDIF
2	     CONTINUE
	   NG1=NG2+1
3	CONTINUE
 
C Mean quadratic error:
C It is a weighted error, and the sum of the weights is 1, so:
	QMOY = SQRT(Q2)
C	WRITE(6,*) ' Mean weighted error on the bispectrum :',QMOY
 
	RETURN
	END
C*******************************************************************
C ATRANS_A compute X(.,N2) = [AT A] X(.,N1)
C*******************************************************************
	SUBROUTINE ATRANS_A(NBETA,NGAMMA,N1,N2)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL X,WEIGHT,RO,BMASK
        COMPLEX XC,YCE
	INTEGER KLM,NGT,IFERMAX
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C4/X(NULL:NBMAX,NULL:3)
	COMMON /C5/NGT(NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C11/IFERMAX
 
C Initialization to zero:
	DO 1 NB=1,NBETA
	  X(NB,N2)=0.
1	CONTINUE
 
	NG1=1
	DO 3 NB = 3,NBETA
	  NG2=NGT(NB)
  	  INDEX=0
	    DO 2 NG=NG1,NG2
	      INDEX=INDEX+1
	      IF (INDEX.LE.IFERMAX.AND.BMASK(NB).NE.0.)THEN
	        K = KLM(1,NG)
	        L = KLM(2,NG)
	        M = KLM(3,NG)
C JLP93: Put WEIGHT=g**2 for AT A
	        YY = (X(K,N1)+X(L,N1)-X(M,N1))*WEIGHT(NG)
C Computing AT [YY] by adding the contribution to X(.,2):
	        X(K,N2) = X(K,N2) + YY
	        X(L,N2) = X(L,N2) + YY
	        X(M,N2) = X(M,N2) - YY
	      ENDIF
2	    CONTINUE
	  NG1=NG2+1
3	CONTINUE
 
	RETURN
	END
C*******************************************************************
C CGRADIENT 
C Resolution of the linear system: 
C                  [ATwA] X(.,0)  = X(.,2)
C with conjugate gradients
C X(.,0) is the unknown (and then solution) PHI 
C X(.,1) is the direction D_n
C X(.,2) is the second member of the equation, and then the residual
C X(.,3) is ATwA D_n   (called Z_n)
C*******************************************************************
	SUBROUTINE CGRADIENT(NBETA,NGAMMA,IT,QMOY)
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL*8 SS,R2_RN,R2_RNPLUS1,R2_PHIPLUS1
	REAL X,OMEGA_N,SUP,WORK
	REAL GAMMA_N,QMOY,RO,BMASK
        COMPLEX YCE,XC
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C4/X(NULL:NBMAX,NULL:3)
C Good problem, with 1.E-08 implies 4 iterations,
C Badly conditionned problem with Nc=25 implies around 40 iterations, so:
	ITMAX=50
 
C STEP 0
C The starting solution PHI_0 is null : X(.,0)=0.
 
C Compute the initial residual R_0 and store it in X(.,2)
	CALL ATRANS(NBETA,NGAMMA,QMOY)
C Note that X(NB,2) is null when BMASK(NB) is null.

C Compute R2_RN : square of the norm of R_0
C  R2_RN is the square norm of R_N
	R2_RN=0.
	DO 1 NB=1,NBETA
C The initial solution is set to 0.	
	     X(NB,0) = 0.
C The first direction D_0 = R_0 is copied to X(.,1)
	     X(NB,1) = X(NB,2)
 	     R2_RN = R2_RN + X(NB,2)*X(NB,2)
1	CONTINUE
          IF(R2_RN.LT.1.E-15)THEN
            WRITE(6,39)
39          FORMAT(' CGRADIENT/Error: Square norm of Residual_{N}=0 !')
            RETURN
          ENDIF
 
	DO 2 IT=1,ITMAX
 
C STEP 1
C Compute Z_N = [AT A] D_N and store it in X(.,3)
	  CALL ATRANS_A(NBETA,NGAMMA,1,3)

C Compute OMEGA_N = R2_RN / (D_N scalar Z_N)
C Warning: SS is very small
	  SS=0.
	  DO 3 NB=1,NBETA
	      SS = SS + X(NB,1)*X(NB,3)*BMASK(NB)
3	  CONTINUE
	  OMEGA_N=R2_RN/SS

C Compute the residual and the next value of PHI
C  R_{N+1} put to X(.,2) and  PHI_{N+1} put to X(.,0) ;
C  R2_RNPLUS1 is the square norm of R_{N+1}
	  R2_RNPLUS1=0.
	  R2_PHIPLUS1=0.
	  DO 4 NB=1,NBETA
C R_[N+1] = R_N - OMEGA_N * Z_N
	     X(NB,2) = X(NB,2) - OMEGA_N*X(NB,3)
C PHI_[N+1] = PHI_N + OMEGA_N * D_N
 	     X(NB,0) = X(NB,0) + OMEGA_N*X(NB,1)
	     R2_RNPLUS1 = R2_RNPLUS1+X(NB,2)*X(NB,2)*BMASK(NB)
	     R2_PHIPLUS1=R2_PHIPLUS1+X(NB,0)*X(NB,0)*BMASK(NB)
4	  CONTINUE
 
C STEP2 : Exit test
          IF(R2_PHIPLUS1.LT.1.E-15)THEN
            WRITE(6,29)
29          FORMAT(' CGRADIENT/Error: Square norm of Delta_{PHI+1}=0 !')
            RETURN
          ENDIF
          WORK=R2_RNPLUS1/R2_PHIPLUS1
	  IF (WORK.LT.1.E-5) THEN
	    RETURN
	  ENDIF
 
C STEP 3
C GAMMA_N = || R_[N+1] ||**2 / || R_N ||**2
	  GAMMA_N=R2_RNPLUS1/R2_RN
C For next step || R_[N+1] ||**2 becomes || R_N ||**2
	  R2_RN=R2_RNPLUS1
C Compute next conjugate direction D_{N+1} which is stored in X(.,1)
CJLP91	  SUM=0
	  DO 5 NB=1,NBETA
	    X(NB,1) = X(NB,2) + GAMMA_N*X(NB,1)
C JLP91 Test if D_N and D_N+1 are orthogonal directions:
CJLP91	    SUM = SUM + X(NB,1)*X(NB,3)
5	  CONTINUE
C JLP91 Test if ATwA D_N and D_N+1 are orthogonal directions:
CJLP91     WRITE(6,26) SUM
C26	  FORMAT(' D_n+1 scalar D_n is :',G12.5)

C End of main loop
2	CONTINUE
 
C STEP 4
	WRITE(6,65) ITMAX,WORK
	WRITE(2,65) ITMAX,WORK
65	FORMAT(' Internal loop: not accurate enough, IT =',I4,
     1  ' Test =',G12.5)
	RETURN
	END
C*******************************************************************
C ERROR_SIMU
C Compute the absolute errors for simulations
C (since model is known and stored in XCR)
C
C E is the mean of the errors of the phase terms 
C       weighted by the value of MODSQ**2:
C   or the mean of the errors of the full term (amplitude and phase)
C E = sqrt( Sum of (RO**2 * |XCR-XC|**2) / Sum of RO**2)
C whereas EP is raw:
C EP = sqrt( Sum of |XCR-XC|**2 / NBETA-2)
C*******************************************************************
	SUBROUTINE ERROR_SIMU(NBETA,E,EP,ERRNAME)
 
C  POUR IR=25
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE,XCR,CC
	REAL*8 SR2
	REAL RO,BMASK
        INTEGER*4 INDX
	CHARACTER ERRNAME*40
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C6/XCR(NULL:NBMAX)
 
C ERREURS GLOBALES
C E  : ERREUR QUADRATIQUE MOYENNE
C EP : ERREUR QUADRATIQUE MOYENNE SUR LA PHASE
	E = 0.
	SR2 = 0.
	EP = 0.
        INDX = 0.
 
	OPEN(9,FILE=ERRNAME,STATUS='UNKNOWN')
	WRITE(9,13)
13      FORMAT(' 0 0. 0. 0. Spectral list, Phase error,',
     1  ' Full quad. error, Cumul. mean phase error')

	DO 1 NB=3,NBETA
         IF(BMASK(NB).NE.0.)THEN
          INDX=INDX+1
          R2 = RO(NB)**2
	  SR2 = SR2 + R2
	  CC = XCR(NB) - XC(NB)
	  RR2 = REAL(CC)**2 + IMAG(CC)**2
	  E = E + R2*RR2
	  EP = EP + RR2
	  WRITE(9,*) NB,SQRT(RR2),SQRT(R2*RR2),SQRT(EP/FLOAT(NB-2))
         ENDIF
1	CONTINUE
 
	IF(SR2.EQ.0)THEN
	 E = 10000.
	 WRITE(6,59)
	 WRITE(2,59)
59	 FORMAT(' Warning: Problem in ERROR, SR2=0 !!!')
	ELSE
C E = sqrt( Sum of (RO**2 * |XCR-XC|**2) / Sum of RO**2)
	 E = SQRT(E/SR2)	
	ENDIF
C EP = sqrt( Sum of |XCR-XC|**2 / NBETA-2)
	EP = SQRT(EP/FLOAT(INDX))
 
	CLOSE(9)

	RETURN
	END
 
C*******************************************************************
C SORTIE.FOR : CRE LES FICHIERS DE SORTIE RES ET IMS
C*******************************************************************
	SUBROUTINE SORTIE(NBETA,ISORTIE,FNAME)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
	COMPLEX XC,YCE,CC
	REAL RE,IM,RO,BMASK
	INTEGER NBCOUV,IXY,NX,NY
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
 
C Central position: 
	IXC = (NX/2) + 1
	IYC = (NY/2) + 1
 
C MISE A ZERO DE RE ET IM
	DO 1 J=1,NY
	  DO 2 I=1,NX
	    RE(I,J) = 0.
	    IM(I,J) = 0.
2	  CONTINUE
1	CONTINUE
 
C FORMATION DES SORTIES DE RE ET IM
	DO 3 NB = 0,NBETA
	  XM = BMASK(NB)*RO(NB)
	  CC = XC(NB)
	  XR = XM*REAL(CC)
	  XI = XM*IMAG(CC)
	  IX=IXY(1,NB)
	  IY=IXY(2,NB)
	  IIX = IXC + IX
	  IIY = IYC + IY
	  RE(IIX,IIY) = XR
	  IM(IIX,IIY) = XI
	  IIX = IXC - IX
	  IIY = IYC - IY
	  RE(IIX,IIY) = XR
	  IM(IIX,IIY) = -XI
3	CONTINUE
 
C ECRITURE DE RE ET IM
C In rei and imi if isortie=1
C In ref and imf if isortie=1
 
	IF (ISORTIE.EQ.1)THEN
C SORTIE INITIALE
	  NAME='rei'
	  COMMENTS=' FFT (rei) of : '//FNAME(1:20)
	  PRINT *,COMMENTS(1:80)
	  CALL JLP_WRITEIMAG(RE,NX,NY,IDIM,NAME,COMMENTS)
C
	  NAME='imi'
	  COMMENTS=' FFT (imi) of : '//FNAME(1:20)
	  PRINT *,COMMENTS(1:80)
	  CALL JLP_WRITEIMAG(IM,NX,NY,IDIM,NAME,COMMENTS)
	ELSE
C SORTIE FINALE
	  NAME='ref'
	  COMMENTS=' FFT (ref) of : '//FNAME(1:20)
	  PRINT *,COMMENTS(1:80)
	  CALL JLP_WRITEIMAG(RE,NX,NY,IDIM,NAME,COMMENTS)
C
	  NAME='imf'
	  COMMENTS=' FFT (imf) of : '//FNAME(1:20)
	  PRINT *,COMMENTS(1:80)
	  CALL JLP_WRITEIMAG(IM,NX,NY,IDIM,NAME,COMMENTS)
	ENDIF
C
	RETURN
	END
C*******************************************************************
C CREYCE1
C Reads a real spectrum (Square modulus only) and bispectrum
C*******************************************************************
	SUBROUTINE CREYCE1(IR,CTE,NBETA,NGAMMA,
     1  NGAMMA_MAX,FNAME,LOWER_MODULUS)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
        REAL RO,BMASK
	REAL BISP(NGMAX,3),RE,IM,WEIGHT,LOWER_MODULUS
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	INTEGER NBCOUV,IXY,KLM,NX,NY
	COMPLEX XC,YCE
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
 
C Input of the modulus:
	WRITE(6,37)
37	FORMAT(' Modulus (squared) of the FFT of the image',/,
     1	' (centered in the frame)')
	FNAME=' '
	CALL JLP_READIMAG(RE,NX,NY,IDIM,FNAME,COMMENTS)
	WRITE(2,57) FNAME(1:14),COMMENTS(1:30)
57	FORMAT(' Sq. modulus: ',A14,' Comments: ',A30)
 
C Initiallizing the real and imaginary part of the spectrum:
	DO 9 J=1,NY
          DO 8 I=1,NX
	    WORK=MAX(RE(I,J),0.)
	    RE(I,J)=SQRT(WORK)
	    IM(I,J)=0.
8          CONTINUE
9        CONTINUE
 
C Central position: 
	IXC = (NX/2) + 1
	IYC = (NY/2) + 1
 
C MISE EN PLACE DU MODULE DANS RO
	DO 1 NB = 0,NBETA
	  IIX = IXC + IXY(1,NB)
	  IIY = IYC + IXY(2,NB)
C Modulus (here that is the real part) :
	  RO(NB)=RE(IIX,IIY)
C Phase factor:
	  XC(NB)=(1.,0.)
1	CONTINUE
 
	DO I=0,4
	   WRITE(6,23) I,RO(I)
	   WRITE(2,23) I,RO(I)
23	   FORMAT(' RO(',I2,') = ',F8.5)
	END DO
 
C Troncation:
	IREM=0
C Mask to discard some values of the spectral list:
        DO 91 NB=0,NBETA
          RO(NB)=RO(NB)/RO(0)
C Check that RO(NB) greater than LOWER_MODULUS:
C Remember that LOWER_MODULUS can be negative...
	  IF(RO(NB).LT.LOWER_MODULUS.OR.
     1       RO(NB).LE.0) THEN
             BMASK(NB)=0.
	     IREM=IREM+1
          ELSE
             BMASK(NB)=1.
          ENDIF
91      CONTINUE
 
        WRITE(6,47) NBETA-IREM
        WRITE(2,47) NBETA-IREM
47	FORMAT(' Number of terms of the spectral list',
     1  ' after correction: ',I4)

C FACTEUR DE PHASE DU SPECTRE CALE EN TRANSLATION : XC(.)
	CALL TRANS(IR,NBETA)

C FACTEUR DE PHASE DU BISPECTRE : YCE
	WRITE(6,41)
41	FORMAT(' Bispectrum of the image (phase term)',
     1	' (NOT centered in the frame)')
	NAME=' '
	CALL JLP_READIMAG(BISP,NX1,NY1,NGMAX,NAME,COMMENTS)
	WRITE(2,58) NAME(1:14),COMMENTS(1:30)
58	FORMAT(' Bispectrum: ',A14,' comments: ',A30)
	IF((NY1.LT.2).OR.(NX1.NE.NGAMMA_MAX))THEN
	   WRITE(6,45)
	   WRITE(2,45)
45	   FORMAT(' CREYCE1/FATAL ERROR: Size of bispectrum inconsistent',
     1	' with IRMAX')
	   STOP
	ENDIF
C
	DO 55 NG=1,NGAMMA
	  XR = BISP(NG,1)
	  XI = BISP(NG,2)
C MODULE
	  XM = SQRT(XR*XR + XI*XI)
C Phase factor:
	  IF(XM.EQ.0.)THEN
	    YCE(NG) = (1.,0.)
	  ELSE
	    YCE(NG) = CMPLX(XR/XM,XI/XM)
	  ENDIF
55	CONTINUE
 
	DO I=1,5
	   WRITE(2,*) ' YCE(',I,') =',YCE(I)
	   WRITE(6,*) ' YCE(',I,') =',YCE(I)
	END DO
 
C Computing the weights:
	IF(CTE.GT.0)THEN
	   CALL BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSEIF(CTE.LT.0)THEN
C Version with SNR stored in 3rd line of bispectrum:
	  CALL BISP_WEIGHT2(BISP(1,3),NBETA,NGAMMA,IR,CTE)
C	  CALL BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE)
	  CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
        ELSE 
	   WRITE(6,34)
	   WRITE(2,34)
34         FORMAT(' Weights set to unity, and then normalized')
	   DO 35 NG=1,NGAMMA
	     WEIGHT(NG)=1.
	     DO 36 KK=1,3
	       K = KLM(KK,NG)
	       IF(BMASK(K).EQ.0.)WEIGHT(NG)=0.
36           CONTINUE
35        CONTINUE
	  CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ENDIF
 
	RETURN
	END
C*******************************************************************
C CREYCE2
C Reads the modulus (squared) and bispectrum derived from simulations,
C re and im (available since it is a simulation)
C*******************************************************************
	SUBROUTINE CREYCE2(IR,CTE,NBETA,NGAMMA,
     1  NGAMMA_MAX,FNAME,LOWER_MODULUS)
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
	REAL*8 ERRMOD,ERRBISPW
        REAL RO,BMASK
	REAL MODSQ(IDIM,IDIM)
	REAL BISP(NGMAX,3),RE,IM,WEIGHT,LOWER_MODULUS
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	INTEGER NBCOUV,IXY,KLM,NX,NY
	COMPLEX XC,YCE,XCR,YY1
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C6/XCR(NULL:NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
 
C Input of the modulus:
	WRITE(6,42)
42	FORMAT(' Modulus (squared) of the FFT of the image',/,
     1	' (centered in the frame)')
	FNAME=' '
	CALL JLP_READIMAG(MODSQ,NX,NY,IDIM,FNAME,COMMENTS)
	WRITE(2,57) FNAME(1:14),COMMENTS(1:30)
57	FORMAT(' Sq. modulus: ',A14,' Comments: ',A30)
 
C Reading RE and IM :
	WRITE(6,37)
37	FORMAT(' REAL PART OF THE FFT OF THE IMAGE',
     1	' (centered in the frame)')
	NAME=' '
	CALL JLP_READIMAG(RE,NX,NY,IDIM,NAME,COMMENTS)
	WRITE(2,54) NAME(1:14),COMMENTS(1:30)
54	FORMAT(' Real part of the fft: ',A14,' comments: ',A30)
 
	WRITE(6,38)
38	FORMAT(' IMAGINARY PART OF THE FFT OF THE IMAGE',
     1	' (centered in the frame)')
	NAME=' '
	CALL JLP_READIMAG(IM,NX,NY,IDIM,NAME,COMMENTS)
	WRITE(2,56) NAME(1:14),COMMENTS(1:30)
56	FORMAT(' Imag. part of the fft: ',A14,' comments: ',A30)
 
C Assume that zero frequency is at IXC,IYC: 
	IXC = (NX/2) + 1
	IYC = (NY/2) + 1
	XW1=SQRT(RE(IXC,IYC)**2+IM(IXC,IYC)**2)
	XW2=SQRT(MODSQ(IXC,IYC))
	WRITE(2,77) XW2,XW1
	WRITE(6,77) XW2,XW1
77      FORMAT(' Central value of input and computed modulus:',
     1  2(1X,G12.5))
 
	OPEN(9,FILE='err_mod.dat',STATUS='UNKNOWN')
        WRITE(9,13)
13      FORMAT(' 0 0. 1. 1. Spectral list, modulus error,'
     1   ' modulus (simulated and theoretical)') 
C MISE EN PLACE DU MODULE ET DE LA PHASE DANS RO ET XC
	INULL=0
	ERRMOD=0.
	DO 1 NB = 0,NBETA
	  IIX = IXC + IXY(1,NB)
	  IIY = IYC + IXY(2,NB)
	  XR = RE(IIX,IIY)
	  XI = IM(IIX,IIY)
	  XM = SQRT(XR*XR + XI*XI)
C Phase factor:
	  IF(XM.EQ.0.)THEN
	    XC(NB) = (1.,0.)
            INULL=INULL+1
C	    WRITE(6,78) NB
78	    FORMAT(' Warning: XC(',I5,') is null!') 
	  ELSE
	    XC(NB) = CMPLX(XR/XM,XI/XM)
	  ENDIF
C Modulus:
	  W1=MAX(0.,MODSQ(IIX,IIY))
	  RO(NB) = SQRT(W1)
	  IF(W1.EQ.0)THEN
            INULL=INULL+1
C	    WRITE(6,75) NB
75	    FORMAT(' CREYCE2/Warning: RO(',I5,') is null!') 
	  ENDIF
	  WORK=(XM/XW1-RO(NB)/XW2)**2
	  ERRMOD=ERRMOD+WORK
          WRITE(9,*) NB,SQRT(WORK),RO(NB)/XW2,XM/XW1
1	CONTINUE
	ERRMOD=SQRT(ERRMOD/FLOAT(NBETA))
	CLOSE(9)
	IF(INULL.GT.0)THEN
	  WRITE(6,79) INULL
	  WRITE(2,79) INULL
79	  FORMAT(' CREYCE2/Modulus null for',I5,' values')
	ENDIF
	WRITE(2,43) ERRMOD
	WRITE(6,43) ERRMOD
43	FORMAT(' rms error of the modulus :',1PE11.3)
 
C Output of some errors of the modulus:
	DO I=0,4
	   IIX = IXC + IXY(1,I)
	   IIY = IYC + IXY(2,I)
	   XR = RE(IIX,IIY)
	   XI = IM(IIX,IIY)
	   XM = SQRT(XR*XR + XI*XI)/XW1
	   W1 = RO(I)/XW2
	   WRITE(6,67) I,W1,XM
	   WRITE(2,67) I,W1,XM
67	   FORMAT(' Normalised RO and computed modulus:',I5,
     1     2(1X,F8.5))
	END DO
 
C Troncation:
	IREM=0
C Mask to discard some values of the spectral list:
        DO 91 NB=0,NBETA
          RO(NB)=RO(NB)/RO(0)
C Check that RO(NB) greater than LOWER_MODULUS:
C Remember that LOWER_MODULUS can be negative...
	  IF(RO(NB).LT.LOWER_MODULUS.OR.
     1       RO(NB).LE.0) THEN
             BMASK(NB)=0.
	     IREM=IREM+1
          ELSE
             BMASK(NB)=1.
          ENDIF
91      CONTINUE
 
        WRITE(6,47) NBETA-IREM
        WRITE(2,47) NBETA-IREM
47	FORMAT(' Number of terms of the spectral list',
     1  ' after correction: ',I4)

C FACTEUR DE PHASE DU SPECTRE CALE EN TRANSLATION : XC(.)
	CALL TRANS(IR,NBETA)

C Reference (since simulation)
	DO 3 NB=1,NBETA
	  XCR(NB)=XC(NB)
C Set initial solution XC to zero:
C	  XC(NB)=(1.,0.)
3	CONTINUE
 
C FACTEUR DE PHASE DU BISPECTRE : YCE
	WRITE(6,41)
41	FORMAT(' Bispectrum of the image (phase term)',
     1	' (NOT centered in the frame)')
	NAME=' '
	CALL JLP_READIMAG(BISP,NX1,NY1,NGMAX,NAME,COMMENTS)
	WRITE(2,58) NAME(1:14),COMMENTS(1:30)
58	FORMAT(' Bispectrum: ',A14,' comments: ',A30)
	IF((NY1.LT.2).OR.(NX1.NE.NGAMMA_MAX))THEN
	   WRITE(6,45)
	   WRITE(2,45)
45	   FORMAT(' CREYCE2/FATAL ERROR: Size of bispectrum inconsistent',
     1	' with IRMAX')
	   STOP
	ENDIF

C Computing the weights:
	IF(CTE.GT.0)THEN
	   CALL BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSEIF(CTE.LT.0)THEN
C Version with SNR stored in 3rd line of bispectrum:
	  CALL BISP_WEIGHT2(BISP(1,3),NBETA,NGAMMA,IR,CTE)
C	  CALL BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE)
	  CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
        ELSE 
	   WRITE(6,34)
	   WRITE(2,34)
34         FORMAT(' Weights set to unity, and then normalized')
	   DO NG=1,NGAMMA
	     WEIGHT(NG)=1.
	     DO KK=1,3
	       K = KLM(KK,NG)
	       IF(BMASK(K).EQ.0.)WEIGHT(NG)=0.
	     END DO
	   END DO
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ENDIF
 
C
	ERRBISP=0.
	ERRBISPW=0.
	INULL=0
	DO 55 NG=1,NGAMMA
	  XR = BISP(NG,1)
	  XI = BISP(NG,2)
C Modulus
	  XM = SQRT(XR*XR + XI*XI)
C Phase factor 
	  IF(XM.EQ.0.)THEN
	    YCE(NG) = (1.,0.)
            INULL=INULL+1
C	    WRITE(2,73) NG
C	    WRITE(6,73) NG
73	    FORMAT(' CREYCE2/Warning: BISPEC(',I5,') is null!') 
	  ELSE
	    YCE(NG) = CMPLX(XR/XM,XI/XM)
	  ENDIF
	  YY1=XC(KLM(1,NG))*XC(KLM(2,NG))*CONJG(XC(KLM(3,NG)))
	  WORK=(REAL(YY1)-REAL(YCE(NG)))**2
     1	+(IMAG(YY1)-IMAG(YCE(NG)))**2
	  ERRBISP=ERRBISP+WORK
C          IF(WORK.GT.3.99)THEN
C	       WRITE(6,688) WORK,KLM(1,NG),KLM(2,NG),KLM(3,NG),
C     1          NG,REAL(YCE(NG)),IMAG(YCE(NG)),
C     1         REAL(YY1),IMAG(YY1),BISP(NG,3)
C	       WRITE(2,688) WORK,KLM(1,NG),KLM(2,NG),KLM(3,NG),
C     1          NG,REAL(YCE(NG)),IMAG(YCE(NG)),
C     1         REAL(YY1),IMAG(YY1),BISP(NG,3)
C688	   FORMAT(' Big error: WORK=',F8.5,' K,L,M=',3(1X,I5),/,
C     1     ' Input, comp bisp & snr:',I5,
C     1     2(1X,F8.5),1X,2(1X,F8.5),1X,1PG10.3)
C	END DO
C          ENDIF
	  ERRBISPW=ERRBISPW+WORK*WEIGHT(NG)
55	CONTINUE
        IF(INULL.GT.0)THEN
          WRITE(2,81) INULL
          WRITE(6,81) INULL
81	  FORMAT(' CREYCE2/Warning: null bispectrum for',I5,' values')
	ENDIF
	ERRBISP=SQRT(ERRBISP/FLOAT(NGAMMA))
	ERRBISPW=SQRT(ERRBISPW)
	WRITE(2,44) ERRBISP,ERRBISPW
	WRITE(6,44) ERRBISP,ERRBISPW
44	FORMAT(' rms error of the bispectrum :',1PE12.5,/,
     1  ' weighted rms error of the bispectrum :',1PE12.5)
 
	DO I=1,5
	  YY1=XC(KLM(1,I))*XC(KLM(2,I))*CONJG(XC(KLM(3,I)))
	   WRITE(6,68) I,REAL(YCE(I)),IMAG(YCE(I)),
     1         REAL(YY1),IMAG(YY1),BISP(I,3)
	   WRITE(2,68) I,REAL(YCE(I)),IMAG(YCE(I)),
     1         REAL(YY1),IMAG(YY1),BISP(I,3)
68	   FORMAT(' Input, computed bisp & snr:',I2,
     1     2(1X,F8.5),1X,2(1X,F8.5),1X,1PG10.3)
	END DO
 
	RETURN
	END
C*****************************************************************
C BISP_WEIGHT0
C Computing the weights:
C JLP Version, not very good...
C*****************************************************************
	SUBROUTINE BISP_WEIGHT0(NGAMMA,IR,CTE)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL WEIGHT,CTE
	REAL DEGRAD,SIGB,SIG2,S22,S33,WORK
	INTEGER NBCOUV,IXY,KLM
	INTEGER NGAMMA,IR
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C8/WEIGHT(NGMAX)

	DEGRAD = 3.141592654/180.
	CTE1 = CTE*DEGRAD/FLOAT(IR**2)
 
	DO 4 NG=1,NGAMMA
	  S22=IR*IR
          S33=0.
	     DO 6 KK=1,3
	       K = KLM(KK,NG)
 	       WORK = IXY(1,K)**2 + IXY(2,K)**2
               S22 = MIN(S22,WORK)
               S33 = MAX(S33,WORK)
6	     CONTINUE

C To make the axes sensitive to the other dimension: 
C	S22=MAX(S22,0.1)
C        S33=S22*S33

C SIGB = ECART TYPE (OU STANDARD DEVIATION) DE GAMMA
	  SIGB = CTE1*S33
 
C VALEUR MOYENNE DU CARRE DE [2 SIN( (GAMMA - GAMMA BARRE)/2 ) ]
C CETTE VALEUR MOYENNE MAJORE LA PRECEDENTE
          WORK=-0.5*(SIGB**2)
	  SIG2 = 1.-EXP(WORK)
 
C WEIGHT(NG) : (do not limit the weight values ...)
	    IF(SIG2.GT.1.E-9)THEN
	      WEIGHT(NG) = 1./SIG2
	    ELSE
	      WRITE(6,46) NG
	      WRITE(2,46) NG
46	      FORMAT(' BISP_WEIGHT0/Warning: Null sigma for NG=',I8)
	      WEIGHT(NG) = 1.E9
	    ENDIF
 
4	CONTINUE
 
	WRITE(6,234)
	WRITE(2,234)
C234	FORMAT(' WEIGHT0/New version. Both axes...')
234	FORMAT(' WEIGHT0/New version. Maxi...')

	RETURN
	END
C*****************************************************************
C BISP_WEIGHT11
C Computing the weights:
C Old version (before 25-07-91) 
C
C*****************************************************************
	SUBROUTINE BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL WEIGHT,CTE,RO,BMASK
	REAL DEGRAD,SIGB,SIG2
	INTEGER NBCOUV,IXY,KLM
	INTEGER NGAMMA,IS2,IR
	COMPLEX XC,YCE
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)

	DEGRAD = 3.141592654/180.
	CTE1 = CTE*DEGRAD/FLOAT(IR**2)
 
	DO 4 NG=1,NGAMMA
	  IS2=0
	  WORK=RO(KLM(1,NG))*RO(KLM(2,NG))*RO(KLM(3,NG))
          WORK=WORK*BMASK(KLM(1,NG))*BMASK(KLM(2,NG))*BMASK(KLM(3,NG))
          IF(WORK.EQ.0)THEN
            WEIGHT(NG)=0.
          ELSE
	     DO 6 KK=1,3
	       K = KLM(KK,NG)
 	       IS2 = IS2 + IXY(1,K)**2 + IXY(2,K)**2
6	     CONTINUE
 
C SIGB = ECART TYPE (OU STANDARD DEVIATION) DE GAMMA
	  SIGB = CTE1*FLOAT(IS2)
 
C VALEUR MOYENNE DU CARRE DE [2 SIN( (GAMMA - GAMMA BARRE)/2 ) ]
C CETTE VALEUR MOYENNE MAJORE LA PRECEDENTE
          WORK=-0.5*(SIGB**2)
	  SIG2 = 1.-EXP(WORK)
 
C WEIGHT(NG) : (do not limit the weight values ...)
	    IF(SIG2.GT.1.E-9)THEN
	      WEIGHT(NG) = 1./SIG2
	    ELSE
	      WRITE(6,46) NG
	      WRITE(2,46) NG
46	      FORMAT(' BISP_WEIGHT11/Warning: Null sigma for NG=',I8)
	      WEIGHT(NG) = 1.E9
	    ENDIF
          ENDIF
 
4	CONTINUE
 
	WRITE(6,234)
	WRITE(2,234)
234	FORMAT(' WEIGHT11/Old version. ')

	RETURN
	END
C*****************************************************************
C BISP_WEIGHT2
C Computing the weights with SNR  stored in bispectrum file (line 3)
C
C*****************************************************************
	SUBROUTINE BISP_WEIGHT2(BISP_SNR,NBETA,NGAMMA,IR,CTE)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL WEIGHT,CTE,SUM,BISP_SNR(*),MEAN_SNR
	INTEGER NGAMMA,IR,NGT
	COMMON /C8/WEIGHT(NGMAX)

	SUM=0.
	DO 4 NG=1,NGAMMA
	  WEIGHT(NG) = MAX(0.,BISP_SNR(NG))
C JLP93 : keep BISP_SNR, do not put any **0.5 **1.5 or **2, or anything else
C          WEIGHT(NG)=WEIGHT(NG)
          SUM=SUM+WEIGHT(NG)
4	CONTINUE
 
        MEAN_SNR = SUM/FLOAT(NGAMMA)

	WRITE(6,234) SUM,MEAN_SNR
	WRITE(2,234) SUM,MEAN_SNR
234	FORMAT(' BISP_WEIGHT2/SNR weight Initial sum:',1PE12.5,
     1  ' Mean bisp. SNR:',1PE12.5)

	RETURN
	END
C*****************************************************************
C BISP_WEIGHT22
C Computing the weights with SNR  stored in bispectrum file (line 3)
C Same as WEIGHT2 but takes modulus into account
C (Not very good... but best results are obtained without calling
C   SIGM2 afterwards, i.e. without dividing by SIGMA or Nber_of_terms)
C
C*****************************************************************
	SUBROUTINE BISP_WEIGHT22(BISP_SNR,NBETA,NGAMMA,IR,CTE)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL WEIGHT,CTE,SUM,BISP_SNR(*),MEAN_SNR
	INTEGER NGAMMA,IR,NGT
        COMPLEX YCE,XC
        REAL RO,BMASK
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)

C JLP93: add modulus 
C First get the relative weights according to the modulus: 
	WMAX=0.
	WMIN=RO(0)
        IG1=1
	DO NB=3,NBETA
          IG2=NGT(NB)
	  DO NG=IG1,IG2
            WEIGHT(NG)=BMASK(NB)*RO(NB)
            IF(WEIGHT(NG).NE.0.)WMIN=MIN(WMIN,WEIGHT(NG))
            WMAX=MAX(WMAX,WEIGHT(NG))
          ENDDO
          IG1=IG2+1
        ENDDO

	WRANGE=WMAX/WMIN
	WRITE(2,24) WRANGE,WMIN,WMAX
	WRITE(6,24) WRANGE,WMIN,WMAX
24	FORMAT(' WEIGHT22/Modulus: Initial range ',G11.4,
     1         ' MIN, MAX:',2(1X,G11.4))

C Troncation (max range to 100):
        WMIN=WMAX/100.
	DO NG=1,NGAMMA
           WEIGHT(NG)=MAX(WEIGHT(NG),WMIN)
C           IF(MOD(NG,100).EQ.1.AND.NG.LT.1000.) PRINT *,' Weight1',WEIGHT(NG) 
        ENDDO

C Old version:
	SUM=0.
	DO 4 NG=1,NGAMMA
C Old Version:
	  WEIGHT(NG) = MAX(0.,BISP_SNR(NG))
C JLP93:
	  WEIGHT(NG) = WEIGHT(NG)*MAX(0.,BISP_SNR(NG))
          IF(MOD(NG,100).EQ.1.AND.NG.LT.1000.) PRINT *,' Weight2',WEIGHT(NG) 
C JLP93 : keep BISP_SNR, do not put any **0.5 **1.5 or **2, or anything else
C          WEIGHT(NG)=WEIGHT(NG)
          SUM=SUM+WEIGHT(NG)
4	CONTINUE
 
        MEAN_SNR = SUM/FLOAT(NGAMMA)

	WRITE(6,234) SUM,MEAN_SNR
	WRITE(2,234) SUM,MEAN_SNR
234	FORMAT(' BISP_WEIGHT22/SNR weight Initial sum:',1PG11.4,
     1  ' Mean bisp. SNR:',1PG11.4)

	RETURN
	END
C**********************************************************************
C Output bispectrum errors:
C**********************************************************************
	SUBROUTINE ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,ARRAY,NBX,NBY)

C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	INTEGER*4 KLM,NGT,IFERMAX,NBX,NBY
	REAL BISP(NGMAX,3),RO,BMASK,WEIGHT,ARRAY(NBX,*)
        REAL*8 ERROR1
	CHARACTER ERRNAME*(*) 
	COMPLEX XC,YCE,C1,C2
        REAL MOD1,MOD2,RE1,IM1,ERROR0,ERROR2
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C5/NGT(NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C11/IFERMAX
 
C Erase input array:
        DO J=1,NBY
          DO I=1,NBX
            ARRAY(I,J)=0.
          END DO
        END DO

	WRITE(6,42) ERRNAME
	WRITE(2,42) ERRNAME
42	FORMAT(' Output of bispectrum errors (unweighted) in: ',A)
 
C Opening output ASCII file:
C        OPEN(3,FILE=ERRNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL')

C        WRITE(3,52) NG,ERROR1,ERROR2,RE1,IM1
52      FORMAT('0 0. 0. 0. 0. Bisp.list, error,',
     1         ' cumul.error, error_re, error_im') 

C 
        ERROR1=0
	NG1=1
	DO 3 NB=3,NBETA
	  NG2=NGT(NB)
	  INDEX=0
	    DO 2 NG=NG1,NG2
	      INDEX=INDEX+1
	        IF(INDEX.LE.IFERMAX)THEN
	          K=KLM(1,NG)
	          L=KLM(2,NG)
	          M=KLM(3,NG)
C Computed bispectrum (from the spectrum of the previous iteration):
C minus experimental bispectrum:
	          C1=XC(K)*XC(L)*CONJG(XC(M))
	          C2=YCE(NG)
                  MOD1=REAL(C1)**2+IMAG(C1)**2
                  MOD2=REAL(C2)**2+IMAG(C2)**2
	          RE1=REAL(C1)/MOD1-REAL(C2)/MOD2
	          IM1=IMAG(C1)/MOD1-IMAG(C2)/MOD2
                  ERROR0=RE1*RE1+RE2*RE2
                  ARRAY(NB,NG-NG1+1)=SQRT(ERROR0)
                  ERROR1=ERROR1+ERROR0
                  IF(NG.GT.2)THEN
                    ERROR2 = SQRT(ERROR1/FLOAT(NG-2))
                  ELSE
                    ERROR2 = 0.
                  ENDIF
C                  WRITE(3,53) NG,SQRT(ERROR0),ERROR2,RE1,IM1
53                FORMAT(I5,1X,4(E12.4,1X))
	        ENDIF
2	     CONTINUE
	   NG1=NG2+1
3	CONTINUE

        ERROR2 = SQRT(ERROR1/FLOAT(NG-2))
        WRITE(6,68) ERROR2
        WRITE(2,68) ERROR2
68      FORMAT(/,' Final bispectrum error (rms unweighted):',G10.4,/)

C        CLOSE(3)

        RETURN
	END
C*****************************************************************
C BISP_WEIGHT1
C Computing the weights:
C Old: Takes the modulus into account 
C New: Takes a simulated modulus into acount
C
C*****************************************************************
	SUBROUTINE BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE
        REAL RO,BMASK
	REAL WEIGHT,CTE,W1,WMIN,WMAX,RANGE_MAX1,RANGE_MAX2
	INTEGER NGT,NGAMMA,IR,KLM,NBCOUV,IXY
	COMMON /C2/KLM(3,NGMAX)
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C8/WEIGHT(NGMAX)
	COMMON /C5/NGT(NBMAX)

C 20 is too small, even 40 is a bit small...
	RANGE_MAX1=800.
	RANGE_MAX2=100000.
	WRITE(2,23) RANGE_MAX1,RANGE_MAX2
	WRITE(6,23) RANGE_MAX1,RANGE_MAX2
23	FORMAT(' Weight3/Modulus, R1*R2*R3/sum*radius**2',/,
     1  ' Maximum range per column, and global: ',2(1X,F8.2))

C First get the relative weights according to the modulus: 
	WMAX=0.
	WMIN=RO(0)*RO(0)*RO(0)
	DO 203 NG=1,NGAMMA
          WORK=BMASK(KLM(1,NG))*BMASK(KLM(2,NG))*BMASK(KLM(3,NG))
	  WEIGHT(NG)=WORK*RO(KLM(1,NG))*RO(KLM(2,NG))*RO(KLM(3,NG))
          IF(WEIGHT(NG).NE.0.)WMIN=MIN(WMIN,WEIGHT(NG))
          WMAX=MAX(WMAX,WEIGHT(NG))
C          WEIGHT(NG)=SQRT(WEIGHT(NG))
C	  WEIGHT(NG)=MIN(RO(KLM(1,NG)),RO(KLM(2,NG)))
C	  WEIGHT(NG)=MIN(WEIGHT(NG),RO(KLM(3,NG)))
203     CONTINUE 

	WRANGE=WMAX/WMIN
	WRITE(2,24) WRANGE
	WRITE(6,24) WRANGE
24	FORMAT(' Initial range of the weights:',G12.5)

C Now modulation according to the column (ie NB)
	WMAX=0.
	WMIN=RO(0)*RO(0)*RO(0)
	IG1=1
	DO NB=3,NBETA
          IG2=NGT(NB)

C Troncation in the column to regularize the problem:
	  WWMIN=RO(0)*RO(0)*RO(0)
	  WWMAX=0.
	    DO NG=IG1,IG2
              IF(WEIGHT(NG).NE.0.)WWMIN=MIN(WWMIN,WEIGHT(NG))
              WWMAX=MAX(WWMAX,WEIGHT(NG))
	    END DO

	  IF(WWMAX.NE.0.)THEN

	  WWMAX=WWMIN*RANGE_MAX1
          SUM0=0.
	    DO NG=IG1,IG2
C To check (JLP91)
              IF(WEIGHT(NG).GT.WWMAX)THEN
	       WRITE(2,69) NB,WEIGHT(NG),WWMAX
	       WRITE(6,69) NB,WEIGHT(NG),WWMAX
69             FORMAT(' Column troncation, NB, wold,new',I5,2(1X,1PG12.5))
	      ENDIF
	      WEIGHT(NG)=MIN(WEIGHT(NG),WWMAX) 
	      SUM0=SUM0+WEIGHT(NG)
            END DO

C Normalizes the weights to obtain a 1/||frequency||**2 law 
C for the spectral list: 
C **4 is good, but weight range is too large for conditionning...
	  XW=4 + IXY(1,NB)**2+IXY(2,NB)**2
	  XW=SUM0*XW*XW
	    DO NG=IG1,IG2
	      WEIGHT(NG)=WEIGHT(NG)/XW
              IF(WEIGHT(NG).NE.0.)WMIN=MIN(WMIN,WEIGHT(NG))
              WMAX=MAX(WMAX,WEIGHT(NG))
            ENDDO

	  ENDIF

	  IG1=IG2+1
	ENDDO

	WRANGE=WMAX/WMIN
	WRITE(2,25) WRANGE
	WRITE(6,25) WRANGE
25	FORMAT(' Range of the weights before troncation: ',G12.5)
        PRINT *,' WMAX,WMIN',WMAX,WMIN

C Reajustment of the range:
	WMAX=RANGE_MAX2*WMIN
        PRINT *,' WMAX',WMAX
	DO NG=1,NGAMMA
	 IF(WEIGHT(NG).GT.WMAX)PRINT *,' NG,WEIGHT(NG)',WEIGHT(NG),NG
	 WEIGHT(NG)=MIN(WEIGHT(NG),WMAX)
        ENDDO

        RETURN
	END
C**********************************************************************
	SUBROUTINE NORMALIZE_L1(ARRAY,NPTS,NORM)
	REAL ARRAY(*),NORM
	REAL*8 SSUM
	INTEGER NPTS

C	WRITE(6,23)
23	FORMAT(' NORMALIZE_L1/Version 26-09-91')

	 SSUM=0.
	 DO I=1,NPTS
	   SSUM=SSUM+ABS(ARRAY(I))
	 END DO

C Normalisation of the weigths:
	IF(SSUM.EQ.0)THEN
	  WRITE(2,24)
	  WRITE(6,24)
24	  FORMAT(' NORMALIZE_L1/Fatal error: norm is null')
	  STOP
	ENDIF

	DO I=1,NPTS
	   ARRAY(I) = ARRAY(I)/SSUM
	END DO

	NORM=SQRT(SSUM)

        RETURN
	END
C**********************************************************
	SUBROUTINE NORMALIZE_L2(ARRAY,NPTS,NORM)
	REAL ARRAY(*),NORM
	REAL*8 SSUM
	INTEGER NPTS

	 SSUM=0.
	 DO I=1,NPTS
	   SSUM=SSUM+ARRAY(I)*ARRAY(I)
	 END DO

C Normalisation of the weigths:
	IF(SSUM.EQ.0)THEN
	  WRITE(2,24)
	  WRITE(6,24)
24	  FORMAT(' NORMALIZE_L2/Fatal error: norm is null')
	  STOP
	ENDIF

	NORM=SQRT(SSUM)

	DO I=1,NPTS
	   ARRAY(I) = ARRAY(I)/NORM
	END DO


        RETURN
	END
C**********************************************************
C Largest eigen value:
C**********************************************************
	SUBROUTINE EIGEN_VALUE1(NBETA,NGAMMA,XLAMBDA1)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL X,XNORM,XLAMBDA1,RO,BMASK
	COMPLEX XC,YCE
	INTEGER NBETA,NGAMMA,NB,I
	COMMON /C4/X(NULL:NBMAX,NULL:3)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)

C Initialization:
	DO NB=1,NBETA
	  X(NB,3)=BMASK(NB)
	END DO

	DO I=1,20

C Projection of this error onto E^+
	  CALL PROJ_EPLUS(NBETA,3,1)

C Normalization: 
	  CALL NORMALIZE_L2(X(1,1),NBETA,XNORM)

C          WRITE(6,45) I-1,XNORM
C45	  FORMAT(' Eigen_value1/ Iteration ',I3,' Norm: ',G12.5)

C Computing ATwA of X(.,1). 
C Output in X(.,3)
          CALL ATRANS_A(NBETA,NGAMMA,1,3)

	ENDDO

	XLAMBDA1 = XNORM

	RETURN
	END
C**********************************************************
C Smallest eigen value:
C Power method applied to I - AtwA/lambda1
C**********************************************************
	SUBROUTINE EIGEN_VALUE2(NBETA,NGAMMA,XLAMBDA1,XLAMBDA2)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL X,XNORM,XLAMBDA2,RO,BMASK
	COMPLEX XC,YCE
	INTEGER NBETA,NGAMMA,NB,I
	COMMON /C4/X(NULL:NBMAX,NULL:3)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)

C Initialization:
	DO NB=1,NBETA
	  X(NB,1)=BMASK(NB)
	END DO

C Main loop:
	DO I=1,20

C Projection of this error onto E^+
C        X(1,1)=0.
C        X(2,1)=0.
	CALL PROJ_EPLUS(NBETA,1,1)

C Normalization: 
	  CALL NORMALIZE_L2(X(1,1),NBETA,XNORM)

C          WRITE(6,45) I-1,XNORM
C45	  FORMAT(' Eigen_value2/ Iteration ',I3,' Norm: ',G12.5)

C Computing ATwA of X(.,1). 
C Output in X(.,3)
	  CALL ATRANS_A(NBETA,NGAMMA,1,3)

C Now computing I - ATwA/lambda1 :
	  DO NB=1,NBETA
	    X(NB,1) = X(NB,1) - X(NB,3)/XLAMBDA1
	  END DO

	ENDDO

C mu = 1 - lambda2/lambda1
C Thus lambda2=lambda1*(1-mu)
	XLAMBDA2=XLAMBDA1*(1.-XNORM)


C********************************
C Computing the norm:
	I=0
	IF(I.EQ.1)THEN
C Normalization of the eigen vector: 
	  CALL NORMALIZE_L2(X(1,1),NBETA,XNORM)

C Projection to E^+
	  CALL PROJ_EPLUS(NBETA,1,1)
C Computing the norm:
	  CALL NORMALIZE_L2(X(1,1),NBETA,XNORM)
	  WRITE(6,71) XNORM
	  WRITE(2,71) XNORM
71        FORMAT(' XNORM(proj): ',1PG12.5)
        ENDIF

	RETURN
	END
C**********************************************************
	SUBROUTINE PROJ_IRMAXI(ARRAY,NBETA,IR_MAXI)
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL ARRAY(*),W1
	INTEGER NBETA,IR_MAXI,NBCOUV,IXY
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)

C Troncates to IR_MAXI
	 INULL=0
	 DO NB=1,NBETA
	    W1=IXY(1,NB)**2+IXY(2,NB)**2
	    W1= SQRT(W1)
	    IF(W1.GT.FLOAT(IR_MAXI))THEN
	      ARRAY(NB)=0.
	      INULL=INULL+1
	    ENDIF
	 END DO

	IF(INULL.GT.0)THEN
	  WRITE(2,79) INULL
	  WRITE(6,79) INULL
79	  FORMAT(' PROJ_IRMAXI/Warning: INULL=',I5)
	ENDIF

	RETURN
	END
C**********************************************************************
C Output SNR map according to the adopted weights
C Generating SNR and sigma maps (needed by DIANE)
C*******************************************************************
        SUBROUTINE OUTPUT_SNR1(SIGM,NBETA,FNAME)
	PARAMETER(IDIM=256)

C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	COMPLEX XC,YCE,CC
	REAL RE,IM,RO,BMASK,SIGM(*)
	INTEGER NBCOUV,IXY,NX,NY
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
	COMMON /C3/XC(NULL:NBMAX),YCE(NGMAX),RO(NULL:NBMAX),
     1             BMASK(NULL:NBMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
 
C Central position: 
	IXC = (NX/2) + 1
	IYC = (NY/2) + 1
 
C Initializing RE and IM:
C 1/SIGM will be in RE and SIGM in IM:
	DO 1 J=1,NY
	  DO 2 I=1,NX
	    RE(I,J) = 0.
	    IM(I,J) = 0.
2	  CONTINUE
1	CONTINUE
 
C Generating SNR and sigma maps (needed by DIANE):
	DO 3 NB = 0,NBETA
	  IX=IXY(1,NB)
	  IY=IXY(2,NB)
	  IIX = IXC + IX
	  IIY = IYC + IY
          SIGG=MAX(5.E-2,SIGM(NB))
	  RE(IIX,IIY) = 1./SIGG
	  IM(IIX,IIY) = SIGG*RO(NB)
	  IIX = IXC - IX
	  IIY = IYC - IY
	  RE(IIX,IIY) = 1./SIGG
	  IM(IIX,IIY) = SIGG*RO(NB)
3	CONTINUE
 
C Output of SNR image: 
	  NAME='snr1'
	  COMMENTS=' 1/SIGMA_PHASE of : '//FNAME(1:20)
	  WRITE(6,54) COMMENTS(1:80)
54        FORMAT(' Output of SNR=1/SIGMA_PHASE in image file: snr1')
	  CALL JLP_WRITEIMAG(RE,NX,NY,IDIM,NAME,COMMENTS)

C Output of SIGMA
	  NAME='sigma'
	  COMMENTS=' MODULUS * SIGMA_PHASE of : '//FNAME(1:20)
	  WRITE(6,55) COMMENTS(1:80)
55        FORMAT(' Output of MODULUS * SIGMA_PHASE in image file: sigma')
	  CALL JLP_WRITEIMAG(IM,NX,NY,IDIM,NAME,COMMENTS)

	RETURN
	END
C**********************************************************************
	INCLUDE 'jlp_cover_30.for'
C**********************************************************************
