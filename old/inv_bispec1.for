C****************************************************************
C Same as inv_bispec_30
C But use C routines for uv coverage (jlp_cover_mask.c)
C Warning: This prog does not use yet BMASK (---> check if this is right later!)
C
C PROGRAM INV_BISPEC.FOR : 2DPA (2D PHASE A)
C To compute the phase of the spectrum from the bispectrum,
C assuming a full pupil.
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
C Version: 04-02-94
C From A. Lannes P.FOR (January 1988)
C****************************************************************
	PROGRAM INV_BISPEC
 
C  POUR IR=25
C        PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
C        PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)
 
C  POUR IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	PARAMETER(IDIM=256)
 
C XC: Spectrum (list): phase factor
C RO: Spectrum (list): modulus
C YCE: Bispectrum (list)
C XCR: Reference spectrum (list), phase factor  (not known for real observ.)
C COMPLEX XC(NBMAX+1),YCE(NGMAX),XCR(NBMAX+1)
        INTEGER PNTR_XC,PNTR_YCE,PNTR_XCR
C REAL SIGM(NBMAX+1),RO(NBMAX+1),X(NBMAX+1,4),BMASK(NBMAX+1)
        INTEGER PNTR_SIGM,PNTR_RO,PNTR_X,PNTR_BMASK
C Kernel vectors: VECXY
C REAL VECXY(NBMAX+1,2),WEIGHT(NGMAX)
        INTEGER PNTR_VECXY,PNTR_WEIGHT
	COMPLEX CC
	INTEGER IFERMAX,NX,NY,NBX,NBY,IDIM_X
        INTEGER NX_MASK,NY_MASK,PNTR_MASK,NCLOSURE_MAX
	REAL RE,IM,LOWER_MODULUS,DEGRAD,SIG_MAX
        REAL EXIT_TOLERANCE
	CHARACTER FNAME*40,COMMENTS*80,DATE*24,ERRNAME*40
 
        INTEGER*4 MEMSIZE,PNTR_ARRAY,MADRID(1)

C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A  NBMAX)
C COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
C COMMON /C2/KLM(3,NGMAX)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
C Simply to give a limit to the number of closure relations to
C be taken into account (when this number would be really too big...)
C For virtual memory:
        COMMON /VMR/MADRID
 
	CALL JLP_BEGIN
        CALL JLP_RANDOM_INIT(1)
	OPEN(2,FILE='inv_bispec.log',STATUS='UNKNOWN',ERR=998)
C Max. number of closure phase relations allowed for each pixel of
C the uv-coverage:
	IFERMAX = 10000

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
C IFERMAX is an internal upper limit, whereas NCLOSURE_MAX is set 
C when computing the data:
	WRITE(2,52) DATE,IOPT,IFERMAX
	WRITE(6,52) DATE,IOPT,IFERMAX
52	FORMAT(' Program inv_bispec1 Version 10-02-94 ',/,
     1	A,/,' Option = ',I3,/,
     1  ' Maximum number of closure relations to be processed here:',I8)
 
C Angular reference for bispectral perturbation simulations (degrees) 
C Used for weight computation and weight selection... 
C Cette reference est associee au carre de la norme de la frequence IR
	CTE = 20.
 
C Test for the main iteration (in degrees)
C This means that when the largest angular correction
C is larger than EXIT_TOLERANCE, we exit from the main loop:
	EXIT_TOLERANCE = 0.1
        NCLOSURE_MAX=10000
 
	WRITE(6,92)
92	FORMAT(' Radius (IR) of uv-coverage, phase error ',
     1    '(0 if all weights=1, <0 if modulus weighted), exit_tolerance (degrees),'
     1    /,' nbeta, lower_modulus, max sigma, sigma_value_when_modulus_is_null',
     1    ' Max nclosure (used for input data computation)',/,
     1    '  (Ex: 12,20,0.1,220,0.,0.8,0.7,10000)')
	READ(5,*) IR,CTE,EXIT_TOLERANCE,NBETA,LOWER_MODULUS,SIG_MAX,
     1    SIGMA_NULL,NCLOSURE_MAX
	WRITE(2,59) IR,CTE,EXIT_TOLERANCE,LOWER_MODULUS,SIG_MAX,SIGMA_NULL
59	FORMAT(' IR = ',I4,' CTE = ',F8.2,
     1   ' EXIT_TOLERANCE = ',G12.4,' Lower_modulus: ',G12.5,/,
     1   ' Sig_max:',G10.3,' Sig_null:',G10.3)

C UV coverage:
C CALL COVERA(IR,NBETA_MAX,NGAMMA_MAX)
        NX_MASK=2*IR+4
        NY_MASK=NX_MASK
        MEMSIZE=NX_MASK*NY_MASK*4
        CALL JLP_GETVM(PNTR_MASK,MEMSIZE)
        CALL FILL_MASK(MADRID(PNTR_MASK),NX_MASK,NY_MASK)
C IFERMAX is an internal upper limit, whereas NCLOSURE_MAX is set 
C when computing the data:
        CALL COVERA_MASK(MADRID(PNTR_MASK),NX_MASK,NY_MASK,
     1                   IR,NCLOSURE_MAX,NBETA_MAX,NGAMMA_MAX)
        CALL JLP_FREEVM(PNTR_MASK,MEMSIZE)
 
        IF(NBETA.GT.NBETA_MAX)THEN
	  WRITE(6,229) NBETA,NBETA_MAX
	  WRITE(2,229) NBETA,NBETA_MAX
229	  FORMAT(' JLP_BISPEC/Error: NBETA (wanted) =',I4,
     1         ' whereas NBETA_MAX =',I5,/,' I correct it to BETA_MAX')
          NBETA = NBETA_MAX
        ENDIF

C NGAMMA=NGT(NBETA)
        CALL COVER_NGT(NGAMMA,NBETA)
	WRITE(6,29) NBETA_MAX,NGAMMA_MAX,NCLOSURE_MAX,NBETA,NGAMMA
	WRITE(2,29) NBETA_MAX,NGAMMA_MAX,NCLOSURE_MAX,NBETA,NGAMMA
29	FORMAT(' NBETA_MAX (input data):',I5,' NGAMMA_MAX (input data):',I8,/,
     1         ' Maximum number of closure relations (input data): ',I8,/,
     1         ' NBETA (used here):',I5,' Corresponding NGAMMA:',I8)

C Max number of iterations in the main loop:
	ITTMAX = 15
 
	WRITE(6,53) CTE,EXIT_TOLERANCE
53	FORMAT(' Angular constant for bispectral noise estimation',
     1 F8.3,' (degrees)',/,' EXIT_TOLERANCE (for exit check) :',E12.3)
 
C Allocation of memory
C COMPLEX XC(NBMAX+1),XCR(NBMAX+1)
	MEMSIZE=(NBETA+1)*8
	CALL JLP_GETVM(PNTR_XC,MEMSIZE)
	CALL JLP_GETVM(PNTR_XCR,MEMSIZE)
C COMPLEX YCE(NGMAX)
	MEMSIZE=(NGAMMA+1)*8
	CALL JLP_GETVM(PNTR_YCE,MEMSIZE)
C REAL WEIGHT(NGMAX)
	MEMSIZE=(NGAMMA+1)*4
	CALL JLP_GETVM(PNTR_WEIGHT,MEMSIZE)
C REAL RO(NBMAX+1),SIGM(NBMAX+1),BMASK(NBMAX+1)
	MEMSIZE=(NBETA+1)*4
	CALL JLP_GETVM(PNTR_RO,MEMSIZE)
	CALL JLP_GETVM(PNTR_SIGM,MEMSIZE)
	CALL JLP_GETVM(PNTR_BMASK,MEMSIZE)
C REAL VECXY(NBMAX+1,2)
	MEMSIZE=(NBETA+1)*2*4
	CALL JLP_GETVM(PNTR_VECXY,MEMSIZE)
C REAL X(NBMAX+1,4)
	MEMSIZE=(NBETA+1)*4*4
	CALL JLP_GETVM(PNTR_X,MEMSIZE)
C First dimension of array X
        IDIM_X=NBETA+1

C Creates YCE : bispectrum phasor,
C and initial set of weights 
	IF(IOPT.EQ.0)THEN
 	  CALL CREYCE_SIMU(IR,CTE,NBETA,NGAMMA,FNAME,MADRID(PNTR_BMASK),
     1             MADRID(PNTR_RO),MADRID(PNTR_YCE),MADRID(PNTR_WEIGHT),
     1             MADRID(PNTR_XCR),MADRID(PNTR_XC))
	ELSEIF(IOPT.EQ.1) THEN
 	  CALL CREYCE1(IR,CTE,NBETA,NGAMMA,NGAMMA_MAX,FNAME,
     1       LOWER_MODULUS,MADRID(PNTR_BMASK),MADRID(PNTR_RO),
     1       MADRID(PNTR_YCE),MADRID(PNTR_WEIGHT),MADRID(PNTR_XC))
	ELSE
 	  CALL CREYCE2(IR,CTE,NBETA,NGAMMA,NGAMMA_MAX,FNAME,
     1       LOWER_MODULUS,MADRID(PNTR_BMASK),MADRID(PNTR_RO),MADRID(PNTR_YCE),
     1       MADRID(PNTR_WEIGHT),MADRID(PNTR_XCR),MADRID(PNTR_XC))
	ENDIF
 
C******************************************************
C To initialize the solution,
C we solve the problem with the recursive method (Weigelt,...) 
C The initial spectral phase is stored in XC(.)
	CALL RECURSIVE(MADRID(PNTR_SIGM),NBETA,NGAMMA,SIGMA_NULL,IFERMAX,
     1            MADRID(PNTR_BMASK),MADRID(PNTR_YCE),MADRID(PNTR_WEIGHT),
     1            MADRID(PNTR_XC),MADRID(PNTR_RO))
 
C New version of the weights:
C JLP93:
        IF(CTE.EQ.0.)THEN
	  CALL WEIGHTS_SIGM1(MADRID(PNTR_SIGM),NBETA,NGAMMA,SIG_MAX,
     1         LOWER_MODULUS,IFERMAX,MADRID(PNTR_BMASK),MADRID(PNTR_WEIGHT))
        ELSE
	  CALL WEIGHTS_SIGM2(MADRID(PNTR_SIGM),NBETA,NGAMMA,SIG_MAX,
     1         LOWER_MODULUS,IFERMAX,MADRID(PNTR_BMASK),MADRID(PNTR_WEIGHT))
        ENDIF

C Output of SNR map:
        CALL OUTPUT_SNR1(MADRID(PNTR_SIGM),NBETA,FNAME,MADRID(PNTR_RO))

C Output the errors of initial bispectrum:
        IF(IOPT.LT.-1045)THEN
	  ERRNAME='bisp_error1.dat'
          NBX=NBETA
          NBY=NBETA/2
          MEMSIZE=NBX*NBY*4
          CALL JLP_GETVM(PNTR_ARRAY,MEMSIZE)
	  CALL ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,
     1                       MADRID(PNTR_ARRAY),NBX,NBY,
     1                       IFERMAX,MADRID(PNTR_YCE),MADRID(PNTR_XC))
	   ERRNAME='bisperr1'
	   COMMENTS='Quadratic errors'
           CALL JLP_WRITEIMAG(MADRID(PNTR_ARRAY),NBX,NBY,NBX,
     1  ERRNAME,COMMENTS)
           CALL JLP_FREEVM(PNTR_ARRAY,MEMSIZE)
        ENDIF

C Initial error estimation:
	IF(IOPT.NE.1)THEN
	  ERRNAME='error1.dat'
	  CALL ERROR_SIMU(NBETA,EI,EPI,ERRNAME,
     1                    MADRID(PNTR_BMASK),MADRID(PNTR_RO),
     1                    MADRID(PNTR_XCR),MADRID(PNTR_XC))
	  WRITE(6,64) EI,EPI
	  WRITE(2,64) EI,EPI
64	  FORMAT(' Comparison with the model (since simulation)',/,
     1  ' Initial rms error of the spectrum',1PG12.5,/,
     1	' Initial rms error of the phase factor of the spect.',1PG12.5)
	ENDIF
 
C Projection of the initial solution onto E^+
	CALL EPLUS(NBETA,MADRID(PNTR_BMASK),MADRID(PNTR_XC),
     1             MADRID(PNTR_X),IDIM_X,MADRID(PNTR_VECXY))
 
C Other alternative (Null phases)
C        CALL ZERO_PHASE(MADRID(PNTR_XC),NBETA)
 
C SORTIE DES PARTIES REELLES ET IMAGINAIRES DE XC
C In rei and imi (i for initial)
	CALL SORTIE(NBETA,1,FNAME,MADRID(PNTR_BMASK),MADRID(PNTR_RO),
     1               MADRID(PNTR_XC))

C Main loop: least square non linear fit. 
	DEGRAD = 3.141592/180.
	EXIT_TOLERANCE = EXIT_TOLERANCE*DEGRAD
        CALL LSQUARES1(EXIT_TOLERANCE,NBETA,NGAMMA,IFERMAX,
     1           MADRID(PNTR_BMASK),MADRID(PNTR_YCE),MADRID(PNTR_WEIGHT),
     1           MADRID(PNTR_XC),MADRID(PNTR_X),IDIM_X,ITTMAX)
 
C Files ref et imf (f for final)
	CALL SORTIE(NBETA,2,FNAME,MADRID(PNTR_BMASK),MADRID(PNTR_RO),
     1               MADRID(PNTR_XC))
 
C Calage en translation
	CALL TRANS(IR,NBETA,MADRID(PNTR_BMASK),MADRID(PNTR_XC))

C ERREURS FINALES GLOBALES ET BILAN
	IF(IOPT.NE.1)THEN
	   ERRNAME='error2.dat'
           CALL FINAL_ERROR(IR,NBETA,EF,EPF,ERRNAME,
     1                   MADRID(PNTR_BMASK),MADRID(PNTR_RO),
     1                   MADRID(PNTR_XCR),MADRID(PNTR_XC),
     1                   MADRID(PNTR_X),IDIM_X,MADRID(PNTR_VECXY),SIGM)
 
	   WRITE(2,65) EI,EF,EPI,EPF
	   WRITE(6,65) EI,EF,EPI,EPF
65	  FORMAT(' Comparison with the model (since simulation)',/,
     1  ' Initial & final rms error of the spectrum:',2(1X,F8.3),/,
     1  ' Initial & final rms error of the phase factor (spect):',
     1   2(1X,F8.3))
 
	ENDIF
 
999	WRITE(6,*) ' Log File in "inv_bispec.log"'
C
C Output the errors of final bispectrum:
        I=0
        IF(I.EQ.1234)THEN
	  ERRNAME='bisp_error2.dat'
          NBX=NBETA
          NBY=NBETA/2
          MEMSIZE=NBX*NBY*4
          CALL JLP_GETVM(PNTR_ARRAY,MEMSIZE)
	  CALL ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,
     1                       MADRID(PNTR_ARRAY),NBX,NBY,
     1                       IFERMAX,MADRID(PNTR_YCE),MADRID(PNTR_XC))
	  ERRNAME='bisperr2'
	  COMMENTS='Quadratic errors'
          CALL JLP_WRITEIMAG(MADRID(PNTR_ARRAY),NBX,NBY,NBX,
     1  ERRNAME,COMMENTS)
          CALL JLP_FREEVM(PNTR_ARRAY,MEMSIZE)
        ENDIF

C Eigenvalues:
        I=0
        IF(I.EQ.1234)THEN
	CALL EIGEN_VALUE1(NBETA,NGAMMA,XLAMBDA1,IFERMAX,
     1                    MADRID(PNTR_BMASK),MADRID(PNTR_WEIGHT),
     1                    MADRID(PNTR_X),IDIM_X,MADRID(PNTR_VECXY))
	CALL EIGEN_VALUE2(NBETA,NGAMMA,XLAMBDA1,XLAMBDA2,IFERMAX,
     1                    MADRID(PNTR_BMASK),MADRID(PNTR_WEIGHT),
     1                    MADRID(PNTR_X),IDIM_X,MADRID(PNTR_VECXY))
	XLAMBDA2=MAX(XLAMBDA2,1.E-12)
        WRITE(6,48) XLAMBDA1,XLAMBDA2,XLAMBDA1/XLAMBDA2
        WRITE(2,48) XLAMBDA1,XLAMBDA2,XLAMBDA1/XLAMBDA2
48	FORMAT(' Largest eigen value: ',G12.5,
     1  /,' Smallest eigen value: ',G12.5,
     1  /,' Conditionning number: ',G12.5)
        ENDIF

	CLOSE(2)
C Free memory:
	MEMSIZE=(NBETA+1)*8
	CALL JLP_FREEVM(PNTR_XC,MEMSIZE)
	CALL JLP_FREEVM(PNTR_XCR,MEMSIZE)
	MEMSIZE=(NGAMMA+1)*8
	CALL JLP_FREEVM(PNTR_YCE,MEMSIZE)
	MEMSIZE=(NGAMMA+1)*4
	CALL JLP_GETVM(WEIGHT,MEMSIZE)
	MEMSIZE=(NBETA+1)*4
	CALL JLP_FREEVM(PNTR_RO,MEMSIZE)
	CALL JLP_FREEVM(PNTR_SIGM,MEMSIZE)
	CALL JLP_FREEVM(PNTR_BMASK,MEMSIZE)
	MEMSIZE=(NBETA+1)*4*4
	CALL JLP_GETVM(PNTR_X,MEMSIZE)
	MEMSIZE=(NBETA+1)*2*4
	CALL JLP_GETVM(PNTR_VECXY,MEMSIZE)

	CALL JLP_END
	STOP
998	WRITE(6,*) ' Fatal error opening inv_bispec.log'
	CALL JLP_END
	STOP
	END
C*******************************************************************
C CREYCE_SIMU : CRE LES TERMES DE PHASE DU SPECTRE ET DU BISPECTRE
C Output in real/imaginary (from artificial bispectrum)
C Centers the Fourier transform
C*******************************************************************
	SUBROUTINE CREYCE_SIMU(IR,CTE,NBETA,NGAMMA,FNAME,BMASK,
     1                         RO,YCE,WEIGHT,XCR,XC)
	PARAMETER(IDIM=256)
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	REAL RE,IM,WEIGHT(*),RO(*),BMASK(*)
        INTEGER IR,NBETA,NGAMMA,IXC,IYC,NX,NY,NB,INULL,IIX,IIY
        REAL XR,XI,XM 
	COMPLEX XCR(*),XC(*),YCE(*)
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
C IIX = IXC + IXY(1,NB)
C IIY = IYC + IXY(2,NB)
          CALL COVER_IXY(IIX,IIY,NB)
          IIX = IXC + IIX
          IIY = IYC + IIY
	  XR = RE(IIX,IIY)
	  XI = IM(IIX,IIY)
C Modulus:
	  XM = SQRT(XR*XR + XI*XI)
	  RO(NB+1) = XM
C Phase factor:
	  IF(XM.EQ.0.)THEN
	    XC(NB+1) = (1.,0.)
            INULL=INULL+1
C	    WRITE(6,78) NB
78	    FORMAT(' CREYCE_SIMU/Warning: XC(',I5,') is null!') 
	  ELSE
	    XC(NB+1) = CMPLX(XR/XM,XI/XM)
	  ENDIF
1	CONTINUE
	    IF(INULL.GT.0)THEN
	      WRITE(6,79) INULL
	      WRITE(2,79) INULL 
79	      FORMAT(' CREYCE_SIMU/Warning: Modulus is null for'
     1               ,I5,' values') 
	    ENDIF

	DO I=0,4
	   WRITE(2,*) ' RO(',I,') =',RO(I+1)
	   WRITE(6,*) ' RO(',I,') =',RO(I+1)
	   WRITE(2,*) ' XC(',I,') =',XC(I+1)
	   WRITE(6,*) ' XC(',I,') =',XC(I+1)
	END DO
 
C Mask to discard some values of the spectral list:
        DO 91 NB=0,NBETA
          BMASK(NB+1)=1.
C Check that RO(NB+1) greater than LOWER_MODULUS:
C Remember that LOWER_MODULUS can be negative...
          RO(NB+1)=RO(NB+1)/RO(0+1)
	  IF(RO(NB+1).LT.LOWER_MODULUS.OR.
     1       RO(NB+1).LE.0)BMASK(NB+1)=0.
91      CONTINUE
 
C FACTEUR DE PHASE DU SPECTRE CALE EN TRANSLATION : XC(.)
	CALL TRANS(IR,NBETA,BMASK,XC)

C Reference (since simulation)
	DO 3 NB=0,NBETA
	  XCR(NB+1)=XC(NB+1)
3	CONTINUE
 
C FACTEUR DE PHASE DU BISPECTRE : YCE
	DO 2 NG=1,NGAMMA
C  YCE(NG)=XC(KLM(1,NG))*XC(KLM(2,NG))*CONJG(XC(KLM(3,NG)))
          CALL COVER_KLM(KLM1,1,NG)
          CALL COVER_KLM(KLM2,2,NG)
          CALL COVER_KLM(KLM3,3,NG)
	  YCE(NG)=XC(KLM1+1)*XC(KLM2+1)*CONJG(XC(KLM3+1))
2	CONTINUE
 
	DO I=1,5
	   WRITE(2,*) ' YCE(',I,') =',YCE(I)
	   WRITE(6,*) ' YCE(',I,') =',YCE(I)
	END DO
 
C Computing the weights:
	IF(CTE.GT.0)THEN
	   CALL BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSEIF(CTE.LT.0)THEN
	   CALL BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSE
	   WRITE(6,34)
	   WRITE(2,34)
34         FORMAT(' Weights set to unity, and then normalized')
	   DO 96 NG=1,NGAMMA
	     WEIGHT(NG)=1.
	     DO 95 KK=1,3
C  K = KLM(KK,NG)
               CALL COVER_KLM(K,KK,NG)
	       IF(BMASK(K+1).EQ.0.)WEIGHT(NG)=0.
95           CONTINUE
96         CONTINUE
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ENDIF
 
C Perturbation of YCE
	DEGRAD = 3.141592/180.
	CTE1 = CTE * DEGRAD /FLOAT(IR**2)
 
	DO 4 NG=1,NGAMMA
	  IS2=0
	    DO 6 KK=1,3
C  K = KLM(KK,NG)
               CALL COVER_KLM(K,KK,NG)
C IS2 = IS2 + IXY(1,K)**2 + IXY(2,K)**2
               CALL COVER_IXY(IXY1,IXY2,K)
 	       IS2 = IS2 + IXY1**2 + IXY2**2
6	    CONTINUE
 
C Random generation (Gaussian law, (1.,0.)) of DGAMMA = DELTA GAMMA
          CALL JLP_RANDOM_GAUSS(WORK)
	  DGAMMA = WORK*CTE1*FLOAT(IS2)
 	  YCE(NG) = YCE(NG)*CMPLX(COS(DGAMMA),SIN(DGAMMA))
4	CONTINUE
 
	RETURN
	END
C*******************************************************************
C TRANS.FOR  
C fait le calage en translation du facteur de phase spectral XC(.)
C  XC(1+1) = exp i Beta(1)
C
C Array TX is defined as:
C  TX(1) = exp -i Beta(1)
C  TX(-1)= exp i Beta(1)
C  TX(2) = exp -i 2 Beta(1)
C   ....
C  TX(N) = exp -i N Beta(1)
C
C OUTPUT:
C  If coordinates of frequency #NB are: (IXY1,IXY2)
C  XC(NB+1) = XC(NB+1) * TX(IXY1) * TY(IXY2)
C  i.e.:
C  exp i Beta_new = exp i Beta_old * exp -i IXY1 Beta(1) * exp -i IXY2 Beta(2)
C*******************************************************************
	SUBROUTINE TRANS(IR,NBETA,BMASK,XC)
 
	PARAMETER(IRMAX=100)
	REAL BMASK(*)
	COMPLEX XC(*)
	COMPLEX TX(IRMAX+1),TY(IRMAX+1),CX,CCX,CY,CCY
 
        IF(IR.GT.IRMAX)THEN
           PRINT *,' TRANS/Fatal error, maximum IR=',IRMAX
           STOP
        ENDIF

	CX=CONJG(XC(1+1))
	CY=CONJG(XC(2+1))
 
C TABLEAU DE LA FORME LINEAIRE CONCERNE DU NOYAU
	CCX=(1.,0.)
	CCY=(1.,0.)
	TX(0+1)=(1.,0.)
	TY(0+1)=(1.,0.)
 
	DO 1 I=1,IR
C
	 CCX=CCX*CX
	 TX(I+1)=CCX
C
	 CCY=CCY*CY
	 TY(I+1)=CCY
1	CONTINUE
 
C CALAGE EN TRANSLATION
	DO 2 NB=0,NBETA
C XC(NB) = XC(NB) * TX(IXY(1,NB)) * TY(IXY(2,NB))
            CALL COVER_IXY(IXY1,IXY2,NB)
              IF(IXY1.GE.0)THEN
                CX = TX(IXY1+1)
              ELSE
                CX = CONJG(TX(1-IXY1))
              ENDIF
            CY = TY(IXY2+1)
	    XC(NB+1) = XC(NB+1)*CX*CY
            XC(NB+1) = XC(NB+1)*BMASK(NB+1)
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
	SUBROUTINE RECURSIVE(SIGM,NBETA,NGAMMA,SIGMA_NULL,IFERMAX,
     1                        BMASK,YCE,WEIGHT,XC,RO)
	COMPLEX XC(*),YCE(*),CC1
	REAL BMASK(*),WEIGHT(*),SIGM(*),RO(*)
C Double precision is necessary for large sums !!!
        REAL*8 SUMXR1,SUMXI1,SUMSQR,SUMSQI,SIGR,SIGI,SUMWEIGHTS
        REAL*8 RNORM
        REAL*4 XR1,XI1,XW,SIGMA_NULL
	INTEGER K,L,M,NG1,INDX,IFERMAX,NVAL
 
	XC(0+1)=(1.,0.)
	XC(1+1)=(1.,0.)
	XC(2+1)=(1.,0.)
	NG1=1
	INULL=0
 
	DO 1 NB=3,NBETA
          SUMXI1=0.
          SUMXR1=0.
          SUMSQR=0.
          SUMSQI=0.
          SUMWEIGHTS=0.
C JLP94: sigma has a meaning only if more than 2 values are involved in the
C mean... 
          NVAL=0
C  NG2=NGT(NB)
          CALL COVER_NGT(NG2,NB)
	  INDX=0
C
	  DO 2 NG=NG1,NG2
	    INDX=INDX+1
	     IF (INDX.LE.IFERMAX)THEN
C  K=KLM(1,NG)
C  L=KLM(2,NG)
C  M=KLM(3,NG)
               CALL COVER_KLM(K,1,NG)
               CALL COVER_KLM(L,2,NG)
               CALL COVER_KLM(M,3,NG)
	       XW=SQRT(WEIGHT(NG))
               IF(XW.GT.0.)THEN
                 NVAL=NVAL+1
                 SUMWEIGHTS=SUMWEIGHTS+XW
	           IF (M.EQ.NB) THEN
C CAS 1 (M = NB)
	             CC1=XC(K+1)*XC(L+1)*CONJG(YCE(NG))
	           ELSE
C CAS 2 (L = NB)
	             CC1=YCE(NG)*CONJG(XC(K+1))*XC(M+1)
	           ENDIF
                 XR1=REAL(CC1)
                 XI1=IMAG(CC1)
                 SUMSQR=SUMSQR+XR1*XR1*XW
                 SUMSQI=SUMSQI+XI1*XI1*XW
                 SUMXR1=SUMXR1+XR1*XW
                 SUMXI1=SUMXI1+XI1*XW
C End of test on XW 
               ENDIF
C End of test on IFERMAX
	     ENDIF
2	  CONTINUE
 
C NORMALISATION DU TERME DE PHASE
	  RNORM=SQRT(SUMXR1*SUMXR1+SUMXI1*SUMXI1)
	  IF(RNORM.GT.0.)THEN
	    XC(NB+1)=CMPLX(SUMXR1/RNORM,SUMXI1/RNORM)

C JLP94: sigma has a meaning only if more than 2 values are involved in the
C mean... 
            IF(NVAL.GT.1)THEN
C SIGR**2 = Sum of (weight_i *(X_i - mean)**2) / Sum of (weight_i) 
C or Sum of (weight_i X_i**2)/Sum of weight_i - mean**2
C Here the weighted mean of real(CC) is SUMXR1/SUMWEIGHTS:
               SIGR=SUMSQR/SUMWEIGHTS-(SUMXR1/SUMWEIGHTS)**2
               SIGI=SUMSQI/SUMWEIGHTS-(SUMXI1/SUMWEIGHTS)**2
C Note that double precision variables are needed for large sums !!!
                IF(SIGR.LT.0..OR.SIGI.LT.0.)THEN
C Just in case of a (round-off??) problem:
                  PRINT *,' FATAL problem in RECURSIVE... for NB=',NB
                  PRINT *,' SIGM(NB):',SIGM(NB+1),'SIGI:',SIGI,'SIGR',SIGR
                  PRINT *,' RNORM:',RNORM,'SUMWEIGHTS:',SUMWEIGHTS
                  PRINT *,' SUMSQR:',SUMSQR,' SUMSQI:',SUMSQI
                  PRINT *,' SUMXR1:',SUMXR1,' SUMXI1:',SUMXI1
                  PRINT *,' NG1, NG2',NG1,NG2
                  PRINT *,' WEIGHT(NG1,NG2)',(WEIGHT(III),III=NG1,NG2)
                  STOP
                ELSE
                  SIGM(NB+1)=SQRT(SIGI+SIGR)
                ENDIF
C If only one value has been taken into account for the mean computation:
            ELSE
               SIGM(NB+1)=SIGMA_NULL
            ENDIF
	  ELSE
C The recursive process has been interrupted:
	    INULL=INULL+1
	    XC(NB+1)=(1.,0.)
C Maximum sigma is one, but I set it to SIGMA_NULL since we assume that
C phase is not important if modulus is too small
            SIGM(NB+1)=SIGMA_NULL
C Allows breaks only if in the list of null modulus:
	    IF(BMASK(NB+1).NE.0.)THEN
C JLP94: Then add this frequency to BMASK, i.e., will not recover the
C phase of this frequency
              BMASK(NB+1)=0.
	      WRITE(6,78) NB
78	      FORMAT(/,' RECURSIVE/Warning: too many null modulii: '
     1     ' XC(',I5,') will remain undetermined!') 
C Debug 94...
              I=0
              IF(I.EQ.1234)THEN
	      WRITE(2,82) NB,SUMXR1,SUMXI1,SUMWEIGHTS,SUMSQR,SUMSQI
	      WRITE(6,82) NB,SUMXR1,SUMXI1,SUMWEIGHTS,SUMSQR,SUMSQI
82	   FORMAT(' RECURSIVE/Warning: Recursive process interrupts at NB=',I5,
     1      /,' SUMXR1,SUMXI1,SUMWEIGHTS,SUMSQR,SUMSQI: ',5(G12.5,1X))
              CALL COVER_NGT(IG1,NB-1)
              CALL COVER_NGT(IG2,NB)
              CALL COVER_NGT(IG3,NB+1)
              WRITE(6,83) IG1,IG2,IG3,IFERMAX,NG1,NG2,INDX,
     1                    BMASK(NB),BMASK(NB+1),BMASK(NB+2),WEIGHT(IG2)
83            FORMAT(' NGT(NB-1), NGT(NB), NGT(NB+1): ',3(I8,1X),
     1               ' (IFERMAX=',I8,')',/,' NG1=',I5,' NG2=',I5,' INDX=',I5,
     1               ' BMASK(NB), BMASK(NB+1), BMASK(NB+2) = ',3(G12.5,1X),/
     1               ' WEIGHT(NGT(NB)) = ',G12.5)
	     DO NG=NG1,NG2
               CALL COVER_KLM(K,1,NG)
               CALL COVER_KLM(L,2,NG)
               CALL COVER_KLM(M,3,NG)
               PRINT *,' NG,K,L,M',NG,K,L,M
               PRINT *,'RO(K),RO(L),RO(M)',RO(K+1),RO(L+1),RO(M+1)
               PRINT *,' WEIGHT(NG)',WEIGHT(NG)
              END DO
C End of I=1234...
              ENDIF
C	      STOP
	    ENDIF
	  ENDIF
	  NG1=NG2+1
1	CONTINUE
 
	  IF(INULL.GT.0)THEN
	    WRITE(2,79) INULL,SIGMA_NULL,1./SIGMA_NULL 
	    WRITE(6,79) INULL,SIGMA_NULL,1./SIGMA_NULL 
79	    FORMAT(' RECURSIVE/Warning: Null modulus (i.e., spectral'
     1     ,' list is reduced) for',I5,' values',/,
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
     1      LOWER_MODULUS,IFERMAX,BMASK,WEIGHT)
        REAL*8 SUMSQR,SUMSQI
        REAL*4 LOWER_MODULUS
	REAL BMASK(*),WEIGHT(*),SIGR,SIGI,SIGM(*),XNUMB,SIG_MAX
	INTEGER K,L,M,NG1,N_OUT,IFERMAX
 
	NG1=1
        N_OUT=0
        SUM_WEIGHTS=0.
        SIGM(4)=5.E-2
        SIGM(5)=5.E-2
 
        OPEN(3,FILE='sigma.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
	DO 1 NB=3,NBETA
C  NG2=NGT(NB)
          CALL COVER_NGT(NG2,NB)
          WRITE(3,*) NB,SIGM(NB+1)
C Caution to avoid SIGM = 0, and problems when computing 1/SIGM in SNR map:
          SIGM(NB+1)=MAX(5.E-2,SIGM(NB+1))
C Troncation to reduce the range of weights:
          IF(SIGM(NB+1).GE.SIG_MAX.OR.BMASK(NB+1).EQ.0.) THEN
              N_OUT=N_OUT+1
              BMASK(NB+1)=0.
	      DO 3 NG=NG1,NG2
                WEIGHT(NG)=0.
3             CONTINUE
            ELSE
	      DO 2 NG=NG1,NG2
C Previous weight divided by this sigma: 
C I have tried without dividing: gives worse results 
C I have tried with **1, **0.8, **0.2 instead: gives worse results 
C JLP93
               WEIGHT(NG)=WEIGHT(NG)/SQRT(SIGM(NB+1))
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
     1      LOWER_MODULUS,IFERMAX,BMASK,WEIGHT)
	REAL BMASK(*),WEIGHT(*),SIGM(*)
        REAL*4 LOWER_MODULUS
	REAL SIGR,SIGI,XNUMB,SIG_MAX
        REAL*8 SUM_WEIGHTS
	INTEGER K,L,M,NG1,N_OUT,IFERMAX,NBETA,NGAMMA
 
	NG1=1
        N_OUT=0
        SUM_WEIGHTS=0.
        SIGM(3+1)=5.E-2
        SIGM(4+1)=5.E-2
 
        OPEN(3,FILE='sigma.dat',STATUS='UNKNOWN',ACCESS='SEQUENTIAL')
	DO 1 NB=3,NBETA
C  NG2=NGT(NB)
          CALL COVER_NGT(NG2,NB)
          WRITE(3,*) NB,SIGM(NB+1)
C Caution to avoid SIGM = 0, and problems when computing 1/SIGM in SNR map:
          SIGM(NB+1)=MAX(5.E-2,SIGM(NB+1))
C Troncation to reduce the range of weights:
          IF(SIGM(NB+1).GE.SIG_MAX.OR.BMASK(NB+1).EQ.0.) THEN
              N_OUT=N_OUT+1
              BMASK(NB+1)=0.
	      DO 3 NG=NG1,NG2
                WEIGHT(NG)=0.
3             CONTINUE
          ELSE
	      DO 2 NG=NG1,NG2
C JLP93
C (Unity weights: gives worse results when not dividing by SIGM)
               WEIGHT(NG)=WEIGHT(NG)/SIGM(NB+1)
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
	SUBROUTINE EPLUS(NBETA,BMASK,XC,X,IDIM_X,VECXY)
	COMPLEX XC(*),CC
	REAL X(IDIM_X,4),VECXY(IDIM_X,2),BMASK(*),COS1,SIN1
 
C VECTEURS DE BASE DU NOYAU
	CALL NOYAU(NBETA,BMASK,VECXY,IDIM_X)
 
C ALPHA_0 EN X(.,1)
	DO 1 NB=0,NBETA
	    CC=XC(NB+1)
	    X(NB+1,1) = BMASK(NB+1)*ATAN2( IMAG(CC), REAL(CC) )
1	CONTINUE
 
C ALPHA_0^+  EN X(.,1)
	CALL PROJ_EPLUS(NBETA,1,1,X,IDIM_X,VECXY)
 
C FACTEUR DE PHASE INITIAL CORRESPONDANT
	DO 2 NB=0,NBETA
	  ALPHA=X(NB+1,1)
          COS1=BMASK(NB+1)*COS(ALPHA)
          SIN1=BMASK(NB+1)*SIN(ALPHA)
	  XC(NB+1) = CMPLX(COS1,SIN1)
2	CONTINUE
 
	RETURN
	END
C*******************************************************************
C NOYAU:
C DEFINIT LES VECTEURS DE BASE DU NOYAU
C*******************************************************************
	SUBROUTINE NOYAU(NBETA,BMASK,VECXY,IDIM_X)
	REAL VECXY(IDIM_X,2),XNORM,BMASK(*)
	REAL*8 SUM
 
C The base vectors are formed with the coordinates of the uv vectors. 
	DO 1 NB=0,NBETA
C  VECXY(NB,1)=BMASK(NB)*IXY(1,NB)
C  VECXY(NB,2)=BMASK(NB)*IXY(2,NB)
            CALL COVER_IXY(IXY1,IXY2,NB)
	    VECXY(NB+1,1)=BMASK(NB+1)*FLOAT(IXY1)
	    VECXY(NB+1,2)=BMASK(NB+1)*FLOAT(IXY2)
1	CONTINUE

C Normalization of the first vector:
	CALL NORMALIZE_L2(VECXY(1+1,1),NBETA,XNORM)

C Orthogonalization:
	SUM=0.
	DO NB=0,NBETA
        SUM=SUM+VECXY(NB+1,1)*VECXY(NB+1,2)
	END DO

	DO 2 NB=0,NBETA
	    VECXY(NB+1,2)=VECXY(NB+1,2)-SUM*VECXY(NB+1,1)
2	CONTINUE

C Normalization of the second vector:
	CALL NORMALIZE_L2(VECXY(1+1,2),NBETA,XNORM)
	
	RETURN
	END
C*******************************************************************
C PROJ_EPLUS
C Projects X(.,N1) onto E^+ and stores result in X(.,N2)
C*******************************************************************
	SUBROUTINE PROJ_EPLUS(NBETA,N1,N2,X,IDIM_X,VECXY)
	INTEGER NBETA,N1,N2
	REAL X(IDIM_X,4),VECXY(IDIM_X,2),DELTA
	REAL*8 XX,YY
 
C Scalar product of X(.,N1) with the base vectors of the kernel
	XX=0.
	YY=0.
 
	DO 1 NB=1,NBETA
	  XX=XX+X(NB+1,N1)*VECXY(NB+1,1)
	  YY=YY+X(NB+1,N1)*VECXY(NB+1,2)
1	CONTINUE
 
	DELTA = SQRT(XX**2 + YY**2)
C	WRITE(6,*) ' Distance to E+: DELTA =',DELTA
C	WRITE(2,*) ' Distance to E+: DELTA =',DELTA
 
C Projection onto E^+ and stored in X(.,N2)
	DO NB=1,NBETA
	    X(NB+1,N2) = X(NB+1,N1) - XX*VECXY(NB+1,1) - YY*VECXY(NB+1,2)
        END DO
 
	RETURN
	END
C*******************************************************************
C ATRANS compute the second member of the system to be solved.
C This member is stored in X(.,3) which is the initial
C residual R_0
C
C  R_0 = AT ( PSI - A PHI_0)
C Here the residual is IMAG(exper.bispectrum * CONJ(estimated bispectrum))) :
C*******************************************************************
	SUBROUTINE ATRANS(NBETA,NGAMMA,QMOY,IFERMAX,BMASK,YCE,WEIGHT,XC,
     1                    X,IDIM_X)
	COMPLEX XC(*),YCE(*),C1,C2,C3
	INTEGER NBETA,NGAMMA
	REAL*8 Q2
	REAL QMOY,X(IDIM_X,4),WEIGHT(*),BMASK(*)
 
C Initialization to zero:
	DO 1 NB=1,NBETA
	  X(NB+1,3)=0.
1	CONTINUE
 
C Weighted quadratic measure:
	Q2=0.
 
	NG1=1
	DO 3 NB=3,NBETA
C  NG2=NGT(NB)
          CALL COVER_NGT(NG2,NB)
	  INDEX=0
	    DO 2 NG=NG1,NG2
	      INDEX=INDEX+1
	        IF(INDEX.LE.IFERMAX.AND.BMASK(NB+1).NE.0.)THEN
C  K=KLM(1,NG)
C  L=KLM(2,NG)
C  M=KLM(3,NG)
                  CALL COVER_KLM(K,1,NG)
                  CALL COVER_KLM(L,2,NG)
                  CALL COVER_KLM(M,3,NG)
C Computed bispectrum (from the spectrum of the previous iteration):
	          C1=XC(K+1)*XC(L+1)*CONJG(XC(M+1))
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
C Computing AT [YY] by adding the contribution to X(.,3):
	          X(K+1,3)=X(K+1,3)+YY
	          X(L+1,3)=X(L+1,3)+YY
	          X(M+1,3)=X(M+1,3)-YY
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
	SUBROUTINE ATRANS_A(NBETA,NGAMMA,N1,N2,IFERMAX,BMASK,WEIGHT,
     1                      X,IDIM_X)
	REAL X(IDIM_X,4),WEIGHT(*),BMASK(*)
	INTEGER IFERMAX
 
C Initialization to zero:
	DO 1 NB=1,NBETA
	  X(NB+1,N2)=0.
1	CONTINUE
 
	NG1=1
	DO 3 NB = 3,NBETA
C  NG2=NGT(NB)
          CALL COVER_NGT(NG2,NB)
  	  INDEX=0
	    DO 2 NG=NG1,NG2
	      INDEX=INDEX+1
	      IF (INDEX.LE.IFERMAX.AND.BMASK(NB+1).NE.0.)THEN
C  K=KLM(1,NG)
C  L=KLM(2,NG)
C  M=KLM(3,NG)
                CALL COVER_KLM(K,1,NG)
                CALL COVER_KLM(L,2,NG)
                CALL COVER_KLM(M,3,NG)
C JLP93: Put WEIGHT=g**2 for AT A
	        YY = (X(K+1,N1)+X(L+1,N1)-X(M+1,N1))*WEIGHT(NG)
C Computing AT [YY] by adding the contribution to X(.,2):
	        X(K+1,N2) = X(K+1,N2) + YY
	        X(L+1,N2) = X(L+1,N2) + YY
	        X(M+1,N2) = X(M+1,N2) - YY
	      ENDIF
2	    CONTINUE
	  NG1=NG2+1
3	CONTINUE
 
	RETURN
	END
C*******************************************************************
C CGRADIENT 
C Resolution of the linear system: 
C                  [ATwA] X(.,1)  = X(.,3)
C with conjugate gradients
C X(.,1) is the unknown (and then solution) PHI 
C X(.,2) is the direction D_n
C X(.,3) is the second member of the equation, and then the residual
C X(.,4) is ATwA D_n   (called Z_n)
C*******************************************************************
	SUBROUTINE CGRADIENT(NBETA,NGAMMA,IT,QMOY,IFERMAX,BMASK,
     1                       YCE,WEIGHT,XC,X,IDIM_X)
	REAL*8 SS,R2_RN,R2_RNPLUS1,R2_PHIPLUS1
	REAL X(IDIM_X,4),OMEGA_N,SUP,WORK
	REAL GAMMA_N,QMOY,BMASK(*),WEIGHT(*)
        INTEGER IFERMAX,IT
        COMPLEX YCE(*),XC(*)

C Good problem, with 1.E-08 implies 4 iterations,
C Badly conditionned problem with Nc=25 implies around 40 iterations, so:
	ITMAX=40
 
C STEP 0
C The starting solution PHI_0 is null : X(.,1)=0.
 
C Compute the initial residual R_0 and store it in X(.,3)
	CALL ATRANS(NBETA,NGAMMA,QMOY,IFERMAX,BMASK,YCE,WEIGHT,XC,
     1              X,IDIM_X)
C Note that X(NB+1,3) is null when BMASK(NB+1) is null.

C Compute R2_RN : square of the norm of R_0
C  R2_RN is the square norm of R_N
	R2_RN=0.
	DO 1 NB=1,NBETA
C The initial solution is set to 0.	
	     X(NB+1,1) = 0.
C The first direction D_0 = R_0 is copied to X(.,2)
	     X(NB+1,2) = X(NB+1,3)
 	     R2_RN = R2_RN + X(NB+1,3)*X(NB+1,3)
1	CONTINUE
          IF(R2_RN.LT.1.E-15)THEN
            WRITE(6,39)
39          FORMAT(' CGRADIENT/Error: Square norm of Residual_{N}=0 !')
            RETURN
          ENDIF
 
	DO 2 IT=1,ITMAX
 
C STEP 1
C Compute Z_N = [AT A] D_N and store it in X(.,4)
	  CALL ATRANS_A(NBETA,NGAMMA,2,4,IFERMAX,BMASK,WEIGHT,X,IDIM_X)

C Compute OMEGA_N = R2_RN / (D_N scalar Z_N)
C Warning: SS is very small
	  SS=0.
	  DO 3 NB=1,NBETA
	      SS = SS + X(NB+1,2)*X(NB+1,4)*BMASK(NB+1)
3	  CONTINUE
	  OMEGA_N=R2_RN/SS

C Compute the residual and the next value of PHI
C  R_{N+1} put to X(.,3) and  PHI_{N+1} put to X(.,1) ;
C  R2_RNPLUS1 is the square norm of R_{N+1}
	  R2_RNPLUS1=0.
	  R2_PHIPLUS1=0.
	  DO 4 NB=1,NBETA
C R_[N+1] = R_N - OMEGA_N * Z_N
	     X(NB+1,3) = X(NB+1,3) - OMEGA_N*X(NB+1,4)
C PHI_[N+1] = PHI_N + OMEGA_N * D_N
 	     X(NB+1,1) = X(NB+1,1) + OMEGA_N*X(NB+1,2)
	     R2_RNPLUS1 = R2_RNPLUS1+X(NB+1,3)*X(NB+1,3)*BMASK(NB+1)
	     R2_PHIPLUS1=R2_PHIPLUS1+X(NB+1,1)*X(NB+1,1)*BMASK(NB+1)
4	  CONTINUE
 
C STEP2 : Exit test
          IF(R2_PHIPLUS1.LT.1.E-15)THEN
            WRITE(6,29)
            WRITE(2,29)
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
C Compute next conjugate direction D_{N+1} which is stored in X(.,2)
CJLP91	  SUM=0
	  DO 5 NB=1,NBETA
	    X(NB+1,2) = X(NB+1,3) + GAMMA_N*X(NB+1,2)
C JLP91 Test if D_N and D_N+1 are orthogonal directions:
CJLP91	    SUM = SUM + X(NB+1,2)*X(NB+1,4)
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
	SUBROUTINE ERROR_SIMU(NBETA,E,EP,ERRNAME,BMASK,RO,XCR,XC)
	COMPLEX XC(*),XCR(*),CC
	REAL RO(*),BMASK(*)
        REAL E,EP,R2,RR2
	REAL*8 SR2
        INTEGER*4 INDX,NB,NBETA
	CHARACTER ERRNAME*40
 
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
         IF(BMASK(NB+1).NE.0.)THEN
          INDX=INDX+1
          R2 = RO(NB+1)**2
	  SR2 = SR2 + R2
	  CC = XCR(NB+1) - XC(NB+1)
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
C SORTIE : CREE LES FICHIERS DE SORTIE RES ET IMS
C*******************************************************************
	SUBROUTINE SORTIE(NBETA,ISORTIE,FNAME,BMASK,RO,XC)
	PARAMETER(IDIM=256)
	COMPLEX XC(*),CC
	REAL RE,IM,RO(*),BMASK(*)
	INTEGER NX,NY
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
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
C JLP94, try...
C  XM = BMASK(NB+1)*RO(NB+1)
	  XM = RO(NB+1)
	  CC = XC(NB+1)
	  XR = XM*REAL(CC)
	  XI = XM*IMAG(CC)
C IX=IXY(1,NB)
C IY=IXY(2,NB)
          CALL COVER_IXY(IX,IY,NB)
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
     1  NGAMMA_MAX,FNAME,LOWER_MODULUS,BMASK,RO,YCE,WEIGHT,XC)
	PARAMETER(IDIM=256)
        REAL RO(*),BMASK(*)
	REAL RE,IM,WEIGHT(*),LOWER_MODULUS
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	INTEGER MADRID(1),NX,NY,PNTR_BISP,MEMSIZE
	COMPLEX XC(*),YCE(*)
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
C For virtual memory:
        COMMON /VMR/MADRID
 
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
C IIX = IXC + IXY(1,NB)
C IIY = IYC + IXY(2,NB)
          CALL COVER_IXY(IIX,IIY,NB)
	  IIX = IXC + IIX
	  IIY = IYC + IIY
C Modulus (here that is the real part) :
	  RO(NB+1)=RE(IIX,IIY)
C Phase factor:
	  XC(NB+1)=(1.,0.)
1	CONTINUE
 
	DO I=0,4
	   WRITE(6,23) I,RO(I+1)
	   WRITE(2,23) I,RO(I+1)
23	   FORMAT(' RO(',I2,') = ',F8.5)
	END DO
 
C Troncation:
	IREM=0
C Mask to discard some values of the spectral list:
        DO 91 NB=0,NBETA
          RO(NB+1)=RO(NB+1)/RO(0+1)
C Check that RO(NB+1) greater than LOWER_MODULUS:
C Remember that LOWER_MODULUS can be negative...
	  IF(RO(NB+1).LT.LOWER_MODULUS.OR.
     1       RO(NB+1).LE.0) THEN
             BMASK(NB+1)=0.
	     IREM=IREM+1
          ELSE
             BMASK(NB+1)=1.
          ENDIF
91      CONTINUE
 
        WRITE(6,47) NBETA-IREM
        WRITE(2,47) NBETA-IREM
47	FORMAT(' Number of terms of the spectral list',
     1  ' after correction: ',I4)

C FACTEUR DE PHASE DU SPECTRE CALE EN TRANSLATION : XC(.)
	CALL TRANS(IR,NBETA,BMASK,XC)

C FACTEUR DE PHASE DU BISPECTRE : YCE
	WRITE(6,41)
41	FORMAT(' Bispectrum of the image (phase term)',
     1	' (NOT centered in the frame)')
C Three lines (real, imag, snr) for the input bispectrum:
	NAME=' '
	CALL JLP_VM_READIMAG(PNTR_BISP,NX1,NY1,NAME,COMMENTS)
	WRITE(2,58) NAME(1:14),COMMENTS(1:30)
58	FORMAT(' Bispectrum: ',A14,' comments: ',A30)
	IF((NY1.LT.2).OR.(NX1.NE.NGAMMA_MAX))THEN
	   WRITE(6,45)
	   WRITE(2,45)
45	   FORMAT(' CREYCE1/FATAL ERROR: Size of bispectrum inconsistent',
     1	' with IRMAX')
	   STOP
	ENDIF

C Load BISP array to YCE array:
	CALL CREYCE_LOAD_BISP(MADRID(PNTR_BISP),NGAMMA,NGAMMA_MAX,YCE)
 
C Computing the weights:
	IF(CTE.GT.0)THEN
	   CALL BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSEIF(CTE.LT.0)THEN
C Version with SNR stored in 3rd line of bispectrum:
	  CALL BISP_WEIGHT2(MADRID(PNTR_BISP),NBETA,NGAMMA,
     1                      NGAMMA_MAX,IR,CTE,WEIGHT)
C	  CALL BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	  CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
        ELSE 
	   WRITE(6,34)
	   WRITE(2,34)
34         FORMAT(' Weights set to unity, and then normalized')
	   DO 35 NG=1,NGAMMA
	     WEIGHT(NG)=1.
	     DO 36 KK=1,3
C  K = KLM(KK,NG)
               CALL COVER_KLM(K,KK,NG)
	       IF(BMASK(K+1).EQ.0.)WEIGHT(NG)=0.
36           CONTINUE
35        CONTINUE
	  CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ENDIF
 
        MEMSIZE=NGAMMA_MAX*3*4
        CALL JLP_FREEVM(PNTR_BISP,MEMSIZE)
	RETURN
	END
C*******************************************************************
C CREYCE_LOAD_BISP
C Load input data from file to YCE array
C (Needed to be called by CREYCE1 because of dynamical allocation MADRID(...))
C*******************************************************************
	SUBROUTINE CREYCE_LOAD_BISP(BISP,NGAMMA,NGAMMA_MAX,YCE)
        REAL BISP(NGAMMA_MAX,3)
        REAL XR,XI,XM
	INTEGER NG,NGAMMA,INULL
	COMPLEX YCE(*)

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
55	CONTINUE
 
        IF(INULL.GT.0)THEN
          WRITE(2,81) INULL
          WRITE(6,81) INULL
81	  FORMAT(' CREYCE_LOAD_BISP/Warning: null bispectrum for',I5,' values')
	ENDIF

	DO I=1,5
	   WRITE(2,*) ' YCE(',I,') =',YCE(I)
	   WRITE(6,*) ' YCE(',I,') =',YCE(I)
	END DO
 
	RETURN
	END
C*******************************************************************
C CREYCE2
C Reads the modulus (squared) and bispectrum derived from simulations,
C re and im (available since it is a simulation)
C*******************************************************************
	SUBROUTINE CREYCE2(IR,CTE,NBETA,NGAMMA,NGAMMA_MAX,FNAME,
     1                     LOWER_MODULUS,BMASK,RO,YCE,WEIGHT,XCR,XC)
	PARAMETER(IDIM=256)
	REAL*8 ERRMOD,ERRBISPW
        REAL RO(*),BMASK(*)
	REAL MODSQ(IDIM,IDIM)
	REAL PNTR_BISP,RE,IM,WEIGHT(*),LOWER_MODULUS
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
	INTEGER NX,NY,MADRID(1),MEMSIZE
	COMPLEX XC(*),YCE(*),XCR(*),YY1
	COMMON /C10/RE(IDIM,IDIM),IM(IDIM,IDIM),NX,NY
        COMMON /VMR/MADRID
 
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
C IIX = IXC + IXY(1,NB)
C IIY = IYC + IXY(2,NB)
          CALL COVER_IXY(IIX,IIY,NB)
	  IIX = IXC + IIX
	  IIY = IYC + IIY
	  XR = RE(IIX,IIY)
	  XI = IM(IIX,IIY)
	  XM = SQRT(XR*XR + XI*XI)
C Phase factor:
	  IF(XM.EQ.0.)THEN
	    XC(NB+1) = (1.,0.)
            INULL=INULL+1
C	    WRITE(6,78) NB
78	    FORMAT(' Warning: XC(',I5,') is null!') 
	  ELSE
	    XC(NB+1) = CMPLX(XR/XM,XI/XM)
	  ENDIF
C Modulus:
	  W1=MAX(0.,MODSQ(IIX,IIY))
	  RO(NB+1) = SQRT(W1)
	  IF(W1.EQ.0)THEN
            INULL=INULL+1
C	    WRITE(6,75) NB
75	    FORMAT(' CREYCE2/Warning: RO(',I5,') is null!') 
	  ENDIF
	  WORK=(XM/XW1-RO(NB+1)/XW2)**2
	  ERRMOD=ERRMOD+WORK
          WRITE(9,*) NB,SQRT(WORK),RO(NB+1)/XW2,XM/XW1
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
C IIX = IXC + IXY(1,NB)
C IIY = IYC + IXY(2,NB)
           CALL COVER_IXY(IIX,IIY,NB)
           IIX = IXC + IIX
           IIY = IYC + IIY
	   XR = RE(IIX,IIY)
	   XI = IM(IIX,IIY)
	   XM = SQRT(XR*XR + XI*XI)/XW1
	   W1 = RO(I+1)/XW2
	   WRITE(6,67) I,W1,XM
	   WRITE(2,67) I,W1,XM
67	   FORMAT(' Normalised RO and computed modulus:',I5,
     1     2(1X,F8.5))
	END DO
 
C Troncation:
	IREM=0
C Mask to discard some values of the spectral list:
        DO 91 NB=0,NBETA
          RO(NB+1)=RO(NB+1)/RO(0+1)
C Check that RO(NB+1) greater than LOWER_MODULUS:
C Remember that LOWER_MODULUS can be negative...
	  IF(RO(NB+1).LT.LOWER_MODULUS.OR.
     1       RO(NB+1).LE.0) THEN
             BMASK(NB+1)=0.
	     IREM=IREM+1
          ELSE
             BMASK(NB+1)=1.
          ENDIF
91      CONTINUE
 
        WRITE(6,47) NBETA-IREM
        WRITE(2,47) NBETA-IREM
47	FORMAT(' Number of terms of the spectral list',
     1  ' after correction: ',I4)

C FACTEUR DE PHASE DU SPECTRE CALE EN TRANSLATION : XC(.)
	CALL TRANS(IR,NBETA,BMASK,XC)

C Reference (since simulation)
	DO 3 NB=1,NBETA
	  XCR(NB+1)=XC(NB+1)
C Set initial solution XC to zero:
C	  XC(NB+1)=(1.,0.)
3	CONTINUE
 
C***********************************************************************
C FACTEUR DE PHASE DU BISPECTRE : YCE
	WRITE(6,41)
41	FORMAT(' Bispectrum of the image (phase term)',
     1	' (NOT centered in the frame)')
	NAME=' '
	CALL JLP_VM_READIMAG(PNTR_BISP,NX1,NY1,NAME,COMMENTS)
	WRITE(2,58) NAME(1:14),COMMENTS(1:30)
58	FORMAT(' Bispectrum: ',A14,' comments: ',A30)
	IF((NY1.LT.2).OR.(NX1.NE.NGAMMA_MAX))THEN
	   WRITE(6,45)
	   WRITE(2,45)
45	   FORMAT(' CREYCE2/FATAL ERROR: Size of bispectrum inconsistent',
     1	' with IRMAX')
	   STOP
	ENDIF

C Load BISP array to YCE array:
	CALL CREYCE_LOAD_BISP(MADRID(PNTR_BISP),NGAMMA,NGAMMA_MAX,YCE)

C Computing the weights:
	IF(CTE.GT.0)THEN
	   CALL BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ELSEIF(CTE.LT.0)THEN
C Version with SNR stored in 3rd line of bispectrum:
	  CALL BISP_WEIGHT2(MADRID(PNTR_BISP),NBETA,NGAMMA,
     1                      NGAMMA_MAX,IR,CTE,WEIGHT)
C	  CALL BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	  CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
        ELSE 
	   WRITE(6,34)
	   WRITE(2,34)
34         FORMAT(' Weights set to unity, and then normalized')
	   DO NG=1,NGAMMA
	     WEIGHT(NG)=1.
	     DO KK=1,3
C  K = KLM(KK,NG)
               CALL COVER_KLM(K,KK,NG)
	       IF(BMASK(K+1).EQ.0.)WEIGHT(NG)=0.
	     END DO
	   END DO
	   CALL NORMALIZE_L1(WEIGHT,NGAMMA,XNORM)
	ENDIF
 
C Compute the errors (since XCR is available):
	ERRBISP=0.
	ERRBISPW=0.
	DO 55 NG=1,NGAMMA
C  YY1=XCR(KLM(1,NG))*XCR(KLM(2,NG))*CONJG(XCR(KLM(3,NG)))
          CALL COVER_KLM(KLM1,1,NG)
          CALL COVER_KLM(KLM2,2,NG)
          CALL COVER_KLM(KLM3,3,NG)
C Theoretical bispectrum:
	  YY1=XCR(KLM1+1)*XCR(KLM2+1)*CONJG(XCR(KLM3+1))
C True error:
	  WORK=(REAL(YY1)-REAL(YCE(NG)))**2
     1	+(IMAG(YY1)-IMAG(YCE(NG)))**2
	  ERRBISP=ERRBISP+WORK
C          IF(WORK.GT.3.99)THEN
C	       WRITE(6,688) WORK,KLM1,KLM2,KLM3,
C     1          NG,REAL(YCE(NG)),IMAG(YCE(NG)),
C     1         REAL(YY1),IMAG(YY1),BISP(NG,3)
C	       WRITE(2,688) WORK,KLM1,KLM2,KLM3,
C     1          NG,REAL(YCE(NG)),IMAG(YCE(NG)),
C     1         REAL(YY1),IMAG(YY1),BISP(NG,3)
C688	   FORMAT(' Big error: WORK=',F8.5,' K,L,M=',3(1X,I5),/,
C     1     ' Input, comp bisp & snr:',I5,
C     1     2(1X,F8.5),1X,2(1X,F8.5),1X,1PG10.3)
C	END DO
C          ENDIF
	  ERRBISPW=ERRBISPW+WORK*WEIGHT(NG)
55	CONTINUE
	ERRBISP=SQRT(ERRBISP/FLOAT(NGAMMA))
	ERRBISPW=SQRT(ERRBISPW)
	WRITE(2,44) ERRBISP,ERRBISPW
	WRITE(6,44) ERRBISP,ERRBISPW
44	FORMAT(' rms error of the bispectrum :',1PE11.3,/,
     1  ' weighted rms error of the bispectrum :',1PE11.3)
 
C Just debug mode:
	DO I=1,5
C YY1=XCR(KLM(1,I))*XCR(KLM(2,I))*CONJG(XCR(KLM(3,I)))
          CALL COVER_KLM(KLM1,1,NG)
          CALL COVER_KLM(KLM2,2,NG)
          CALL COVER_KLM(KLM3,3,NG)
	  YY1=XCR(KLM1+1)*XCR(KLM2+1)*CONJG(XCR(KLM3+1))
	   WRITE(6,68) I,REAL(YCE(I)),IMAG(YCE(I)),
     1         REAL(YY1),IMAG(YY1)
	   WRITE(2,68) I,REAL(YCE(I)),IMAG(YCE(I)),
     1         REAL(YY1),IMAG(YY1)
68	   FORMAT(' Input and computed bisp :',I2,
     1     2(1X,F8.5),1X,2(1X,F8.5))
	END DO
 
        MEMSIZE=NGAMMA_MAX*3*4
        CALL JLP_FREEVM(PNTR_BISP,MEMSIZE)
	RETURN
	END
C*****************************************************************
C BISP_WEIGHT0
C Computing the weights:
C JLP Version, not very good...
C*****************************************************************
	SUBROUTINE BISP_WEIGHT0(NGAMMA,IR,CTE,WEIGHT)
	REAL WEIGHT(*),CTE
	REAL DEGRAD,SIGB,SIG2,S22,S33,WORK
	INTEGER NGAMMA,IR

	DEGRAD = 3.141592/180.
	CTE1 = CTE*DEGRAD/FLOAT(IR**2)
 
	DO 4 NG=1,NGAMMA
	  S22=IR*IR
          S33=0.
	     DO 6 KK=1,3
C  K = KLM(KK,NG)
               CALL COVER_KLM(K,KK,NG)
C  WORK = IXY(1,K)**2 + IXY(2,K)**2
               CALL COVER_IXY(IXY1,IXY2,K)
 	       WORK = IXY1**2 + IXY2**2
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
	SUBROUTINE BISP_WEIGHT11(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
	REAL WEIGHT(*),CTE,RO(*),BMASK(*)
	REAL DEGRAD,SIGB,SIG2
	INTEGER NGAMMA,IS2,IR

	DEGRAD = 3.141592/180.
	CTE1 = CTE*DEGRAD/FLOAT(IR**2)
 
	DO 4 NG=1,NGAMMA
	  IS2=0
C  WORK=RO(KLM(1,NG))*RO(KLM(2,NG))*RO(KLM(3,NG))
          CALL COVER_KLM(KLM1,1,NG)
          CALL COVER_KLM(KLM2,2,NG)
          CALL COVER_KLM(KLM3,3,NG)
	  WORK=RO(KLM1+1)*RO(KLM2+1)*RO(KLM3+1)
          WORK=WORK*BMASK(KLM1+1)*BMASK(KLM2+1)*BMASK(KLM3+1)
          IF(WORK.EQ.0)THEN
            WEIGHT(NG)=0.
          ELSE
	     DO 6 KK=1,3
C  K = KLM(KK,NG)
               CALL COVER_KLM(K,KK,NG)
C IS2 = IS2 + IXY(1,K)**2 + IXY(2,K)**2
               CALL COVER_IXY(IXY1,IXY2,K)
 	       IS2 = IS2 + IXY1**2 + IXY2**2
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
	SUBROUTINE BISP_WEIGHT2(BISP,NBETA,NGAMMA,NGAMMA_MAX,IR,CTE,WEIGHT)
	REAL WEIGHT(*),CTE,SUM,BISP(NGAMMA_MAX,3),MEAN_SNR
	INTEGER NGAMMA,IR

	SUM=0.
	DO 4 NG=1,NGAMMA
	  WEIGHT(NG) = MAX(0.,BISP(NG,3))
C JLP93 : keep BISP, do not put any **0.5 **1.5 or **2, or anything else
C          WEIGHT(NG)=WEIGHT(NG)
          SUM=SUM+WEIGHT(NG)
4	CONTINUE
 
        MEAN_SNR = SUM/FLOAT(NGAMMA)

	WRITE(6,234) SUM,MEAN_SNR
	WRITE(2,234) SUM,MEAN_SNR
234	FORMAT(' BISP_WEIGHT2/SNR weight Initial sum:',1PG11.4,
     1  ' Mean bisp. SNR:',1PG11.4)

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
	SUBROUTINE BISP_WEIGHT22(BISP,NBETA,NGAMMA,NGAMMA_MAX,IR,CTE,BMASK,
     1                           RO,WEIGHT)
	REAL WEIGHT(*),CTE,SUM,BISP(NGAMMA_MAX,3),MEAN_SNR
	INTEGER NGAMMA,IR,NBETA
        REAL RO(*),BMASK(*)

C JLP93: add modulus 
C First get the relative weights according to the modulus: 
	WMAX=0.
	WMIN=RO(0+1)
        IG1=1
	DO NB=3,NBETA
C  IG2=NGT(NB)
          CALL COVER_NGT(IG2,NB)
	  DO NG=IG1,IG2
            WEIGHT(NG)=BMASK(NB+1)*RO(NB+1)
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
	  WEIGHT(NG) = MAX(0.,BISP(NG,3))
C JLP93:
	  WEIGHT(NG) = WEIGHT(NG)*MAX(0.,BISP(NG,3))
          IF(MOD(NG,100).EQ.1.AND.NG.LT.1000.) PRINT *,' Weight2',WEIGHT(NG) 
C JLP93 : keep BISP, do not put any **0.5 **1.5 or **2, or anything else
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
C Warning: This prog does not use BMASK (---> check if this is right later!)
C**********************************************************************
	SUBROUTINE ERROR_BISPECT(NBETA,NGAMMA,ERRNAME,
     1                           ARRAY,NBX,NBY,IFERMAX,YCE,XC)
	INTEGER*4 IFERMAX,NBX,NBY
	REAL ARRAY(NBX,*)
        REAL*8 ERROR1
	CHARACTER ERRNAME*(*) 
	COMPLEX XC(*),YCE(*),C1,C2
        REAL MOD1,MOD2,RE1,IM1,ERROR0,ERROR2
 
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
C  NG2=NGT(NB)
          CALL COVER_NGT(NG2,NB)
	  INDEX=0
	    DO 2 NG=NG1,NG2
	      INDEX=INDEX+1
	        IF(INDEX.LE.IFERMAX)THEN
C  K=KLM(1,NG)
C  L=KLM(2,NG)
C  M=KLM(3,NG)
                  CALL COVER_KLM(K,1,NG)
                  CALL COVER_KLM(L,2,NG)
                  CALL COVER_KLM(M,3,NG)
C Computed bispectrum (from the spectrum of the previous iteration):
C minus experimental bispectrum:
	          C1=XC(K+1)*XC(L+1)*CONJG(XC(M+1))
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
	SUBROUTINE BISP_WEIGHT1(NBETA,NGAMMA,IR,CTE,BMASK,RO,WEIGHT)
        REAL RO(*),BMASK(*)
	REAL WEIGHT(*),CTE,W1,WMIN,WMAX,RANGE_MAX1,RANGE_MAX2
	INTEGER NGAMMA,IR

C 20 is too small, even 40 is a bit small...
	RANGE_MAX1=800.
	RANGE_MAX2=100000.
	WRITE(2,23) RANGE_MAX1,RANGE_MAX2
	WRITE(6,23) RANGE_MAX1,RANGE_MAX2
23	FORMAT(' Weight3/Modulus, R1*R2*R3/sum*radius**2',/,
     1  ' Maximum range per column, and global: ',2(1X,F8.2))

C First get the relative weights according to the modulus: 
	WMAX=0.
	WMIN=RO(0+1)*RO(0+1)*RO(0+1)
	DO 203 NG=1,NGAMMA
C  WORK=BMASK(KLM(1,NG))*BMASK(KLM(2,NG))*BMASK(KLM(3,NG))
          CALL COVER_KLM(KLM1,1,NG)
          CALL COVER_KLM(KLM2,2,NG)
          CALL COVER_KLM(KLM3,3,NG)
          WORK=BMASK(KLM1+1)*BMASK(KLM2+1)*BMASK(KLM3+1)
	  WEIGHT(NG)=WORK*RO(KLM1+1)*RO(KLM2+1)*RO(KLM3+1)
          IF(WEIGHT(NG).NE.0.)WMIN=MIN(WMIN,WEIGHT(NG))
          WMAX=MAX(WMAX,WEIGHT(NG))
C          WEIGHT(NG)=SQRT(WEIGHT(NG))
C	  WEIGHT(NG)=MIN(RO(KLM1+1),RO(KLM2+1))
C	  WEIGHT(NG)=MIN(WEIGHT(NG),RO(KLM3+1))
203     CONTINUE 

	WRANGE=WMAX/WMIN
	WRITE(2,24) WRANGE
	WRITE(6,24) WRANGE
24	FORMAT(' Initial range of the weights:',G12.5)

C Now modulation according to the column (ie NB)
	WMAX=0.
	WMIN=RO(0+1)*RO(0+1)*RO(0+1)
	IG1=1
	DO NB=3,NBETA
C  IG2=NGT(NB)
          CALL COVER_NGT(IG2,NB)

C Troncation in the column to regularize the problem:
	  WWMIN=RO(0+1)*RO(0+1)*RO(0+1)
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
C  XW=4 + IXY(1,NB)**2+IXY(2,NB)**2
          CALL COVER_IXY(IXY1,IXY2,NB)
	  XW=4 + IXY1**2 + IXY2**2
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

	 SSUM=0.
	 DO I=1,NPTS
	   SSUM=SSUM+ABS(ARRAY(I))
	 END DO

C Normalisation of the weigths:
	IF(SSUM.EQ.0.)THEN
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
C To normalize a vector with L2 norm
C Called to normalize the two vectors of the Kernel
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
	SUBROUTINE EIGEN_VALUE1(NBETA,NGAMMA,XLAMBDA1,IFERMAX,
     1                          BMASK,WEIGHT,X,IDIM_X,VECXY)
	REAL X(IDIM_X,4),VECXY(IDIM_X,2),XNORM,XLAMBDA1,BMASK(*),WEIGHT(*)
	INTEGER NBETA,NGAMMA,NB,I,IFERMAX

C Initialization:
	DO NB=1,NBETA
	  X(NB+1,4)=BMASK(NB+1)
	END DO

	DO I=1,20

C Projection of this error onto E^+
	  CALL PROJ_EPLUS(NBETA,3,1,X,IDIM_X,VECXY)

C Normalization: 
	  CALL NORMALIZE_L2(X(1+1,2),NBETA,XNORM)

C          WRITE(6,45) I-1,XNORM
C45	  FORMAT(' Eigen_value1/ Iteration ',I3,' Norm: ',G12.5)

C Computing ATwA of X(.,2). 
C Output in X(.,4)
          CALL ATRANS_A(NBETA,NGAMMA,2,4,IFERMAX,BMASK,WEIGHT,X,IDIM_X)

	ENDDO

	XLAMBDA1 = XNORM

	RETURN
	END
C**********************************************************
C Smallest eigen value:
C Power method applied to I - AtwA/lambda1
C**********************************************************
	SUBROUTINE EIGEN_VALUE2(NBETA,NGAMMA,XLAMBDA1,XLAMBDA2,
     1                          IFERMAX,BMASK,WEIGHT,X,IDIM_X,VECXY)
	REAL X(IDIM_X,4),VECXY(IDIM_X,2),XNORM,XLAMBDA2,BMASK(*),WEIGHT(*)
	INTEGER NBETA,NGAMMA,NB,I,IFERMAX

C Initialization:
	DO NB=1,NBETA
	  X(NB+1,2)=BMASK(NB+1)
	END DO

C Main loop:
	DO I=1,20

C Projection of this error onto E^+
C        X(1+1,2)=0.
C        X(2+1,2)=0.
	CALL PROJ_EPLUS(NBETA,2,2,X,IDIM_X,VECXY)

C Normalization: 
	  CALL NORMALIZE_L2(X(1+1,2),NBETA,XNORM)

C          WRITE(6,45) I-1,XNORM
C45	  FORMAT(' Eigen_value2/ Iteration ',I3,' Norm: ',G12.5)

C Computing ATwA of X(.,2). 
C Output in X(.,4)
	  CALL ATRANS_A(NBETA,NGAMMA,2,4,IFERMAX,BMASK,WEIGHT,X,IDIM_X)

C Now computing I - ATwA/lambda1 :
	  DO NB=1,NBETA
	    X(NB+1,2) = X(NB+1,2) - X(NB+1,4)/XLAMBDA1
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
	  CALL NORMALIZE_L2(X(1+1,2),NBETA,XNORM)

C Projection to E^+
	  CALL PROJ_EPLUS(NBETA,2,2,X,IDIM_X,VECXY)
C Computing the norm:
	  CALL NORMALIZE_L2(X(1+1,2),NBETA,XNORM)
	  WRITE(6,71) XNORM
	  WRITE(2,71) XNORM
71        FORMAT(' XNORM(proj): ',1PG12.5)
        ENDIF

	RETURN
	END
C**********************************************************************
C Output SNR map according to the adopted weights
C Generating SNR and sigma maps (needed by DIANE)
C*******************************************************************
        SUBROUTINE OUTPUT_SNR1(SIGM,NBETA,FNAME,RO)
	PARAMETER(IDIM=256)
	REAL RE,IM,RO(*),SIGM(*)
	INTEGER NX,NY
	CHARACTER NAME*40,COMMENTS*80,FNAME*40
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
C IX=IXY(1,NB)
C IY=IXY(2,NB)
          CALL COVER_IXY(IX,IY,NB)
	  IIX = IXC + IX
	  IIY = IYC + IY
          SIGG=MAX(5.E-2,SIGM(NB+1))
	  RE(IIX,IIY) = 1./SIGG
	  IM(IIX,IIY) = SIGG*RO(NB+1)
	  IIX = IXC - IX
	  IIY = IYC - IY
	  RE(IIX,IIY) = 1./SIGG
	  IM(IIX,IIY) = SIGG*RO(NB+1)
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
C Fills the mask with ones
C**********************************************************************
        SUBROUTINE FILL_MASK(MASK,NX,NY)
        INTEGER I,NXY,NX,NY
        REAL*4 MASK(*)
        NXY=NX*NY 
          DO I=1,NXY
            MASK(I)=1.
          END DO
        RETURN
        END
C**********************************************************************
C Main loop
C Least square minimization, non-linear least square fit
C
C**********************************************************************
        SUBROUTINE LSQUARES1(EXIT_TOLERANCE,NBETA,NGAMMA,IFERMAX,BMASK,YCE,
     1                       WEIGHT,XC,X,IDIM_X,ITTMAX)
        REAL WEIGHT(*),X(IDIM_X,4),BMASK(*)
	COMPLEX XC(*),YCE(*)
        REAL COSR,SINR,R,SUP,QMOY,QMOY0,EXIT_TOLERANCE,DELQ
	INTEGER IFERMAX,NX,NY,NB,NBETA,NGAMMA,ITT,ITTMAX,IT
 
	QMOY0 = 0.
 
	DO 3 ITT=1,ITTMAX
C Computes phase error term X(.,1) for PHI^+
C Iterative solution by conjugate gradients.
 
	  CALL CGRADIENT(NBETA,NGAMMA,IT,QMOY,IFERMAX,BMASK,YCE,
     1                   WEIGHT,XC,X,IDIM_X)
	  WRITE(6,77) QMOY,ITT,IT
	  WRITE(2,77) QMOY,ITT,IT
77	  FORMAT(' rms bisp error',1PG12.5,
     1    ' ITT =',I3,' Internal IT =',I4,' done')
 
C SUP : Maximum of the solution PHI^+ X(.,1) (L1 norm).
	  SUP=0.
	  DO 4 NB=1,NBETA
              R = X(NB+1,1)*BMASK(NB+1)
	      R=ABS(R)
	      SUP=AMAX1(R,SUP)
4	  CONTINUE
 
C Test of convergence: 
	  IF (SUP.LT.EXIT_TOLERANCE) THEN
	   WRITE(6,*) ' Normal exit: SUP .LT. EXIT_TOLERANCE '
	   WRITE(2,*) ' Normal exit: SUP .LT. EXIT_TOLERANCE '
	   GO TO 6
	  ENDIF
 
C Generating the new value of the spectral term XC(.)
	  DO 5 NB=0,NBETA
            R=X(NB+1,1)*BMASK(NB+1)
	    COSR=COS(R)
	    SINR=SIN(R)
	    XC(NB+1)=XC(NB+1)*CMPLX(COSR,SINR)
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

	RETURN
	END
C**********************************************************************
C FINAL_ERROR 
C Compute errors when object is available (simulations)
C
C**********************************************************************
        SUBROUTINE FINAL_ERROR(IR,NBETA,EF,EPF,ERRNAME,BMASK,RO,XCR,XC,
     1                   X,IDIM_X,VECXY,SIGM)
	COMPLEX XC(*),XCR(*),CC
	REAL X(IDIM_X,4),VECXY(IDIM_X,2)
	REAL RO(*),BMASK(*),SIGM(*)
        REAL EF,EPF,RADDEG,RR
        INTEGER NBETA,NB,IRS,IR
        CHARACTER ERRNAME*(*)

	  CALL ERROR_SIMU(NBETA,EF,EPF,ERRNAME,BMASK,RO,XCR,XC)
 
C ECART DE PHASE ANGULAIRE POINT PAR POINT EXPRIME EN RADIAN
	  DO 10 NB = 0,NBETA
	    CC = XCR(NB+1)
	    XX = ATAN2(IMAG(CC),REAL(CC))
	    CC = XC(NB+1)
	    X(NB+1,1) = BMASK(NB+1)*(XX - ATAN2(IMAG(CC),REAL(CC)))
10	  CONTINUE
 
C Projection of this error X(.,1) onto E^+
	  CALL PROJ_EPLUS(NBETA,1,1,X,IDIM_X,VECXY)
 
C ECART ANGULAIRE CORRESPONDANT AUX U CROISSANTS EN NORME LE LONG
C D'UN RAYON; FINALEMENT EXPRIME EN DEGRE
	  RADDEG = 180./3.141592
 
	  DO 11 IRS = 1,IR
C NBS = NBCOUV(IRS,0)
            CALL COVER_NBCOUV(NBS,IRS,0,IR)
	    RR = BMASK(NBS+1)*RADDEG*X(NBS+1,1)
	    WRITE(6,72) IRS,NBS,RR,SIGM(NBS+1)
	    WRITE(2,72) IRS,NBS,RR,SIGM(NBS+1)
72	    FORMAT(' #',I2,' Spect. index ',I4,' Phase err (deg) ',
     1      F8.3,' Sigma ',F8.3)
11	  CONTINUE
 
	 RETURN
         END
C*****************************************************************
C Set complex array XC to (1.,0.)
C*****************************************************************
        SUBROUTINE ZERO_PHASE(XC,NBETA)
        COMPLEX XC(*)
        INTEGER NB,NBETA

	DO 2 NB=0,NBETA
	  XC(NB+1)=(1.,0.)
2	CONTINUE

	RETURN
	END
