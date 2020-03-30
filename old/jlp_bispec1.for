C++********************************************************************
C Set of routines:
C BISPEC2, BISPEC1, OUT_BISP
C
C JLP
C Version 16-06-92
C--*******************************************************************
C*******************************************************************
C BISPEC2
C NEW VERSION: No photon noise corection, and sum of phasors only!
C
C Reads a spectrum from an real image (from observations)
C
C Computes the bispectrum and spectrum lists from RE and IM (input arrays)
C
C Integrates the squared modulus of the spectrum and
C the phase factor of the bispectrum. Corrects from photon noise effects.
C
C Photon noise correction (cf JOSA 2, 14, Wirnitzer):
C (When spectrum normalized to one in the center)
C <i(u)>**2 = E(D...)/N**2 - 1/N
C
C <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
C
C Input:
C RE(NX,NY)
C IM(NX,NY)
C
C Output:
C MODSQ: Sum of the modulus squared
C YCE1(.,1) and YCE1(.,2): phase factor of the bispectrum (real and imaginary)
C YCE1(.,3): sum of square real parts of the bispectrum
C YCE1(.,4): sum of square imaginary parts of the bispectrum
C
C*******************************************************************
	SUBROUTINE BISPEC2(RE,IM,MODSQ,SNRM,NX,NY,IDIM,
     1	YCE1,IR,NBETA,NGAMMA)
 
C  POUR IR=25
	PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL MODSQ(IDIM,*),YCE1(NGMAX,4)
	REAL RE(IDIM,*),IM(IDIM,*),SNRM(IDIM,*)
	REAL W2,W1,WORK1,WORK2,WORK
	INTEGER K1,K2,K3
	REAL WR1,WR2,WR3,WI1,WI2,WI3,AAR,AAI
	REAL XR(NULL:NBMAX),XI(NULL:NBMAX)
C	 COMPLEX CC1,XC(NULL:NBMAX)
 
C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A NBMAX )
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
	COMMON /C2/KLM(3,NGMAX)
	
C Check the FFT (zero frequency is at (1,1)):
	XPHOT1=RE(1,1)**2+IM(1,1)**2
	XPHOT1=SQRT(XPHOT1)
C	PRINT *,' BISPEC2/ XPHOT1 = ',XPHOT1
 
	DO 9 NB=0,NBETA
C IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
C element NB of the spectral list. 
C As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
C we do the following transformation:
	  IIX=MOD(IXY(1,NB)+NX,NX)+1
	  IIY=MOD(IXY(2,NB)+NY,NY)+1
	  XR(NB)=RE(IIX,IIY)
	  XI(NB)=IM(IIX,IIY)
9	CONTINUE
 
C Phase factor of the bispectrum (with bispectral list):
C and correction from photon noise effects:
	W2=2./(XPHOT1*XPHOT1)
	DO 2 NG=1,NGAMMA
C This part takes a long time:
C 	  CC1=XC(KLM(1,NG))*XC(KLM(2,NG))*CONJG(XC(KLM(3,NG)))
C     1	-XC(KLM(1,NG))*CONJG(XC(KLM(1,NG)))
C     1	-XC(KLM(2,NG))*CONJG(XC(KLM(2,NG)))
C     1	-XC(KLM(3,NG))*CONJG(XC(KLM(3,NG)))
C	  YCE1(NG,1)=REAL(CC1)+YCE1(NG,1)+2*XPHOT1
CC	  IF(NG.NE.-1234)GOTO 34
CC	  CC1=XC(KLM(1,NG))*XC(KLM(2,NG))*CONJG(XC(KLM(3,NG)))
CC	  YCE1(NG,1)=REAL(CC1)+YCE1(NG,1)
CC	  YCE1(NG,2)=AIMAG(CC1)+YCE1(NG,2)
CC         YCE1(NG,3)=1.
C New version (30% shorter):
34	  K1=KLM(1,NG)
	  K2=KLM(2,NG)
	  K3=KLM(3,NG)
	  WR1=XR(K1)
	  WR2=XR(K2)
	  WR3=XR(K3)
	  WI1=XI(K1)
	  WI2=XI(K2)
	  WI3=XI(K3)
C  CC1=XC1*XC2*CONJG(XC3): AAR*WR3+AAI*WI3
	  AAR=WR1*WR2-WI1*WI2
	  AAI=WR2*WI1+WR1*WI2
	  WORK1=AAR*WR3+AAI*WI3
	  WORK2=AAI*WR3-AAR*WI3
C CC1=XC1*XC2*CONJG(XC3)-XC1*CONJG(XC1)-XC2*CONJG(XC2)-XC3*CONJG(XC3)
C Try without noise correction:
C	  WORK1=WORK1
C     1	-(WR1*WR1+WR2*WR2+WR3*WR3
C     1	+WI1*WI1+WI2*WI2+WI3*WI3)/XPHOT1+W2
	  YCE1(NG,1)=YCE1(NG,1)+WORK1
	  YCE1(NG,2)=YCE1(NG,2)+WORK2

C Estimation of the noise with the sum of squares:
	  YCE1(NG,3)=YCE1(NG,3)+WORK1*WORK1
	  YCE1(NG,4)=YCE1(NG,4)+WORK2*WORK2
2	CONTINUE
 
C Squared modulus for the output (corrected from photon noise)
C	W1=1./XPHOT1
	DO 79 J=1,NY
	  DO 78 I=1,NX
            WORK=RE(I,J)*RE(I,J)+IM(I,J)*IM(I,J)
C Photon noise correction:
C	    MODSQ(I,J)=MODSQ(I,J)+WORK-W1
C No correction:
	    MODSQ(I,J)=MODSQ(I,J)+WORK
	    SNRM(I,J)=SNRM(I,J)+WORK*WORK
78        CONTINUE
79	CONTINUE
 
	RETURN
	END
C Estimate of the noise (valid only if sum of phasors):
C To have an approximation of the mean,
C waits until 30 frames have been processed:
C	  IF(IFRAME.LT.-30)THEN
C	   WORK3=YCE1(NG,1)/FLOAT(IFRAME)
C	   WORK4=YCE1(NG,2)/FLOAT(IFRAME)
C	   WORK5=SQRT(WORK3**2+WORK4**2)
C	   WORK1=WORK1-WORK3/WORK5
C	   WORK2=WORK2-WORK4/WORK5
C           YCE1(NG,3)=YCE1(NG,3)+WORK1*WORK1+WORK2*WORK2
C           IF(NG.EQ.1)THEN
C	     PRINT *,'MEAN YCE1,2, MEAN ERR, ERR',WORK3/WORK5,
C     1      WORK4/WORK5,YCE1(NG,3)/FLOAT(IFRAME-30),WORK1**2+WORK2**2
C	     ENDIF
C          ENDIF
C*******************************************************************
C BISPEC1
C Same as BISPEC2, but older version, and sum of full bispectrum terms.
C (amplitude and phase)
C Reads a spectrum from an real image (from observations)
C
C Computes the bispectrum and spectrum lists from RE and IM
C already computed before calling this routine.
C
C Integrates the squared modulus of the spectrum and
C the phase factor of the bispectrum. Corrects from photon noise effects.
C
C Photon noise correction (cf JOSA 2, 14, Wirnitzer):
C (When spectrum normalized to one in the center)
C <i(u)>**2 = E(D...)/N**2 - 1/N
C
C <i(3)(u,v)> = E(D(3)(u,v)/N**3 - E(D(2)(u)/N**2)/N - E(D(2)(v)... +2/N**3)
C
C Input:
C RE(NX,NY)
C IM(NX,NY)
C
C Output:
C MODSQ: Sum of the modulus squared
C YCE1: phase factor of the bispectrum (real and imaginary)
C YCE1(.,3): sum of square real parts of the bispectrum
C YCE1(.,4): sum of square imaginary parts of the bispectrum
C
C*******************************************************************
	SUBROUTINE BISPEC1(RE,IM,MODSQ,SNRM,NX,NY,IDIM,
     1	YCE1,IR,NBETA,NGAMMA)
 
C  POUR IR=25
	PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
C	PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL MODSQ(IDIM,*),YCE1(NGMAX,4)
	REAL RE(IDIM,*),IM(IDIM,*),SNRM(IDIM,*)
	REAL W2,W1,WORK1,WORK2,WORK
	INTEGER K1,K2,K3
	REAL WR1,WR2,WR3,WI1,WI2,WI3,AAR,AAI
	REAL XR(NULL:NBMAX),XI(NULL:NBMAX)
C@	COMPLEX CC1,XC(NULL:NBMAX)
 
C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A NBMAX )
	COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
	COMMON /C2/KLM(3,NGMAX)
	
C Check the FFT (zero frequency is at (1,1)):
	XPHOT1=RE(1,1)**2+IM(1,1)**2
	XPHOT1=SQRT(XPHOT1)
C	PRINT *,' BISPEC1/ XPHOT1 = ',XPHOT1
 
C Normalization (necessary for photon noise correction):
	DO J=1,NY
	  DO I=1,NX
	    RE(I,J)=RE(I,J)/XPHOT1
	    IM(I,J)=IM(I,J)/XPHOT1
	  END DO
	END DO
 
	DO 9 NB=0,NBETA
C IXY(1,NB) and IXY(2,NB) are the X,Y coordinates of the
C element NB of the spectral list. 
C As they (i.e. IXY) can be negative, and that zero frequency is at (0,0),
C we do the following transformation:
	  IIX=MOD(IXY(1,NB)+NX,NX)+1
	  IIY=MOD(IXY(2,NB)+NY,NY)+1
	  XR(NB)=RE(IIX,IIY)
	  XI(NB)=IM(IIX,IIY)
C@	  XM=XR(NB)*XR(NB)+XI(NB)*XI(NB)
C Spectrum:
C@	    XC(NB)=CMPLX(XR(NB)/XM,XI(NB)/XM)
 
9	CONTINUE
 
C Phase factor of the bispectrum (with bispectral list):
C and correction from photon noise effects:
	W2=2./(XPHOT1*XPHOT1)
	DO 2 NG=1,NGAMMA
C This part takes a long time:
C 	  CC1=XC(KLM(1,NG))*XC(KLM(2,NG))*CONJG(XC(KLM(3,NG)))
c     1	-XC(KLM(1,NG))*CONJG(XC(KLM(1,NG)))
c     1	-XC(KLM(2,NG))*CONJG(XC(KLM(2,NG)))
c     1	-XC(KLM(3,NG))*CONJG(XC(KLM(3,NG)))
c	  YCE1(NG,1)=REAL(CC1)+YCE1(NG,1)+2*XPHOT1
c	  YCE1(NG,2)=AIMAG(CC1)+YCE1(NG,2)
C New version (30% shorter):
	  K1=KLM(1,NG)
	  K2=KLM(2,NG)
	  K3=KLM(3,NG)
	  WR1=XR(K1)
	  WR2=XR(K2)
	  WR3=XR(K3)
	  WI1=XI(K1)
	  WI2=XI(K2)
	  WI3=XI(K3)
C  CC1=XC1*XC2*CONJG(XC3): AAR*WR3+AAI*WI3
	  AAR=WR1*WR2-WI1*WI2
	  AAI=WR2*WI1+WR1*WI2
	  WORK1=AAR*WR3+AAI*WI3
	  WORK2=AAI*WR3-AAR*WI3
C CC1=XC1*XC2*CONJG(XC3)-XC1*CONJG(XC1)-XC2*CONJG(XC2)-XC3*CONJG(XC3)
C With photon noise correction:
	  WORK1=WORK1
     1	-(WR1*WR1+WR2*WR2+WR3*WR3
     1	+WI1*WI1+WI2*WI2+WI3*WI3)/XPHOT1+W2
	  YCE1(NG,1)=YCE1(NG,1)+WORK1
	  YCE1(NG,2)=YCE1(NG,2)+WORK2

C Estimation of the noise with the sum of squares:
	  YCE1(NG,3)=YCE1(NG,3)+WORK1*WORK1
	  YCE1(NG,4)=YCE1(NG,4)+WORK2*WORK2
2	CONTINUE
 
C Squared modulus for the output (corrected from photon noise)
	W1=1./XPHOT1
	DO 79 J=1,NY
	  DO 78 I=1,NX
            WORK=RE(I,J)*RE(I,J)+IM(I,J)*IM(I,J)
C Photon noise correction:
            WORK=WORK-W1
	    MODSQ(I,J)=MODSQ(I,J)+WORK
	    SNRM(I,J)=SNRM(I,J)+WORK*WORK
78        CONTINUE
79	CONTINUE
 
	RETURN
	END
 
C Example how to retrieve the phase:
c	    XR=RE(I,J)
c	    XI=IM(I,J)
c	      IF(XR.EQ.0) THEN
c	        XP=PI/2.
c	        IF(XI.LT.0.)XP=-PI/2.
c	      ELSE
c	        XP=ATAN(XI/XR)
c	        IF(XR.LT.0.)XP=XP+PI
c	      ENDIF
c	     PHA(I,J)=PHA(I,J)+XP
C
C*******************************************************************
C OUT_BISP
C To output a bispectrum
C
C Input:
C YCE1
C Output:
C PHA: phase
C MOD: modulus
C*******************************************************************
	SUBROUTINE OUT_BISP(YCE1,NGAMMA,MOD,PHA,NX,NY,IDIM)
 
C  POUR IR=25
	PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
C	PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)

	REAL MOD(IDIM,*),PHA(IDIM,*),YCE1(NGMAX,2)
	REAL PI
	INTEGER NX,NY,NGAMMA
 
C KLM(1,.) POUR K ; KLM(2,.) POUR L ; KLM(3,.) POUR M ;
	COMMON /C2/KLM(3,NGMAX)
 
	PI=3.14159265358979323846
 
C Erasing the arrays:
	DO J=1,NY
	  DO I=1,NX
	    PHA(I,J)=0.
	    MOD(I,J)=0.
	  END DO
	END DO
 
C Phase of the bispectrum:
	DO 2 NG=1,NGAMMA
	    XR=YCE1(NG,1)
	    XI=YCE1(NG,2)
	    XM=SQRT(XR*XR+XI*XI)
	      IF(XR.EQ.0) THEN
	        XP=PI/2.
	        IF(XI.LT.0.)XP=-PI/2.
	      ELSE
	        XP=ATAN(XI/XR)
	        IF(XR.LT.0.)XP=XP+PI
	      ENDIF
	  IIX=KLM(1,NG)+1
	  IIY=KLM(2,NG)+1
	  IF((IIX.LE.NX).AND.(IIY.LE.NY))THEN
	    PHA(IIX,IIY)=XP
	    MOD(IIX,IIY)=XM
	  ENDIF
C BISP(U,V)=BISP(V,U):
c	  IIX=KLM(2,NG)+1
c	  IIY=KLM(1,NG)+1
c	  PHA(IIX,IIY)=XP
c	  MOD(IIX,IIY)=XM
2	CONTINUE
 
	RETURN
	END
