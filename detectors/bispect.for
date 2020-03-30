C****************************************************************
C BISPECT.FOR
C Program to compute the bispectrum of a single frame
C and test the routines used in DECODE_MAMA
C JLP
C Version: 10-06-92
C****************************************************************
	PROGRAM BISPECT
 
C  POUR IR=25
	PARAMETER(IDIM=256,NGMAX=187566)

C  POUR IR=30
C	PARAMETER(IDIM=256,NGMAX=388400)

	REAL IMAGE(IDIM,IDIM),IM(IDIM,IDIM),SNRM(IDIM,IDIM)
	REAL YCE1(NGMAX,4),PSF(IDIM,IDIM),MODSQ(IDIM,IDIM)
	INTEGER IR,NBETA,NGAMMA
	CHARACTER NAME1*40,NAME2*40,NAME*40,COMMENTS*80
C
	CALL JLP_BEGIN
 
	OPEN(2,FILE='bispect.log',STATUS='UNKNOWN',
     1	ACCESS='SEQUENTIAL')
 
C UV coverage:
	WRITE(6,*) ' Radius of uv-coverage (IR) ?'
	READ(5,*) IR
	CALL COVERA(IR,NBETA,NGAMMA)
 
C FORMAT FICHIER
	CALL JLP_INQUIFMT
 
C Input image:
	WRITE(6,37) IDIM,IDIM
37	FORMAT(' Input object : (max size :',I4,'x',I4,')')
	NAME1=' '
	CALL JLP_READIMAG(IMAGE,NX,NY,IDIM,NAME1,COMMENTS)
 
C Simulates the effects of the telescope on the modulus of the FFT:
C@	WRITE(6,38)
C@ 38	FORMAT(' Input psf (modulus squared): (same size as object)')
C@	NAME2=' '
C2	CALL JLP_READIMAG(PSF,NX,NY,IDIM,NAME2,COMMENTS)
 
C@	CALL SIMU_PSF(IMAGE,MIM,PSF,NX,NY,IDIM)
 
 
C Null imaginary part:
	CALL ERASE(IM,NX,NY,IDIM)

C Spectrum by FFT:
        CALL FFT_2D(IMAGE,IM,NX,NY,IDIM,1)

C Preparing the next step:	
	CALL ERASE(MODSQ,NX,NY,IDIM)
	CALL ERASE(SNRM,NX,NY,IDIM)
	CALL ERASE(YCE1,NGAMMA,4,NGMAX)

C Bispectrum:
C Processing this image now:  bispec1 is with photon noise correction
C bispec2 is without
C	CALL BISPEC1(IMAGE,IM,MODSQ,SNRM,NX,NY,IDIM,
C     1	YCE1,IR,NBETA,NGAMMA)
	IFRAME=1
	CALL BISPEC2(IMAGE,IM,MODSQ,SNRM,NX,NY,IDIM,
     1	YCE1,IR,NBETA,NGAMMA,IFRAME)
 
C Output:
	WRITE(6,42)
	WRITE(2,42)
42	FORMAT(' Squared modulus: output in b_modsq ',/,
     1	' Real part:      b_real',/,
     1	' Imaginary part: b_imag')
 
	CALL RECENTRE(MODSQ,MODSQ,NX,NY,IDIM)
	CALL RECENTRE(IMAGE,IMAGE,NX,NY,IDIM)
	CALL RECENTRE(IM,IM,NX,NY,IDIM)

	WRITE(COMMENTS,43) NAME1(1:20)
43	FORMAT(' MODSQ of ',A20)
	NAME='b_modsq'
	CALL JLP_WRITEIMAG(MODSQ,NX,NY,IDIM,NAME,COMMENTS)
 
	WRITE(COMMENTS,44) NAME1(1:20)
C@	WRITE(COMMENTS,44) NAME1(1:20),NAME2(1:14)
44	FORMAT(' Real part of ',A20)
C244	FORMAT(' Real part of ',A20,' PSF: ',A14)
	NAME='b_real'
	CALL JLP_WRITEIMAG(IMAGE,NX,NY,IDIM,NAME,COMMENTS)
 
	WRITE(COMMENTS,45) NAME1(1:20)
C@45	FORMAT(' Imag. part of ',A20,' PSF: ',A14)
45	FORMAT(' Imag. part of ',A20)
	NAME='b_imag'
	CALL JLP_WRITEIMAG(IM,NX,NY,IDIM,NAME,COMMENTS)
 
C	WRITE(COMMENTS,46) NAME1(1:20)
C46	FORMAT(' Transformation of: ',A20)
C@46	FORMAT(' Transformation of: ',A20,' PSF: ',A14)
C	NAME='b_object'
C	CALL JLP_WRITEIMAG(IMAGE,NX,NY,IDIM,NAME,COMMENTS)
 
C Bispectrum : YCE1
	DO I=1,NGAMMA
	  YCE1(I,3)=1.
	  YCE1(I,4)=1.
	END DO
	WRITE(2,41)
	WRITE(6,41)
41	FORMAT(' BISPECTRUM (List) in b_bisp ')
	WRITE(2,51) NAME1(1:20)
	WRITE(COMMENTS,51) NAME1(1:20)
C@51	FORMAT(' Bispectrum of ',A20,' PSF: ',A14)
51	FORMAT(' Bispectrum of ',A20)
	NAME='b_bisp'
	CALL JLP_WRITEIMAG(YCE1,NGAMMA,4,NGMAX,NAME,COMMENTS)
 
	WRITE(6,*) ' Log File in "bispect.log"'
	CLOSE(2)
	CALL JLP_END
	STOP
	END
C******************************************************************
	SUBROUTINE ERASE(ARRAY,NX,NY,IDIM)
	REAL ARRAY(IDIM,*)
	INTEGER NX,NY
 
	DO J=1,NY
	  DO I=1,NX
	    ARRAY(I,J)=0.
	  END DO
	END DO
 
	RETURN
	END
C******************************************************************
C Subroutine SIMU_PSF
C To take into account the effects of the PSF of the telescope
C******************************************************************
	SUBROUTINE SIMU_PSF(IMAGE,MIM,PSF,NX,NY,IDIM)
	REAL IMAGE(IDIM,*),MIM(IDIM,*),PSF(IDIM,*)
	INTEGER NX,NY
 
C FFT of the image:
	CALL FFT_2D(IMAGE,MIM,NX,NY,IDIM,1)
 
C Effect of PSF (PSF normalized to 1. in the center):
	IXC=NX/2+1
	IYC=NY/2+1
	XPSF=SQRT(PSF(IXC,IYC))
	DO J=1,NY
	  DO I=1,NX
	      RR=MAX(0.,PSF(I,J))
	      RR=SQRT(RR)/XPSF
	      IMAGE(I,J)=IMAGE(I,J)*RR
	      MIM(I,J)=MIM(I,J)*RR
	  END DO
	END DO
 
C Inverse FFT:
	CALL FFT_2D(IMAGE,MIM,NX,NY,IDIM,-1)
 
	RETURN
	END
C******************************************************************
	include 'jlp_cover.for'
	include 'fft_jlp.for'
	include 'jlp_bispec1.for'
