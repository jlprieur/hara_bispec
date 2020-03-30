C++****************************************************************
C Program SNR_BISPEC
C To compute the SNR of the bispectrum (phase factor)
C (not finished yet!)
C Using formula from Dainty and Greenaway (1979, JOSA, 69,5,786)
C
C JLP
C Version: 30-03-90
C--****************************************************************
	PROGRAM SNR_BISPEC
	PARAMETER(IDIM=256)
	CHARACTER NAME*40,COMMENTS*80
	REAL MODSQ(IDIM,IDIM),SNRM(IDIM,IDIM)
 
C Format of the files:
	PRINT *,' WARNING: all the files should not have odd',
     1	' numbers for NX and NY (size in X and Y)'
	CALL JLP_INQUIFMT
 
C Input of the bispectrum:
	WRITE(6,37)
37	FORMAT(' Bispectrum list :')
	NAME=' '
	CALL JLP_READIMAG(BISP1,NX,NY,IDIM,NAME,COMMENTS)
 
	CALL SNR_BISP(BISP1,NX,NY,IDIM,SNRM)
 
	WRITE(COMMENTS,11) NAME(1:10)
11	FORMAT(' SNR of the spectrum:',A10)
	NAME=' '
	CALL JLP_WRITEIMAG(SNRM,NX,NY,IDIM,NAME,COMMENTS)
 
	CALL JLP_END
	END
 
C*******************************************************************
C SNR_MODSQ
C Computes the SNR of the modulus of the spectrum (squared)
C Uses the formula from Dainty and Greenaway (1979, JOSA, 69,5,786)
C*******************************************************************
	SUBROUTINE SNR_MODSQ(MODSQ,NX,NY,IDIM,SNRM)
	INTEGER NX,NY,IDIM,IXC,IYC,W1
	REAL MODSQ(IDIM,*),SNRM(IDIM,*)
 
C Center:
	IXC=(NX/2)+1
	IYC=(NY/2)+1
C Number of photons:
	IF(SNRM(IXC,IYC).LE.0)THEN
	  WRITE(6,93)
93	  FORMAT(' Fatal error: the mean number of photons',
     1	'per frame is null!')
	  CALL EXIT
	ENDIF
	PHOT=SQRT(SNRM(IXC,IYC))
 
C Main loop:
C Formula from Dainty and Greenaway (1979, JOSA, 69,5,786)
	DO 91 J=0,NY/2
	  J1=IYC+J
	  J2=IYC+2*J
	  DO 92 I=0,NX/2
	    I1=IXC+I
	    I2=IXC+2*I
	    W1=MODSQ(I1,J1)/PHOT
	    IF((I2.GT.NX).OR.(J2.GT.NY))THEN
	      W2=W1/(1.+W1)
	    ELSE
	      W2=W1/SQRT((1.+W1)*(1.+W1)+MODSQ(I2,J2)/PHOT**2)
	    ENDIF
	    SNRM(IXC+I,IYC+J)=W2
	    SNRM(IXC+I,IYC-J)=W2
	    SNRM(IXC-I,IYC+J)=W2
	    SNRM(IXC-I,IYC-J)=W2
92	  CONTINUE
91	CONTINUE
	
	RETURN
	END
