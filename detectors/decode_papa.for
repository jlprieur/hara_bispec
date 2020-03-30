C++**************************************************************************
C Program DECODE_PAPA
C To decode speckle files generated by PAPA type camera
C (X,Y,X,Y,etc) in 8 bits for X and 8 bits for Y.
C Output: bispectrum phase factor, modulus (squared) of the spectrum,
C         and long integration.
C
C Output file for the details of the input parameters: DECODE_PAPA.DAT
C
C Tests: fft_nag.for or fft_jlp.for give exactly the same FFT,
C        and JLP version goes as fast as the nag version (maybe a bit
C        faster even)
C
C JLP
C Version: 14-06-92
C--**************************************************************************
C Nota: C@ points to debugging instructions
C
	PROGRAM DECODE_PAPA

C  POUR IR=25
	PARAMETER (NGMAX=187566,IDIM=256)

C  POUR IR=30
C	PARAMETER (NGMAX=388400,IDIM=256)

C MAX number of photons/frame:
        PARAMETER (MAXPH=2000,MAXPH2=4000)

	REAL IMAGE(IDIM,IDIM),LONG_INTEG(IDIM,IDIM)
C Squared modulus, and bispectral phase factor:
	REAL MODSQ(IDIM,IDIM),IMAGINARY(IDIM,IDIM),FFIELD(IDIM,IDIM)
	REAL YCE1(NGMAX,4),SNRM(IDIM,IDIM)
C !! integer*2 !!!!
	INTEGER*2 IBUF(MAXPH),IX,IY
	INTEGER*4 NX1,NY1,ISKIP
	INTEGER*4 NCARA,IFRAME,NFRAME,NPHOTONS
	LOGICAL*1 IBUFT(MAXPH2),IX1(2),IY1(2)
	CHARACTER LONGNAME*40
	CHARACTER FILE1*40,FFNAME*40,COMMENTS*80,NAME*40
	EQUIVALENCE (IBUF(1),IBUFT(1))
	EQUIVALENCE (IX1(1),IX)
	EQUIVALENCE (IY1(1),IY)
	COMMON/MFRAMES/LONG_INTEG
	COMMON/BISP1/MODSQ,SNRM,IMAGINARY,YCE1,IR,NBETA,NGAMMA
 

10	FORMAT(A)
 
	CALL JLP_BEGIN
	CALL JLP_INQUIFMT
 
	OPEN(2,STATUS='UNKNOWN',FILE='decode_papa.log',ERR=999)
	WRITE(2,58)
	WRITE(6,58)
58	FORMAT(' Program DECODE_PAPA Version 14-06-92',/,
     1	' Compression of images to 128x128',/)
        IFACT=2
	NX1=256/IFACT
	NY1=256/IFACT
 
C UV Coverage:
	WRITE(6,*) ' Radius (IR) of uv-coverage ?'
	READ(5,*) IR
	CALL COVERA(IR,NBETA,NGAMMA)
 
C Speckle data:
95	WRITE(6,*) ' INPUT FILE (SPECKLE DATA) :'
	READ(5,10) FILE1
C	OPEN(10,FILE=FILE1,
C     1	FORM='UNFORMATTED',STATUS='OLD',ERR=95)
C Imode=0 (read only)
	IMODE=0
	CALL JLP_OSDOPEN(FILE1,40,IMODE,IFID,ISTAT)
	 IF(ISTAT.NE.0)THEN
	    WRITE(6,*) ' Fatal error opening the data input file'
	    GOTO 998
	 ENDIF
 
	WRITE(6,98)
98	FORMAT(' Number or photons per frame (max=2000)',
     1	' number of frames, and number of bytes to skip at the beginning')
	READ(5,*) NPHOTONS,NFRAME,ISKIP
        NPHOTONS=MIN(2000,NPHOTONS)
	IF(ISKIP.GE.1)THEN
         NCARA=1
	 DO I=1,ISKIP
	   CALL JLP_OSDREAD(IFID,IBUFT,NCARA,ISTAT)
	 END DO
	 PRINT 86,ISKIP,ISTAT
86	 FORMAT(' Skipping ',I6,' bytes at the beginning, ISTAT=',I5)
	ENDIF
 
C Enter the flat field:
        WRITE(6,73)
73      FORMAT(/,' ********ENTER THE FLAT FIELD:*******')
        FFNAME=' '
        CALL JLP_READIMAG(FFIELD,NXFF,NYFF,IDIM,FFNAME,COMMENTS)
          IF(NXFF.NE.256.OR.NYFF.NE.256)THEN
            WRITE(2,76)
76          FORMAT(' Fatal error: Inconsistent size of the flat field')
            GOTO 998
          ENDIF
C Normalization of the flat field:
        CALL NORMALIZE(FFIELD,256,256,IDIM)

C Header for the logfile:
	WRITE(2,92) FILE1,NPHOTONS,NFRAME
92	FORMAT(' Input file:',A,/,
     1	' Number of photons per frame :',I8,/,
     1	' Number of frames (wanted):',I8,/)
 
	 XPHOT=0.
	 CALL ERASE_IMAGE(IMAGINARY,NX1,NY1,IDIM)
	 CALL ERASE_IMAGE(MODSQ,NX1,NY1,IDIM)
	 CALL ERASE_IMAGE(SNRM,NX1,NY1,IDIM)

C Read data as a continuous flow, with sets of 1024 points
	 NCARA=2*NPHOTONS
	 DO IFRAME=1,NFRAME
	   CALL JLP_OSDREAD(IFID,IBUFT,NCARA,ISTAT)
	   NBUF=NCARA/2
	   IF(ISTAT.NE.0.OR.NBUF.NE.NPHOTONS)THEN
	      IF(K.EQ.1)THEN
	      WRITE(6,*) ' Fatal error reading the data in the input file'
	      CALL JLP_OSDCLOSE(IFID,ISTAT)
	      GOTO 998
	      ELSE
C Normal exit from the loop:
	      GOTO 9999
	      ENDIF
	   ENDIF
C To know how it goes:
           IF(MOD(IFRAME,10).EQ.1) THEN
	     WRITE(2,111)IFRAME
	     WRITE(6,111)IFRAME
111	     FORMAT(' FRAME #',I8)
           ENDIF
C Swaping (since data coming from the Vax):	
C	   CALL SWAP(IBUF,NBUF)
C Decoding now, integrating the individual frames,
C mean spectrum, bispectrum :
	    DO I=1,NCARA,2
C For Sun: IX1(1) MSB, IX1(2) LSB.
	      IX=0
	      IY=0
	      IX1(2)=IBUFT(I)
	      IY1(2)=IBUFT(I+1)
	      IX=IX+1
	      IY=IY+1
C Add the new photon now:
C Flat field correction, and integration of the photons on IMAGE:
             IXX=IX/IFACT
             IYY=IY/IFACT
	     IMAGE(IXX,IYY)=IMAGE(IXX,IYY)+FFIELD(IX,IY)
	    END DO

	  CALL PAPA_SPIMAGE(IMAGE,NX1,NY1)
	 END DO
 
C End:
C9999	CLOSE(10)
9999	CALL JLP_OSDCLOSE(IFID,ISTAT)
	IF(ISTAT.NE.0)THEN
	  WRITE(6,*) ' Fatal error closing the data input file'
	ENDIF
 
C Recentre Fourier transforms:
        CALL RECENTRE(MODSQ,MODSQ,NX1,NY1,IDIM)
        CALL RECENTRE(SNRM,SNRM,NX1,NY1,IDIM)

C Diagnostic:
	NFRAME=IFRAME-1
	IXC=NX1/2+1
	IYC=NY1/2+1
	W1=MODSQ(IXC,IYC)/FLOAT(NFRAME)
	XPHOT=NPHOTONS
	WRITE(6,88) XPHOT,NFRAME,W1
	WRITE(2,88) XPHOT,NFRAME,W1
88	FORMAT(' Mean Nphotons/frame:',F14.2,'  Nframes :',I8,/,
     1	' Maximum of MODSQ : ',G12.4,2X)
	W1=MODSQ(IXC,IYC)
	CALL MEAN_FRAMES(MODSQ,NX1,NY1,IDIM,W1)
	CALL MEAN_FRAMES(SNRM,NX1,NY1,IDIM,W1)
 
C Simply divide the integrated frames by the number of frames:
	W1=REAL(NFRAME)
	CALL MEAN_FRAMES(YCE1,NGAMMA,4,NGMAX,W1)
	CALL MEAN_FRAMES(LONG_INTEG,NX1,NY1,IDIM,W1)
 
C Squared modulus:
	WRITE(6,39)
39	FORMAT(/,' ********OUTPUT OF THE MEAN SQUARED MODULUS:*******')
	   WRITE(COMMENTS,38) FILE1(1:12),NFRAME,XPHOT
38	   FORMAT('MODSQ ',A12,' Nf,Np:',X,I10,X,F9.1)
	   WRITE(6,10) COMMENTS
	   WRITE(2,10) COMMENTS
	   NAME='modsq'
	   CALL JLP_WRITEIMAG(MODSQ,NX1,NY1,IDIM,NAME,COMMENTS)
	   WRITE(2,36)
	   WRITE(6,36)
36	   FORMAT(' Mean spectrum (modulus squared) in: modsq')
 
C Bipsectrum:
	   NAME='bisp1'
	   CALL JLP_WRITEIMAG(YCE1,NGAMMA,4,NGMAX,NAME,COMMENTS)
	   WRITE(2,37)
37	   FORMAT(' Mean bispectrum in: bisp1 (list')
 
	WRITE(6,49)
49	FORMAT(/,' ****OUTPUT OF THE LONG INTEGRATION FRAME:*****')
	   WRITE(COMMENTS,48) FILE1(1:14),NFRAME,XPHOT
48	   FORMAT(' ',A14,' Nf,Np:',X,I10,X,F9.1)
	   WRITE(6,10) COMMENTS
	   WRITE(2,10) COMMENTS
	   LONGNAME='long'
	   CALL JLP_WRITEIMAG(LONG_INTEG,NX1,NY1,IDIM,LONGNAME,
     1	COMMENTS)
	   WRITE(2,47)
	   WRITE(6,47)
47	   FORMAT(' Long integ. (mean) of the frames in: long')
 
C End:
998	CLOSE(2)
	PRINT *,' Logfile in: decode_papa.log'
	CALL JLP_END
	STOP
999	PRINT *,' Fatal error opening: decode_papa.log'
	CALL JLP_END
	STOP
	END
C********************************************************************
C Subroutine PAPA_SPIMAGE
C Integrates the counts on speckle images
C********************************************************************
	SUBROUTINE PAPA_SPIMAGE(IMAGE,NX1,NY1)
C  POUR IR=25
	PARAMETER (NGMAX=187566,IDIM=256)

C  POUR IR=30
C	PARAMETER (NGMAX=388400,IDIM=256)

	REAL*4 MODSQ(IDIM,IDIM),YCE1(NGMAX,4)
	INTEGER*4 IR,NBETA,NGAMMA
	REAL*4 IMAGE(IDIM,IDIM)
	REAL*4 IMAGINARY(IDIM,IDIM)
	INTEGER*4 NX1,NY1
 
	COMMON/BISP1/MODSQ,SNRM,IMAGINARY,YCE1,IR,NBETA,NGAMMA
	
C Long integration:
	CALL SUM_FRAMES(IMAGE,NX1,NY1)
 
C Fourier Transform:
C        WRITE(6,*) ' Start FFT'
        CALL MYFOURN2(IMAGE,IMAGINARY,NX1,NY1,IDIM)

C Computes the mean spectrum and bispectrum (correction from photon noise)
C        PRINT *,' Calling BISPEC1' 
C Nota: 02-12-91: BISPEC1 seems better than BISPEC2
	CALL BISPEC1(IMAGE,IMAGINARY,MODSQ,SNRM,NX1,NY1,IDIM,
     1	YCE1,IR,NBETA,NGAMMA)
C	CALL BISPEC2(IMAGE,IMAGINARY,MODSQ,SNRM,NX1,NY1,IDIM,
C     1	YCE1,IR,NBETA,NGAMMA)
 
C Prepares the next step:
C Cleans the previous image:
	CALL ERASE_IMAGE(IMAGE,NX1,NY1,IDIM)
 
	RETURN
	END
C********************************************************************
C Subroutine ERASE
C Erases the image
C********************************************************************
	SUBROUTINE ERASE_IMAGE(IMAGE,NX,NY,IDIM)
	REAL*4 IMAGE(IDIM,*)
	INTEGER*4 NX,NY
	DO J=1,NY
	  DO I=1,NX
	     IMAGE(I,J)=0.
	  END DO
	END DO
	RETURN
	END
C********************************************************************
C Subroutine SUM_FRAMES
C To add up all the frames
C********************************************************************
	SUBROUTINE SUM_FRAMES(IMAGE,NX,NY)
	PARAMETER (IDIM=256)
	REAL*4 LONG_INTEG(IDIM,IDIM)
	REAL*4 IMAGE(IDIM,*)
	COMMON/MFRAMES/LONG_INTEG
	
C Add up all the frames:
	DO J=1,NY
	  DO I=1,NX
	    LONG_INTEG(I,J)=LONG_INTEG(I,J)+IMAGE(I,J)
	  END DO
	END DO
 
	RETURN
	END
C********************************************************************
C Subroutine MEAN_FRAMES
C Computes the mean frame
C********************************************************************
	SUBROUTINE MEAN_FRAMES(A,NX,NY,IDIM,XQUO)
	REAL*4 A(IDIM,*),XNUMB,XQUO
	INTEGER*4 NX,NY
	
	DO J=1,NY
	   DO I=1,NX
	      A(I,J)=A(I,J)/XQUO
	   END DO
	END DO
	
	RETURN
	END
 
C********************************************************************
C Swap the two bytes of integer*2 values.
C
C IN   : Input/Output array
C NCAR : Number of INTEGER*2 to swap
c************************************************************
        SUBROUTINE SWAP (IN, NCAR)
	INTEGER*2 IN(*), AUXI2
	INTEGER*4 NCAR
	LOGICAL*1 LOW, AUXL1(2)
	EQUIVALENCE (AUXI2, AUXL1(1))
 
        DO I=1, NCAR
	 AUXI2=IN(I)
	 LOW  =AUXL1(1)
C Swap
	 AUXL1(1)=AUXL1(2)
	 AUXL1(2)=LOW
	 IN(I)  =AUXI2
	 END DO
	 RETURN
 
	 END
C********************************************************************
C Subroutine NORMALIZE
C To normalize the flat field on the working area
C
C JLP91: Please note that we take 1/FFIELD in order to multiply with
C the flat field... (easier to handle with padding filter and zeroes...)
C********************************************************************
        SUBROUTINE NORMALIZE(FFIELD,NX,NY,IDIM)
        REAL*4 FFIELD(IDIM,*)
        REAL*8 XMEAN
        REAL*4 MIN_VALUE,MAX_VALUE
        CHARACTER NAME*40,COMMENTS*80

C Computing the mean on the good area:
        XMEAN=0.
        DO J=1,NY
          DO I=1,NX
             FFIELD(I,J)=MAX(0.,FFIELD(I,J))
             XMEAN=XMEAN+FFIELD(I,J)
          END DO
        END DO
        XMEAN=XMEAN/FLOAT(NX*NY)
        WRITE(6,58) XMEAN
        WRITE(2,58) XMEAN
58      FORMAT(' First mean value of the Flat Field: ',G12.5)

C Then normalization to 1:
C This value should not be too small, otherwise it distorts the spectrum!!!
C (0.01 is too small, 0.333 not too bad,...)
        MIN_VALUE=0.1
        MAX_VALUE=10.
        SUM=0.
        DO J=1,NY
          DO I=1,NX
             FFIELD(I,J)=FFIELD(I,J)/XMEAN
             FFIELD(I,J)=MAX(FFIELD(I,J),MIN_VALUE)
             FFIELD(I,J)=MIN(FFIELD(I,J),MAX_VALUE)
             SUM=SUM+FFIELD(I,J)
          END DO
        END DO
        XMEAN=SUM/FLOAT(NX*NY)
        WRITE(6,59) XMEAN
        WRITE(2,59) XMEAN
59      FORMAT(' Second mean of the Flat Field (limits 0.1, 10.): ',G12.5)

C Normalization to 1 and inversion of the flat field:
        DO J=1,NY
          DO I=1,NX
            FFIELD(I,J)=XMEAN/FFIELD(I,J)
          END DO
        END DO

C Test:
        COMMENTS='Inverse flat field'
        NAME='ffield_test'
        PRINT *,' Output of the synthetic Flat Field: ',NAME
        CALL JLP_WRITEIMAG(FFIELD,NX,NY,IDIM,NAME,
     1       COMMENTS)

        RETURN
        END
C*************************************************************
	INCLUDE 'jlp_bispec1.for'
	include 'jlp_cover.for'
	include 'fft_jlp.for'
C	include 'fft_nag.for'
C********************************************************************
        SUBROUTINE MYFOURN2(TR,TI,NX,NY,IDIM)
	REAL TR(IDIM,*),TI(IDIM,*)
        REAL DATA(600*600*2)
        INTEGER NN(2)

C Rearranging the data for MYFOURN1
C (real and complex value ordered in a standard fortran array...)
        K=1
        DO J=1,NY
          DO I=1,NX
          DATA(K)=TR(I,J)
          DATA(K+1)=0.
          K = K+2
          END DO
        END DO

C Actual FFT:
        NN(1)=NX
        NN(2)=NY
        CALL MYFOURN1(DATA,NN)

C Rearranging the data
        K=1
        DO J=1,NY
          DO I=1,NX
          TR(I,J)=DATA(K)
          TI(I,J)=DATA(K+1)
          K = K+2
          END DO
        END DO

        RETURN
        END
C++*******************************************************************
C FOURN1
C From "Numerical Recipees" p 451
C
C Replaces DATA by its NDIM-dimensional disccrete Fourier transform.
C NN is an integer array of length NDIM, containing the lengths of each
C dimension (number of complex values), which MUST all be powers of 2.
C DATA is a real array of length twice the product of these lengths, in
C which the data are stored as in a multidimensional complex Fortran array.
C
C Restriction of fourn1 to ndim=2, and only direct FFT
C
C JLP
C Version 03-11-91
C--*********************************************************************
	SUBROUTINE MYFOURN1(DATA,NN)
 
C Double precision for trigonometric recurrences:
	real*8 wr,wi,wpr,wpi,wtemp,theta
 
	integer*4 ntot,nprev,n,i1,i2,i3,i2rev,i3rev
	integer*4 ibit,ifp1,ifp2,ip1,ip2,k1,k2,idim,nrem
	integer*4 nn(*)
	real*4 data(*),tempr,tempi
 
C Compute total number of complex values:
	ntot=1
	do 11 idim=1,2
	  ntot=ntot*nn(idim)
11	continue
 
	nprev=1
 
C Main loop over the dimensions:
	do 18 idim=1,2
	  n=nn(idim)
	  nrem=ntot/(n*nprev)
	  ip1=2*nprev
	  ip2=ip1*n
	  ip3=ip2*nrem
	  i2rev=1
 
C This is the bit reversal routine
C (to set up the book keeping for the next step):
	  do 14 i2=1,ip2,ip1
	      if(i2.lt.i2rev)then
	         do 13 i1=i2,i2+ip1-2,2
	            do 12 i3=i1,ip3,ip2
	               i3rev=i2rev+i3-i2
C Swap real array  (i3 <-> i3rev) and imaginary array (i3+1 <-> i3rev+1)
	               tempr=data(i3)
	               tempi=data(i3+1)
	               data(i3)=data(i3rev)
	               data(i3+1)=data(i3rev+1)
	               data(i3rev)=tempr
	               data(i3rev+1)=tempi
12	            continue
13	         continue
              endif
              ibit=ip2/2
1	      if((ibit.ge.ip1).and.(i2rev.gt.ibit))then
	          i2rev=i2rev-ibit
	          ibit=ibit/2
	      go to 1
	      endif
	      i2rev=i2rev+ibit
14	   continue
 
C Here begins the Danielson-Lanczos section of the routine:
C (i.e. iterative processing of the data)
	ifp1=ip1
2	if(ifp1.lt.ip2)then
	  ifp2=2*ifp1
C Initialize for the trig. recurrence:
	  theta=-6.28318530717959/(ifp2/ip1)
	  wpr=-2.d0*dsin(0.5d0*theta)**2
	  wpi=dsin(theta)
	  wr=1.d0
	  wi=0.d0
	  do 17 i3=1,ifp1,ip1
	    do 16 i1=i3,i3+ip1-2,2
	      do 15 i2=i1,ip3,ifp2
	        k1=i2
	        k2=k1+ifp1
C Danielson-Lanczos formula:
	        tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
	        tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
	        data(k2)=data(k1)-tempr
	        data(k2+1)=data(k1+1)-tempi
	        data(k1)=data(k1)+tempr
	        data(k1+1)=data(k1+1)+tempi
15	      continue
16	    continue
 
C Trigonometric recurrence:
	    wtemp=wr
	    wr=wr*wpr-wi*wpi+wr
	    wi=wi*wpr+wtemp*wpi+wi
17	    continue
	  ifp1=ifp2
	  go to 2
	  endif
	  nprev=n*nprev
18	continue

	return
	end