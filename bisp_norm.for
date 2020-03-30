C++***************************************************************
C Program to normalize the bispectrum
C
C JLP
C Version of 23-10-92
C--***************************************************************
	PROGRAM NORM_BISPEC 
	INTEGER*4 IN1,OUT2,MADRID(1)
	CHARACTER NAME*40,COMMENTS*80
	COMMON /VMR/MADRID
 
C To get the possibility of command line
	CALL JLP_BEGIN
 
C Inquire the format (input/output) :
	CALL JLP_INQUIFMT
 
        PRINT 81
81	FORMAT(' Program NORM_BISPEC to normalize a bispectrum ')

C Input :
	NAME=' '
	CALL JLP_VM_READIMAG(IN1,NX1,NY1,NAME,COMMENTS)
 
C Get memory space:
	ISIZE=NX1*NY1*4
	CALL JLP_GETVM(OUT2,ISIZE)

C Processing
	CALL NORMALIZE(MADRID(IN1),MADRID(OUT2),NX1,NY1)

C Output :
	WRITE(COMMENTS,84) NAME(1:40)
84	FORMAT(' Normalized version of ',A40)
	NAME=' '
	  CALL JLP_WRITEIMAG(MADRID(OUT2),NX1,NY1,NX1,
     1	NAME,COMMENTS)
 
99	CALL JLP_END
	STOP
	END
C*********************************************************************
	SUBROUTINE NORMALIZE(IN1,OUT2,NX1,NY1)
	INTEGER NX1,NY1,NNULL
	INTEGER I,J
	REAL*4 IN1(NX1,*),OUT2(NX1,*),XNORM

        NNULL=0
	DO I=1,NX1
            XNORM=IN1(I,1)*IN1(I,1)+IN1(I,2)*IN1(I,2)
            IF(XNORM.NE.0)THEN
                 XNORM=SQRT(XNORM)
                 OUT2(I,1)=IN1(I,1)/XNORM
                 OUT2(I,2)=IN1(I,2)/XNORM
            ELSE
              NNULL=NNULL+1
            ENDIF
        ENDDO

C Copy third line (SNR) if present
        IF(NY1.GE.3)THEN
           OUT2(I,3)=IN1(I,3)
        ENDIF

        IF(NNULL.GT.0)THEN
          PRINT *,'BISP_NORMALIZE/Warning: ',NNULL,' null complex terms'
        ENDIF

	RETURN
	END
