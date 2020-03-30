C++***************************************************************
C Program to divide two bispectra 
C
C JLP
C Version of 23-10-92
C--***************************************************************
	PROGRAM BISP_DIVIDE
	INTEGER*4 IN1,IN2,OUT,MADRID(1)
	CHARACTER NAME1*40,NAME2*40,NAME_OUT*40,COMMENTS*80
	COMMON /VMR/MADRID
 
C To get the possibility of command line
	CALL JLP_BEGIN
 
C Inquire the format (input/output) :
	CALL JLP_INQUIFMT
 
        PRINT 81
81	FORMAT(' Program BISP_DIVIDE to divide one bispectrum by another')

C Input :
        PRINT *,' First bispectrum:'
	NAME1=' '
	CALL JLP_VM_READIMAG(IN1,NX1,NY1,NAME1,COMMENTS)
C
        PRINT *,' Second bispectrum:'
	NAME2=' '
	CALL JLP_VM_READIMAG(IN2,NX2,NY2,NAME2,COMMENTS)
 
C Error handling:
        IF(NX1.NE.NX2)THEN
          PRINT *,' BISP_DIVIDE/Fatal error: incompatible sizes'
          CALL JLP_END
          STOP
        ENDIF

C Get memory space:
	ISIZE=NX1*NY1*4
	CALL JLP_GETVM(OUT,ISIZE)

C Processing
	CALL DIVIDE(MADRID(IN1),MADRID(IN2),MADRID(OUT),NX1,NY1)

C Output :
	WRITE(COMMENTS,84) NAME1(1:20),NAME2(1:20)
84	FORMAT(A20,' divided by ',A20)
	NAME_OUT=' '
	  CALL JLP_WRITEIMAG(MADRID(OUT),NX1,NY1,NX1,
     1	NAME_OUT,COMMENTS)
 
99	CALL JLP_END
	STOP
	END
C*********************************************************************
C To divide two complex:
C Real = Real1 * Real2 + Ima1 * Ima2
C Ima  = Ima1  * Real2 - Real1 * Ima2 
C********************************************************************
	SUBROUTINE DIVIDE(IN1,IN2,OUT,NX1,NY1)
	INTEGER NX1,NY1
	INTEGER I,J
	REAL*4 IN1(NX1,*),IN2(NX1,*),OUT(NX1,*)

	DO I=1,NX1
          OUT(I,1)=IN1(I,1)*IN2(I,1)+IN1(I,2)*IN2(I,2)
          OUT(I,1)=IN1(I,2)*IN2(I,1)-IN1(I,1)*IN2(I,2)
        ENDDO

C Copy third line (SNR) if present
        IF(NY1.GE.3)THEN
           OUT(I,3)=IN1(I,3)
        ENDIF

	RETURN
	END
