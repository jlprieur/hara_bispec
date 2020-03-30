C*****************************************************************
C Main program to plot UV coverage, and numbering
C
C JLP 
C Version of 20-05-92
C---------------------------------------------------------------------
      PROGRAM PLOT_UV
 
C  POUR IR=25
      PARAMETER(IRMAX=25,MIRMAX=-IRMAX)
      PARAMETER(NBMAX=980,NGMAX=187566,NULL=0)

C  POUR IR=30
C      PARAMETER(IRMAX=30,MIRMAX=-IRMAX)
C      PARAMETER(NBMAX=1410,NGMAX=388400,NULL=0)
 
      PARAMETER (IDIM=1500,KCUR=2)
      REAL*4 X(IDIM,KCUR),Y(IDIM,KCUR)
      INTEGER*4 NBCOUV,IXY,NPTS(KCUR),IOPTION
      CHARACTER CHAR1*20,CHAR2*20,TITLE*40,DATE*24
      CHARACTER PLOTDEV*32,FILENAME*40,ANS*1
      CHARACTER NCHAR(KCUR)*4,LABEL(IDIM,KCUR)*4
 
C NBCOUV( X DE -IRMAX A IRMAX,  Y DE 0 A IRMAX )
C IXY( 1 POUR X ; 2 POUR Y,  NB DE 0 A NBMAX )
      COMMON /C1/NBCOUV(MIRMAX:IRMAX,NULL:IRMAX),IXY(2,NULL:NBMAX)
 
10      FORMAT(A)
 
      WRITE(6,56)
56    FORMAT(' Program PLOT_UV  / Version 03-06-92',
     1       ' Menu:',/,
     1       ' 1= uv coverage of a full disc-like aperture',/,
     1       ' 2= uv coverage contained in a file',/,
     1       ' 3=Network with telescope coordinates in a file',/
     1       ' Enter your choice: ')
      READ(5,*) IOPTION

C******************************************************************
C OPTION 1: u-v coverage of a full disc-like aperture
      IF(IOPTION.EQ.1) THEN
        PRINT *,' Radius (IR) of uv-coverage ?'
        READ(5,*) IR
        CALL COVERA(IR,NBETA,NGAMMA)
 
        WRITE(6,46) 
46      FORMAT(' Maximum NBETA to troncate the spectral list :')
        READ(5,*) I
        IF(I.GT.0.AND.I.LT.NBETA) NBETA=I

C Last value is 0:
        I=0
        X(NBETA+1,1)=0.
        Y(NBETA+1,1)=0.
        WRITE(LABEL(NBETA+1,1),8) 0
8        FORMAT(I3)
 
      DO I=1,NBETA
        X(I,1)=IXY(1,I)
        Y(I,1)=IXY(2,I)
        WRITE(LABEL(I,1),8) I
C Write also the symmetric relative to 0 (curve #2):
        X(I,2)=-X(I,1)
        Y(I,2)=-Y(I,1)
        WRITE(LABEL(I,2),8) -I
      END DO
 
      NPTS(1)=NBETA+1
      NPTS(2)=NBETA
 
C******************************************************************
C OPTION=2: uv-coverage contained in a file: 
      ELSEIF(IOPTION.EQ.2) THEN
        WRITE(6,45)
45      FORMAT(' Input data file ( u-coord, v-coord, frequency_number):')
        CALL RREADFILE(X(1,1),Y(1,1),NPTS(1),X(1,2),Y(1,2),
     1  NPTS(2),1)

        NBETA=NPTS(1)

        DO I=1,NBETA
          WRITE(LABEL(I,1),8) INT(Y(I,2)) 
          PRINT*,' X,Y,LABEL',X(I,1),Y(I,1),LABEL(I,1)
        END DO

        DO I=1,NBETA
C First write label before writing over Y(I,2) !!!
          WRITE(LABEL(I,2),8) -INT(Y(I,2)) 
C Write also the symmetric relative to 0 (curve #2):
          X(I,2)=-X(I,1)
          Y(I,2)=-Y(I,1)
      END DO
 
C Last value is 0:
        I=0
        X(NBETA+1,1)=0.
        Y(NBETA+1,1)=0.
        WRITE(LABEL(NBETA+1,1),8) 0

        NPTS(1)=NBETA+1
        NPTS(2)=NBETA
 
C******************************************************************
C OPTION=3: Telescope network:
C******************************************************************
      ELSE
        WRITE(6,47)
47      FORMAT(' Input data file ( Xcoord, Ycoord, Telescope_number):')
        CALL RREADFILE(X(1,1),Y(1,1),NPTS(1),X(1,2),Y(1,2),
     1  NPTS(2),1)

        DO I=1,NPTS(1)
          WRITE(LABEL(I,1),8) INT(Y(I,2)) 
          PRINT*,' X,Y,LABEL',X(I,1),Y(I,1),LABEL(I,1)
        END DO

      ENDIF
 
C******************************************************************
C Input of symbol codes:
      IF(IOPTION.LE.2)THEN
        PRINT *,'Symbols for positive and negative Y? ("93,83" or "92,82")'
        READ(5,*) I1,I2
        WRITE(NCHAR(1),'(I2)') I1
        WRITE(NCHAR(2),'(I2)') I2
        KCURVE=2
      ELSE
        PRINT *,'Symbol? ("93" or "83")'
        READ(5,*) I1
        WRITE(NCHAR(1),'(I2)') I1
        KCURVE=1
      ENDIF

C******************************************************************
88      WRITE(6,*) ' Output graphic device (&xterm, &square?)'
      READ(5,10) PLOTDEV
 
      PRINT *,'Title  ?'
      READ(5,10) TITLE
      CHAR1=' '
      CHAR2=' '
 
      CALL PLOT_UV1(X,Y,LABEL,NPTS,IDIM,KCURVE,
     1      CHAR1,CHAR2,TITLE,NCHAR,PLOTDEV)
 
      PRINT *,' Do you want a copy of this graph ? (Y)'
      READ(5,10) ANS
      IF(ANS.NE.'n'.AND.ANS.NE.'N')GOTO 88
 
      END
C*****************************************************************
C Subroutine PLOT_UV1 
C
C Possibility of drawing labels for each point (LABEL(.,.))
C
C---------------------------------------------------------------------
      SUBROUTINE PLOT_UV1(X,Y,LABEL,NPTS,NMAX,KCURVE,
     1      CHAR1,CHAR2,TITLE,NCHAR,PLOTDEV)
      REAL*4 X(NMAX,*),Y(NMAX,*)
      REAL*4 EXPAND,X1,X2,Y1,Y2
      INTEGER*4 NPTS(*),ISYMB,ISIZE,ICODE,MAX_LENGTH
      CHARACTER CHAR1*20,CHAR2*20,TITLE*40,DATE*24
      CHARACTER PLOTDEV*32,PLOTDEV1*32
      CHARACTER NCHAR(*)*4,WORD*4,LABEL(NMAX,*)*4
      CHARACTER BUFFER*80
      LOGICAL CONNECTED,FILE,COPY
 
C Comon block with "NEWPLOT" and "WINDOW_LIMITS"
      COMMON/PARAMETERS/OFFX,OFFY,AXLEN,AYLEN,XMIN,YMIN,
     1      XMAX,YMAX,TDX,TDY
 
C******************************************************
C Compute the parameters of the frame:      
 
C Compute the minimum and maximum :
      CALL NEWSCALE(X,NMAX,KCURVE,NPTS,XMIN,XMAX)
      CALL NEWSCALE(Y,NMAX,KCURVE,NPTS,YMIN,YMAX)
      
      WRITE(6,45) XMIN,XMAX,YMIN,YMAX 
45    FORMAT(' XMIN,XMAX,YMIN,YMAX :', 4(G10.4,1X),/,
     1   ' Enter new values if wanted ')
      READ(5,10) BUFFER
10    FORMAT(A)
      X1=0
      X2=0
      READ(BUFFER,*,ERR=98,END=98) X1,X2,Y1,Y2
      IF(X1.NE.X2)THEN
       XMIN=X1
       XMAX=X2
       YMIN=Y1
       YMAX=Y2
      ENDIF
98    CONTINUE 
C Get the window parameters for the graph:  (reads the file GRAPHIC.KER)
C Get AXLEN,AYLEN,OFFX,OFFY:
      COPY=.FALSE.
      CALL WINDOW_LIMITS(PLOTDEV,ICODE,FILE,COPY)
 
C Set same scale in X and Y:
        IPLAN=1
        CALL JLP_SETUP_PLOT(OFFX,OFFY,
     1    AXLEN,AYLEN,XMIN,XMAX,YMIN,YMAX,IPLAN)


C Select the device
        TEX_FLAG=0
        I=200
        J=300
        PLOTDEV1=PLOTDEV(2:)
        CALL JLP_SPDEVICE(PLOTDEV1,I,J,TEX_FLAG)

C********************************************************
C Decode the plotting device: (terminal or laser printer)
        PLOTDEV1=PLOTDEV(2:32)
 
C********************************************************
C***************************************************
C Drawing the frame ....
        CALL JLP_SPBOX(CHAR1,CHAR2,TITLE)

C********************************************************
C Drawing a symbol at the location of each point :
        EXPAND=1.4
        MAX_LENGTH=3
        ANGLE=0.
        IDRAW=1
 
      DO K=1,KCURVE
C First decode the symbol :
          READ(NCHAR(K),500) ISYMB,ISIZE
500       FORMAT(I1,I1)
 
C Computing shift in window coordinates between symbol and label:
C (Window coordinate unit is 200*ISIZE for the symbol)
        IF(FILE)THEN
C Good for Postscript:
          OFX=-1000*ISIZE
          OFY=200*ISIZE
        ELSE
C Good for xterm:
          OFX=-300*ISIZE
          OFY=200*ISIZE
        ENDIF
        DO I=1,NPTS(K)
C Drawing the chosen symbol at each point :
              IX=(X(I,K)-XMIN)*TDX+OFFX
              IY=(Y(I,K)-YMIN)*TDY+OFFY
              CALL JLP_SYMBOL(IX,IY,ISIZE,ISYMB)
C Drawing the label: 
              WORD=LABEL(I,K)
              IX=IX+OFX
              IY=IY+OFY
              CALL JLP_SPLABEL(WORD,MAX_LENGTH,IX,IY,ANGLE,EXPAND,
     1        IDRAW,LENGTH)
        END DO
      END DO
 
C********************************************************
C Send graph to screen: 
        CALL JLP_GFLUSH
 
C Write the date if laser print :
      IF(FILE)THEN
        CALL JLP_DATE_TIME(DATE)
        EXPAND=0.35
C        CALL MGOSETEXPAND(EXPAND)
C        CALL MGORELOCATE(0.8,1.2)
C        CALL MGOLABEL(24,DATE)
      ELSE
         CALL JLP_WHERE(X,Y,IN_FRAME)
      ENDIF

 
C Closes plot:
        CALL JLP_SPCLOSE
 
      RETURN
      END
C********************************************************************
      INCLUDE 'jlp_cover.for'
 
