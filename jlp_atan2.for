C+*****************************************
C Set of routines to make C and Fortran programs run identically
C
C JLP
C Version 14-09-95
C--**********************************************
C ALPHA in radians
C**********************************************
       SUBROUTINE JLP_COS(COS1,ALPHA)
       REAL COS1, ALPHA
       COS1=COS(ALPHA)
       RETURN
       END
C**********************************************
       SUBROUTINE JLP_SIN(SIN1,ALPHA)
       REAL SIN1, ALPHA
       SIN1=SIN(ALPHA)
       RETURN
       END
C**********************************************
       SUBROUTINE JLP_ATAN2(ALPHA,IM1,RE1)
       REAL ALPHA, RE1, IM1
       ALPHA = ATAN2(IM1,RE1) 
       RETURN
       END
