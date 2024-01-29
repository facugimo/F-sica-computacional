      PROGRAM INTEGRA
        IMPLICIT NONE

        REAL*8 X1,X2,A,A1,A2,PI,B,H
        INTEGER IOPT,M

        COMMON/OPTION/IOPT

C       INICIALITZEM LES VARIABES QUE NECESITAREM

        A=57.91D6
        B=57.67D6
        PI=4.D0*ATAN(1.D0)

        X1=2*A
        X2=4*A
        IOPT=1

        OPEN (1,FILE="resultsP3.dat")

        DO M=1,14
           H=(X2-X1)/(2**M)
           CALL MYINTEGRATOR(X1,X2,M,1,A1)
           CALL MYINTEGRATOR(X1,X2,M,2,A2)
           WRITE (1,10) H,A1,A2
           ENDDO

        CLOSE(1)

        END


      SUBROUTINE MYFUNCI(X,F,IOPT)
        REAL*8 X,F,B,A
        INTEGER IOPT

        PARAMETER(B=56.67D6)
        PARAMETER(A=57.91D6)
C       DEFINIM LA FUNCIO QUE INTEGRAREM PER TROBAR L'ÀREA
        IF (IOPT.EQ.1) THEN
           F=2*B*SQRT(1-((X-3*A)**2)/(A**2))
C       DEFINIM LA FUNCIO QUE INTEGRAREM PER TROBAR L'ÀREA AMB EL CANVI DE VARIABLE
        ELSE IF (IOPT.EQ.2) THEN
           F=2*B*A*COS(X)**2
           ENDIF
C	EVITEM ELS INFINITS ASIGNANT EL VALOR 0
        IF (1/F.EQ.0) THEN
           F=0
           ENDIF

        END SUBROUTINE


      SUBROUTINE MYINTEGRATOR(X1,X2,M,IM,VAL)
        REAL*8 X1,X2,VAL,F
        INTEGER M,IM,IOPT

        COMMON/OPTION/IOPT

        VAL=0
        H=(X2-X1)/(2.D0**M)

C       METODE DEL TRAPEZI COMPOSTA
        IF (IM.EQ.1) THEN
           CALL MYFUNCI(X1,F,IOPT)
           VAL=VAL+F*H/2.D0

           CALL MYFUNCI(X2,F,IOPT)
           VAL=VAL+F*H/2.D0

           DO I=1,(2**M)-1
              CALL MYFUNCI(X1+H*I,F,IOPT)
              VAL=VAL+H*F
              ENDDO

C       METODE DE SIMPSON
        ELSE IF (IM.EQ.2) THEN
           CALL MYFUNCI(X1,F,IOPT)
           VAL=VAL+F*H/3.D0

           CALL MYFUNCI(X2,F,IOPT)
           VAL=VAL+F*H/3.D0

           DO I=1,(2**M)-1
              CALL MYFUNCI(X1+I*H,F,IOPT)
              IF (MOD(I,2).EQ.0) THEN
                 VAL=VAL+F*2*H/3.D0
              ELSE
                 VAL=VAL+F*4*H/3.D0
                 ENDIF
              ENDDO
        ENDIF
        END SUBROUTINE

