      PROGRAM ZEROS

        IMPLICIT NONE

        DOUBLE PRECISION POLI,P,V,PX1,PX2,PUNTS,X,XX,T,PP,DERIVADA,TT(5)
     *  ,ISO(100,5),ISO1(100),ISO2(100),ISO3(100),ISO4(100),ISO5(100)
        INTEGER N,I,MI,J

        DIMENSION P(4),PUNTS(6),X(30),XX(10,30)

        COMMON/POLIS/P,N

C     INICIALITZEM TOTES LES VARIABLES QUE HAUREM DE FER SERVIR
        N=3
        P(1)=-1.D0
        P(2)=6.D0
        P(3)=-9.D0

        PUNTS(1)=0.45D0
        PUNTS(2)=0.51D0
        PUNTS(3)=0.9D0
        PUNTS(4)=1.D0
        PUNTS(5)=1.01D0
        PUNTS(6)=1.34D0


C     CREEM L'ARXIU QUE ENS SERVIRÀ PER FER LA FIGURA "poliP4.dat"
        OPEN(1,FILE="poliP4.dat")
        DO I=130,600
          V=I/300.D0
          T=0.93
          P(4)=4.D0*T
          PX1=POLI(N,P,V)

          T=0.98
          P(4)=4*T
          PX2=POLI(N,P,V)
          WRITE (1,*) V,PX1,PX2
          ENDDO
        CLOSE(1)

C     BUSQUEM PEL METODE SE LA BISECCIÓ (I ESCRIBIM EN UN ARXIU PER GENERAR
C     LA FIGURA) UN DELS ZEROS DE LA FUNCIÓ
        T=0.91
        P(4)=4*T
        CALL BISECCIO(1.D0/3.D0+0.1D0,2.D0,1.D-12)

C     ARA BUSQUEM ZEROS PEL MÈTODE DE NEWTON-RAPHSON PARTINT DE DIFERENTS PUNTS
        MI=10

        DO I=1,6
          CALL NR(PUNTS(I),MI,X)
          DO J=1,MI
            XX(I,J)=X(J)
            ENDDO
          ENDDO

        OPEN(1,FILE="nrP4.dat")
        DO J=1,MI
          WRITE (1,"(I4,2X,10(F16.12,2X))") J,(XX(I,J),I=1,6)
          ENDDO
        CLOSE(1)

C
        TT(1)=0.86D0
        TT(2)=0.88D0
        TT(3)=0.92D0
        TT(4)=0.94D0
        TT(5)=1.07D0

        DO I=1,100
          V=1.D0/3.D0+0.1D0+(I-1)*(4.D0-1.D0/3.d0-0.1D0)/99.D0
          DO J=1,5
            ISO(I,J)=PP(V,TT(J))
            ENDDO
          ENDDO

C        DO I=1,100

        END



      DOUBLE PRECISION FUNCTION POLI(N,P,V)
        DOUBLE PRECISION P,V
        INTEGER N,I

        DIMENSION P(N+1)
        POLI=0

        DO I=0,N
          POLI=POLI+P(I+1)*V**(I)
          ENDDO
        END FUNCTION

      DOUBLE PRECISION FUNCTION DPOLI(N,P,V)
        DOUBLE PRECISION P,V
        INTEGER N,I

        DIMENSION P(N+1)
        DPOLI=0

        DO I=1,N
          DPOLI=DPOLI+P(I+1)*V**(I-1)*I
        ENDDO
        END FUNCTION


      SUBROUTINE BISECCIO(X1,X2,MAXERROR)
        IMPLICIT NONE

        DOUBLE PRECISION X1,X2,MAXERROR,A,B,ZERO,POLI,P,ERROR
        INTEGER N,I

        DIMENSION P(4)

        COMMON/POLIS/P,N

        A=X1
        B=X2
        I=1
        ERROR=B-A

        OPEN (1,FILE="bisecP4.dat")
        DO WHILE (ERROR.GT.MAXERROR)
          ERROR=B-A
          ZERO=(A+B)/2.D0

          IF (POLI(N,P,A)*POLI(N,P,ZERO).LT.0.D0) THEN
            B=ZERO

          ELSE IF (POLI(N,P,B)*POLI(N,P,ZERO).LT.0) THEN
            A=ZERO

          ELSE IF (POLI(N,P,ZERO).EQ.0) THEN
            EXIT

          ELSE IF (POLI(N,P,A)*POLI(N,P,B).GT.0) THEN
            WRITE (1,*) "La funcion no cambia de signo en ",A,B
            EXIT
            ENDIF
          WRITE (1,*) I,ZERO,ERROR

          I=I+1

          ENDDO
        CLOSE(1)
        END SUBROUTINE



      SUBROUTINE NR(X0,MAXITER,X)
        DOUBLE PRECISION X0,X,X1,P,POLI,DPOLI,DP
        INTEGER N,MAXITER

        DIMENSION P(4)
        DIMENSION X(MAXITER)

        COMMON/POLIS/P,N

        DP=DPOLI(N,P,X0)
        IF (DP.EQ.0) THEN
          DP=-20.D0
          ENDIF

        X1=X0

        DO I=1,MAXITER
          X(I)=X1
          DP=DPOLI(N,P,X(I))

          IF (DP.EQ.0) THEN
            DP=-20.D0
            ENDIF

          X1=X(I)-POLI(N,P,X(I))/DP
          ENDDO
        END SUBROUTINE



      DOUBLE PRECISION FUNCTION PP(V,T)
        DOUBLE PRECISION V,T

        PP=8.D0*T/(3.D0*V-1)-3.D0/V**2
        END FUNCTION



      SUBROUTINE DERIVADA (NDAT,X,FUNCI,DFUNCI)
        INTEGER NDAT,I
        DOUBLE PRECISION X(NDAT),FUNCI(NDAT),DFUNCI(NDAT),H

        H=X(2)-X(1)

        DFUNCI(1)=(FUNCI(2)-FUNCI(1))/H
        DFUNCI(NDAT)=(FUNCI(NDAT)-FUNCI(NDAT-1))/H

        DO I=2,NDAT-1
          DFUNCI(I)=(FUNCI(I+1)-FUNCI(I-1))/(2*H)
          ENDDO
        END SUBROUTINE




