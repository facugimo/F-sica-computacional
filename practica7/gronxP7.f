      PROGRAM GRONXP7
        IMPLICIT NONE
        REAL*8 PI,T0,TF,TN,WN,G,L0
        INTEGER N

        PI=4.D0*ATAN(1.D0)
        G=9.8D0
        L0=1.98D0
        WN=SQRT(G/L0)
        TN=2.D0*PI/WN
        T0=0.D0
        TF=20.D0*TN
        N=10000

        CALL EULER_MILLORAT(T0,TF,0.03D0,0.D0,N,"aaa.dat")

        END

      SUBROUTINE EULER_MILLORAT(T0,TF,THETA0,P0,N,ARXIU)
        IMPLICIT NONE
        INTEGER N,I,J
        REAL*8 T0,TF,THETA0,P0,T1,T2,THETA1,THETA2,THETA3,P1,P2,P3,H,G
     *         ,L,L0,M,ALPHA,THETA_MAX
        CHARACTER*7 ARXIU

C  Inicialitzem les constants que haurem de fer servir
        H=(TF-T0)/REAL(N-1)
        G=9.8D0
        L0=1.98D0
        M=1.3D0


C  Guardem a l'arxiu les condicions inicials (t=0s)
        OPEN(1,FILE=ARXIU)

        DO J=1,50
          ALPHA=REAL(J)/10.D0*SQRT(G/L0)
          T1=T0
          THETA1=THETA0
          P1=P0
          L=L0

          THETA_MAX=ABS(THETA1)
C  Fem un pas amb el mètode d'Euler per trobar un segon punt (t=h)
          L=L0*(1+0.05*SIN(ALPHA*H)
          THETA2=THETA0+H*P0
          P2=P0-G/L*H*SIN(THETA0)

          IF (ABS(THETA2).GT.THETAMAX) THEN
            THETA_MAX=THETA2
            ENDIF

C  I a continuació realitzem els calculs mitjançant el mètode d'Euler
C  millorat amb els dos valors immediatament anteriors
          DO I=2,N-1
            T1=T0+I*H
            L=L0*(1+0.05*SIN(ALPHA*T1)
            THETA3=THETA1+2*H*P2
            P3=P1-2*H*G/L*SIN(THETA2)

            IF (ABS(THETA3).GT.THETA_MAX) THEN
              THETA_MAX=THETA3
              ENDIF

C  Reanomenem les variables per tal que representin els dos ultims valors calculats
            P1=P2
            P2=P3
            THETA1=THETA2
            THETA2=THETA3
            ENDDO

          WRITE(1,*) ALPHA,THETA_MAX
          ENDDO
        CLOSE(1)
        END SUBROUTINE

C  PER FALTA DE TEMPS NOMES HE POGUT IMPLEMENTA LA VARIACIÓ DE L AMB EL TEMPS
