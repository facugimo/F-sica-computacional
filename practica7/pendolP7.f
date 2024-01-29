      PROGRAM PENDOLP7
        IMPLICIT NONE
        REAL*8 PI,T0,TF,TN,WN,G,L
        INTEGER N

        PI=4.D0*ATAN(1.D0)
        G=9.8D0
        L=1.35D0
        WN=SQRT(G/L)
        TN=2.D0*PI/WN
        T0=0.D0
        TF=22.D0*TN
        N=10000

        CALL EULER_MILLORAT(T0,TF,0.03D0,0.D0,N,"e11.dat")
        CALL EULER_MILLORAT(T0,TF,0.12D0,0.D0,N,"e12.dat")
        CALL EULER_MILLORAT(T0,TF,0.25D0,0.D0,N,"e13.dat")

        CALL EULER_MILLORAT(T0,TF,1.5D0,0.D0,N,"e21.dat")
        CALL EULER_MILLORAT(T0,TF,2.5D0,0.D0,N,"e22.dat")
        CALL EULER_MILLORAT(T0,TF,PI-0.01D0,0.D0,N,"e23.dat")
        CALL EULER_MILLORAT(T0,TF,PI-0.001D0,0.D0,N,"e24.dat")

        CALL EULER_MILLORAT(T0,TF,1.2D0,0.D0,N,"en1.dat")
        CALL EULER_MILLORAT(T0,TF,PI-0.001D0,0.D0,N,"en2.dat")

        CALL EULER_MILLORAT(T0,TF,0.D0,2*WN-0.01D0,N,"tr1.dat")
        CALL EULER_MILLORAT(T0,TF,0.D0,2*WN+0.01D0,N,"tr2.dat")

        CALL EULER_MILLORAT(T0,TF,1.86D0,0.D0,500,"en1.dat")
        CALL EULER_MILLORAT(T0,TF,1.86D0,0.D0,4000,"en2.dat")
        CALL EULER_MILLORAT(T0,TF,1.86D0,0.D0,10000,"en3.dat")
        CALL EULER_MILLORAT(T0,TF,1.86D0,0.D0,25000,"en4.dat")

        END

      SUBROUTINE EULER_MILLORAT(T0,TF,THETA0,P0,N,ARXIU)
        IMPLICIT NONE
        INTEGER N,I
        REAL*8 T0,TF,THETA0,P0,T1,T2,THETA1,THETA2,THETA3,P1,P2,P3,H,G
     *         ,L,E_CIN,E_POT,E_TOT,M
        CHARACTER*7 ARXIU

C  Inicialitzem les constants que haurem de fer servir
        H=(TF-T0)/REAL(N-1)
        G=9.8D0
        L=1.35D0
        M=1.3D0

C  Guardem a l'arxiu les condicions inicials (t=0s)
        OPEN(1,FILE=ARXIU)
        T1=T0
        THETA1=THETA0
        P1=P0
        E_CIN=M*P1**2*L**2/2.D0
        E_POT=-M*G*L*COS(THETA1)
        E_TOT=E_CIN+E_POT

        WRITE(1,*) T1,THETA1,P1,E_CIN,E_POT,E_TOT

C  Fem un pas amb el mètode d'Euler per trobar un segon punt (t=h)
        THETA2=THETA0+H*P0
        P2=P0-G/L*H*SIN(THETA0)
        E_CIN=M*P2**2*L**2/2.D0
        E_POT=-M*G*L*COS(THETA2)
        E_TOT=E_CIN+E_POT

        WRITE(1,*) T1+H,THETA2,P2,E_CIN,E_POT,E_TOT

C  I a continuació realitzem els calculs mitjançant el mètode d'Euler
C  millorat amb els dos valors immediatament anteriors
        DO I=2,N-1
          T1=T0+I*H
          THETA3=THETA1+2*H*P2
          P3=P1-2*H*G/L*SIN(THETA2)
          E_CIN=M*(P3**2)*(L**2)/2.D0
          E_POT=-M*G*L*COS(THETA3)
          E_TOT=E_CIN+E_POT

C  Reanomenem les variables per tal que representin els dos ultims valors calculats
          P1=P2
          P2=P3
          THETA1=THETA2
          THETA2=THETA3
          WRITE(1,*) T1,THETA3,P3,E_CIN,E_POT,E_TOT
          ENDDO
        CLOSE(1)
        END SUBROUTINE
