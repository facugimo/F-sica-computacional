      PROGRAM MIRUNGEKUTTA
        IMPLICIT NONE
        INTEGER NEQ,N,I,K
        PARAMETER (NEQ=2)
        REAL*8 T,H,Y0(NEQ),Y(NEQ),T0,TF,B,BETA(7),W,G,YE,Y00,W1,NN(4)

        COMMON/CONDICIONS/B,W

C  Definim les constants que farem servir
        N=10000
        W=1.6D0

C  Apartat 2a)

        BETA(1)=0.D0
        BETA(2)=0.96D0
        BETA(3)=3.136D0
        BETA(4)=4.8D0

        OPEN(10,FILE="11.dat")
        OPEN(20,FILE="12.dat")
        OPEN(30,FILE="13.dat")
        OPEN(40,FILE="14.dat")

        DO K=1,4
C  Definim les condicions inicials
          B=BETA(K)
          Y0(1)=1.D0
          Y0(2)=0.D0
          T0=0.D0
          TF=15.D0/W
          T=T0
          Y=Y0

          H=(TF-T0)/(N-1)

C  I fem servir el mètode de Runge-Kutta 4 per calcular la 
C  posició a tot el rang de temps
          DO I=1,N
            WRITE (10*k,*) T,Y(1)
            CALL MIRK4(T,H,Y0,NEQ,Y)
            Y0=Y
            T=T0+H*REAL(I)
            ENDDO
          CLOSE(10*K)
        ENDDO

C-------------------------------------------------------------------

C  Apartat 2B)

        BETA(5)=1.28D0
        BETA(6)=2.24D0
        BETA(7)=3.104D0

        OPEN(50,FILE="21.dat")
        OPEN(60,FILE="22.dat")
        OPEN(70,FILE="23.dat")

        DO K=5,7
C  Definim les condicions inicials
          N=40
          B=BETA(K)
          W1=W*SQRT(1.D0-B**2/(4.D0*W**2))
          G=B/2
          Y0(1)=1.D0
          Y00=Y0(1)
          Y0(2)=-G
          T0=0.D0
          TF=15.D0/W
          T=T0
          Y=Y0

          H=(TF-T0)/(N-1)

C  I fem servir el mètode de Runge-Kutta 4 per calcular la 
C  posició a tot el rang de temps
          DO I=1,N
            YE=Y00*EXP(-G*T)*COS(W1*T)
            WRITE (10*k,*) T,Y(1),YE
            CALL MIRK4(T,H,Y0,NEQ,Y)
            T=T0+H*REAL(I)
            Y0=Y
            ENDDO
          CLOSE(10*K)
        ENDDO
C------------------------------------------------------------------------------

C  Apartat 3)

        BETA(5)=1.28D0

        OPEN(10,FILE="31.dat")
        OPEN(20,FILE="32.dat")
        OPEN(30,FILE="33.dat")
        OPEN(40,FILE="34.dat")

        NN(1)=100
        NN(2)=200
        NN(3)=500
        NN(4)=20000


        DO k=1,4
C  Definim les condicions inicials
          N=NN(K)
          B=BETA(5)
          Y0(1)=1.D0
          Y0(2)=0.D0
          T0=0.D0
          TF=15.D0/W
          T=T0
          Y=Y0
          YE=Y0(1)

          H=(TF-T0)/(N-1)

C  I fem servir el mètode de Runge-Kutta 4 i el de Euler 
C  per calcular la posició a tot el rang de temps
          DO I=1,N
            WRITE (10*K,*) T,Y(1),YE

            CALL EULER(T,H,Y0,NEQ,YE)
            CALL MIRK4(T,H,Y0,NEQ,Y)
            T=T0+H*REAL(I)
            Y0=Y
            ENDDO
          CLOSE(50)
        ENDDO


        END




      SUBROUTINE MIRK4(X,DX,YY,NEQU,YOUT)
        IMPLICIT NONE
        INTEGER NEQU
        REAL*8 X,DX,YY(NEQU),YOUT(NEQU),K1(NEQU),K2(NEQU),K3(NEQU),
     *         K4(NEQU)

C  Calculem els coeficients k1,k2,k3 i k4 a partir de l'equació diferencial
        CALL DERIVADES(X,YY,K1,NEQU)
        CALL DERIVADES(X+DX/2.D0,YY+K1*DX/2.D0,K2,NEQU)
        CALL DERIVADES(X+DX/2.D0,YY+K2*DX/2.D0,K3,NEQU)
        CALL DERIVADES(X+DX,YY+K3*DX,K4,NEQU)

C  Fem servir els coeficients calculats per trobar el seguent pas y(t+h)
        YOUT=YY+DX/6.D0*(K1+2.D0*K2+2.D0*K3+K4)
        END SUBROUTINE



      SUBROUTINE DERIVADES (X,Y,DY,NEQU)
        INTEGER NEQU,I
        REAL*8 X,Y(NEQU),DY(NEQU),W,B

        COMMON/CONDICIONS/B,W

        W=1.6D0

C  Calculem les derivades a partir de l'equació diferencial
        DO I=1,NEQU-1
          DY(I)=Y(I+1)
          ENDDO

        DY(NEQU)=-W**2*Y(1)-B*Y(2)

        END SUBROUTINE


      SUBROUTINE EULER (X,H,Y,NEQU,YOUT)
        INTEGER NEQU
        REAL*8 X,H,Y(NEQU),YOUT

        YOUT=Y(1)+H*Y(2)
        END SUBROUTINE
