      PROGRAM FERMIONS
        IMPLICIT NONE

        CALL SUBGAUSS()
        CALL MCARLOMD()

        END

C  Apartat 1
      SUBROUTINE SUBGAUSS()
        IMPLICIT NONE
        REAL*8 X,R,PHI,Y,ERRORS(100),P,INTEGRAL,HISTO(100),ERROR,
     *         I2
        INTEGER ISEED,I,J,N

        COMMON/RANDOM/Y
        DIMENSION Y(10000)

        N=10000
        ISEED=16484786
        INTEGRAL=0
        I2=0

        CALL SRAND(ISEED)

        DO J=1,100
          HISTO(J)=0
          ERRORS(J)=0.D0
          ENDDO

        DO I=1,N
          R=SQRT(-2.D0*LOG(RAND()))
          PHI=8.D0*ATAN(1.D0)*RAND()
          Y(I)=R*COS(PHI)
          
          DO J=1,100
            IF (Y(I).LE.(-5.D0+REAL(J)/10.D0)) THEN
              HISTO(J)=HISTO(J)+1.D0/REAL(N)
              EXIT
              ENDIF
            P=HISTO(J)
            ERRORS(J)=(P*(1.D0-P))**2
            ENDDO
          ENDDO

        OPEN(1,FILE="histogaus.dat")

        DO J=1,100
          WRITE(1,*) -5.D0+REAL(J)/10.D0,HISTO(J),ERRORS(J)
          ENDDO
        CLOSE(1)

        DO I=1,N
          INTEGRAL=INTEGRAL+COSH(Y(I))**2
          I2=I2+COSH(Y(I))**4
          ENDDO
        INTEGRAL=INTEGRAL/REAL(N)
        I2=I2/REAL(N)**2
        ERROR=SQRT(ABS(I2-INTEGRAL**2)/REAL(N))

        OPEN(1,FILE="valorsI3.dat")

        WRITE(1,*) INTEGRAL,ERROR

        CLOSE(1)
        END SUBROUTINE

      SUBROUTINE MCARLOMD()
        REAL*8 Y,X(5,2000),PROD,SUMA
        INTEGER K,N,NMAX,I,J

        COMMON/RANDOM/Y

        DIMENSION Y(10000)

        DO I=1,5
          DO J=1,2000
            X(I,J)=Y(-2000+2000*I+J)
            ENDDO
          ENDDO

        OPEN(1,FILE="multiMD.dat")

        DO K=1,10
          INTEGRAL=0
          NMAX=200*K
          DO N=1,NMAX
            PROD=1.D0
            DO I=1,4
              DO J=I+1,5
               PROD=PROD*(X(I,N)-X(J,N))**2
                ENDDO
              ENDDO
            SUMA=0
            DO I=1,5
              SUMA=SUMA+X(I,N)
              ENDDO
            INTEGRAL=INTEGRAL+SUMA**2*PROD*256.D0/675.D0
            ENDDO
          WRITE(1,*) INTEGRAL/REAL(NMAX)
          ENDDO
        CLOSE(1)

        END SUBROUTINE























