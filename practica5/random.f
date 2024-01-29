      PROGRAM RANDOM
      IMPLICIT NONE

      INTEGER A,B

      A=12000
      B=50
      CALL SUBUNI(A,B)

      A=8000
      B=80
      CALL SUBAIR(A,B)
      END



      SUBROUTINE SUBUNI(NDIM,NCAIXES)
      IMPLICIT NONE

      INTEGER NDIM,NCAIXES,ISEED,I,J,K,X(NCAIXES)
      DOUBLE PRECISION XX(NDIM),MITJA,MOM(9),B,DESVEST,
     *                 RMITJA,RVAR,RDESVEST,SIG(NCAIXES)

      B=2.D0
      ISEED=16484786
      MITJA=0

      RMITJA=1.D0
      RVAR=1.D0/3.D0
      RDESVEST=SQRT(RVAR)

      CALL SRAND(ISEED)

      OPEN(1,FILE="mynumsP5.dat")
      DO I=1,NDIM
        XX(I)=B*RAND()
        MITJA=MITJA+XX(I)/NDIM
        WRITE (1,*) XX(I)
        ENDDO
      CLOSE(1)

      DO I=1,9
        MOM(I)=0
      ENDDO
      
      OPEN (1,FILE="histo1P5.dat")

      DO I=1,NCAIXES
        X(I)=0
        DO J=1,NDIM
          DO K=1,9
            MOM(K)=MOM(K)+(XX(J)-MITJA)**(K+1)/(NDIM)
            ENDDO
          IF(XX(J).LT.REAL(I)/REAL(NCAIXES).AND.
     *       XX(J).GE.REAL(I-1)/REAL(NCAIXES)) THEN
            X(I)=X(I)+1
            ENDIF
          ENDDO
          SIG(I)=REAL(NCAIXES)/(B*SQRT(REAL(NDIM)))*
     *           (REAL(X(I))/REAL(NDIM)*(1-REAL(X(I))/REAL(NDIM)))
          WRITE (1,100) REAL(I)*B/REAL(NCAIXES),REAL(X(I))/REAL(NDIM),
     *                  SIG(I)
        ENDDO
 100  FORMAT(1P,E12.1,2X,E12.6,2X,E12.6)




c      DESVEST=SQRT(MOM(1))
c      WRITE(1,*) "         MITJA       VARIANCIA   DESVIACIÓ ESTANDARD"
c      WRITE(1,10) MITJA,MOM(1),DESVEST
c      WRITE(1,20) RMITJA,RVAR,RDESVEST
c 10   FORMAT("CALCULAT",E12.4,E12.4,2X,E12.4)
c 20   FORMAT("REAL    ",E12.4,E12.4,2X,E12.4)
      
      END SUBROUTINE



      SUBROUTINE SUBAIR(NACCEP,NCAIXES)
        IMPLICIT NONE
        INTEGER NACCEP,NCAIXES,ISEED,I,J
        DOUBLE PRECISION Y,P,PI,A,B,M,PY,X(NACCEP),XX(NACCEP),
     *            MITJA,VAR,DESVEST,SIG(NCAIXES)

        PI=4.D0*ATAN(1.D0)
        A=-PI
        B=PI
        M=1.D0/4.D0
        I=0
        MITJA=0
        VAR=0
        DESVEST=0

        ISEED=16484786
        CALL SRAND(ISEED)
        DO WHILE (I.NE.NACCEP)
          Y=A+(B-A)*RAND()
          P=M*RAND()
          PY=M*SIN(ABS(Y))
          IF (P.LE.PY) THEN
            I=I+1
            XX(I)=Y
            MITJA=MITJA+Y
            ENDIF
          ENDDO

        MITJA=MITJA/NACCEP

        OPEN (1,FILE="histo2P5.dat")
        DO I=1,NCAIXES
          DO J=1,NACCEP
            IF (XX(J).LT.A+(B-A)*REAL(I)/REAL(NCAIXES).AND.
     *         XX(J).GE.A+(B-A)*REAL(I-1)/REAL(NCAIXES)) THEN
              X(I)=X(I)+1
              ENDIF
            ENDDO
            SIG(I)=REAL(NCAIXES)/((B-A)*SQRT(REAL(NACCEP)))*
     *           (REAL(X(I))/REAL(NACCEP)*(1-REAL(X(I))/REAL(NACCEP)))
          WRITE (1,*) A+REAL(I)*(B-A)/REAL(NCAIXES),X(I)/NACCEP,SIG(I)
          ENDDO
        END SUBROUTINE




















