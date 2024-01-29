      PROGRAM LAPLACE
        INTEGER IOPT
        REAL*8 TI

        IOPT = 1

        open (10, file="1011.dat")
        TI = 10.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="1012.dat")
        TI = 20.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="1013.dat")
        TI = 10000.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        IOPT = 2

        open (10, file="1021.dat")
        TI = 10.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="1022.dat")
        TI = 20.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="1023.dat")
        TI = 10000.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        IOPT = 3

        open (10, file="1031.dat")
        TI = 10.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="1032.dat")
        TI = 20.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="1033.dat")
        TI = 10000.D0
        CALL DIFUSSIO (TI, IOPT)
        CLOSE (10)

        open (10, file="101.dat")
        TI = 20.D0
        CALL DIFUSSIO2 (TI)
        CLOSE (10)


        END


      SUBROUTINE DIFUSSIO (TI, IOPT)
        INTEGER IOPT, I, J, LX, LY, NX, NY, K
        PARAMETER (LX=35.D0, LY = 25.D0)
        REAL*8 H, TOLD, TNEW, TI, ERROR, ERRORMAX, W, P, P0, R1, R2,
     *         X, Y
        PARAMETER (H = 0.5D0, P0 = 10.D0, W = 1.5D0, ERRORMAX = 0.01D0)
        PARAMETER (NX = LX/H, NY = LY/H)

        DIMENSION TOLD (0:NX, 0:NY), TNEW(0:NX, 0:NY)

        ERROR = 10.D0

C  Definim les condicions de contorn
        DO I = 0, NX
C  Temperatura a (x,Ly)
          TOLD (I, NY) = 20.D0
          TNEW (I, NY) = 20.D0
C  Temperatura a (x,0)
          TOLD (I, 0)  = 20.D0
          TNEW (I, 0)  = 20.D0
          ENDDO

        DO J = 0, NY
C  Temperatura a (0,y)
          TOLD (0, J)  = 20.D0
          TNEW (0, J)  = 20.D0
C  Temperatura a (lx,y)
          TOLD (NX, J) = 20.D0
          TNEW (NX, J) = 20.D0
          ENDDO

C  Temperatura interior
        DO I=1, NX-1
          DO J=1, NY-1
            TOLD(I,J) = TI
            TNEW(I,J) = TI
            ENDDO
          ENDDO

        K = 0

C  Gauss-Seidel
        IF (IOPT.EQ.1) THEN

          DO WHILE (ERROR.GT.ERRORMAX)
            k = k+1
            ERROR = 0.D0
C  Calculem la nova temperatura amb el mètode de Gauss-Seidel i l'error
            DO I=1, NX-1
              DO J=1, NY-1
                X = I * H
                Y = J * H
                R1 = SQRT((X-7)**2  + (Y-20)**2)
                R2 = SQRT((X-25)**2 + (Y-5)**2)
                P = P0 * EXP(-(R1-5.D0)**2 / 0.5**2) +
     *              P0 * EXP(-(R2-3.D0)**2 / 0.5**2)
                TNEW(I,J) = (TOLD(I+1,J) + TOLD(I-1,J) +
     *                      TOLD(I,J+1) + TOLD(I,J-1) + P*H*H) / 4.D0
                ERROR = ERROR + ABS((TOLD(I,J) - TNEW(I, J)) /
     *                  (TOLD(I,J) + TNEW(I, J)))
                ENDDO
              ENDDO

C  Reasignem els valors de Told per els que acabem de calcular i calculem la temperatura al punt (8,10)
            DO I=1, NX
              DO J=1, NY
                TOLD(I,J) = TNEW(I,J)
                ENDDO
              ENDDO
            TI = TNEW(16, 20)
            WRITE (10,*) K, TI
            ENDDO

C  Jacobi
        ELSE IF (IOPT.EQ.2) THEN

          DO WHILE (ERROR.GT.ERRORMAX)
            k = k+1
            ERROR = 0.D0
C  Calculem la nova temperatura amb el mètode de Jacobi i l'error
            DO I=1, NX-1
              DO J=1, NY-1
                X = I * H
                Y = J * H
                R1 = SQRT((X-7)**2  + (Y-20)**2)
                R2 = SQRT((X-25)**2 + (Y-5)**2)
                P = P0 * EXP(-(R1-5.D0)**2 / 0.5**2) +
     *              P0 * EXP(-(R2-3.D0)**2 / 0.5**2)
                TNEW(I,J) = (TNEW(I+1,J) + TNEW(I-1,J) + 
     *                      TNEW(I,J+1) + TNEW(I,J-1) + P*H*H) / 4.D0
                ERROR = ERROR + ABS((TOLD(I,J) - TNEW(I, J)) /
     *                  (TOLD(I,J) + TNEW(I, J)))
                ENDDO
              ENDDO

C  Reasignem els valors de Told per els que acabem de calcular i calculem la temperatura interior
            DO I=1, NX
              DO J=1, NY
                TOLD(I,J) = TNEW(I,J)
                ENDDO
              ENDDO
            TI = TNEW(16, 20)
            WRITE (10,*) K, TI
            ENDDO

C  Sobre-relaxació
        ELSE IF (IOPT.EQ.3) THEN
          DO WHILE (ERROR.GT.ERRORMAX)
            k = k+1
            ERROR = 0.D0
C  Calculem la nova temperatura amb el mètode de sobre-relaxació i l'error
            DO I=1, NX-1
              DO J=1, NY-1
                X = I * H
                Y = J * H
                R1 = SQRT((X-7)**2  + (Y-20)**2)
                R2 = SQRT((X-25)**2 + (Y-5)**2)
                P = P0 * EXP(-(R1-5.D0)**2 / 0.5**2) +
     *              P0 * EXP(-(R2-3.D0)**2 / 0.5**2)
                TNEW(I,J) = TOLD(I,J) + W * 
     *                      (TOLD(I+1,J) + TNEW(I-1,J) + 
     *                      TOLD(I,J+1) + TNEW(I,J-1) -
     *                      4.D0 * TOLD(I,J) + P*H*H) / 4.D0
                ERROR = ERROR + ABS((TOLD(I,J) - TNEW(I, J)) /
     *                  (TOLD(I,J) + TNEW(I, J)))
                ENDDO
              ENDDO

C  Reasignem els valors de Told per els que acabem de calcular i calculem la temperatura interior
            DO I=1, NX
              DO J=1, NY
                TOLD(I,J) = TNEW(I,J)
                ENDDO
              ENDDO
            TI = TNEW(16, 20)
            WRITE (10,*) K, TI
            ENDDO

          ENDIF
        END SUBROUTINE



      SUBROUTINE DIFUSSIO2 (TI)
        INTEGER I, J, LX, LY, NX, NY
        PARAMETER (LX=35.D0, LY = 25.D0)
        REAL*8 H, TOLD, TNEW, TI, ERROR, ERRORMAX, W, P, P0, R1, R2,
     *         X, Y
        PARAMETER (H = 0.1D0, ERRORMAX = 0.01D0, W = 1.5D0, P0 = 10.D0)
        PARAMETER (NX = LX/H, NY = LY/H)

        DIMENSION TOLD (0:NX, 0:NY), TNEW(0:NX, 0:NY)

        ERROR = 10.D0

C  Definim les condicions de contorn
        DO I = 0, NX
C  Temperatura a (x,Ly)
          TOLD (I, NY) = 20.D0
          TNEW (I, NY) = 20.D0
C  Temperatura a (x,0)
          TOLD (I, 0)  = 20.D0
          TNEW (I, 0)  = 20.D0
          ENDDO

        DO J = 0, NY
C  Temperatura a (0,y)
          TOLD (0, J)  = 20.D0
          TNEW (0, J)  = 20.D0
C  Temperatura a (lx,y)
          TOLD (NX, J) = 20.D0
          TNEW (NX, J) = 20.D0
          ENDDO

C  Temperatura interior
        DO I=1, NX-1
          DO J=1, NY-1
            TOLD(I,J) = TI
            TNEW(I,J) = TI
            ENDDO
          ENDDO

        DO WHILE (ERROR.GT.ERRORMAX)
          ERROR = 0.D0
C  Calculem la nova temperatura amb el mètode de sobre-relaxació i l'error
          DO I=1, NX-1
            DO J=1, NY-1
              X = I * H
              Y = J * H
              R1 = SQRT((X-7)**2  + (Y-20)**2)
              R2 = SQRT((X-25)**2 + (Y-5)**2)
              P = P0 * EXP(-(R1-5.D0)**2 / 0.5**2) +
     *            P0 * EXP(-(R2-3.D0)**2 / 0.5**2)

              TNEW(I,J) = TOLD(I,J) + W * 
     *                    (TOLD(I+1,J) + TNEW(I-1,J) + 
     *                    TOLD(I,J+1) + TNEW(I,J-1) -
     *                    4.D0 * TOLD(I,J) + P*H*H) / 4.D0

              ERROR = ERROR + ABS((TOLD(I,J) - TNEW(I, J)) /
     *                (TOLD(I,J) + TNEW(I, J)))
              ENDDO
            ENDDO

C  Reasignem els valors de Told per els que acabem de calcular
          DO I=1, NX
            DO J=1, NY
              TOLD(I,J) = TNEW(I,J)
              ENDDO
            ENDDO
          ENDDO

C  Finalment escribim la temperatura a cada punt
        DO I=1, NX
          DO J=1, NY
            WRITE (10,*) I, J, TNEW(I,J)
            ENDDO
            WRITE (10,*) ""
          ENDDO

        END SUBROUTINE
















