      PROGRAM MIRUNGEKUTTA
        IMPLICIT NONE

C  Y(1) = theta1
C  Y(2) = theta2
C  Y(3) = thetapunt1
C  Y(4) = thetapunt2

c           E_POT = -(M1 + M2) * G * L1 * COS(Y0(1)) -
c     *             M2 * G * L2 * COS(Y0(2))
c           E_CIN = M1 * L1**2 * Y0(3)**2 / 2 + M2 * 
c     *             (L1**2 * Y0(3)**2 + L2**2 * Y0(4)**2 + 
c     *             2 * L1 * L2 * Y0(3) * Y0(4) * COS(Y0(1) - Y0(2))) / 2
c           E_TOT = E_POT + E_CIN



        INTEGER NEQU, I
        PARAMETER(NEQU=4)
        REAL*8 M1, M2, L1, L2, G, PI, T0, TF, DT, T, Y0(NEQU),
     *         Y(NEQU), DTS(5), Y00(2), YDELTA0, DELTA, K, K0

        COMMON/CONST/L1, M1, L2, M2, G

        G=9.8D0
        PI=4.D0 * ATAN(1.D0)


C  Cas límit: M1 >> M2, L1 = L2
        M1 = 1.D0
        M2 = 0.1D0
        L1 = 1.D0
        L2 = 1.D0

C  Condicions de contorn
        T0 = 0.D0
        TF = 20.D0
        DT = 0.01D0
        T  = T0

C  Petites oscilacions
        Y0(1) = 0.1D0
        Y0(2) = 0.2D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

C  Creem el fitxer de dades per el primer cas

        OPEN(10, FILE="11P9.dat")

        DO WHILE (T.LE.TF)
           WRITE (10,*) T,Y0(1), Y0(2)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO

        WRITE (10,*) T,Y0(1), Y0(2)
        CLOSE(10)


C  Condicions de contorn
        T  = T0

C  Grans oscilacions
        Y0(1) = 2.D0
        Y0(2) = 3.D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

C  Creem el fitxer de dades per el segon cas
        OPEN(10, FILE="12P9.dat")

        DO WHILE (T.LE.TF)
           WRITE (10,*) T,Y0(1), Y0(2)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO

        WRITE (10,*) T,Y0(1), Y0(2)
        CLOSE(10)


C  Cas límit: M1 = M2, L1 >> L2
        M1 = 1.D0
        M2 = 1.D0
        L1 = 1.D0
        L2 = 0.1D0

C  Condicions de contorn
        T  = T0

C  Petites oscilacions
        Y0(1) = 0.1D0
        Y0(2) = 0.2D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

C  Creem el fitxer de dades per el tercer cas

        OPEN(10, FILE="13P9.dat")

        DO WHILE (T.LE.TF)
           WRITE (10,*) T,Y0(1), Y0(2)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO

        WRITE (10,*) T,Y0(1), Y0(2)
        CLOSE(10)


C  Condicions de contorn
        T  = T0

C  Grans oscilacions
        Y0(1) = 2.D0
        Y0(2) = 3.D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

C  Creem el fitxer de dades per el quart cas

        OPEN(10, FILE="14P9.dat")

        DO WHILE (T.LE.TF)
           WRITE (10,*) T,Y0(1), Y0(2)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO

        WRITE (10,*) T,Y0(1), Y0(2)
        CLOSE(10)


C------------------------------------------------------------------------
C------------------------------------------------------------------------


C  Estudi dels modes normals: M1 = M2, L1 = L2
        M1 = 1.D0
        M2 = 1.D0
        L1 = 1.D0
        L2 = 1.D0

C  Condicions de contorn
        T0 = 0.D0
        TF = 20.D0
        DT = 0.01D0
        T  = T0

C  Condicions inicials del primer cas
        Y0(1) = 0.1D0
        Y0(2) = 0.15D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

C  Creem el fitxer de dades per el primer cas

        OPEN(10, FILE="21P9.dat")

        DO WHILE (T.LE.TF)
           WRITE (10,*) T,Y0(1), Y0(2)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO

        WRITE (10,*) T,Y0(1), Y0(2)
        CLOSE(10)


C  Condicions de contorn
        T  = T0
        DT = 0.01

C  Condicions inicials del segon cas
        Y0(1) = 3.1D0
        Y0(2) = -3.1D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

C  Creem el fitxer de dades per el segon cas

        OPEN(10, FILE="22P9.dat")

        DO WHILE (T.LE.TF)
           WRITE (10,*) T,Y0(1), Y0(2)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO

        WRITE (10,*) T,Y0(1), Y0(2)
        CLOSE(10)


C------------------------------------------------------------------------
C------------------------------------------------------------------------


C  Estudi de la convergència del mètode amb l'increment de t
        DTS(1) = 0.5D0
        DTS(2) = 0.1D0
        DTS(3) = 0.05D0
        DTS(4) = 0.01D0
        DTS(5) = 0.001D0

        OPEN(10, FILE="31P9.dat")

C  Condicions de contorn
        T0 = 0.D0
        TF = 20.D0
        DT = 0.01D0
        T  = T0

C  Condicions inicials del primer cas

        DO I=1, 5
           Y0(1)  = 0.1D0
           Y0(2)  = -0.15D0
           Y0(3)  = 0.D0
           Y0(4)  = 0.D0
           Y00(1) = Y0(1)
           Y00(2) = Y0(2)

           T  = T0
           DT = DTS(I)
           
           DO WHILE (T.LE.TF)
              CALL MIRK4(T,DT,Y0,NEQU,Y)
              T = T + DT
              Y0 = Y
              ENDDO

           DO WHILE (T.GE.T0)
              CALL MIRK4(T,DT,Y0,NEQU,Y)
              T = T - DT
              Y0 = Y
              ENDDO


           WRITE (10,*) DT,ABS(Y0(1)-Y00(1)), ABS(Y0(2)-Y00(2))

           ENDDO

        CLOSE(10)


C  Estudi de la convergència del mètode amb l'increment de t (2)

        OPEN(10, FILE="32P9.dat")

C  Condicions inicials del segon cas

        DO I=1, 5
           Y0(1)  = PI-0.5D0
           Y0(2)  = PI-0.1D0
           Y0(3)  = 0.D0
           Y0(4)  = 0.D0
           Y00(1) = Y0(1)
           Y00(2) = Y0(2)

           T  = T0
           DT = DTS(I)

           
           DO WHILE (T.LE.TF)
              CALL MIRK4(T,DT,Y0,NEQU,Y)
              T = T + DT
              Y0 = Y
              ENDDO

           DO WHILE (T.GE.T0)
              CALL MIRK4(T,DT,Y0,NEQU,Y)
              T = T - DT
              Y0 = Y
              ENDDO

           WRITE (10,*) DT,ABS(Y0(1)-Y00(1)), ABS(Y0(2)-Y00(2))
           ENDDO

        CLOSE(10)



C------------------------------------------------------------------------
C------------------------------------------------------------------------


C  M1 = M2, L1 = L2
        M1 = 1.D0
        M2 = 1.D0
        L1 = 1.D0
        L2 = 1.D0

C  Condicions de contorn
        T0 = 0.D0
        TF = 6.D0
        DT = 0.01D0

C  Condicions inicials del cas DELTA = 0
        Y0(1) = PI - 0.5D0
        Y0(2) = PI - 0.2D0
        Y0(3) = 0.D0
        Y0(4) = 0.D0

        DO WHILE (T.LE.TF)
           CALL MIRK4(T,DT,Y0,NEQU,Y)
           T = T + DT
           Y0 = Y
           ENDDO
        K = M1 * L1**2 * Y0(3)**2 / 2 + M2 * 
     *      (L1**2 * Y0(3)**2 + L2**2 * Y0(4)**2 + 
     *      2 * L1 * L2 * Y0(3) * Y0(4) * COS(Y0(1) - Y0(2))) / 2

C  Condicions de contorn
        YDELTA0 = Y(1)
        K0 = K
        DT = 0.01D0

C  Creem el fitxer de dades

        OPEN(10, FILE="41P9.dat")

C  Farem servir 10 deltes diferents de 0
        DO I=1, 10

C  Condicions inicials
           DELTA = I * 0.02
           Y0(1) = PI - 0.5D0
           Y0(2) = PI - 0.2D0 + DELTA
           Y0(3) = 0.D0
           Y0(4) = 0.D0
           T  = T0

           DO WHILE (T.LE.TF)
              CALL MIRK4(T,DT,Y0,NEQU,Y)
              T = T + DT
              Y0 = Y
              ENDDO
           K = M1 * L1**2 * Y0(3)**2 / 2 + M2 * 
     *         (L1**2 * Y0(3)**2 + L2**2 * Y0(4)**2 + 
     *         2 * L1 * L2 * Y0(3) * Y0(4) * COS(Y0(1) - Y0(2))) / 2

           WRITE (10,*) DELTA,Y0(1) - YDELTA0, K - K0
           ENDDO
        CLOSE(10)


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
        IMPLICIT NONE
        INTEGER NEQU
        DOUBLE PRECISION X,Y(NEQU),DY(NEQU)
        DOUBLE PRECISION L1,M1,L2,M2,G,PT1,PT2,DD
        COMMON/CONST/L1,M1,L2,M2,G


 
        DD=Y(2)-Y(1)

        DY(1)=Y(3)
        DY(2)=Y(4)

        DY(3) = (M2 * L1 * Y(3)**2 * SIN(DD) * COS(DD) +
     *          M2 * G * SIN(Y(2)) * COS(DD) +
     *          M2 * L2 * Y(4)**2 * SIN(DD) -
     *          (M1 + M2) * G* SIN(Y(1))) /
     *          ((M1 + M2) * L1 - M2 * L1 * COS(DD)**2)

        DY(4) = (-M2 * L2 * Y(4)**2 * SIN(DD) * COS(DD) +
     *          (M1 + M2) * (G * SIN(Y(1)) * COS(DD) -
     *          L1 * Y(3)**2 * SIN(DD) - G * SIN(Y(2)))) /
     *          ((M1 + M2) * L2 - M2 * L2 * COS(DD)**2)

        END
