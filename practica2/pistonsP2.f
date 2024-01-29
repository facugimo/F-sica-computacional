      PROGRAM PISTONSP2
	IMPLICIT NONE

	INTEGER T
	REAL W,L,T_,X,XI,YI,YZ,YL
	DIMENSION X(5)
	DIMENSION XI(41)
	DIMENSION YI(41)
c	Declarem i dimensionem les variables que farem servir
	
	COMMON/POSICIONS/XI,YI	
C	Els vectors XI i YI es trobaran en un bloc common
	
	W=3.3
	L=11.1
C	Inicialitzem les variables amb els valors donats
	
	
	OPEN (1,FILE="posiP2.dat")

	DO T=0,40
	   T_=T/10.
	   CALL POSI(W,L,T_,X)
	   WRITE (1,10) T_,X(1),X(2),X(3),X(4),X(5)
 10		FORMAT (F12.1,F12.05,F12.05,F12.05,F12.05,F12.5)
	   ENDDO
C	Calculem la posició de cada pistó en l rang de temps donat
C	i guardem el resultat en un arxiu
	CLOSE (1)

	OPEN (2,FILE="posiP2.dat")

	DO T=1,41
	   READ (2,*) XI(T),T_,T_,YI(T)
	   ENDDO
	CLOSE(1)
C	Tornem a obrir l'arxiu i llegim la columna del temps i la del
C	pisto escollit i guardem els valors als vectors del bloc common

	OPEN (1,FILE="interP2.dat")

	DO T=1,2000
	   T_=REAL(T)*4./2000.
	   CALL INTERPO(T_,YZ,YL)
	   WRITE (1,20) T_,YZ,YL
 20		FORMAT (F12.3,F12.5,F12.5)
	   ENDDO
C	Ara calculem la interpolació amb 2000 valors de temps i ho
C	guardem en un arxiu

	CLOSE(1)

	END


      REAL FUNCTION RADI(I)
	INTEGER I

	RADI = EXP(REAL(I)/3.)
	END FUNCTION
C	calculem Radi(i)

      REAL FUNCTION PHI(N)
	REAL PI
	INTEGER N
	PI=2*ASIN(1.)
	PHI = REAL(N)*PI/3.
	END FUNCTION
C	Calculem phi(n)

      SUBROUTINE POSI(W,L,T,X)
	INTEGER J
	REAL W,L,T,X,PHI,RADI
	DIMENSION X(5)

	DO J=1,5
	   R=RADI(J)
	   X(J) = RADI(J)*COS(W*T+PHI(J))+
     *   	  SQRT(L**2-(RADI(J)**2)*SIN(W*T+PHI(J)))
	   ENDDO
	
      END SUBROUTINE
C	Calculem la posició dels 5 pistons en funcio del radi i del temps

      SUBROUTINE INTERPO (X,YZERO,YLIN)
	REAL X,YZERO,YLIN
	REAL XI,YI
	DIMENSION XI(41),YI(41)
	COMMON/POSICIONS/XI,YI

	DO I=1,41
	   IF (X.EQ.XI(I)) THEN
		YLIN = YI(I)
		YZERO=YI(I)
		EXIT

	   ELSE IF (X.LT.XI(I)) THEN
		YLIN=((YI(I)-YI(I-1))/(XI(I)-XI(I-1)))*X
     * 		  +YI(I-1)-(YI(I)-YI(I-1))/(XI(I)-XI(I-1))*XI(I-1)

		YZERO=YI(I-1)
		
C	y=mx+n	m=(y2-y1)/(x2-x1)	n=y2-mx2
		EXIT
		ENDIF
	   ENDDO
	END SUBROUTINE
C	Ara fem una interpolació a partir dels vectors al bloc common
