      PROGRAM PRACTICA1

	INTEGER IESTADO
	REAL N, PN, SN
C Declarem les variables

	SN=0.
	
	WRITE (*,*) "ESCRIU K"
	READ (*,*,IOSTAT=IESTADO) N
C Demanem el valor de k per calcular P_k

	DO WHILE ((IESTADO.NE.0).OR.(N.LT.50).OR.(N.GT.100))
	   WRITE(*,*) "SI US PLAU, INTRODUEIX UN NOMBRE ENTER ENTRE", 
     A		       " 1 I 200:"
	   READ (*,*,IOSTAT=IESTADO) N
	   ENDDO
C En cas de no entregar un valor vàlid per a k es torna a demanar
C especificant el rang de valors acceptat [50,100]

	CALL P_K(N,PN)
	WRITE (*,*) PN
C Calculem P_k per a la k introduida amb la subrutina P_K(K,PK)

	CALL S_N(43,SN)
	WRITE(*,*) SN
C Ara calculem quan val el sumatori de P_k (k=1,...,43) amb la 
C subrutina S_N (N,SN) 

	OPEN (1,FILE="resulP1.dat")
	DO I=3,41,2
	   CALL S_N(I,PN)
	   WRITE(1,*) I,PN
	   ENDDO
C I per últim calculem els sumatoris de P_k (k=1,...,N (N=3,5,...,41))

	END



      SUBROUTINE P_K(K,PK)
	REAL K, PK
	PK = (REAL(K)**2/3.+4*K-5)
	END SUBROUTINE
C La subrutina P_K(K,PK) calcula P_K per a una k donada

      SUBROUTINE S_N(N,SN)
	INTEGER N, I
	REAL SN, PN
	SN=0
	
	DO I=1,N
	   CALL P_K(REAL(I),PN)
	   SN=SN+PN
	   ENDDO
	END SUBROUTINE
C La subrutina S_N(S,SN) calcula el sumatori de P_k (k=1,...,N)
C per a una N donada cridant recursivament la subrutina P_K (K,PK)
