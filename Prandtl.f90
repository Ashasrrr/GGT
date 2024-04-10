
!************************************************************************************************
SUBROUTINE Prandtl(U,V,P,ro,NI,NJ,Su_max,EPS,dx,dy,NU,U0,H, gama, ro0,P0, eps_G, mu,sp_max)
    IMPLICIT NONE
	INTEGER NI,NJ,S,I,J,K,su_max, sp_max,sp
	REAL*8 dx,dy,NU,U0,H, gama, p0, mu
	REAL*8 EPS,ERRO,SU,SV, eps_G, G, G0,Ut(NJ), Vt(NJ), Pt
    REAL*8 MaxU, MaxV, ro0, ro(NI)
        REAL*8,DIMENSION(NI,NJ):: U,V,P
    REAL*8, DIMENSION(NJ) :: A, B, C, D, QU

    QU(:) = 0.0


    G0 = ro0 * U0 * h
    Ut = 0.
    Vt = 0.
    pt = 0.
       DO I=2,NI
                   u(i,:) = u(i-1,:)
           v(i,:) = v(i-1,:)
         p(i,:)=p(i-1,1)
         ro(i)=ro(i-1)
         sp = 0
         do !Sp
           u(i,:) = u(i-1,:)
           v(i,:) = v(i-1,:)
           ERRO = 1000.0
           k=0
           DO WHILE (K .LE. Su_max .AND. ERRO .GT. EPS) !Su
             DO J=2,NJ-1
             	A(j) = -ro(i) * v(i, j-1) / (2. * dy) - mu / dy**2
             	B(j) = ro(i) * u(i, j) / dx + 2. * mu / dy**2
            	C(j) = ro(i) * v(i, j+1) / (2. * dy) - mu / dy**2
               	D(j) = ro(i-1) * u(i-1, j)**2 / dx - (p(i,j) - p(i-1, j)) / dx
             ENDDO
             A(1) = 0.0
             B(1) = 1.0
             C(1) = 0.0
             D(1) = 0.0

             A(NJ) = -1.0
             B(NJ) = 1.0
             C(NJ) = 0.0
             D(NJ) = 0.0

             call Solve_Tridiagonal_Matrix(NJ, A, B, C, D, QU)

             DO J=1,NJ
                U(I,J) = QU(J)
             ENDDO

             DO J=2,NJ
		v(i,j) = v(i, j-1) - dy / (2. * dx) * (u(i, j) + u(i, j-1) - &
                        ro(i-1) / ro(i) * (u(i-1,j) + u(i-1, j-1)))
             ENDDO

             ERRO = 0.0
             MaxU = 0.0
             MaxV = 0.0
             Su=0
             Sv=0
                DO J=1,NJ
                    MaxU = max( abs(U(I,J)), MaxU)
                    MaxV = max( abs(V(I,J)), MaxV)
                    SU = max( abs( U(I,J) - Ut(J) ), SU )
                    SV = max( abs( V(I,J) - Vt(J) ), SV )
                ENDDO
                SU = SU / MaxU
                SV = SV / MaxV
                ERRO = MAX(SU,SV)


                Vt = v(i,:)
                Ut = u(i,:)
             K = K+1

           ENDDO !Su

           print*, k
           !расчёт расхода
           G = 0.
           DO j = 2, nj
               G = G + u(i,j) * ro(i) * dy
           END DO
            !print*, sp, g
           !расчёт давления
           IF (ABS(G0 - G) / G0 < EPS_G) THEN
                print*, abs(p(i,1) - pt) / abs(p(i,1)*100)
                EXIT
           END IF
           pt = p(i,1)
           p(i,:) = p(i,1) + 0.1 * G0 / h**2 * (G - G0)
           ro(i) = ro0 * ( p(i,1)/P0)**(1./gama)

           sp = sp + 1

         ENDDO !Sp
       Enddo !i

	END SUBROUTINE
