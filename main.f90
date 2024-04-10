        Program Pr
        Implicit none

        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER NI, NJ, NITER, Su_MAX, PressureDensityCoupling
        INTEGER I,J,S, sp_max
        REAL*8 L,H,U0,MU,Nu,Ro0,P0, gama
        REAL*8 dx,dy,CFL,EPS, eps_G
        REAL*8,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL*8,ALLOCATABLE :: U_n(:,:),V_n(:,:),P_n(:,:),R_n(:,:)
        REAL*8,ALLOCATABLE :: U(:,:),V(:,:),P(:,:),R(:,:), tw(:), delt(:), Re_x(:), tw_b(:), ro(:,:)

        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) EPS
        read(IO,*) EPS_G
        read(IO,*) Su_MAX
        read(IO,*) Sp_max


        read(IO,*) U0
        read(IO,*) MU
        read(IO,*) RO0
        read(IO,*) P0
        read(IO,*) gama

        CLOSE(IO)

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

!----------------- Node variables -----------------------------
        allocate(U_n(NI,NJ))  ! Velocity U
        allocate(V_n(NI,NJ))  ! Velocity V
        allocate(P_n(NI,NJ))  ! Pressure
        allocate(U(NI,NJ))
        allocate(V(NI,NJ))
        allocate(P(NI,NJ))
        allocate(tw(NI))
        allocate(delt(NI))
	allocate(Re_x(NI))
	allocate(tw_b(NI))
	allocate(ro(NI,NJ))

!----------------- Coordinate of nodes ------------------------
        dx=L/(NI-1)
        dy=H/(NJ-1)

        DO I=1,NI
          DO J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          END DO
        END DO

!----------------- Parameters ------------------------

        NU=MU/Ro0

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'Reh= ', U0*H/NU*2
        write(*,*)'S= ', Su_Max
        print*, 'gamma= ', gama

!----------------- Initial fields -----------------------------

    ro(1,:) = ro0

    p(1,:) = P0

    u(1,:) = u0

    v(1,:) = 0D0


!---------------- Solve Prandtl equations ---------------------



        write(*,*) 'Solve Prandtl equations'
        call  Prandtl(U,V,P,ro,NI,NJ,Su_max,EPS,dx,dy,NU,U0,H, gama, ro0,P0, eps_G, mu,sp_max)
 !----------------- Output data ------------------------------

        write(*,*) 'Operate data'

 	do i=1, Ni
     tw(i) = MU*(U(i,2)-U(i,1))/dy/Ro0/U0*2
!	 Re_x(i) = U0*dx*(i-1)/nu
!	 tw_b(i) = 0.664/sqrt(re_x(i))
	enddo
        write(*,*) 'Output data'
        Open(IO,FILE='Results.plt')
        Call Output_Fields(IO,NI,NJ,X_Node,Y_Node,U,V,P,tw,delt, Re_x, tw_b)
        Close(IO)

        END PROGRAM
