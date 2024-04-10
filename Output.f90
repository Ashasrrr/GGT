!************************************************************************************************
       SUBROUTINE Output_Fields(IO,NI,NJ,X,Y,U,V,P,cf,delt,Re_x,cf_b)
         IMPLICIT NONE

         INTEGER NI,NJ,IO, I, j
         REAL*8,DIMENSION(NI,NJ):: X,Y
         REAL*8,DIMENSION(NI,NJ):: U,V,P
	 real*8 cf(Ni),delt(ni), Re_X(ni), cf_b(NI)


         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "dP"'
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ)
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI,1:NJ)
         Write(IO,'(100E25.16)') V(1:NI,1:NJ)
         Write(IO,'(100E25.16)') 1e5 - P(1:NI,1:NJ)
         CLOSE(1)
open(2,FILE="2.plt")
write(2,*) 'VARIABLES = "Re_x", "Cf", "delt", "cf_b"'
write(2,*)  'ZONE I=',NI-1
do i=2, Ni
write(2,*) x(i,1), Cf(i), delt(i), cf_b(i)
enddo

       END  SUBROUTINE

!************************************************************************************************
