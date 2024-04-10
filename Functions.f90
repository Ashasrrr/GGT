!*********************************************************************************************
       SUBROUTINE Solve_Tridiagonal_Matrix(NX, ka, kb, kc, kd, P)
       IMPLICIT NONE

       integer :: i, NX
       real(8) :: P(NX), ka(NX), kb(NX), kc(NX), kd(NX), alf(NX), bet(NX)

       alf(2) = - kc(1)/kb(1)
       bet(2) = kd(1)/kb(1)

       do i = 2, NX - 1
         alf(i+1) = - kc(i)/(kb(i) + ka(i) * alf(i))
         bet(i+1) = (kd(i) - ka(i) * bet(i)) / (kb(i) + ka(i) * alf(i))
       end do
       P(NX) = (kd(NX) - ka(NX)*bet(NX))/(kb(NX) + ka(NX)*alf(NX))

       do i = NX - 1, 1, -1
         P(i) = alf(i+1) * P(i + 1) + bet(i+1)
       end do

       END SUBROUTINE
