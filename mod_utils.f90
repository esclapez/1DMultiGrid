!###################################################################################
!#       file: mod_utils.f90
!#       desc: Module contaning several generic fonctionality 
!#    author : various
!#      date : 09/2018
!###################################################################################
MODULE mod_utils

   USE mod_prec_defs

   IMPLICIT NONE

   CONTAINS

!###################################################################################
!# subroutine: tridiag
!#       desc: Tridiagonal solve used for implicit velocity, tempe, spec solves.
!#    author : various
!#      date : 09/2018
!###################################################################################
SUBROUTINE tridiag(a,b,c,r,u,gam,n)

   IMPLICIT NONE
   
   INTEGER :: n
   REAL(pr), DIMENSION(n) :: a, b, c, r, u, gam
   INTEGER :: j
   REAL(pr) :: bet 
   
   IF (b(1) == 0) print *,'CANT HAVE B(1) = ZERO'
   
   bet = b(1)
   u(1) = r(1)/bet
   
   DO j = 2, n
      gam(j) = c(j-1) / bet
      bet = b(j) - a(j)*gam(j)
      IF (bet == 0) THEN
        print *,'TRIDIAG FAILED '
        stop
      END IF
      u(j) = (r(j)-a(j)*u(j-1))/bet
   END DO
   
   DO j = n-1, 1, -1
     u(j) = u(j) - gam(j+1)*u(j+1)
   END DO

END SUBROUTINE tridiag

END MODULE mod_utils
