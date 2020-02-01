PROGRAM mgDriver

   USE mod_prec_defs
   USE class_LinOp, ONLY : LinOp_t
   USE class_ABecLap, ONLY : ABecLap_t
   USE class_ABecCecLap, ONLY : ABecCecLap_t
   USE class_mg1d, ONLY : mg1d_t

   IMPLICIT NONE

   INTEGER, PARAMETER :: npts = 32
   REAL(pr), PARAMETER :: lx = 0.01_pr
   REAL(pr) :: dx
   REAL(pr) :: leftBCval = 0.0_pr
   REAL(pr) :: rightBCval = 5.0_pr
   REAL(pr) :: alpha
   REAL(pr) :: beta
   REAL(pr) :: eta
   REAL(pr), DIMENSION(1:npts) :: Avec
   REAL(pr), DIMENSION(1:npts+1) :: Bvec
   REAL(pr), DIMENSION(1:npts+1) :: Cvec
   REAL(pr), DIMENSION(1:npts) :: init
   REAL(pr), DIMENSION(1:npts) :: sol
   REAL(pr), DIMENSION(1:npts) :: rhs

   REAL(pr) :: r_tol, a_tol

   INTEGER :: n

   TYPE(ABecLap_t), TARGET :: ABec
   TYPE(ABecCecLap_t), TARGET :: ABecCec
   CLASS(LinOp_t), POINTER :: GenericLinOpPtr
   TYPE(mg1d_t) :: mg1d

   WRITE(*,*) " Testing Fortran OO MG class " 

   dx = lx / npts

   init(:) = 0.0_pr
   sol(:) = 0.0_pr
   rhs(:) = 1.0_pr
   alpha = 0.0_pr
   beta = -10.0_pr
   eta = 0.0_pr
   Avec(:) = 0.0_pr
   Bvec(:) = 1.0_pr
   Cvec(:) = 0.0_pr

!  Set ABec Linear operator
   CALL ABec%define(npts,dx)
   CALL ABec%setScalars(alpha,beta)
   CALL ABec%setA(Avec(:))
   CALL ABec%setB(Bvec(:))
   CALL ABec%setBCType(1,1)
   CALL ABec%setBCVals(leftBCval, rightBCval)

!  Since we are in 1D, we can perform of tridiag solve to get the exact solution
   sol(:) = ABec%tridiagSolve(npts, rhs(:))
   DO n = 1, npts
      WRITE(*,*) sol(n)
   END DO
   
!  Initialize the multigrid solver: pass a pointer to the linear operator
!  and solve to a certain tolerance
   sol(:) = 0.0_pr
   CALL mg1d%define(ABec)
   CALL mg1d%setVerbose(1)
   r_tol = 1.0e-14
   a_tol = 1.0e-14
   CALL mg1d%solve(npts, sol(:), rhs(:), r_tol, a_tol)
   DO n = 1, npts
      WRITE(*,*) sol(n)
   END DO

!  Update the linear operator and solve again
   sol(:) = 0.0_pr
   sol(:) = 0.0_pr
   Bvec(17:33) = 100.0_pr
   CALL ABec%setB(Bvec(:))
   CALL mg1d%solve(npts, sol(:), rhs(:), r_tol, a_tol)
   DO n = 1, npts
      WRITE(*,*) sol(n)
   END DO

END PROGRAM mgDriver
