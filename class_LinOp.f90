! Definition of an abstract linear operator class
! Some functions are deferred to the extended derived class
MODULE class_LinOp

   USE mod_prec_defs

   IMPLICIT NONE

!  Linear operator
   TYPE, ABSTRACT :: LinOp_t
      INTEGER  :: nb_lvls                                   ! Effective number of levels
      INTEGER  :: nb_lvls_max = 0                           ! User defined max number of levels
      INTEGER, DIMENSION(:), ALLOCATABLE :: nx_lvls         ! Grid size per level. Will help avoid code duplication.
      REAL(pr) :: m_alpha, m_beta, m_eta                    ! Linear operator scalar values
      LOGICAL :: need_update = .TRUE.                       ! Need to update the MG levels coeffs ?
      LOGICAL :: is_allocated = .FALSE.                     ! LinOp is allocated ?
      CONTAINS
          PROCEDURE( LinOp_printType ), DEFERRED         :: MyType
          PROCEDURE( LinOp_alloc ), DEFERRED, PASS       :: define
          PROCEDURE( LinOp_dealloc ), DEFERRED           :: destroy
          PROCEDURE( LinOp_setScalars ), DEFERRED, PASS  :: setScalars
          PROCEDURE( LinOp_setA ), DEFERRED, PASS        :: setA
          PROCEDURE( LinOp_setB ), DEFERRED, PASS        :: setB
          PROCEDURE( LinOp_setC ), DEFERRED, PASS        :: setC
          PROCEDURE( LinOp_setBCType ), DEFERRED, PASS   :: setBCType
          PROCEDURE( LinOp_setBCVals ), DEFERRED, PASS   :: setBCVals
          PROCEDURE( LinOp_show ), DEFERRED,PASS         :: show
          PROCEDURE :: setMaxMGLevels => LinOp_setMaxMGLevels
          PROCEDURE( LinOp_apply ), DEFERRED, PASS       :: apply
          PROCEDURE( LinOp_tridiagSolve ), DEFERRED, PASS:: tridiagSolve
          PROCEDURE :: smooth => LinOp_smooth
          PROCEDURE( LinOp_Fsmooth ), DEFERRED, PASS     :: Fsmooth
          PROCEDURE :: restriction => LinOp_restriction
          PROCEDURE :: prolongation => LinOp_prolongation
          PROCEDURE( LinOp_initMGLevels ), DEFERRED      :: initMGLevels
          PROCEDURE :: check_alloc => LinOp_checkAlloc
          PROCEDURE( LinOp_fillBndy ), DEFERRED, PASS    :: applyBC
!         PROCEDURE :: ABecLap_copy
!         GENERIC :: ASSIGNMENT(=) => ABecLap_copy 
   END TYPE LinOp_t

   ABSTRACT INTERFACE
      SUBROUTINE LinOp_printType(this)
         IMPORT LinOp_t
         CLASS(LinOp_t), INTENT(INOUT) :: this
      END SUBROUTINE LinOp_printType
      SUBROUTINE LinOp_alloc(this,nx,dx)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: nx
         REAL(pr), INTENT(IN) :: dx
      END SUBROUTINE LinOp_alloc
      SUBROUTINE LinOp_dealloc(this)
         IMPORT LinOp_t
         CLASS(LinOp_t), INTENT(INOUT) :: this
      END SUBROUTINE LinOp_dealloc
      SUBROUTINE LinOp_setBCType(this,l_type,r_type)
         IMPORT LinOp_t
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: l_type, r_type
      END SUBROUTINE LinOp_setBCType
      SUBROUTINE LinOp_setBCVals(this,l_val,r_val)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         REAL(pr), INTENT(IN) :: l_val, r_val
      END SUBROUTINE LinOp_setBCVals
      SUBROUTINE LinOp_setScalars(this,alpha,beta,eta)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         REAL(pr), INTENT(IN) :: alpha, beta
         REAL(pr), OPTIONAL :: eta
      END SUBROUTINE LinOp_setScalars
      SUBROUTINE LinOp_setA(this,A_in)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         REAL(pr), DIMENSION(:), INTENT(IN) :: A_in
      END SUBROUTINE LinOp_setA
      SUBROUTINE LinOp_setB(this,B_in)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         REAL(pr), DIMENSION(:), INTENT(IN) :: B_in
      END SUBROUTINE LinOp_setB
      SUBROUTINE LinOp_setC(this,C_in)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         REAL(pr), DIMENSION(:), INTENT(IN) :: C_in
      END SUBROUTINE LinOp_setC
      SUBROUTINE LinOp_fillBndy(this,lvl,nx,BCmode,x_in,x_gc)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: lvl
         INTEGER, INTENT(IN) :: nx
         INTEGER, INTENT(IN) :: BCmode          ! on lvl 1, applyc phys BC only if BCmode == 1
         REAL(pr), DIMENSION(1:nx), INTENT(IN) :: x_in
         REAL(pr), DIMENSION(0:nx+1), INTENT(OUT) :: x_gc
      END SUBROUTINE LinOp_fillBndy
      SUBROUTINE LinOp_show(this,lvl)
         IMPORT LinOp_t
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: lvl
      END SUBROUTINE LinOp_show
      SUBROUTINE LinOp_apply(this,lvl,BCmode,x_in,y_out)  ! y = ABec(x)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: lvl
         INTEGER, INTENT(IN) :: BCmode
         REAL(pr), DIMENSION(:), INTENT(IN) :: x_in
         REAL(pr), DIMENSION(:), INTENT(OUT) :: y_out
      END SUBROUTINE LinOp_apply
      SUBROUTINE LinOp_Fsmooth(this,lvl,nx, redblack,sol,rhs)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: lvl
         INTEGER, INTENT(IN) :: nx
         INTEGER, INTENT(IN) :: redblack
         REAL(pr), DIMENSION(0:nx+1), INTENT(INOUT) :: sol
         REAL(pr), DIMENSION(1:nx), INTENT(IN) :: rhs
      END SUBROUTINE LinOp_Fsmooth
      SUBROUTINE LinOp_initMGLevels(this)
         IMPORT LinOp_t
         CLASS(LinOp_t), INTENT(INOUT) :: this
      END SUBROUTINE LinOp_initMGLevels
      FUNCTION LinOp_tridiagSolve(this,nx_io,rhs) RESULT(x_inout) ! x_inout = ABec^-1(rhs)
         IMPORT LinOp_t
         IMPORT pr
         CLASS(LinOp_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: nx_io
         REAL(pr), DIMENSION(:), INTENT(IN) :: rhs
         REAL(pr), DIMENSION(1:nx_io) :: x_inout
      END FUNCTION LinOp_tridiagSolve
   END INTERFACE

   CONTAINS
      
SUBROUTINE LinOp_setMaxMGLevels(this,mg_max)
   CLASS(LinOp_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: mg_max
   this%nb_lvls_max = mg_max
END SUBROUTINE LinOp_setMaxMGLevels

SUBROUTINE LinOp_checkAlloc(this,err_string)
   CLASS(LinOp_t), INTENT(IN) :: this
   CHARACTER(LEN=*), INTENT(IN) :: err_string
   IF(.NOT. this%is_allocated)THEN
      WRITE(*,*) "LinOp not allocated "//TRIM(err_string)
      STOP
   END IF
END SUBROUTINE LinOp_checkAlloc

SUBROUTINE LinOp_smooth(this,lvl,sol,rhs)
   CLASS(LinOp_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   REAL(pr), DIMENSION(1:this%nx_lvls(lvl)), INTENT(INOUT) :: sol
   REAL(pr), DIMENSION(1:this%nx_lvls(lvl)), INTENT(IN) :: rhs
   REAL(pr), DIMENSION(0:this%nx_lvls(lvl)+1) :: sol_gc
   INTEGER :: rb, BCmode
   CALL this%check_alloc("at smooth")
   BCmode = 0              ! Use phys bc vals for smooth. BC val set to 0 for mg_lvl > 1.
   DO rb = 1, 2 
      CALL this%applyBC(lvl, this%nx_lvls(lvl), BCmode, sol, sol_gc)
      CALL this%Fsmooth(lvl, this%nx_lvls(lvl), rb, sol_gc, rhs)
      sol(1:this%nx_lvls(lvl)) = sol_gc(1:this%nx_lvls(lvl))
   END DO
END SUBROUTINE LinOp_smooth

SUBROUTINE LinOp_restriction(this, lvl_c, coarse, fine)
   CLASS(LinOp_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl_c
   REAL(pr), DIMENSION(1:this%nx_lvls(lvl_c)), INTENT(OUT) :: coarse
   REAL(pr), DIMENSION(1:this%nx_lvls(lvl_c-1)), INTENT(IN) :: fine
   INTEGER :: n
   DO n = 1, this%nx_lvls(lvl_c)
      coarse(n) = 0.5_pr * (fine(2*n-1) + fine(2*n))
   END DO
END SUBROUTINE LinOp_restriction

SUBROUTINE LinOp_prolongation(this, lvl_c, fine, coarse)
   CLASS(LinOp_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl_c
   REAL(pr), DIMENSION(1:this%nx_lvls(lvl_c-1)), INTENT(INOUT) :: fine
   REAL(pr), DIMENSION(1:this%nx_lvls(lvl_c)), INTENT(IN) :: coarse
   INTEGER :: n
   DO n = 1, this%nx_lvls(lvl_c)
      fine(2*n) = fine(2*n) + coarse(n)
      fine(2*n-1) = fine(2*n-1) + coarse(n)
   END DO
END SUBROUTINE LinOp_prolongation

END MODULE class_LinOp
