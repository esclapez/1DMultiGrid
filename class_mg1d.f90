MODULE class_mg1d

   USE mod_prec_defs
   USE class_LinOp, ONLY: LinOp_t

   IMPLICIT NONE

   TYPE mgLvlData_t
      INTEGER :: nx
      REAL(pr), DIMENSION(:), ALLOCATABLE :: sol
      REAL(pr), DIMENSION(:), ALLOCATABLE :: res
      REAL(pr), DIMENSION(:), ALLOCATABLE :: rhs
      REAL(pr), DIMENSION(:), ALLOCATABLE :: cor
      REAL(pr), DIMENSION(:), ALLOCATABLE :: rescor
      CONTAINS
         PROCEDURE :: define => mgLvlData_init
         PROCEDURE :: destroy => mgLvlData_clear
   END TYPE mgLvlData_t

   TYPE mg1d_t
      INTEGER :: nlvl = -1                ! Number of MG levels (from LinOp)
      INTEGER :: fixed_nvcycle = -1
      INTEGER :: max_nvcycle = 200
      INTEGER :: n_smooth_up = 2
      INTEGER :: n_smooth_down = 2
      INTEGER :: n_smooth_bottom = 8
      INTEGER :: verbose = 0
      INTEGER :: finalIter = 0
      REAL(pr) :: abs_tol, rel_tol
      TYPE(mgLvlData_t), DIMENSION(:), ALLOCATABLE :: lvls
      CLASS(LinOp_t), POINTER :: LinOp
      CONTAINS
         PROCEDURE :: define => mg_init
         PROCEDURE :: destroy => mg_clear
         PROCEDURE :: setPreSmooth => mg_setPreSmooth
         PROCEDURE :: setPostSmooth => mg_setPostSmooth
         PROCEDURE :: setBottomSmooth => mg_setBottomSmooth
         PROCEDURE :: setFixedIter => mg_setFixedIter
         PROCEDURE :: setMaxIter => mg_setMaxIter
         PROCEDURE :: setVerbose => mg_setVerbose
         PROCEDURE :: solve => mg_solve
         PROCEDURE :: solutionResidual => mg_solutionResidual
         PROCEDURE :: correctionResidual => mg_correctionResidual
         PROCEDURE :: oneIter => mg_oneIter
         PROCEDURE :: getMGIterCount => mg_getIterCount
   END TYPE mg1d_t

   CONTAINS

SUBROUTINE mgLvlData_init(this, nn)
   CLASS(mgLvlData_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: nn
   this%nx = nn
   ALLOCATE(this%sol(1:nn))
   ALLOCATE(this%res(1:nn))
   ALLOCATE(this%rhs(1:nn))
   ALLOCATE(this%cor(1:nn))
   ALLOCATE(this%rescor(1:nn))
END SUBROUTINE mgLvlData_init

SUBROUTINE mgLvlData_clear(this)
   CLASS(mgLvlData_t), INTENT(INOUT) :: this
   this%nx = 0
   DEALLOCATE(this%sol)
   DEALLOCATE(this%res)
   DEALLOCATE(this%rhs)
   DEALLOCATE(this%cor)
   DEALLOCATE(this%rescor)
END SUBROUTINE mgLvlData_clear

SUBROUTINE mg_init(this, LinOp)      
   CLASS(mg1d_t), INTENT(INOUT) :: this
   CLASS(LinOp_t), INTENT(IN), TARGET :: LinOp
   INTEGER :: mg_lvl
   this%LinOp => LinOp
   this%nlvl = LinOp%nb_lvls
   ALLOCATE(this%lvls(1:this%nlvl))
   CALL this%lvls(1)%define(this%LinOp%nx_lvls(1))
   DO mg_lvl = 2, this%nlvl
      CALL this%lvls(mg_lvl)%define(this%lvls(mg_lvl-1)%nx/2)
   END DO
END SUBROUTINE mg_init

SUBROUTINE mg_clear(this)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER :: mg_lvl
   IF(this%nlvl==-1) RETURN
   this%nlvl = -1
   this%LinOp => NULL()
   DO mg_lvl = 1, this%nlvl
      CALL this%lvls(mg_lvl)%destroy()
   END DO
   DEALLOCATE(this%lvls)
END SUBROUTINE mg_clear

SUBROUTINE mg_setPreSmooth(this, n_smooth)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: n_smooth
   this%n_smooth_down = n_smooth
END SUBROUTINE mg_setPreSmooth

SUBROUTINE mg_setPostSmooth(this, n_smooth)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: n_smooth
   this%n_smooth_up = n_smooth
END SUBROUTINE mg_setPostSmooth

SUBROUTINE mg_setBottomSmooth(this, n_smooth)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: n_smooth
   this%n_smooth_bottom = n_smooth
END SUBROUTINE mg_setBottomSmooth

SUBROUTINE mg_setMaxIter(this, niter)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: niter
   this%max_nvcycle = niter
END SUBROUTINE mg_setMaxIter

SUBROUTINE mg_setFixedIter(this, niter)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: niter
   this%fixed_nvcycle = niter
END SUBROUTINE mg_setFixedIter

SUBROUTINE mg_setVerbose(this, verbose)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: verbose
   this%verbose = verbose
END SUBROUTINE mg_setVerbose

SUBROUTINE mg_solve(this, nn, x, rhs, r_tol, a_tol)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: nn
   REAL(pr), DIMENSION(1:nn), INTENT(INOUT) :: x            ! Overwritten with the solution
   REAL(pr), DIMENSION(1:nn), INTENT(IN) :: rhs
   REAL(pr), INTENT(IN) :: r_tol, a_tol
   REAL(pr) :: resnorm0, rhsnorm0, norminf
   REAL(pr) :: res_target
   INTEGER :: iter, niters
   LOGICAL :: converged
   IF (nn/=this%LinOp%nx_lvls(1)) THEN
      WRITE(*,*) "mg1d::solve : size of input vector incompatible with LinOp"  
      STOP
   END IF
   this%lvls(1)%sol(:) = x(:)
   this%lvls(1)%rhs(:) = rhs(:)
   IF (this%LinOp%need_update) THEN
      CALL this%LinOp%initMGLevels()
   END IF
   CALL this%solutionResidual(1)
   resnorm0 = MAXVAL(ABS(this%lvls(1)%res(:)))
   rhsnorm0 = MAXVAL(ABS(this%lvls(1)%rhs(:)))
   IF (this%verbose>=1) THEN
      WRITE(*,*) "mg1d: Initial rhs               = ", rhsnorm0
      WRITE(*,*) "mg1d: Initial residual (resid0) = ", resnorm0
   END IF
   res_target = MAX(a_tol, MAX(r_tol,1.d-16)*resnorm0);
   IF ( resnorm0 <= res_target ) THEN
      IF (this%verbose>=1) THEN
         WRITE(*,*) "mg1d: No iterations needed"
      END IF
      RETURN
   ELSE
      niters = this%max_nvcycle
      IF ( this%fixed_nvcycle>0 ) niters = this%fixed_nvcycle
      DO iter = 1, niters
         CALL this%oneIter(iter)
         converged = .FALSE.
         CALL this%solutionResidual(1)
         norminf = MAXVAL(ABS(this%lvls(1)%res(:)))
!         write(*,*) this%lvls(1)%res(:)
!         write(*,*) this%lvls(1)%sol(:)
         IF (norminf <= res_target) converged = .TRUE.
         IF ( this%fixed_nvcycle<0 ) THEN
            IF (converged) THEN
               IF (this%verbose>=1) THEN
                  WRITE(*,*) "mg1d: Final Iter. ", iter, " resid, resid/resid0 ", &
                              norminf, " , ", norminf/resnorm0
               END IF
               this%finalIter = iter
               EXIT
            ELSE
               IF (this%verbose>=2) THEN
                  WRITE(*,*) "mg1d: Iter. ", iter, " resid, resid/resid0 ", &
                              norminf, " , ", norminf/resnorm0
               END IF
            END IF
         ELSE
            IF (this%verbose>=2 .AND. iter /= niters) THEN
               WRITE(*,*) "mg1d: Iter. ", iter, " resid, resid/resid0 ", &
                           norminf, " , ", norminf/resnorm0
            END IF
            IF (this%verbose>=1 .AND. iter == niters) THEN
               WRITE(*,*) "mg1d: Final Iter. ", iter, " resid, resid/resid0 ", &
                           norminf, " , ", norminf/resnorm0
            END IF
            this%finalIter = iter
         END IF
      END DO
      IF (.NOT. converged .AND. this%fixed_nvcycle<0 ) THEN
         WRITE(*,*) "mg1d: Failed to converge after ", niters, " iterations. resid, resid/resid0 ",&
                     norminf, " , ", norminf/resnorm0
         STOP
      END IF
      x(:) = this%lvls(1)%sol(:)
   END IF
END SUBROUTINE mg_solve

SUBROUTINE mg_solutionResidual(this,lvl)    ! res = rhs - LinOp(sol)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   CALL this%LinOp%apply( lvl, &
                          1, &
                          this%lvls(lvl)%sol(:), & 
                          this%lvls(lvl)%res(:) )
   this%lvls(lvl)%res(:) = this%lvls(lvl)%rhs(:) - this%lvls(lvl)%res(:)              
END SUBROUTINE mg_solutionResidual

SUBROUTINE mg_correctionResidual(this,lvl,BCmode) ! rescor = res - LinOp(cor)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   INTEGER, INTENT(IN) :: BCmode
   CALL this%LinOp%apply( lvl, &
                          BCmode, & 
                          this%lvls(lvl)%cor(:), & 
                          this%lvls(lvl)%rescor(:) )
   this%lvls(lvl)%rescor(:) = this%lvls(lvl)%res(:) - this%lvls(lvl)%rescor(:)              
END SUBROUTINE mg_correctionResidual

SUBROUTINE mg_oneIter(this, iter)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: iter
   REAL(pr) :: norminf
   INTEGER :: mg_lvl, i
   ! Down
   DO mg_lvl = 1, this%nlvl-1
      this%lvls(mg_lvl)%cor(:) = 0.0_pr
      IF (this%verbose>=3) THEN
         norminf = MAXVAL(ABS(this%lvls(mg_lvl)%res(:)))
         WRITE(*,*) "At level ", mg_lvl, " Down: Before smooth ", norminf
      END IF
      DO i = 1, this%n_smooth_down
         CALL this%LinOp%smooth(mg_lvl, this%lvls(mg_lvl)%cor(:),& 
                                        this%lvls(mg_lvl)%res(:) )
      END DO
      CALL this%correctionResidual(mg_lvl,0)
      IF (this%verbose>=3) THEN
         norminf = MAXVAL(ABS(this%lvls(mg_lvl)%rescor(:)))
         WRITE(*,*) "At level ", mg_lvl, " Down: After smooth ", norminf
      END IF
      CALL this%LinOp%restriction(mg_lvl+1, this%lvls(mg_lvl+1)%res(:), &
                                            this%lvls(mg_lvl)%rescor(:))
   END DO
   ! Bottom
   this%lvls(this%nlvl)%cor(:) = 0.0_pr
   IF (this%verbose>=3) THEN
      norminf = MAXVAL(ABS(this%lvls(this%nlvl)%res(:)))
      WRITE(*,*) "Bottom level Down", norminf
   END IF
   DO i = 1, this%n_smooth_bottom
      CALL this%LinOp%smooth(this%nlvl, this%lvls(this%nlvl)%cor(:),& 
                                        this%lvls(this%nlvl)%res(:) )
   END DO
   IF (this%verbose>=3) THEN
      CALL this%correctionResidual(this%nlvl,0)
      norminf = MAXVAL(ABS(this%lvls(this%nlvl)%rescor(:)))
      WRITE(*,*) "Bottom level Up ", norminf
   END IF
   ! Up
   DO mg_lvl = this%nlvl-1, 1, -1
      CALL this%LinOp%prolongation(mg_lvl+1, this%lvls(mg_lvl)%cor(:), & 
                                             this%lvls(mg_lvl+1)%cor(:))
      CALL this%correctionResidual(mg_lvl,0)
      IF (this%verbose>=3) THEN
         norminf = MAXVAL(ABS(this%lvls(mg_lvl)%rescor(:)))
         WRITE(*,*) "At level ", mg_lvl, " Up: Before smooth ", norminf
      END IF
      DO i = 1, this%n_smooth_up
         CALL this%LinOp%smooth(mg_lvl, this%lvls(mg_lvl)%cor(:),&
                                        this%lvls(mg_lvl)%res(:) )
      END DO
      CALL this%correctionResidual(mg_lvl,0)
      IF (this%verbose>=3) THEN
         norminf = MAXVAL(ABS(this%lvls(mg_lvl)%rescor(:)))
         WRITE(*,*) "At level ", mg_lvl, " Up: After smooth ", norminf
      END IF
   END DO
   ! Update solution on fine
   this%lvls(1)%sol(:) = this%lvls(1)%sol(:) + this%lvls(1)%cor(:)
END SUBROUTINE mg_oneIter

FUNCTION mg_getIterCount(this) RESULT(itercount)
   CLASS(mg1d_t), INTENT(INOUT) :: this
   INTEGER :: itercount
   itercount = this%finalIter
END FUNCTION mg_getIterCount

END MODULE class_mg1d
