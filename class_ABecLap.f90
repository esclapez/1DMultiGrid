MODULE class_ABecLap

   USE mod_prec_defs
   USE class_LinOp, ONLY : LinOp_t
   USE class_ABecLvl, ONLY : ABecLvl_t

   IMPLICIT NONE

   TYPE, EXTENDS(LinOp_t) :: ABecLap_t
      TYPE(ABecLvl_t), DIMENSION(:), ALLOCATABLE :: lvls    ! MG levels tower
      CONTAINS
          PROCEDURE :: MyType
          PROCEDURE :: define
          PROCEDURE :: destroy
          PROCEDURE :: setBCType
          PROCEDURE :: setBCVals
          PROCEDURE :: setScalars
          PROCEDURE :: setA
          PROCEDURE :: setB
          PROCEDURE :: setC
          PROCEDURE :: applyBC
          PROCEDURE :: apply
          PROCEDURE :: show
          PROCEDURE :: Fsmooth
          PROCEDURE :: initMGLevels
          PROCEDURE :: tridiagSolve
   END TYPE ABecLap_t

   CONTAINS

SUBROUTINE MyType(this)
   CLASS(ABecLap_t) , INTENT(INOUT) :: this
   write(*,*) " I'm ABecLap type "
END SUBROUTINE MyType

SUBROUTINE define(this,nx,dx)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: nx
   REAL(pr), INTENT(IN) :: dx
   REAL(pr) :: tmp_lvl_res
   IF(this%is_allocated)THEN
      WRITE(*,*) "ABecLap already initialized !"
      STOP
   END IF
   tmp_lvl_res = LOG(REAL(nx))/LOG(2.0_pr)
   this%nb_lvls = NINT(tmp_lvl_res)
   IF(2**this%nb_lvls /= nx) THEN
      WRITE(*,*) "ABecLap only works for 2^n grids!"
      STOP
   END IF
   this%nb_lvls = this%nb_lvls
   IF(this%nb_lvls_max/=0) this%nb_lvls = MIN(this%nb_lvls,this%nb_lvls_max)
   ALLOCATE(this%lvls(1:this%nb_lvls))    ! Allocate MG tower
   ALLOCATE(this%nx_lvls(1:this%nb_lvls))    ! Allocatenx/lvl array
   call this%lvls(1)%define(1,nx,dx)      ! Define the finest level      
   this%nx_lvls(1) = nx
   this%is_allocated = .TRUE.
END SUBROUTINE define

SUBROUTINE destroy(this)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER :: mg_lvl
   IF( .NOT. this%is_allocated)THEN
      RETURN
   END IF
   DO mg_lvl = this%nb_lvls, 1, -1   
      CALL this%lvls(mg_lvl)%destroy()
   END DO
   this%nb_lvls = 0
   this%nb_lvls_max = 0
   DEALLOCATE(this%lvls)
   DEALLOCATE(this%nx_lvls)
   this%is_allocated = .FALSE.
END SUBROUTINE destroy

SUBROUTINE setBCType(this,l_type,r_type)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: l_type, r_type
   call this%lvls(1)%leftBC%setType(l_type)
   call this%lvls(1)%rightBC%setType(r_type)
END SUBROUTINE setBCType

SUBROUTINE setBCVals(this,l_val,r_val)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   REAL(pr), INTENT(IN) :: l_val, r_val
   IF(this%lvls(1)%leftBC%isDirichlet()) call this%lvls(1)%leftBC%setVal(l_val)
   IF(this%lvls(1)%rightBC%isDirichlet()) call this%lvls(1)%rightBC%setVal(r_val)
END SUBROUTINE setBCVals

SUBROUTINE setScalars(this,alpha,beta,eta)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   REAL(pr), INTENT(IN) :: alpha, beta
   REAL(pr), OPTIONAL :: eta
   call this%check_alloc("at setScalars")
   IF (PRESENT(eta)) WRITE(*,*) "setScalars: eta not used for ABecLap"
   this%m_alpha = alpha
   this%m_beta = beta
END SUBROUTINE setScalars

SUBROUTINE setA(this,A_in)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   REAL(pr), DIMENSION(:), INTENT(IN) :: A_in
   CALL this%check_alloc("at setA")
   IF ( SIZE(A_in) /= this%lvls(1)%nx ) THEN
      WRITE(*,*) "ABecLap::setA : wrong input vector size "
      STOP
   ENDIF
   CALL this%lvls(1)%setA(A_in)
   this%need_update = .TRUE.
END SUBROUTINE setA

SUBROUTINE setB(this,B_in)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   REAL(pr), DIMENSION(:), INTENT(IN) :: B_in
   CALL this%check_alloc("at setB")
   IF ( SIZE(B_in) /= this%lvls(1)%nx+1 ) THEN
      WRITE(*,*) "ABecLap::setB : wrong input vector size "
      STOP
   ENDIF
   CALL this%lvls(1)%setB(B_in)
   this%need_update = .TRUE.
END SUBROUTINE setB

SUBROUTINE setC(this,C_in)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   REAL(pr), DIMENSION(:), INTENT(IN) :: C_in
   WRITE(*,*) "ABecLap::setC should not be called !"
   STOP
END SUBROUTINE setC

SUBROUTINE applyBC(this,lvl,nx,BCmode,x_in,x_gc)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   INTEGER, INTENT(IN) :: nx
   INTEGER, INTENT(IN) :: BCmode          ! on lvl 1, applyc phys BC only if BCmode == 1
   REAL(pr), DIMENSION(1:this%lvls(lvl)%nx), INTENT(IN) :: x_in
   REAL(pr), DIMENSION(0:this%lvls(lvl)%nx+1), INTENT(OUT) :: x_gc
   INTEGER :: n
   x_gc(1:this%lvls(lvl)%nx) = x_in(1:this%lvls(lvl)%nx)
   IF ( lvl == 1 ) THEN                         ! Physical BC on finest level
      IF (this%lvls(lvl)%leftBC%isNeumann()) THEN
         x_gc(0) = x_gc(1)
      ELSE
         x_gc(0) = - x_gc(1)
         IF ( BCmode == 1 ) x_gc(0) = x_gc(0) + 2.0_pr * this%lvls(lvl)%leftBC%getVal()
      END IF
      IF (this%lvls(lvl)%rightBC%isNeumann()) THEN
         x_gc(this%lvls(lvl)%nx+1) = x_gc(this%lvls(lvl)%nx)
      ELSE
         x_gc(this%lvls(lvl)%nx+1) = - x_gc(this%lvls(lvl)%nx) 
         IF ( BCmode == 1 ) x_gc(this%lvls(lvl)%nx+1) = x_gc(this%lvls(lvl)%nx+1) &
                                                        + 2.0_pr * this%lvls(lvl)%rightBC%getVal()       
      END IF
   ELSE                                         ! Homogeneous BC on coarse levels
      IF (this%lvls(lvl)%leftBC%isNeumann()) THEN
         x_gc(0) = x_gc(1)
      ELSE
         x_gc(0) = - x_gc(1)
      END IF   
      IF (this%lvls(lvl)%rightBC%isNeumann()) THEN
         x_gc(this%lvls(lvl)%nx+1) = x_gc(this%lvls(lvl)%nx)
      ELSE
         x_gc(this%lvls(lvl)%nx+1) = - x_gc(this%lvls(lvl)%nx)
      END IF
   END IF
END SUBROUTINE applyBC

SUBROUTINE apply(this,lvl,BCmode,x_in,y_out)  ! y = ABec(x)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   INTEGER, INTENT(IN) :: BCmode
   REAL(pr), DIMENSION(:), INTENT(IN) :: x_in
   REAL(pr), DIMENSION(:), INTENT(OUT) :: y_out
   REAL(pr), DIMENSION(0:this%lvls(lvl)%nx+1) :: x_gc
   REAL(pr) :: Bdxsq
   INTEGER :: n
   CALL this%check_alloc("at apply")
   IF ( SIZE(x_in) /= this%lvls(lvl)%nx ) THEN
      WRITE(*,*) "ABecLap::apply : wrong input vector size "
      STOP
   ENDIF
   IF ( SIZE(y_out) /= this%lvls(lvl)%nx ) THEN
      WRITE(*,*) "ABecLap::apply : wrong output vector size "
      STOP
   ENDIF
   CALL this%applyBC(lvl, this%lvls(lvl)%nx, BCmode, x_in, x_gc)
   Bdxsq = this%m_beta / ( this%lvls(lvl)%dx**2.0_pr )
   DO n = 1, this%lvls(lvl)%nx
      y_out(n) =  this%m_alpha * this%lvls(lvl)%m_a(n) * x_gc(n) &
                - Bdxsq * (   this%lvls(lvl)%m_b(n+1) * ( x_gc(n+1) - x_gc(n)  ) & 
                            - this%lvls(lvl)%m_b(n)   * ( x_gc(n)   - x_gc(n-1)) )
   END DO
END SUBROUTINE apply

SUBROUTINE show(this,lvl)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   REAL(pr) :: Bdxsq, diagcoef, leftcoef, rightcoef
   INTEGER :: n
   CALL this%check_alloc("at show")
   Bdxsq = this%m_beta / ( this%lvls(lvl)%dx**2.0_pr )
   DO n = 1, this%lvls(lvl)%nx
      leftcoef = 0.0_pr
      rightcoef = 0.0_pr
      IF ( n/= 1 ) THEN
         leftcoef = - Bdxsq * this%lvls(lvl)%m_b(n)
         WRITE(*,*) " ( Row: ", n, "; Col: ", n-1, " = ", leftcoef, " )"
      END IF
      diagcoef = this%m_alpha * this%lvls(lvl)%m_a(n) + Bdxsq * ( this%lvls(lvl)%m_b(n+1) + this%lvls(lvl)%m_b(n) )
      WRITE(*,*) " ( Row: ", n, "; Col: ", n,   " = ", diagcoef, " )"
      IF ( n/=this%lvls(lvl)%nx ) THEN
         rightcoef = - Bdxsq * this%lvls(lvl)%m_b(n+1)
         WRITE(*,*) " ( Row: ", n, "; Col: ", n+1, " = ", rightcoef, " )"
      END IF
      WRITE(*,*) " ( Row ", n, " diagonal dom  = |diag| - |offdiags|: ", ABS(diagcoef) - ABS(leftcoef+rightcoef)
   END DO
END SUBROUTINE show

SUBROUTINE Fsmooth(this,lvl,nx,redblack,sol,rhs)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl
   INTEGER, INTENT(IN) :: nx
   INTEGER, INTENT(IN) :: redblack
   REAL(pr), DIMENSION(0:nx+1), INTENT(INOUT) :: sol
   REAL(pr), DIMENSION(1:nx), INTENT(IN) :: rhs
   REAL(pr) :: diag, minusoffdiag, Bdxsq, delta
   REAL(pr) :: cf0, cf1
   INTEGER :: n
   Bdxsq = this%m_beta / ( this%lvls(lvl)%dx**2.0_pr )
   DO n = 1, this%lvls(lvl)%nx
      IF ( mod(n+redblack,2) == 0 ) then
       cf0 = 0.0_pr
       cf1 = 0.0_pr
       IF ( n == 1 ) cf0 = -1.0_pr
       IF ( n == this%lvls(lvl)%nx ) cf1 = -1.0_pr
       delta = Bdxsq * ( this%lvls(lvl)%m_b(n) * cf0 + this%lvls(lvl)%m_b(n+1) * cf1 ) 
!       delta = 0.0_pr
         diag =   this%m_alpha * this%lvls(lvl)%m_a(n) &
                + Bdxsq * ( this%lvls(lvl)%m_b(n) + this%lvls(lvl)%m_b(n+1) )
         minusoffdiag = Bdxsq * (  this%lvls(lvl)%m_b(n+1) * sol(n+1) &
                                 + this%lvls(lvl)%m_b(n)   * sol(n-1) )    
         sol(n) = ( rhs(n) + minusoffdiag - delta * sol(n) ) / ( diag - delta )
      END IF
   END DO
END SUBROUTINE Fsmooth

SUBROUTINE initMGLevels(this)
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER :: mg_lvl
   call this%check_alloc("at initMGlevels")
   DO mg_lvl = 2, this%nb_lvls
      call this%lvls(mg_lvl)%define(mg_lvl,this%lvls(mg_lvl-1)%nx/2,this%lvls(mg_lvl-1)%dx*2.0_pr)
      this%nx_lvls(mg_lvl) = this%lvls(mg_lvl)%nx
      call this%lvls(mg_lvl)%leftBC%setType(1)
      call this%lvls(mg_lvl)%rightBC%setType(1)
      call this%lvls(mg_lvl)%leftBC%setVal(0.0_pr)
      call this%lvls(mg_lvl)%rightBC%setVal(0.0_pr)
      call averageDownMGCell(this%lvls(mg_lvl)%nx, this%lvls(mg_lvl-1)%m_a(:),this%lvls(mg_lvl)%m_a(:))
      call averageDownMGFace(this%lvls(mg_lvl)%nx, this%lvls(mg_lvl-1)%m_b(:),this%lvls(mg_lvl)%m_b(:))
   END DO
   this%need_update = .FALSE.
END SUBROUTINE initMGLevels

FUNCTION tridiagSolve(this,nx_io,rhs)  RESULT(x_inout)! x_inout = ABec^-1(rhs)
   USE mod_utils, ONLY : tridiag
   CLASS(ABecLap_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: nx_io
   REAL(pr), DIMENSION(:), INTENT(IN) :: rhs
   REAL(pr), DIMENSION(1:nx_io) :: x_inout
   REAL(pr), DIMENSION(0:this%lvls(1)%nx+1) :: x_gc
   REAL(pr), DIMENSION(1:this%lvls(1)%nx) :: a, b, c, r, u, gam
   REAL(pr) :: Bdxsq
   INTEGER :: n, n_solve
   call this%check_alloc("at tridiagSolve")
   IF ( SIZE(x_inout) /= this%lvls(1)%nx ) THEN
      WRITE(*,*) "ABecLap::tridiagSolve : wrong input vector size "
      STOP
   ENDIF
   IF ( SIZE(rhs) /= this%lvls(1)%nx ) THEN
      WRITE(*,*) "ABecLap::tridiagSolve : wrong rhs vector size "
      STOP
   ENDIF
   x_inout(:) = 0.0_pr
   call this%applyBC(1, this%lvls(1)%nx, 1, x_inout, x_gc)
   Bdxsq = this%m_beta / ( this%lvls(1)%dx**2.0_pr )
   DO n = 1, this%lvls(1)%nx
      u(n) = 0.0_pr
      a(n) = - Bdxsq * this%lvls(1)%m_b(n)
      c(n) = - Bdxsq * this%lvls(1)%m_b(n+1)
      b(n) = this%m_alpha * this%lvls(1)%m_a(n) - ( a(n) + c(n) )
      r(n) = rhs(n)
      IF ( n == 1 ) THEN
         b(n) = b(n) - a(n)
         r(n) = r(n) - a(n) * x_gc(n-1)
         a(n) = 0.0_pr
      ELSE IF ( n == this%lvls(1)%nx ) THEN
         b(n) = b(n) - c(n)
         r(n) = r(n) - c(n) * x_gc(n+1)
         c(n) = 0.0_pr
      END IF   
   END DO
   n_solve = this%lvls(1)%nx
   call tridiag(a,b,c,r,u,gam,n_solve)
   x_inout(:) = u(1:this%lvls(1)%nx)
END FUNCTION tridiagSolve

SUBROUTINE averageDownMGCell(nx_c, fine, coarse)
   INTEGER, INTENT(IN) :: nx_c
   REAL(pr), DIMENSION(1:*), INTENT(IN) :: fine
   REAL(pr), DIMENSION(1:nx_c), INTENT(OUT) :: coarse
   INTEGER :: n
   DO n = 1, nx_c
      coarse(n) = 0.5_pr * (fine(2*n-1) + fine(2*n))
   END DO
END SUBROUTINE averageDownMGCell

SUBROUTINE averageDownMGFace(nx_c, fine, coarse)
   INTEGER, INTENT(IN) :: nx_c
   REAL(pr), DIMENSION(1:*), INTENT(IN) :: fine
   REAL(pr), DIMENSION(1:nx_c+1), INTENT(OUT) :: coarse
   INTEGER :: n
   DO n = 1, nx_c+1
      coarse(n) = fine(2*n-1)
   END DO
END SUBROUTINE averageDownMGFace

END MODULE class_ABecLap
