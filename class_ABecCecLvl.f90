! Definition of a child ABecCec lvl class from LinOpLvl abstract
MODULE class_ABecCecLvl

   USE mod_prec_defs
   USE class_LinOpLvl, ONLY : LinOpLvl_t

   IMPLICIT NONE

   TYPE, EXTENDS(LinOpLvl_t) :: ABecCecLvl_t
      CONTAINS
          PROCEDURE :: MyType
          PROCEDURE :: define
          PROCEDURE :: destroy
          PROCEDURE :: setC
   END TYPE ABecCecLvl_t

   CONTAINS
      
! ##############################################################################################
! LinOpLvl_t definitions
! ##############################################################################################
SUBROUTINE MyType(this)
   CLASS(ABecCecLvl_t) , INTENT(INOUT) :: this
   write(*,*) " I'm ABecCec type "
END SUBROUTINE MyType

SUBROUTINE define(this,lvl_id,nx,dx)
   CLASS(ABecCecLvl_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: lvl_id
   INTEGER, INTENT(IN) :: nx
   REAL(pr), INTENT(IN) :: dx
   IF(ALLOCATED(this%m_a)) CALL this%destroy()
   this%nx = nx
   this%dx = dx
   this%lvl_idx = lvl_id
   ALLOCATE(this%m_a(1:nx))
   ALLOCATE(this%m_b(1:nx+1))
   ALLOCATE(this%m_c(1:nx+1))
   this%is_allocated = .TRUE.
END SUBROUTINE define

SUBROUTINE destroy(this)
   CLASS(ABecCecLvl_t), INTENT(INOUT) :: this
   this%nx = 0
   this%dx = 0.0_pr
   this%lvl_idx = 0
   IF ( ALLOCATED(this%m_a) ) DEALLOCATE(this%m_a)
   IF ( ALLOCATED(this%m_b) ) DEALLOCATE(this%m_b)
   IF ( ALLOCATED(this%m_c) ) DEALLOCATE(this%m_c)
   this%is_allocated = .FALSE.
END SUBROUTINE destroy

SUBROUTINE setC(this,C_in)
   CLASS(ABecCecLvl_t), INTENT(INOUT) :: this
   REAL(pr), DIMENSION(1:this%nx+1), INTENT(IN) :: C_in
   this%m_c(1:this%nx+1) = C_in(1:this%nx+1)
END SUBROUTINE setC

END MODULE class_ABecCecLvl
