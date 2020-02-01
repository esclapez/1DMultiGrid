MODULE class_LinOpBC

   USE mod_prec_defs

   IMPLICIT NONE

!  Linear operator BC   
   TYPE LinOpBC_t
      INTEGER :: bc_type = -1                               ! 0 = Neumann, 1 = Dirichlet, -1 : Bogus
      REAL(pr) :: bc_val                                     ! BC value, only for Dirichlet
      CONTAINS
         PROCEDURE :: setType => LinOpBC_setType
         PROCEDURE :: isNeumann => LinOpBC_isNeumann
         PROCEDURE :: isDirichlet => LinOpBC_isDirichlet
         PROCEDURE :: setVal => LinOpBC_setDirichVal
         PROCEDURE :: getVal => LinOpBC_getDirichVal
         PROCEDURE :: LinOpBC_copy
         GENERIC :: ASSIGNMENT(=) => LinOpBC_copy 
   END TYPE LinOpBC_t

   CONTAINS

! ##############################################################################################
! LinOpBC_t definitions
! ##############################################################################################
SUBROUTINE LinOpBC_setType(this,type)   
   CLASS(LinOpBC_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: type
   IF((type/=0).AND.(type/=1))THEN
      WRITE(*,*) "LinOpBC::setType : type should be either 0 (Neumann) or 1 (Dirichlet)" 
      STOP
   END IF
   this%bc_type = type
END SUBROUTINE LinOpBC_setType

FUNCTION LinOpBC_isNeumann(this) RESULT(logicalCheck)
   CLASS(LinOpBC_t), INTENT(IN) :: this
   LOGICAL :: logicalCheck
   logicalCheck = .FALSE.
   IF(this%bc_type == 0) logicalCheck = .TRUE.
END FUNCTION LinOpBC_isNeumann

FUNCTION LinOpBC_isDirichlet(this) RESULT(logicalCheck)
   CLASS(LinOpBC_t), INTENT(IN) :: this
   LOGICAL :: logicalCheck
   logicalCheck = .FALSE.
   IF(this%bc_type == 1) logicalCheck = .TRUE.
END FUNCTION LinOpBC_isDirichlet

SUBROUTINE LinOpBC_setDirichVal(this,val)   
   CLASS(LinOpBC_t), INTENT(INOUT) :: this
   REAL(pr), INTENT(IN) :: val
   IF(.NOT. this%isDirichlet())THEN
      WRITE(*,*) "LinOpBC::setVal : BC is not Dirichlet"
      STOP
   END IF
   this%bc_val = val
END SUBROUTINE LinOpBC_setDirichVal

FUNCTION LinOpBC_getDirichVal(this) RESULT(val)
   CLASS(LinOpBC_t), INTENT(IN) :: this
   REAL(pr) :: val
   IF(.NOT. this%isDirichlet())THEN
      WRITE(*,*) "LinOpBC::getVal : BC is not Dirichlet"
      STOP
   END IF
   val = this%bc_val
END FUNCTION LinOpBC_getDirichVal

SUBROUTINE LinOpBC_copy(tothat, this)
   TYPE(LinOpBC_t), INTENT(IN) :: this
   CLASS(LinOpBC_t), INTENT(INOUT) :: tothat
   tothat%bc_type = this%bc_type
   tothat%bc_val = this%bc_val
END SUBROUTINE LinOpBC_copy

END MODULE class_LinOpBC
