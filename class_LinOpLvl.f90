! Definition of an abstract linear operator level class
! Some functions are deferred to the extended derived class
MODULE class_LinOpLvl

   USE mod_prec_defs
   USE class_LinOpBC, ONLY : LinOpBC_t

   IMPLICIT NONE

!  Linear operator MG level
   TYPE, ABSTRACT :: LinOpLvl_t
      INTEGER  :: nx                                        ! Size 
      INTEGER  :: lvl_idx                                   ! Position in the MG tower : 1 = finest
      REAL(pr) :: dx                                        ! Grid size
      TYPE(LinOpBC_t) :: leftBC, rightBC                    ! Boundary conditions   
      REAL(pr), DIMENSION(:), ALLOCATABLE :: m_a, m_b, m_c  ! Linear operator coefficients
      LOGICAL :: is_allocated = .FALSE.                     ! Level is allocated ?
      CONTAINS
          PROCEDURE( LinOpLvl_printType ), DEFERRED         :: MyType
          PROCEDURE( LinOpLvl_alloc ), DEFERRED, PASS       :: define
          PROCEDURE( LinOpLvl_dealloc ), DEFERRED           :: destroy
          PROCEDURE :: setBCType => LinOpLvl_setBCType
          PROCEDURE :: setBCVals => LinOpLvl_setBCVals
          PROCEDURE :: setA => LinOpLvl_setA
          PROCEDURE :: setB => LinOpLvl_setB
          PROCEDURE( LinOpLvl_setC ), DEFERRED, PASS        :: setC
          PROCEDURE :: isDefined => LinOpLvl_isDefined
!          PROCEDURE( LinOpLvl_copy ), DEFERRED              :: copy
!          GENERIC :: ASSIGNMENT(=) => copy 
   END TYPE LinOpLvl_t

   ABSTRACT INTERFACE
      SUBROUTINE LinOpLvl_printType(this)
         IMPORT LinOpLvl_t
         CLASS(LinOpLvl_t), INTENT(INOUT) :: this
      END SUBROUTINE LinOpLvl_printType
      SUBROUTINE LinOpLvl_alloc(this,lvl_id,nx,dx)
         IMPORT LinOpLvl_t
         IMPORT pr
         CLASS(LinOpLvl_t), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: lvl_id
         INTEGER, INTENT(IN) :: nx
         REAL(pr), INTENT(IN) :: dx
      END SUBROUTINE LinOpLvl_alloc
      SUBROUTINE LinOpLvl_dealloc(this)
         IMPORT LinOpLvl_t
         CLASS(LinOpLvl_t), INTENT(INOUT) :: this
      END SUBROUTINE LinOpLvl_dealloc
      SUBROUTINE LinOpLvl_setC(this,C_in)
         IMPORT LinOpLvl_t
         IMPORT pr
         CLASS(LinOpLvl_t), INTENT(INOUT) :: this
         REAL(pr), DIMENSION(1:this%nx+1), INTENT(IN) :: C_in
      END SUBROUTINE LinOpLvl_setC
!      SUBROUTINE LinOpLvl_copy(tothat, this)
!         IMPORT LinOpLvl_t
!         TYPE(LinOpLvl_t), INTENT(IN) :: this
!         CLASS(LinOpLvl_t), INTENT(INOUT) :: tothat
!      END SUBROUTINE LinOpLvl_copy
   END INTERFACE

   CONTAINS
      
! ##############################################################################################
! LinOpLvl_t definitions
! ##############################################################################################

SUBROUTINE LinOpLvl_setBCType(this,l_type,r_type)
   CLASS(LinOpLvl_t), INTENT(INOUT) :: this
   INTEGER, INTENT(IN) :: l_type, r_type
   call this%leftBC%setType(l_type)
   call this%rightBC%setType(r_type)
END SUBROUTINE LinOpLvl_setBCType

SUBROUTINE LinOpLvl_setBCVals(this,l_val,r_val)
   CLASS(LinOpLvl_t), INTENT(INOUT) :: this
   REAL(pr), INTENT(IN) :: l_val, r_val
   IF(this%leftBC%isDirichlet()) call this%leftBC%setVal(l_val)
   IF(this%rightBC%isDirichlet()) call this%rightBC%setVal(r_val)
END SUBROUTINE LinOpLvl_setBCVals

SUBROUTINE LinOpLvl_setA(this,A_in)
   CLASS(LinOpLvl_t), INTENT(INOUT) :: this
   REAL(pr), DIMENSION(1:this%nx), INTENT(IN) :: A_in
   this%m_a(1:this%nx) = A_in(1:this%nx)
END SUBROUTINE LinOpLvl_setA

SUBROUTINE LinOpLvl_setB(this,B_in)
   CLASS(LinOpLvl_t), INTENT(INOUT) :: this
   REAL(pr), DIMENSION(1:this%nx+1), INTENT(IN) :: B_in
   this%m_b(1:this%nx+1) = B_in(1:this%nx+1)
END SUBROUTINE LinOpLvl_setB

FUNCTION LinOpLvl_isDefined(this) RESULT(logicalCheck)
   CLASS(LinOpLvl_t), INTENT(IN) :: this
   LOGICAL :: logicalCheck
   logicalCheck = .FALSE.
   IF(this%is_allocated) logicalCheck = .TRUE.
END FUNCTION LinOpLvl_isDefined

!SUBROUTINE LinOpLvl_copy(tothat, this)
!   TYPE(LinOpLvl_t), INTENT(IN) :: this
!   CLASS(LinOpLvl_t), INTENT(INOUT) :: tothat
!   IF ( tothat%nx /= 0 ) CALL tothat%destroy()
!   tothat%nx = this%nx
!   tothat%dx = this%dx
!   tothat%lvl_idx = this%lvl_idx
!   tothat%leftBC = this%leftBC
!   tothat%rightBC = this%rightBC
!   ALLOCATE(tothat%m_a, source = this%m_a)
!   ALLOCATE(tothat%m_b, source = this%m_b)
!   ALLOCATE(tothat%m_c, source = this%m_c)
!   tothat%is_allocated = this%is_allocated 
!END SUBROUTINE LinOpLvl_copy

END MODULE class_LinOpLvl
