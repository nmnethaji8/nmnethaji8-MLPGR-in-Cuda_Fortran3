PROGRAM test
  USE, INTRINSIC :: iso_c_binding
  INTEGER(KIND=4)::NODN,INNZ
  INTEGER(KIND=4),ALLOCATABLE :: IVVCSR(:),JVVCSR(:)
  REAL(KIND=8),ALLOCATABLE:: FBCSR(:),SKKCSR(:), PTCSR(:)

  ! external functions
  INTERFACE
  SUBROUTINE fortran_solve_csr(n,m,nnz,rowoffset,col,val,rhs,x) BIND(c, name='fortran_solve_csr_')
    USE, INTRINSIC :: iso_c_binding
    integer(kind=4), value :: n, m, nnz
    integer(kind=4):: rowoffset(n+1),col(nnz)
    real(kind=8) :: val(nnz),rhs(n),x(n)
  END SUBROUTINE fortran_solve_csr
  END INTERFACE

  NODN=4
  ALLOCATE(IVVCSR(NODN+1),JVVCSR(NODN*NODN),SKKCSR(NODN*NODN),FBCSR(NODN), PTCSR(NODN))

  INNZ=7
  IVVCSR=(/0,2,4,6,7/)
  JVVCSR(1:INNZ)=(/1,2,2,4,3,4,4/)!(/0,1,1,3,2,3,3/)
  JVVCSR(1:INNZ)=JVVCSR(1:INNZ)-1
  SKKCSR(1:INNZ)=(/10,20,30,40,50,60,80/)
  FBCSR=1
  ! 10 20  0  0
  !  0 30  0 40
  !  0  0 50 60
  !  0  0  0 80 !B=1

  PRINT*,NODN, INNZ
  PRINT*,SIZE(IVVCSR),SIZE(JVVCSR(1:INNZ)),SIZE(SKKCSR(1:INNZ)),SIZE(FBCSR),SIZE(PTCSR)
  call fortran_solve_csr(NODN,NODN,INNZ,IVVCSR,JVVCSR(1:INNZ),SKKCSR(1:INNZ),FBCSR,PTCSR)

  PRINT,PTCSR(:)
END PROGRAM test