!      subroutine bi_f77_init()
!      integer ierr
!      call mpi_init(ierr)
!      return
!      end

*
c      SUBROUTINE BI_F77_GET_CONSTANTS(F77COMMWORLD, SETUP, CONST)
      SUBROUTINE bi_f77_get_constants(F77COMMWORLD, SETUP, CONST)
      INCLUDE 'mpif.h'
      INTEGER F77COMMWORLD, SETUP
      INTEGER CONST(*)
*
      F77COMMWORLD = MPI_COMM_WORLD
      IF (SETUP .NE. 0) THEN
         CONST( 1) = MPI_SUCCESS
         CONST( 2) = MPI_ERR_UNKNOWN
         CONST( 3) = MPI_ERR_OTHER
         CONST( 4) = MPI_ERR_INTERN
         CONST( 5) = MPI_ANY_SOURCE
         CONST( 6) = MPI_UNDEFINED
         CONST( 7) = MPI_STATUS_SIZE
         CONST( 8) = MPI_SOURCE
         CONST( 9) = MPI_TAG
         CONST(10) = MPI_INTEGER
         CONST(11) = MPI_REAL
         CONST(12) = MPI_DOUBLE_PRECISION
         CONST(13) = MPI_COMPLEX
         CONST(14) = MPI_DOUBLE_COMPLEX
         CONST(15) = MPI_PACKED
         CONST(16) = MPI_BYTE
         CONST(17) = MPI_COMM_WORLD
         CONST(18) = MPI_COMM_NULL
         CONST(19) = MPI_TAG_UB
         CONST(20) = MPI_MAX
         CONST(21) = MPI_MIN
         CONST(22) = MPI_SUM
         CONST(23) = MPI_REQUEST_NULL
      END IF
*
      RETURN
      END
