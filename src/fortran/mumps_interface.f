      SUBROUTINE mumps_interface(job,parproc,comm,sym,icntl,n,nz,nz_loc,
     *                           irn_loc,jcn_loc,a_loc,b)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_struc.h'
      TYPE (MUMPS_STRUC) mumps_par
      integer i
      integer job,parproc,comm,sym
      integer n,nz,nz_loc
      integer,target :: irn_loc(nz_loc),jcn_loc(nz_loc)
      integer icntl(20)
      real*8 ,target :: a_loc(nz_loc),b(n)
      save mumps_par
      
c----------------------------------------- put values to fortran90 structure 
      mumps_par%JOB      = job
      mumps_par%COMM     = comm
      mumps_par%SYM      = sym
      mumps_par%PAR      = parproc
      mumps_par%N        = n
      mumps_par%NZ       = nz
      mumps_par%NZ_loc   = nz_loc
      mumps_par%IRN_loc => irn_loc(1:nz_loc)
      mumps_par%JCN_loc => jcn_loc(1:nz_loc)
      mumps_par%A_loc   => a_loc(1:nz_loc)
      mumps_par%RHS     => b(1:n)
c-----------------------------------------------------------------------init
      if (job.eq.-1) THEN
         mumps_par%ICNTL(3)  = 0
         icntl(3)            = 0 
         mumps_par%JOB      = -1
c----------------------------------------------------------- call init stage
         call mumps(mumps_par)
c-------------------------------------put default values of options to icntl
         do i=1,20
            icntl(i) = mumps_par%ICNTL(i)
         enddo
c---------------------------------- set type to distributed assembled matrix 
         mumps_par%ICNTL(18) = 3
         icntl(18)           = 3 
c--------------------------------------------------------------- set I/O off
         mumps_par%ICNTL(3)  = 0
         icntl(3)            = 0 
c-------------------------------------------------- call for analysis phase         
         mumps_par%JOB       = 1    
         call mumps(mumps_par)

         goto 100
      endif
c---------------------------------------------------------------end of init





c-----------------------------------------call factorization and solve phase
      if (job.eq.6) then
c--------------------------------------------------------set default options
         do i=1,20
            mumps_par%ICNTL(i) = icntl(i)
         enddo
c---------------------------------- set type to distributed assembled matrix 
         mumps_par%ICNTL(18) = 3
         icntl(18)           = 3 
c--------------------------------------------------------------- set I/O off
         mumps_par%ICNTL(3)  = 0
         icntl(3)            = 0 
c----------------------------------------------------call for  factorization 
         mumps_par%JOB     = 2
         call mumps(mumps_par)
c----------------------------------------------------------call for solution
         mumps_par%JOB     = 3
         call mumps(mumps_par)
c-------------------------------------------------------------- save options
         do i=1,20
            icntl(i) = mumps_par%ICNTL(i)
         enddo
         goto 100
      endif
c----------------------------------------- end of factorization and solution




c-------------------------------------------------------- call solution only      
      if (job.eq.3) then
c------------------------------------------------- set matrix to distributed
         mumps_par%ICNTL(18) = 3
         icntl(18)           = 3 
c--------------------------------------------------------------- set I/O off
         mumps_par%ICNTL(3)  = 0
         icntl(3)            = 0 
c-------------------------------------------------- set job to solution only         
         mumps_par%JOB      = 3
c----------------------------------------------------------------call solver          
         call mumps(mumps_par)
         goto 100
      endif
c-------------------------------------------------------end of solution only  



100   continue
      RETURN
      END
