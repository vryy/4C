      SUBROUTINE mumps_interface(job,parproc,comm,sym,n,nz,nz_loc,
     *                           irn_loc,jcn_loc,a_loc)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_struc.h'
      TYPE (MUMPS_STRUC) mumps_par
      integer job,parproc,comm,sym
      integer n,nz,nz_loc
      integer irn_loc(nz_loc),jcn_loc(nz_loc)
      real*8  a_loc(nz_loc)
      
c----------------------------------------- put values to fortran90 structure 
      mumps_par%JOB  = job
      mumps_par%COMM = comm
      mumps_par%SYM  = sym
      mumps_par%PAR  = parproc
c------------------------------------------------- call init phase of solver 
      if (job.eq.-1) then
         call mumps(mumps_par)
         goto 100
      endif



  







 












c      call mumps()
100   continue
      RETURN
      END
