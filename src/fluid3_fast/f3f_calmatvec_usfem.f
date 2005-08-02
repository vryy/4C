C-----------------------------------------------------------------------
C
C \file
C \brief evaluate 3D fluid coefficient matrix 
C
C <pre>
C Maintainer: Christiane Foerster
C             foerster@statik.uni-stuttgart.de
C             http://www.uni-stuttgart.de/ibs/members/foerster/
C             0711 - 685-6572
C </pre>
C
C-------------------------------------------------------------------------

C-------------------------------------------------------------------------
C \brief evaluate fluid coefficient matrix
C
C <pre>                                                   chfoe 05/05
C
C In this routine the Gauss point contributions to the elemental coefficient
C matrix of a stabilised fluid2 element are calculated. The procedure is 
C based on the Rothe method of first integrating in time. Hence the 
C resulting terms include coefficients containing time integration variables
C such as theta or delta t which are represented by 'timefac'.
C
C The stabilisation is based on the residuum:
C
C R_M = u + timefac u * grad u - timefac * 2 nu div epsilon(u) 
C     + timefac grad p - rhsint
C
C R_C = div u
C
C The corresponding weighting operators are
C L_M = v + timefac u_old * grad v + timefac v * grad u_old 
C     - timefac * 2 nu alpha div epsilon (v) + timefac beta grad q
C
C L_C = div v
C
C where alpha = -1
C       beta  = -1 
C are sign regulating factors and rhsint differs for different time 
C These factores are worked in now and cannot be changed any more.
C
C integration schemes:
C
C One-step-Theta:
C rhsint = u_old + Theta dt f + (1-Theta) acc_old
C
C BDF2:
C
C generalised alpha:
C
C
C The stabilisation by means of the momentum residuum R_M is of the unusual
C type: 
C    Galerkin parts MINUS sum over elements (stabilising parts)
C The stabilisation by means of the continuity equation R_C is done in the
C usual way:
C    Galerkin parts PLUS sum over elements (stabilising parts)
C
C The calculation proceeds as follows.
C 1) obtain single (linearised) operators of R_M, R_C, L_M and L_C
C 2) build Galerkin terms from them
C 3) build stabilising terms from them
C 4) build Galerkin and stabilising terms of RHS
C
C NOTE: u_old represents the last iteration value. (The most recent one 
C       we've got!)
C 
C NOTE: Galerkin and stabilisation matrices are calculated within one 
C       routine.
C
C 
C Notational remarks:
C
C                    /              \
C                   | u_x,x   u_x,y |
C vderxy = grad u = |               |
C                   | u_y,x   u_y,y |
C                   \               /
C
C            /                         \
C           | u_x,xx   u_x,yy   u_x,xy |
C vderxy2 = |                          |
C           | u_y,xx   u_y,yy   u_y,xy |
C           \                          /
C
C for further comments see comment lines within code.
C
C
C </pre>
C \param **estif      DOUBLE        (o)   ele stiffness matrix
C \param  *eforce     DOUBLE        (o)   ele force vector
C \param  *velint     DOUBLE        (i)   vel at INT point
C \param  *histvec    DOUBLE        (i)   rhs at INT point
C \param  *gridvint   DOUBLE        (i)   gridvel at INT point
C \param **vderxy     DOUBLE        (i)   global vel derivatives
C \param  *vderxy2    DOUBLE        (i)   2nd global vel derivatives
C \param  *funct      DOUBLE        (i)   nat. shape funcs
C \param **derxy      DOUBLE        (i)   global coord. deriv
C \param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
C \param  *edeadng    DOUBLE        (i)   dead load at time n+1
C \param   fac        DOUBLE        (i)   weighting factor
C \param   visc       DOUBLE        (i)   fluid viscosity
C \param   iel        INT           (i)   number of nodes of act. ele
C \param  *hasext     INT           (i)   flag, if element has volume load
C \param   isale      INT           (i)   flag, if ALE or EULER
C
C-------------------------------------------------------------------------
C234567890123456789012345678901234567890123456789012345678901234567890
      subroutine f3fcalmat(estif, eforce, velint, histint, gridvint,
     &                     vderxy, vderxy2, pderxy,funct,
     &                     derxy, derxy2, edeadng, fac, tau, hasext,
     &                     paravec, flagvec, sizevec)
C
      implicit none

C     the vectors containing sizes and parameters
      real*8  paravec(3)
      integer flagvec(2)
      integer sizevec(6)

C     the remaining arguments
      real*8 estif(sizevec(4),sizevec(3),sizevec(3))
      real*8 eforce(sizevec(4),sizevec(3))
      real*8 velint(sizevec(4),3)
      real*8 histint(sizevec(4),3)
      real*8 gridvint(sizevec(4),3)
      real*8 vderxy(sizevec(4),3,3)
      real*8 vderxy2(sizevec(4),6,3)
      real*8 pderxy(sizevec(4),3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 edeadng(sizevec(4),3)
      real*8 fac(sizevec(4))
      real*8 tau(sizevec(4),3)
      integer hasext(sizevec(4))

C     locally required variables and arrays
      integer i, j, ri, ci
      integer cix, ciy, ciz, cip
      integer rix, riy, riz, rip
      integer l ! loop counter for THE LOOP
      real*8  timefac    ! One-step-Theta: timefac = theta*dt
C                          BDF2:           timefac = 2/3 * dt
      real*8  aux(sizevec(4))
      real*8  auxmat(sizevec(4),3,3)

      real*8  tau_M(sizevec(4))
      real*8  tau_Mp(sizevec(4))
      real*8  tau_C(sizevec(4))
      real*8  viscs2(sizevec(4),3,3*sizevec(1))
      real*8  viscous(sizevec(4),3,3,3*sizevec(1))
      real*8  conv_c(sizevec(4),sizevec(1))
      real*8  conv_g(sizevec(4),sizevec(1))
      real*8  conv_r(sizevec(4),3,3*sizevec(1))
      real*8  ugradv(sizevec(4),sizevec(1),3*sizevec(1))
      real*8  conv_old(sizevec(4),3), visc_old(sizevec(4),3)
      real*8  rhsint(sizevec(4),3)
      real*8  time2nue
      real*8  timetauM(sizevec(4)), timetauMp(sizevec(4))
      real*8  ttimetauM(sizevec(4)), ttimetauMp(sizevec(4))
      real*8  timefacfac(sizevec(4))


C========================== initialisation =============================
      timefac = paravec(3)

C     evaluate rhs vector at integration point
      do l=1,sizevec(5)
        if (hasext(l).ne.0) then
          rhsint(l,1) = timefac * edeadng(l,1) + histint(l,1)
          rhsint(l,2) = timefac * edeadng(l,2) + histint(l,2)
          rhsint(l,3) = timefac * edeadng(l,3) + histint(l,3)
        else
          rhsint(l,1) = histint(l,1)
          rhsint(l,2) = histint(l,2)
          rhsint(l,3) = histint(l,3)
        endif
      enddo

      do l=1,sizevec(5)
C       integration factors and koefficients of single terms
        tau_M(l)  = tau(l,1)*fac(l)
        tau_Mp(l) = tau(l,2)*fac(l)
        tau_C(l)  = tau(l,3)*fac(l)

        time2nue      = timefac * 2.0 * paravec(2)
        timetauM(l)   = timefac * tau_M(l)
        timetauMp(l)  = timefac * tau_Mp(l)

        ttimetauM(l)  = timefac * timetauM(l)
        ttimetauMp(l) = timefac * timetauMp(l)
        timefacfac(l) = timefac * fac(l)

C       get numerical representation of single operators
C       Convective term  u_old * grad u_old:
        conv_old(l,1) = vderxy(l,1,1) * velint(l,1)
     &                + vderxy(l,2,1) * velint(l,2)
     &                + vderxy(l,3,1) * velint(l,3)
        conv_old(l,2) = vderxy(l,1,2) * velint(l,1)
     &                + vderxy(l,2,2) * velint(l,2)
     &                + vderxy(l,3,2) * velint(l,3)
        conv_old(l,3) = vderxy(l,1,3) * velint(l,1)
     &                + vderxy(l,2,3) * velint(l,2)
     &                + vderxy(l,3,3) * velint(l,3)
C       Viscous term  div epsilon(u_old)
        visc_old(l,1) = vderxy2(l,1,1)
     &                + 0.5*( vderxy2(l,2,1) + vderxy2(l,4,2)
     &                      + vderxy2(l,3,1) + vderxy2(l,5,3) )
        visc_old(l,2) = vderxy2(l,2,2)
     &                + 0.5*( vderxy2(l,1,2) + vderxy2(l,4,1)
     &                      + vderxy2(l,3,2) + vderxy2(l,6,3) )
        visc_old(l,3) = vderxy2(l,3,3)
     &                + 0.5*( vderxy2(l,1,3) + vderxy2(l,5,1)
     &                      + vderxy2(l,2,3) + vderxy2(l,6,2) )

      enddo ! the loop
C
C
      do i=1,sizevec(2)       ! loop over elemental nodes
        do l=1,sizevec(5)     ! the loop
C         Reactive term  u:  funct
C         linearise convective term
C         convective part u_old * grad (funct)
C         u_old_x * N,x  +  u_old_y * N,y  +  u_old_z * N,z
C           with  N .. form function matrix
          conv_c(l,i) = derxy(l,i,1) * velint(l,1)
     &                + derxy(l,i,2) * velint(l,2)
     &                + derxy(l,i,3) * velint(l,3)
C         reactive part funct * grad (u_old)
C         /                                     \
C         |  u_old_x,x   u_old_x,y   u_old x,z  |
C         |                                     |
C         |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
C         |                                     |
C         |  u_old_z,x   u_old_z,y   u_old_z,z  |
C         \                                     /
C            with  N .. form function matrix
          conv_r(l,1,3*i-2) = vderxy(l,1,1) * funct(i)
          conv_r(l,1,3*i-1) = vderxy(l,2,1) * funct(i)
          conv_r(l,1,3*i)   = vderxy(l,3,1) * funct(i)
          conv_r(l,2,3*i-2) = vderxy(l,1,2) * funct(i)
          conv_r(l,2,3*i-1) = vderxy(l,2,2) * funct(i)
          conv_r(l,2,3*i)   = vderxy(l,3,2) * funct(i)
          conv_r(l,3,3*i-2) = vderxy(l,1,3) * funct(i)
          conv_r(l,3,3*i-1) = vderxy(l,2,3) * funct(i)
          conv_r(l,3,3*i)   = vderxy(l,3,3) * funct(i)
C         viscous term  - grad * epsilon(u):
C             /                                                \
C             |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
C           1 |                                                |
C         - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
C           2 |                                                |
C             |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
C             \                                                /
C            with N_x .. x-line of N
C                 N_y .. y-line of N
          viscs2(l,1,3*i-2) = - derxy2(l,i,1)
     &                        - 0.5 * ( derxy2(l,i,2) + derxy2(l,i,3) )
          viscs2(l,1,3*i-1) = - 0.5 *  derxy2(l,i,4)
          viscs2(l,1,3*i)   = - 0.5 *  derxy2(l,i,5)
          viscs2(l,2,3*i-2) = - 0.5 *  derxy2(l,i,4)
          viscs2(l,2,3*i-1) = - 0.5 * ( derxy2(l,i,1) + derxy2(l,i,3) )
     &                        - derxy2(l,i,2)
          viscs2(l,2,3*i)   = - 0.5 *  derxy2(l,i,6)
          viscs2(l,3,3*i-2) = - 0.5 *  derxy2(l,i,5)
          viscs2(l,3,3*i-1) = - 0.5 *  derxy2(l,i,6)
          viscs2(l,3,3*i)   = - 0.5 * ( derxy2(l,i,1) + derxy2(l,i,2) )
     &                        - derxy2(l,i,3)
C         viscous term (after integr. by parts)
C            /                                             \
C            |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
C          1 |                                             |
C          - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
C          2 |                                             |
C            |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
C            \                                             /
C            with N_x .. x-line of N
C                 N_y .. y-line of N
C                 N_z .. z-line of N
          viscous(l,1,1,3*i-2) = derxy(l,i,1)       ! 1st index:
          viscous(l,1,1,3*i-1) = 0.0                !   loop
          viscous(l,1,1,3*i)   = 0.0                ! 2nd index:
          viscous(l,1,2,3*i-2) = 0.5 * derxy(l,i,2) !   line of epsilon
          viscous(l,1,2,3*i-1) = 0.5 * derxy(l,i,1) ! 3rd index:
          viscous(l,1,2,3*i)   = 0.0                !   column of epsilon
          viscous(l,1,3,3*i-2) = 0.5 * derxy(l,i,3) ! 4th index:
          viscous(l,1,3,3*i-1) = 0.0                !   elemental vel dof
          viscous(l,1,3,3*i)   = 0.5 * derxy(l,i,1)
          viscous(l,2,1,3*i-2) = 0.5 * derxy(l,i,2)
          viscous(l,2,1,3*i-1) = 0.5 * derxy(l,i,1)
          viscous(l,2,1,3*i)   = 0.0
          viscous(l,2,2,3*i-2) = 0.0
          viscous(l,2,2,3*i-1) = derxy(l,i,2)
          viscous(l,2,2,3*i)   = 0.0
          viscous(l,2,3,3*i-2) = 0.0
          viscous(l,2,3,3*i-1) = 0.5 * derxy(l,i,3)
          viscous(l,2,3,3*i)   = 0.5 * derxy(l,i,2)
          viscous(l,3,1,3*i-2) = 0.5 * derxy(l,i,3)
          viscous(l,3,1,3*i-1) = 0.0
          viscous(l,3,1,3*i)   = 0.5 * derxy(l,i,1)
          viscous(l,3,2,3*i-2) = 0.0
          viscous(l,3,2,3*i-1) = 0.5 * derxy(l,i,3)
          viscous(l,3,2,3*i)   = 0.5 * derxy(l,i,2)
          viscous(l,3,3,3*i-2) = 0.0
          viscous(l,3,3,3*i-1) = 0.0
          viscous(l,3,3,3*i)   = derxy(l,i,3)
        enddo ! the loop

        if (flagvec(2).ne.0) then
          do l=1,sizevec(5)     ! the loop
C           convective grid part u_G * grad (funct)
C           u_G,old_x * N,x  +  u_G,old_y * N,y  +  u_G,old_z * N,z
C             with  N .. form function matrix */
            conv_g(l,i) = - derxy(l,i,1) * gridvint(l,1)
     &                    - derxy(l,i,2) * gridvint(l,2)
     &                    - derxy(l,i,3) * gridvint(l,3)
          enddo ! the loop
        endif

        do j=1,sizevec(2)
          do l=1,sizevec(5)     ! the loop
C           ugradv-Term
C          /                                                          \
C          |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
C          |                                                          |
C          |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
C          |                                                          |
C          |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
C          |                                           .              |
C          |  . . .                                        .          |
C          |                                                  Ni*Ni,y |
C          \                                                          /
C           remark: vgradu = ugradv^T

            ugradv(l,i,3*j-2) = derxy(l,i,1) * funct(j)
            ugradv(l,i,3*j-1) = derxy(l,i,2) * funct(j)
            ugradv(l,i,3*j)   = derxy(l,i,3) * funct(j)
          enddo ! the loop
        enddo ! j - loop over nodes
      enddo   ! i - loop over nodes
C
C
C     now build single stiffness terms
      do ri=1,sizevec(2)     ! row index
        do ci=1,sizevec(2)   ! column index
          do l=1,sizevec(5)  ! the loop
C         prepare indices
          cix = ci*4-3
          ciy = ci*4-2
          ciz = ci*4-1
          cip = ci*4
          rix = ri*4-3
          riy = ri*4-2
          riz = ri*4-1
          rip = ri*4
C       ************** integrate element coefficient matrix ************
C===================== GALERKIN part of the matrix =====================
C         a concentration of the following terms:
C          'mass matrix' (u,v)
C          N_c (u_old * grad u, v)
C          N_r (u * grad u_old, v)
            aux(l) = funct(ri) * ( funct(ci)*fac(l) 
     &                             + timefacfac(l) * conv_c(l,ci) )
            estif(l,cix,rix) = estif(l,cix,rix)
     &                       + funct(ri) * conv_r(l,1,3*ci-2) 
     &                       * timefacfac(l) + aux(l)
            estif(l,ciy,rix) = estif(l,ciy,rix)
     &                       + funct(ri) * conv_r(l,1,ci*3-1)
     &                       * timefacfac(l)
            estif(l,ciz,rix) = estif(l,ciz,rix)
     &                       + funct(ri) * conv_r(l,1,ci*3)
     &                       * timefacfac(l)
            estif(l,cix,riy) = estif(l,cix,riy)
     &                       + funct(ri) * conv_r(l,2,ci*3-2)
     &                       * timefacfac(l)
            estif(l,ciy,riy) = estif(l,ciy,riy)
     &                       + funct(ri) * conv_r(l,2,ci*3-1)
     &                       * timefacfac(l) + aux(l)
            estif(l,ciz,riy) = estif(l,ciz,riy)
     &                       + funct(ri) * conv_r(l,2,ci*3)
     &                       * timefacfac(l)
            estif(l,cix,riz) = estif(l,cix,riz)
     &                       + funct(ri) * conv_r(l,3,ci*3-2)
     &                       * timefacfac(l)
            estif(l,ciy,riz) = estif(l,ciy,riz)
     &                       + funct(ri) * conv_r(l,3,ci*3-1)
     &                       * timefacfac(l)
            estif(l,ciz,riz) = estif(l,ciz,riz)
     &                       + funct(ri) * conv_r(l,3,ci*3)
     &                       * timefacfac(l) + aux(l)
C           K (2 * nu * epsilon(u), epsilon(v)) */
            auxmat(l,1,1) = viscous(l,1,1,ri*3-2)*viscous(l,1,1,ci*3-2)
     &                    + viscous(l,1,2,ri*3-2)*viscous(l,1,2,ci*3-2)
     &                    + viscous(l,1,3,ri*3-2)*viscous(l,1,3,ci*3-2)
     &                    + viscous(l,2,1,ri*3-2)*viscous(l,2,1,ci*3-2)
     &                    + viscous(l,2,2,ri*3-2)*viscous(l,2,2,ci*3-2)
     &                    + viscous(l,2,3,ri*3-2)*viscous(l,2,3,ci*3-2)
     &                    + viscous(l,3,1,ri*3-2)*viscous(l,3,1,ci*3-2)
     &                    + viscous(l,3,2,ri*3-2)*viscous(l,3,2,ci*3-2)
     &                    + viscous(l,3,3,ri*3-2)*viscous(l,3,3,ci*3-2)
            auxmat(l,1,2) = viscous(l,1,1,ri*3-2)*viscous(l,1,1,ci*3-1)
     &                    + viscous(l,1,2,ri*3-2)*viscous(l,1,2,ci*3-1)
     &                    + viscous(l,1,3,ri*3-2)*viscous(l,1,3,ci*3-1)
     &                    + viscous(l,2,1,ri*3-2)*viscous(l,2,1,ci*3-1)
     &                    + viscous(l,2,2,ri*3-2)*viscous(l,2,2,ci*3-1)
     &                    + viscous(l,2,3,ri*3-2)*viscous(l,2,3,ci*3-1)
     &                    + viscous(l,3,1,ri*3-2)*viscous(l,3,1,ci*3-1)
     &                    + viscous(l,3,2,ri*3-2)*viscous(l,3,2,ci*3-1)
     &                    + viscous(l,3,3,ri*3-2)*viscous(l,3,3,ci*3-1)
            auxmat(l,1,3) = viscous(l,1,1,ri*3-2)*viscous(l,1,1,ci*3)
     &                    + viscous(l,1,2,ri*3-2)*viscous(l,1,2,ci*3)
     &                    + viscous(l,1,3,ri*3-2)*viscous(l,1,3,ci*3)
     &                    + viscous(l,2,1,ri*3-2)*viscous(l,2,1,ci*3)
     &                    + viscous(l,2,2,ri*3-2)*viscous(l,2,2,ci*3)
     &                    + viscous(l,2,3,ri*3-2)*viscous(l,2,3,ci*3)
     &                    + viscous(l,3,1,ri*3-2)*viscous(l,3,1,ci*3)
     &                    + viscous(l,3,2,ri*3-2)*viscous(l,3,2,ci*3)
     &                    + viscous(l,3,3,ri*3-2)*viscous(l,3,3,ci*3)
            auxmat(l,2,1) = viscous(l,1,1,ri*3-1)*viscous(l,1,1,ci*3-2)
     &                    + viscous(l,1,2,ri*3-1)*viscous(l,1,2,ci*3-2)
     &                    + viscous(l,1,3,ri*3-1)*viscous(l,1,3,ci*3-2)
     &                    + viscous(l,2,1,ri*3-1)*viscous(l,2,1,ci*3-2)
     &                    + viscous(l,2,2,ri*3-1)*viscous(l,2,2,ci*3-2)
     &                    + viscous(l,2,3,ri*3-1)*viscous(l,2,3,ci*3-2)
     &                    + viscous(l,3,1,ri*3-1)*viscous(l,3,1,ci*3-2)
     &                    + viscous(l,3,2,ri*3-1)*viscous(l,3,2,ci*3-2)
     &                    + viscous(l,3,3,ri*3-1)*viscous(l,3,3,ci*3-2)
            auxmat(l,2,2) = viscous(l,1,1,ri*3-1)*viscous(l,1,1,ci*3-1)
     &                    + viscous(l,1,2,ri*3-1)*viscous(l,1,2,ci*3-1)
     &                    + viscous(l,1,3,ri*3-1)*viscous(l,1,3,ci*3-1)
     &                    + viscous(l,2,1,ri*3-1)*viscous(l,2,1,ci*3-1)
     &                    + viscous(l,2,2,ri*3-1)*viscous(l,2,2,ci*3-1)
     &                    + viscous(l,2,3,ri*3-1)*viscous(l,2,3,ci*3-1)
     &                    + viscous(l,3,1,ri*3-1)*viscous(l,3,1,ci*3-1)
     &                    + viscous(l,3,2,ri*3-1)*viscous(l,3,2,ci*3-1)
     &                    + viscous(l,3,3,ri*3-1)*viscous(l,3,3,ci*3-1)
            auxmat(l,2,3) = viscous(l,1,1,ri*3-1)*viscous(l,1,1,ci*3)
     &                    + viscous(l,1,2,ri*3-1)*viscous(l,1,2,ci*3)
     &                    + viscous(l,1,3,ri*3-1)*viscous(l,1,3,ci*3)
     &                    + viscous(l,2,1,ri*3-1)*viscous(l,2,1,ci*3)
     &                    + viscous(l,2,2,ri*3-1)*viscous(l,2,2,ci*3)
     &                    + viscous(l,2,3,ri*3-1)*viscous(l,2,3,ci*3)
     &                    + viscous(l,3,1,ri*3-1)*viscous(l,3,1,ci*3)
     &                    + viscous(l,3,2,ri*3-1)*viscous(l,3,2,ci*3)
     &                    + viscous(l,3,3,ri*3-1)*viscous(l,3,3,ci*3)
            auxmat(l,3,1) = viscous(l,1,1,ri*3)*viscous(l,1,1,ci*3-2)
     &                    + viscous(l,1,2,ri*3)*viscous(l,1,2,ci*3-2)
     &                    + viscous(l,1,3,ri*3)*viscous(l,1,3,ci*3-2)
     &                    + viscous(l,2,1,ri*3)*viscous(l,2,1,ci*3-2)
     &                    + viscous(l,2,2,ri*3)*viscous(l,2,2,ci*3-2)
     &                    + viscous(l,2,3,ri*3)*viscous(l,2,3,ci*3-2)
     &                    + viscous(l,3,1,ri*3)*viscous(l,3,1,ci*3-2)
     &                    + viscous(l,3,2,ri*3)*viscous(l,3,2,ci*3-2)
     &                    + viscous(l,3,3,ri*3)*viscous(l,3,3,ci*3-2)
            auxmat(l,3,2) = viscous(l,1,1,ri*3)*viscous(l,1,1,ci*3-1)
     &                    + viscous(l,1,2,ri*3)*viscous(l,1,2,ci*3-1)
     &                    + viscous(l,1,3,ri*3)*viscous(l,1,3,ci*3-1)
     &                    + viscous(l,2,1,ri*3)*viscous(l,2,1,ci*3-1)
     &                    + viscous(l,2,2,ri*3)*viscous(l,2,2,ci*3-1)
     &                    + viscous(l,2,3,ri*3)*viscous(l,2,3,ci*3-1)
     &                    + viscous(l,3,1,ri*3)*viscous(l,3,1,ci*3-1)
     &                    + viscous(l,3,2,ri*3)*viscous(l,3,2,ci*3-1)
     &                    + viscous(l,3,3,ri*3)*viscous(l,3,3,ci*3-1)
            auxmat(l,3,3) = viscous(l,1,1,ri*3)*viscous(l,1,1,ci*3)
     &                    + viscous(l,1,2,ri*3)*viscous(l,1,2,ci*3)
     &                    + viscous(l,1,3,ri*3)*viscous(l,1,3,ci*3)
     &                    + viscous(l,2,1,ri*3)*viscous(l,2,1,ci*3)
     &                    + viscous(l,2,2,ri*3)*viscous(l,2,2,ci*3)
     &                    + viscous(l,2,3,ri*3)*viscous(l,2,3,ci*3)
     &                    + viscous(l,3,1,ri*3)*viscous(l,3,1,ci*3)
     &                    + viscous(l,3,2,ri*3)*viscous(l,3,2,ci*3)
     &                    + viscous(l,3,3,ri*3)*viscous(l,3,3,ci*3)
            aux(l) = time2nue * fac(l)
            estif(l,cix,rix) = estif(l,cix,rix) + auxmat(l,1,1)*aux(l)
            estif(l,ciy,rix) = estif(l,ciy,rix) + auxmat(l,1,2)*aux(l)
            estif(l,ciz,rix) = estif(l,ciz,rix) + auxmat(l,1,3)*aux(l)
            estif(l,cix,riy) = estif(l,cix,riy) + auxmat(l,2,1)*aux(l)
            estif(l,ciy,riy) = estif(l,ciy,riy) + auxmat(l,2,2)*aux(l)
            estif(l,ciz,riy) = estif(l,ciz,riy) + auxmat(l,2,3)*aux(l)
            estif(l,cix,riz) = estif(l,cix,riz) + auxmat(l,3,1)*aux(l)
            estif(l,ciy,riz) = estif(l,ciy,riz) + auxmat(l,3,2)*aux(l)
            estif(l,ciz,riz) = estif(l,ciz,riz) + auxmat(l,3,3)*aux(l)
C           G (- div v, p)
            estif(l,cip,rix) = estif(l,cip,rix)
     &                     - timefacfac(l) * derxy(l,ri,1) * funct(ci)
            estif(l,cip,riy) = estif(l,cip,riy)
     &                     - timefacfac(l) * derxy(l,ri,2) * funct(ci)
            estif(l,cip,riz) = estif(l,cip,riz)
     &                     - timefacfac(l) * derxy(l,ri,3) * funct(ci)
C           G^T ( div u, q)
            estif(l,cix,rip) = estif(l,cix,rip)
     &                   + timefacfac(l) * funct(ri) * derxy(l,ci,1)
            estif(l,ciy,rip) = estif(l,ciy,rip)
     &                   + timefacfac(l) * funct(ri) * derxy(l,ci,2)
            estif(l,ciz,rip) = estif(l,ciz,rip)
     &                   + timefacfac(l) * funct(ri) * derxy(l,ci,3)

C==================== Stabilisation part of the matrix =================
C           CONVECTIVE stabilisation
C           a concentration of the following two terms:
C           tau_M*timefac*(u, u_old * grad v)
C           -tau_M*timefac*timefac*(u_old * grad u, u_old * grad v)
            aux(l) = conv_c(l,ri) * 
     &              ( timetauM(l)*funct(ci)+ttimetauM(l)*conv_c(l,ci) )
            estif(l,cix,rix) = estif(l,cix,rix) + aux(l)
            estif(l,ciy,riy) = estif(l,ciy,riy) + aux(l)
            estif(l,ciz,riz) = estif(l,ciz,riz) + aux(l)
C           a concentration of the following two terms:
C           -tau_M*timefac*timefac*(u * grad u_old, u_old * grad v)
C           tau_M*timefac*timefac*2*nu*(div epsilon(u), u_old * grad v)
            aux(l) = timetauM(l) * time2nue
            estif(l,cix,rix) = estif(l,cix,rix)
     &            + conv_c(l,ri) * ( conv_r(l,1,ci*3-2) * ttimetauM(l)
     &                             + viscs2(l,1,ci*3-2) * aux(l) )
            estif(l,ciy,rix) = estif(l,ciy,rix)
     &            + conv_c(l,ri) * ( conv_r(l,1,ci*3-1) * ttimetauM(l)
     &                             + viscs2(l,1,ci*3-1) * aux(l) )
            estif(l,ciz,rix) = estif(l,ciz,rix)
     &            + conv_c(l,ri) * ( conv_r(l,1,ci*3) * ttimetauM(l)
     &                             + viscs2(l,1,ci*3) * aux(l) )
            estif(l,cix,riy) = estif(l,cix,riy)
     &            + conv_c(l,ri) * ( conv_r(l,2,ci*3-2) * ttimetauM(l)
     &                             + viscs2(l,2,ci*3-2) * aux(l) )
            estif(l,ciy,riy) = estif(l,ciy,riy)
     &            + conv_c(l,ri) * ( conv_r(l,2,ci*3-1) * ttimetauM(l)
     &                             + viscs2(l,2,ci*3-1) * aux(l) )
            estif(l,ciz,riy) = estif(l,ciz,riy)
     &            + conv_c(l,ri) * ( conv_r(l,2,ci*3) * ttimetauM(l)
     &                             + viscs2(l,2,ci*3) * aux(l) )
            estif(l,cix,riz) = estif(l,cix,riz)
     &            + conv_c(l,ri) * ( conv_r(l,3,ci*3-2) * ttimetauM(l)
     &                             + viscs2(l,3,ci*3-2) * aux(l) )
            estif(l,ciy,riz) = estif(l,ciy,riz)
     &            + conv_c(l,ri) * ( conv_r(l,3,ci*3-1) * ttimetauM(l)
     &                             + viscs2(l,3,ci*3-1) * aux(l) )
            estif(l,ciz,riz) = estif(l,ciz,riz)
     &            + conv_c(l,ri) * ( conv_r(l,3,ci*3) * ttimetauM(l)
     &                             + viscs2(l,3,ci*3) * aux(l) )
C           -tau_M*timefac*timefac*(grad p, u_old * grad v)
            estif(l,cip,rix) = estif(l,cip,rix)
     &            + conv_c(l,ri) * derxy(l,ci,1) * ttimetauM(l)
            estif(l,cip,riy) = estif(l,cip,riy)
     &            + conv_c(l,ri) * derxy(l,ci,2) * ttimetauM(l)
            estif(l,cip,riz) = estif(l,cip,riz)
     &            + conv_c(l,ri) * derxy(l,ci,3) * ttimetauM(l)
C           DIFFUSION part of stabilisation
C           a concentration of the following two terms:
C           tau_M*timefac*2*nu*(u, div epsilon(v))
C           tau_M*timefac*timefac*2*nu*(u_old * grad u, div epsilon(v))
            aux(l) = time2nue * ( funct(ci)*tau_Mp(l) + conv_c(l,ci) 
     &                                                * timetauMp(l) )
            estif(l,cix,rix) = estif(l,cix,rix)
     &                       + viscs2(l,1,ri*3-2) * aux(l)
            estif(l,ciy,rix) = estif(l,ciy,rix)
     &                       + viscs2(l,2,ri*3-2) * aux(l)
            estif(l,ciz,rix) = estif(l,ciz,rix)
     &                       + viscs2(l,3,ri*3-2) * aux(l)
            estif(l,cix,riy) = estif(l,cix,riy)
     &                       + viscs2(l,1,ri*3-1) * aux(l)
            estif(l,ciy,riy) = estif(l,ciy,riy)
     &                       + viscs2(l,2,ri*3-1) * aux(l)
            estif(l,ciz,riy) = estif(l,ciz,riy)
     &                       + viscs2(l,3,ri*3-1) * aux(l)
            estif(l,cix,riz) = estif(l,cix,riz)
     &                       + viscs2(l,1,ri*3) * aux(l)
            estif(l,ciy,riz) = estif(l,ciy,riz)
     &                       + viscs2(l,2,ri*3) * aux(l)
            estif(l,ciz,riz) = estif(l,ciz,riz)
     &                       + viscs2(l,3,ri*3) * aux(l)
C           tau_M*timefac*timefac*2*nu*(u * grad u_old, div epsilon(v))
            aux(l) = timetauMp(l) * time2nue
            estif(l,cix,rix) = estif(l,cix,rix) + aux(l) *
     &                      ( viscs2(l,1,ri*3-2) * conv_r(l,1,ci*3-2)
     &                      + viscs2(l,2,ri*3-2) * conv_r(l,2,ci*3-2)
     &                      + viscs2(l,3,ri*3-2) * conv_r(l,3,ci*3-2) )
            estif(l,cix,riy) = estif(l,cix,riy) + aux(l) *
     &                      ( viscs2(l,1,ri*3-1) * conv_r(l,1,ci*3-2)
     &                      + viscs2(l,2,ri*3-1) * conv_r(l,2,ci*3-2)
     &                      + viscs2(l,3,ri*3-1) * conv_r(l,3,ci*3-2) )
            estif(l,cix,riz) = estif(l,cix,riz) + aux(l) *
     &                      ( viscs2(l,1,ri*3) * conv_r(l,1,ci*3-2)
     &                      + viscs2(l,2,ri*3) * conv_r(l,2,ci*3-2)
     &                      + viscs2(l,3,ri*3) * conv_r(l,3,ci*3-2) )
            estif(l,ciy,rix) = estif(l,ciy,rix) + aux(l) *
     &                      ( viscs2(l,1,ri*3-2) * conv_r(l,1,ci*3-1)
     &                      + viscs2(l,2,ri*3-2) * conv_r(l,2,ci*3-1)
     &                      + viscs2(l,3,ri*3-2) * conv_r(l,3,ci*3-1) )
            estif(l,ciy,riy) = estif(l,ciy,riy) + aux(l) *
     &                      ( viscs2(l,1,ri*3-1) * conv_r(l,1,ci*3-1)
     &                      + viscs2(l,2,ri*3-1) * conv_r(l,2,ci*3-1)
     &                      + viscs2(l,3,ri*3-1) * conv_r(l,3,ci*3-1) )
            estif(l,ciy,riz) = estif(l,ciy,riz) + aux(l) *
     &                      ( viscs2(l,1,ri*3) * conv_r(l,1,ci*3-1)
     &                      + viscs2(l,2,ri*3) * conv_r(l,2,ci*3-1)
     &                      + viscs2(l,3,ri*3) * conv_r(l,3,ci*3-1) )
            estif(l,ciz,rix) = estif(l,ciz,rix) + aux(l) *
     &                      ( viscs2(l,1,ri*3-2) * conv_r(l,1,ci*3)
     &                      + viscs2(l,2,ri*3-2) * conv_r(l,2,ci*3)
     &                      + viscs2(l,3,ri*3-2) * conv_r(l,3,ci*3) )
            estif(l,ciz,riy) = estif(l,ciz,riy) + aux(l) *
     &                      ( viscs2(l,1,ri*3-1) * conv_r(l,1,ci*3)
     &                      + viscs2(l,2,ri*3-1) * conv_r(l,2,ci*3)
     &                      + viscs2(l,3,ri*3-1) * conv_r(l,3,ci*3) )
            estif(l,ciz,riz) = estif(l,ciz,riz) + aux(l) *
     &                      ( viscs2(l,1,ri*3) * conv_r(l,1,ci*3)
     &                      + viscs2(l,2,ri*3) * conv_r(l,2,ci*3)
     &                      + viscs2(l,3,ri*3) * conv_r(l,3,ci*3) )
C           -tau_M*timefac*timefac*4*nu^2(div epsilon(u), div epsilon(v))
            aux(l) = time2nue * time2nue * tau_Mp(l)
            estif(l,cix,rix) = estif(l,cix,rix) + aux(l) *
     &                      ( viscs2(l,1,ri*3-2) * viscs2(l,1,ci*3-2)
     &                      + viscs2(l,2,ri*3-2) * viscs2(l,2,ci*3-2)
     &                      + viscs2(l,3,ri*3-2) * viscs2(l,3,ci*3-2) )
            estif(l,cix,riy) = estif(l,cix,riy) + aux(l) *
     &                      ( viscs2(l,1,ri*3-1) * viscs2(l,1,ci*3-2)
     &                      + viscs2(l,2,ri*3-1) * viscs2(l,2,ci*3-2)
     &                      + viscs2(l,3,ri*3-1) * viscs2(l,3,ci*3-2) )
            estif(l,cix,riz) = estif(l,cix,riz) + aux(l) *
     &                      ( viscs2(l,1,ri*3) * viscs2(l,1,ci*3-2)
     &                      + viscs2(l,2,ri*3) * viscs2(l,2,ci*3-2)
     &                      + viscs2(l,3,ri*3) * viscs2(l,3,ci*3-2) )
            estif(l,ciy,rix) = estif(l,ciy,rix) + aux(l) *
     &                      ( viscs2(l,1,ri*3-2) * viscs2(l,1,ci*3-1)
     &                      + viscs2(l,2,ri*3-2) * viscs2(l,2,ci*3-1)
     &                      + viscs2(l,3,ri*3-2) * viscs2(l,3,ci*3-1) )
            estif(l,ciy,riy) = estif(l,ciy,riy) + aux(l) *
     &                      ( viscs2(l,1,ri*3-1) * viscs2(l,1,ci*3-1)
     &                      + viscs2(l,2,ri*3-1) * viscs2(l,2,ci*3-1)
     &                      + viscs2(l,3,ri*3-1) * viscs2(l,3,ci*3-1) )
            estif(l,ciy,riz) = estif(l,ciy,riz) + aux(l) *
     &                      ( viscs2(l,1,ri*3) * viscs2(l,1,ci*3-1)
     &                      + viscs2(l,2,ri*3) * viscs2(l,2,ci*3-1)
     &                      + viscs2(l,3,ri*3) * viscs2(l,3,ci*3-1) )
            estif(l,ciz,rix) = estif(l,ciz,rix) + aux(l) *
     &                      ( viscs2(l,1,ri*3-2) * viscs2(l,1,ci*3)
     &                      + viscs2(l,2,ri*3-2) * viscs2(l,2,ci*3)
     &                      + viscs2(l,3,ri*3-2) * viscs2(l,3,ci*3) )
            estif(l,ciz,riy) = estif(l,ciz,riy) + aux(l) *
     &                      ( viscs2(l,1,ri*3-1) * viscs2(l,1,ci*3)
     &                      + viscs2(l,2,ri*3-1) * viscs2(l,2,ci*3)
     &                      + viscs2(l,3,ri*3-1) * viscs2(l,3,ci*3) )
            estif(l,ciz,riz) = estif(l,ciz,riz) + aux(l) *
     &                      ( viscs2(l,1,ri*3) * viscs2(l,1,ci*3)
     &                      + viscs2(l,2,ri*3) * viscs2(l,2,ci*3)
     &                      + viscs2(l,3,ri*3) * viscs2(l,3,ci*3) )
C           tau_M*timefac*timefac*2*nu*(grad p, div epsilon(v))
            aux(l) = time2nue * timetauMp(l)
            estif(l,cip,rix) = estif(l,cip,rix) + aux(l) *
     &                       ( viscs2(l,1,ri*3-2) * derxy(l,ci,1)
     &                       + viscs2(l,2,ri*3-2) * derxy(l,ci,2)
     &                       + viscs2(l,3,ri*3-2) * derxy(l,ci,3) )
            estif(l,cip,riy) = estif(l,cip,riy) + aux(l) *
     &                       ( viscs2(l,1,ri*3-1) * derxy(l,ci,1)
     &                       + viscs2(l,2,ri*3-1) * derxy(l,ci,2)
     &                       + viscs2(l,3,ri*3-1) * derxy(l,ci,3) )
            estif(l,cip,riz) = estif(l,cip,riz) + aux(l) *
     &                       ( viscs2(l,1,ri*3) * derxy(l,ci,1)
     &                       + viscs2(l,2,ri*3) * derxy(l,ci,2)
     &                       + viscs2(l,3,ri*3) * derxy(l,ci,3) )
C           PRESSURE part of stabilisation
C           a concentration of the following terms:
C           -tau_M*timefac*(u, grad q)
C           -tau_M*timefac*timefac*(u_old * grad u, grad q)
            aux(l) = funct(ci)*timetauMp(l)+conv_c(l,ci)*ttimetauMp(l)
            estif(l,cix,rip) = estif(l,cix,rip)
     &                       + derxy(l,ri,1) * aux(l)
            estif(l,ciy,rip) = estif(l,ciy,rip)
     &                       + derxy(l,ri,2) * aux(l)
            estif(l,ciz,rip) = estif(l,ciz,rip)
     &                       + derxy(l,ri,3) * aux(l)
C           -tau_M*timefac*timefac*(u * grad u_old, grad q)
            estif(l,cix,rip) = estif(l,cix,rip)+ttimetauMp(l) *
     &                       ( derxy(l,ri,1) * conv_r(l,1,ci*3-2)
     &                       + derxy(l,ri,2) * conv_r(l,2,ci*3-2)
     &                       + derxy(l,ri,3) * conv_r(l,3,ci*3-2) )
            estif(l,ciy,rip) = estif(l,ciy,rip)+ttimetauMp(l) *
     &                       ( derxy(l,ri,1) * conv_r(l,1,ci*3-1)
     &                       + derxy(l,ri,2) * conv_r(l,2,ci*3-1)
     &                       + derxy(l,ri,3) * conv_r(l,3,ci*3-1) )
            estif(l,ciz,rip) = estif(l,ciz,rip)+ttimetauMp(l) *
     &                       ( derxy(l,ri,1) * conv_r(l,1,ci*3)
     &                       + derxy(l,ri,2) * conv_r(l,2,ci*3)
     &                       + derxy(l,ri,3) * conv_r(l,3,ci*3) )
C           tau_M*timefac*timefac*2*nu*(div epsilon(u), grad q)
            aux(l) = timetauMp(l) * time2nue
            estif(l,cix,rip) = estif(l,cix,rip) + aux(l) *
     &                       ( derxy(l,ri,1) * viscs2(l,1,ci*3-2)
     &                       + derxy(l,ri,2) * viscs2(l,2,ci*3-2)
     &                       + derxy(l,ri,3) * viscs2(l,3,ci*3-2) )
            estif(l,ciy,rip) = estif(l,ciy,rip) + aux(l) *
     &                       ( derxy(l,ri,1) * viscs2(l,1,ci*3-1)
     &                       + derxy(l,ri,2) * viscs2(l,2,ci*3-1)
     &                       + derxy(l,ri,3) * viscs2(l,3,ci*3-1) )
            estif(l,ciz,rip) = estif(l,ciz,rip) + aux(l) *
     &                       ( derxy(l,ri,1) * viscs2(l,1,ci*3)
     &                       + derxy(l,ri,2) * viscs2(l,2,ci*3)
     &                       + derxy(l,ri,3) * viscs2(l,3,ci*3) )
C           -tau_M*timefac*timefac*(grad p, grad q)
            estif(l,cip,rip) = estif(l,cip,rip) 
     &                       + timefac * timetauMp(l) * (
     &                              derxy(l,ri,1) * derxy(l,ci,1)
     &                            + derxy(l,ri,2) * derxy(l,ci,2)
     &                            + derxy(l,ri,3) * derxy(l,ci,3) )
C           R(u_old) * L_conv STABILISATION
C           a concentration of the following terms:
C           -tau_M*timefac*(u_old, u * grad v)
C           -tau_M*timefac*timefac*(u_old * grad u_old, u * grad v)
C           tau_M*timefac*timefac*2*nu*(div epsilon(u_old), u * grad v)
C           -tau_M*timefac*timefac*(grad p_old, u * grad v)
C           linear part of RHS stabilisation (goes into matrix)
C           tau_M*timefac*(rhsint, u * grad v)
            aux(l) = - timetauM(l) * time2nue
            estif(l,cix,rix) = estif(l,cix,rix)
     &                  + ( ( velint(l,1) - rhsint(l,1) ) * timetauM(l)
     &                  + ( conv_old(l,1) + pderxy(l,1) ) * ttimetauM(l)
     &                    + visc_old(l,1) * aux(l) )*ugradv(l,ri,ci*3-2)
            estif(l,ciy,rix) = estif(l,ciy,rix)
     &                  + ( ( velint(l,1) - rhsint(l,1) ) * timetauM(l)
     &                  + ( conv_old(l,1) + pderxy(l,1) ) * ttimetauM(l)
     &                    + visc_old(l,1) * aux(l) )*ugradv(l,ri,ci*3-1)
            estif(l,ciz,rix) = estif(l,ciz,rix)
     &                  + ( ( velint(l,1) - rhsint(l,1) ) * timetauM(l)
     &                  + ( conv_old(l,1) + pderxy(l,1) ) * ttimetauM(l)
     &                    + visc_old(l,1) * aux(l) )*ugradv(l,ri,ci*3)
            estif(l,cix,riy) = estif(l,cix,riy)
     &                  + ( ( velint(l,2) - rhsint(l,2) ) * timetauM(l)
     &                  + ( conv_old(l,2) + pderxy(l,2) ) * ttimetauM(l)
     &                    + visc_old(l,2) * aux(l) )*ugradv(l,ri,ci*3-2)
            estif(l,ciy,riy) = estif(l,ciy,riy)
     &                  + ( ( velint(l,2) - rhsint(l,2) ) * timetauM(l)
     &                  + ( conv_old(l,2) + pderxy(l,2) ) * ttimetauM(l)
     &                    + visc_old(l,2) * aux(l) )*ugradv(l,ri,ci*3-1)
            estif(l,ciz,riy) = estif(l,ciz,riy)
     &                  + ( ( velint(l,2) - rhsint(l,2) ) * timetauM(l)
     &                  + ( conv_old(l,2) + pderxy(l,2) ) * ttimetauM(l)
     &                    + visc_old(l,2) * aux(l) )*ugradv(l,ri,ci*3)
            estif(l,cix,riz) = estif(l,cix,riz)
     &                  + ( ( velint(l,3) - rhsint(l,3) ) * timetauM(l)
     &                  + ( conv_old(l,3) + pderxy(l,3) ) * ttimetauM(l)
     &                    + visc_old(l,3) * aux(l) )*ugradv(l,ri,ci*3-2)
            estif(l,ciy,riz) = estif(l,ciy,riz)
     &                  + ( ( velint(l,3) - rhsint(l,3) ) * timetauM(l)
     &                  + ( conv_old(l,3) + pderxy(l,3) ) * ttimetauM(l)
     &                    + visc_old(l,3) * aux(l) )*ugradv(l,ri,ci*3-1)
            estif(l,ciz,riz) = estif(l,ciz,riz)
     &                  + ( ( velint(l,3) - rhsint(l,3) ) * timetauM(l)
     &                  + ( conv_old(l,3) + pderxy(l,3) ) * ttimetauM(l)
     &                    + visc_old(l,3) * aux(l) )*ugradv(l,ri,ci*3)
C           CONTINUITY equation stabilisation
C           tau_C*timefac*timefac*(div u, div v)
            aux(l) = timefac * timefac * tau_C(l)
            estif(l,cix,rix) = estif(l,cix,rix)
     &                       + derxy(l,ri,1)*derxy(l,ci,1)*aux(l)
            estif(l,ciy,rix) = estif(l,ciy,rix)
     &                       + derxy(l,ri,1)*derxy(l,ci,2)*aux(l)
            estif(l,ciz,rix) = estif(l,ciz,rix)
     &                       + derxy(l,ri,1)*derxy(l,ci,3)*aux(l)
            estif(l,cix,riy) = estif(l,cix,riy)
     &                       + derxy(l,ri,2)*derxy(l,ci,1)*aux(l)
            estif(l,ciy,riy) = estif(l,ciy,riy)
     &                       + derxy(l,ri,2)*derxy(l,ci,2)*aux(l)
            estif(l,ciz,riy) = estif(l,ciz,riy)
     &                       + derxy(l,ri,2)*derxy(l,ci,3)*aux(l)
            estif(l,cix,riz) = estif(l,cix,riz)
     &                       + derxy(l,ri,3)*derxy(l,ci,1)*aux(l)
            estif(l,ciy,riz) = estif(l,ciy,riz)
     &                       + derxy(l,ri,3)*derxy(l,ci,2)*aux(l)
            estif(l,ciz,riz) = estif(l,ciz,riz)
     &                       + derxy(l,ri,3)*derxy(l,ci,3)*aux(l)
          enddo  ! the loop
        enddo  ! ci - column loop
C    *************** integrate element force vector ********************
C    ================== Galerkin part of the RHS =======================
        do l=1,sizevec(5)  ! the loop
C         'Original' RHS, concentrated
C         (rhsint, v)
C         from Nonlinearity of Galerkin stiffness ---*/
C         timefac*(u_old * grad u_old, v) */
          eforce(l,rix) = eforce(l,4*ri-3) + funct(ri) 
     &           * ( rhsint(l,1)*fac(l) + conv_old(l,1)*timefacfac(l) )
          eforce(l,riy) = eforce(l,4*ri-2) + funct(ri) 
     &           * ( rhsint(l,2)*fac(l) + conv_old(l,2)*timefacfac(l) )
          eforce(l,riz) = eforce(l,riz) + funct(ri) 
     &           * ( rhsint(l,3)*fac(l) + conv_old(l,3)*timefacfac(l) )
C         from 'time integration' of conservation equation
C         (dt-timefac)*(div u_old, q)
C         this is switched off at the moment!
C         aux = - (dt-timefac) * fac;
C         eforce[rip+3] += *divuold * funct[ri] * aux;
C
C     =============== Stabilisation part of the RHS ====================
C         'Original' RHS
C         tau_M*timefac*2*nu*(rhsint, div epsilon(v))
          aux(l) = time2nue * tau_Mp(l)
          eforce(l,rix) = eforce(l,rix)
     &                   + ( rhsint(l,1) * viscs2(l,1,ri*3-2)
     &                     + rhsint(l,2) * viscs2(l,2,ri*3-2)
     &                     + rhsint(l,3) * viscs2(l,3,ri*3-2) ) * aux(l)
          eforce(l,riy) = eforce(l,riy)
     &                   + ( rhsint(l,1) * viscs2(l,1,ri*3-1)
     &                     + rhsint(l,2) * viscs2(l,2,ri*3-1)
     &                     + rhsint(l,3) * viscs2(l,3,ri*3-1) ) * aux(l)
          eforce(l,riz) = eforce(l,riz)
     &                   + ( rhsint(l,1) * viscs2(l,1,ri*3)
     &                     + rhsint(l,2) * viscs2(l,2,ri*3)
     &                     + rhsint(l,3) * viscs2(l,3,ri*3) ) * aux(l)
C         -tau_M*timefac*(rhsint, grad q)
          eforce(l,rip) = eforce(l,rip) +
     &                  (  rhsint(l,1) * derxy(l,ri,1)
     &                   + rhsint(l,2) * derxy(l,ri,2)
     &                   + rhsint(l,3) * derxy(l,ri,3) ) * timetauMp(l)
C         Terms resulting from stiffness linearisation
C         a concentration of the following:
C         -tau_M*timefac*(u_old, u_old * grad v)
C       tau_M*timefac*timefac*2*nu*(div epsilon(u_old), u_old * grad v)
C         -tau_M*timefac*timefac*(grad p_old, u_old * grad v)
          aux(l) = - timetauM(l) * time2nue
          eforce(l,rix) = eforce(l,rix)
     &                + conv_c(l,ri) * ( velint(l,1)*timetauM(l)
     &                                 + visc_old(l,1)*aux(l)
     &                                 + pderxy(l,1)*ttimetauM(l) )
          eforce(l,riy) = eforce(l,riy)
     &                + conv_c(l,ri) * ( velint(l,2)*timetauM(l)
     &                                 + visc_old(l,2)*aux(l)
     &                                 + pderxy(l,2)*ttimetauM(l) )
          eforce(l,riz) = eforce(l,riz)
     &                + conv_c(l,ri) * ( velint(l,3)*timetauM(l)
     &                                 + visc_old(l,3)*aux(l)
     &                                 + pderxy(l,3)*ttimetauM(l) )
C         -tau_M*2*timefac*timefac*(u_old * grad u_old, u_old * grad v)
          aux(l) = ttimetauM(l) * 2.0
          eforce(l,rix) = eforce(l,rix)
     &                  + conv_old(l,1) * conv_c(l,ri) * aux(l)
          eforce(l,riy) = eforce(l,riy)
     &                  + conv_old(l,2) * conv_c(l,ri) * aux(l)
          eforce(l,riz) = eforce(l,riz)
     &                  + conv_old(l,3) * conv_c(l,ri) * aux(l)
C        tau_M*timefac*timefac*2*nu*(u_old * grad u_old, div epsilon(v))
          aux(l) = timetauMp(l) * time2nue
          eforce(l,rix) = eforce(l,rix) + aux(l) *
     &                 ( conv_old(l,1) * viscs2(l,1,ri*3-2)
     &                  +conv_old(l,2) * viscs2(l,2,ri*3-2)
     &                  +conv_old(l,3) * viscs2(l,3,ri*3-2) )
          eforce(l,riy) = eforce(l,riy) + aux(l) *
     &                 ( conv_old(l,1) * viscs2(l,1,ri*3-1)
     &                  +conv_old(l,2) * viscs2(l,2,ri*3-1)
     &                  +conv_old(l,3) * viscs2(l,3,ri*3-1) ) 
          eforce(l,riz) = eforce(l,riz) + aux(l) *
     &                 ( conv_old(l,1) * viscs2(l,1,ri*3)
     &                  +conv_old(l,2) * viscs2(l,2,ri*3)
     &                  +conv_old(l,3) * viscs2(l,3,ri*3) )  
C         -tau_M*timefac*timefac*(u_old * grad u_old, grad q)
          eforce(l,rip) = eforce(l,rip) + ttimetauMp(l) *
     &                  ( conv_old(l,1) * derxy(l,ri,1)
     &                   +conv_old(l,2) * derxy(l,ri,2)
     &                   +conv_old(l,3) * derxy(l,ri,3) )
        enddo ! the loop
      enddo  ! ri - row loop

C     care for ALE terms
      if(flagvec(2).ne.0) then
        do ri=1,sizevec(2)     ! row index
          do ci=1,sizevec(2)   ! column index
            do l=1,sizevec(5)  ! the loop

C             N_c (-u_G * grad u, v)
              aux(l) = timefacfac(l) * funct(ri) * conv_g(l,ci)
              estif(l,cix,rix) = estif(l,cix,rix) + aux(l)
              estif(l,ciy,riy) = estif(l,ciy,riy) + aux(l)
              estif(l,ciz,riz) = estif(l,ciz,riz) + aux(l)
C             -tau_M*timefac*timefac*(-u_G * grad u, u_old * grad v)
              aux(l) = ttimetauM(l) * conv_c(l,ri) * conv_g(l,ci)
              estif(l,cix,rix) = estif(l,cix,rix) + aux(l)
              estif(l,ciy,riy) = estif(l,ciy,riy) + aux(l)
              estif(l,ciz,riz) = estif(l,ciz,riz) + aux(l)
C             CONVECTIVE GRID stabilisation
C             a concentration of the following terms:
C             -tau_M*timefac*(u, -u_G * grad v)
C             -tau_M*timefac*timefac*(u_old * grad u, -u_G * grad v)
C             -tau_M*timefac*timefac*(-u_G * grad u, -u_G * grad v)
              aux(l) = conv_g(l,ri) * ( ttimetauM(l) * ( conv_c(l,ri)
     &            + conv_g(l,ci) ) + timetauM(l) * funct(ci) )
              estif(l,cix,rix) = estif(l,cix,rix) + aux(l)
              estif(l,ciy,riy) = estif(l,ciy,riy) + aux(l)
              estif(l,ciz,riz) = estif(l,ciz,riz) + aux(l)
C             a concentration of the following two terms:
C             -tau_M*timefac*timefac*(u * grad u_old, -u_G * grad v)
C             tau_M*timefac*timefac*2*nu*(div epsilon(u), -u_G * grad v)
              aux(l) = timetauM(l) * time2nue
              estif(l,cix,rix) = estif(l,cix,rix)
     &              + conv_g(l,ri) * ( conv_r(l,1,ci*3-2) * ttimetauM(l)
     &                               + viscs2(l,1,ci*3-2) * aux(l) )
              estif(l,ciy,rix) = estif(l,ciy,rix)
     &              + conv_g(l,ri) * ( conv_r(l,1,ci*3-1) * ttimetauM(l)
     &                               + viscs2(l,1,ci*3-1) * aux(l) )
              estif(l,ciz,rix) = estif(l,ciz,rix)
     &              + conv_g(l,ri) * ( conv_r(l,1,ci*3) * ttimetauM(l)
     &                               + viscs2(l,1,ci*3) * aux(l) )
              estif(l,cix,riy) = estif(l,cix,riy)
     &              + conv_g(l,ri) * ( conv_r(l,2,ci*3-2) * ttimetauM(l)
     &                               + viscs2(l,2,ci*3-2) * aux(l) )
              estif(l,ciy,riy) = estif(l,ciy,riy)
     &              + conv_g(l,ri) * ( conv_r(l,2,ci*3-1) * ttimetauM(l)
     &                               + viscs2(l,2,ci*3-1) * aux(l) )
              estif(l,ciz,riy) = estif(l,ciz,riy)
     &              + conv_g(l,ri) * ( conv_r(l,2,ci*3) * ttimetauM(l)
     &                               + viscs2(l,2,ci*3) * aux(l) )
              estif(l,cix,riz) = estif(l,cix,riz)
     &              + conv_g(l,ri) * ( conv_r(l,3,ci*3-2) * ttimetauM(l)
     &                               + viscs2(l,3,ci*3-2) * aux(l) )
              estif(l,ciy,riz) = estif(l,ciy,riz)
     &              + conv_g(l,ri) * ( conv_r(l,3,ci*3-1) * ttimetauM(l)
     &                               + viscs2(l,3,ci*3-1) * aux(l) )
              estif(l,ciz,riz) = estif(l,ciz,riz)
     &              + conv_g(l,ri) * ( conv_r(l,3,ci*3) * ttimetauM(l)
     &                               + viscs2(l,3,ci*3) * aux(l) )
C             -tau_M*timefac*timefac*(grad p, -u_G * grad v)
              estif(l,cip,rix) = estif(l,cip,rix)
     &                         + conv_g(l,ri) * derxy(l,ci,1)
     &                         * ttimetauM(l)
              estif(l,cip,riy) = estif(l,cip,riy)
     &                         + conv_g(l,ri) * derxy(l,ci,2)
     &                         * ttimetauM(l)
              estif(l,cip,riz) = estif(l,cip,riz)
     &                         + conv_g(l,ri) * derxy(l,ci,3)
     &                         * ttimetauM(l)
C             tau_M*timefac*timefac*2*nu*(-u_G * grad u, div epsilon(v))
              aux(l) = timetauMp(l) * time2nue * conv_g(l,ci)
              estif(l,cix,rix) = estif(l,cix,rix)
     &                         + viscs2(l,1,ri*3-2) * aux(l)
              estif(l,ciy,rix) = estif(l,ciy,rix)
     &                         + viscs2(l,2,ri*3-2) * aux(l)
              estif(l,ciz,rix) = estif(l,ciz,rix)
     &                         + viscs2(l,3,ri*3-2) * aux(l)
              estif(l,cix,riy) = estif(l,cix,riy)
     &                         + viscs2(l,1,ri*3-1) * aux(l)
              estif(l,ciy,riy) = estif(l,ciy,riy)
     &                         + viscs2(l,2,ri*3-1) * aux(l)
              estif(l,ciz,riy) = estif(l,ciz,riy)
     &                         + viscs2(l,3,ri*3-1) * aux(l)
              estif(l,cix,riz) = estif(l,cix,riz)
     &                         + viscs2(l,1,ri*3) * aux(l)
              estif(l,ciy,riz) = estif(l,ciy,riz)
     &                         + viscs2(l,2,ri*3) * aux(l)
              estif(l,ciz,riz) = estif(l,ciz,riz)
     &                         + viscs2(l,3,ri*3) * aux(l)
C             -tau_M*timefac*timefac*(-u_G * grad u, grad q)
              aux(l) = conv_g(l,ci) * ttimetauMp(l)
              estif(l,cix,rip) = estif(l,cix,rip)
     &                         + derxy(l,ri,1) * aux(l)
              estif(l,ciy,rip) = estif(l,ciy,rip)
     &                         + derxy(l,ri,2) * aux(l)
              estif(l,ciz,rip) = estif(l,ciz,rip)
     &                         + derxy(l,ri,3) * aux(l)
            enddo  ! the loop
          enddo  ! ci - column loop
          do l=1,sizevec(5)  ! the loop
C           -tau_M*timefac*(rhsint, -u_G * grad v)
            eforce(l,rix) = eforce(l,rix)
     &                    + rhsint(l,1) * conv_g(l,ri) * timetauM(l)
            eforce(l,riy) = eforce(l,riy)
     &                    + rhsint(l,2) * conv_g(l,ri) * timetauM(l)
            eforce(l,riz) = eforce(l,riz)
     &                    + rhsint(l,3) * conv_g(l,ri) * timetauM(l)
C         ALE:
C         -tau_M*timefac*timefac*(u_old * grad u_old, u_old * grad v)
            eforce(l,rix) = eforce(l,rix)
     &                    + conv_old(l,1) * conv_g(l,ri) * ttimetauM(l)
            eforce(l,riy) = eforce(l,riy)
     &                    + conv_old(l,2) * conv_g(l,ri) * ttimetauM(l)
            eforce(l,riz) = eforce(l,riz)
     &                    + conv_old(l,3) * conv_g(l,ri) * ttimetauM(l)
          enddo ! the loop
        enddo  ! ri - row loop
      endif  ! ALE

      return
      end subroutine




C-----------------------------------------------------------------------
C \brief Gauss point contributions for integration of boundary forces
C
C <pre>                                                     chfoe 05/05
C
C This routine evaluates the Gauss point vaulues of the residual vector 
C of one element taking stabilisation effects into account. Only the 
C residual of the momentum equation R_M is calculated.
C
C R_M = u + timefac u * grad u - timefac * 2 nu div epsilon(u) 
C     + timefac grad p - rhsint
C
C The residual contains stabilisation of the type
C
C Sum_over_k (R_M, tau L_M)_k with
C
C L_M = v + timefac u_old * grad v + timefac v * grad u_old 
C     - timefac * 2 nu alpha div epsilon (v) + timefac beta grad q
C
C where alpha = -1
C       beta  = -1 
C
C timefac depends on the time integration scheme:
C
C One-step theta:
C
C timefac = theta * dt
C
C BDF2:
C
C timefac = 2/3 * dt
C
C NOTE: this works perfectly only when the fluid is solved via usfem
C
C </pre>
C \param  *eforce    DOUBLE    (o)    element force vector (residual)
C \param  *velint    DOUBLE    (i)    (converged) vel. at int.-point
C \param   histvec   DOUBLE    (i)    histroy data
C \param **vderxy    DOUBLE    (i)    velocity gradient at int.-point
C \param **vderxy2   DOUBLE    (i)    second vel. derivatives at int.-point
C \param  *funct     DOUBLE    (i)    natural shape functions
C \param **derxy     DOUBLE    (i)    shape function derivatives
C \param **derxy2    DOUBLE    (i)    second shape funct. derivs
C \param  *edeadng   DOUBLE    (i)    body forces
C \param  *press     DOUBLE    (i)    pressure at Gauss point
C \param   gradp[2]  DOUBLE    (i)    pressure gradient at GP
C \param   fac       DOUBLE    (i)    integration factor
C \param   visc      DOUBLE    (i)    fluid viscosity
C \param   iel       INT       (i)    number of elemental nodes
C \param  *hasext    INT       (i)    flag, if there is body force
C \param   is_ale    INT       (i)    flag, if it's ale or Euler
C \return void
C 
C ----------------------------------------------------------------------
C23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine f3fcalresvec(eforce,velint,histint,vderxy,vderxy2,
     &                        funct,derxy,derxy2,edeadng,aleconvint,
     &                        presint,pderxy,fac,hasext,tau,
     &                        paravec,flagvec,sizevec)
      implicit none

C     the vectors containing sizes and parameters
      real*8  paravec(4)
      integer flagvec(2)
      integer sizevec(6)

C     the remaining arguments
      real*8 eforce(sizevec(4),sizevec(3))
      real*8 velint(sizevec(4),3)
      real*8 histint(sizevec(4),3)
      real*8 vderxy(sizevec(4),3,3)
      real*8 vderxy2(sizevec(4),6,3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 edeadng(sizevec(4),3)
      real*8 aleconvint(sizevec(4),3)
      real*8 presint(sizevec(4))
      real*8 pderxy(sizevec(4),3)
      real*8 fac(sizevec(4))
      integer hasext(sizevec(4))
      real*8 tau(sizevec(4),3)
      
C     locally required variables and arrays
      integer i, ri
      integer rix, riy, riz
      integer l ! loop counter for THE LOOP
      real*8  timefac    ! One-step-Theta: timefac = theta*dt
C                          BDF2:           timefac = 2/3 * dt
      real*8  invtime    !  1 / timefac
      real*8  dt
      real*8  twovisc

      real*8  tau_M(sizevec(4)) 
      real*8  tau_Mp(sizevec(4))
      real*8  tau_C(sizevec(4))
      real*8  viscous(sizevec(4),3,3,3*sizevec(1))
      real*8  viscs2(sizevec(4),3,3*sizevec(1))
      real*8  visc2(sizevec(4),3)
      real*8  eps_u(sizevec(4),3,3)
      real*8  conv_c(sizevec(4),sizevec(1))
      real*8  rhsint(sizevec(4),3)
      real*8  resid(sizevec(4),3)

C========================== initialisation =============================
      timefac = paravec(3)
      dt      = paravec(4)
      invtime = 1 / timefac
      twovisc = 2.0 * paravec(2)

C     evaluate rhs vector at integration point ...
C      ... including actual velocity and convective term
      do l=1,sizevec(5)
        if(flagvec(1).ne.0) then
          if(hasext(l).ne.0) then
            rhsint(l,1) = timefac * ( edeadng(l,1)
     &                  - vderxy(l,1,1) * aleconvint(l,1)
     &                  - vderxy(l,2,1) * aleconvint(l,2)
     &                  - vderxy(l,3,1) * aleconvint(l,3) )
     &                  + histint(l,1) - velint(l,1)
            rhsint(l,2) = timefac * ( edeadng(l,2)
     &                  - vderxy(l,1,2) * aleconvint(l,1)
     &                  - vderxy(l,2,2) * aleconvint(l,2)
     &                  - vderxy(l,3,2) * aleconvint(l,3) )
     &                  + histint(l,2) - velint(l,2)
            rhsint(l,3) = timefac * ( edeadng(l,3)
     &                  - vderxy(l,1,3) * aleconvint(l,1)
     &                  - vderxy(l,2,3) * aleconvint(l,2)
     &                  - vderxy(l,3,3) * aleconvint(l,3) )
     &                  + histint(l,3) - velint(l,3)
          else
            rhsint(l,1) = histint(l,1) - velint(l,1) + timefac * (
     &                  - vderxy(l,1,1) * aleconvint(l,1)
     &                  - vderxy(l,2,1) * aleconvint(l,2)
     &                  - vderxy(l,3,1) * aleconvint(l,3) )
            rhsint(l,2) = histint(l,2) - velint(l,2) + timefac * (
     &                  - vderxy(l,1,2) * aleconvint(l,1)
     &                  - vderxy(l,2,2) * aleconvint(l,2)
     &                  - vderxy(l,3,2) * aleconvint(l,3) )
            rhsint(l,3) = histint(l,3) - velint(l,3) + timefac * (
     &                  - vderxy(l,1,3) * aleconvint(l,1)
     &                  - vderxy(l,2,3) * aleconvint(l,2)
     &                  - vderxy(l,3,3) * aleconvint(l,3) )
          endif
        else ! no ALE
          if(hasext(l).ne.0) then
            rhsint(l,1) = timefac * ( edeadng(l,1)
     &                  - vderxy(l,1,1) * velint(l,1)
     &                  - vderxy(l,2,1) * velint(l,2)
     &                  - vderxy(l,3,1) * velint(l,3) )
     &                  + histint(l,1) - velint(l,1)
            rhsint(l,2) = timefac * ( edeadng(l,2)
     &                  - vderxy(l,1,2) * velint(l,1)
     &                  - vderxy(l,2,2) * velint(l,2)
     &                  - vderxy(l,3,2) * velint(l,3) )
     &                  + histint(l,2) - velint(l,2)
            rhsint(l,3) = timefac * ( edeadng(l,3)
     &                  - vderxy(l,1,3) * velint(l,1)
     &                  - vderxy(l,2,3) * velint(l,2)
     &                  - vderxy(l,3,3) * velint(l,3) )
     &                  + histint(l,3) - velint(l,3)
          else
            rhsint(l,1) = histint(l,1) - velint(l,1) + timefac * (
     &                  - vderxy(l,1,1) * velint(l,1)
     &                  - vderxy(l,2,1) * velint(l,2)
     &                  - vderxy(l,3,1) * velint(l,3) )
            rhsint(l,2) = histint(l,2) - velint(l,2) + timefac * (
     &                  - vderxy(l,1,2) * velint(l,1)
     &                  - vderxy(l,2,2) * velint(l,2)
     &                  - vderxy(l,3,2) * velint(l,3) )
            rhsint(l,3) = histint(l,3) - velint(l,3) + timefac * (
     &                  - vderxy(l,1,3) * velint(l,1)
     &                  - vderxy(l,2,3) * velint(l,2)
     &                  - vderxy(l,3,3) * velint(l,3) )
          endif
        endif
      enddo

      do l=1,sizevec(5)
        tau_M(l)  = tau(l,1)*fac(l)
        tau_Mp(l) = tau(l,2)*fac(l)
        tau_C(l)  = tau(l,3)*fac(l)

C       viscous term (after integr. by parts)
C         /                                             \
C         |  2 u_x,x    u_x,y + u_y,x    u_x,z + u_z,x  |
C       1 |                                             |
C       - |  u_y,x + u_x,y    2 u_y,y    u_y,z + u_z,y  |
C       2 |                                             |
C         |  u_x,z + u_z,x    u_y,z + u_z,y    2 u_z,z  |
C         \                                             /
        eps_u(l,1,1) = vderxy(l,1,1)
        eps_u(l,1,2) = 0.5 * ( vderxy(l,2,1) + vderxy(l,1,2) )
        eps_u(l,1,3) = 0.5 * ( vderxy(l,3,1) + vderxy(l,1,3) )
        eps_u(l,2,1) = eps_u(l,1,2)
        eps_u(l,2,2) = vderxy(l,2,2)
        eps_u(l,2,3) = 0.5 * ( vderxy(l,3,2) + vderxy(l,2,3) )
        eps_u(l,3,1) = eps_u(l,1,3)
        eps_u(l,3,2) = eps_u(l,2,3)
        eps_u(l,3,3) = vderxy(l,3,3)
C 
C       viscous term (without integr. by parts)
C       Viscous term  div epsilon(u_old)
        visc2(l,1) = vderxy2(l,1,1) 
     &             + 0.5 * ( vderxy2(l,2,1) + vderxy2(l,4,2)
     &                     + vderxy2(l,3,1) + vderxy2(l,5,3) )
        visc2(l,2) = vderxy2(l,2,2)
     &             + 0.5 * ( vderxy2(l,1,2) + vderxy2(l,4,1)
     &                     + vderxy2(l,3,2) + vderxy2(l,6,3) )
        visc2(l,3) = vderxy2(l,3,3) 
     &             + 0.5 * ( vderxy2(l,1,3) + vderxy2(l,5,1)
     &                     + vderxy2(l,2,3) + vderxy2(l,6,2) )

        resid(l,1) = rhsint(l,1) + timefac 
     &             * ( twovisc * visc2(l,1) - pderxy(l,1) )
        resid(l,2) = rhsint(l,2) + timefac 
     &             * ( twovisc * visc2(l,2) - pderxy(l,2) )
        resid(l,3) = rhsint(l,3) + timefac 
     &             * ( twovisc * visc2(l,3) - pderxy(l,3) )

      enddo
C     get partially integrated terms
      do i=1,sizevec(2)
        do l=1,sizevec(5)
C         viscous term (after integr. by parts)
C            /                                             \
C            |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
C          1 |                                             |
C          - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
C          2 |                                             |
C            |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
C            \                                             /
C            with N_x .. x-line of N
C                 N_y .. y-line of N
C                 N_z .. z-line of N
          viscous(l,1,1,3*i-2) = derxy(l,i,1)       ! 1st index:
          viscous(l,1,1,3*i-1) = 0.0                !   loop
          viscous(l,1,1,3*i)   = 0.0                ! 2nd index:
          viscous(l,1,2,3*i-2) = 0.5 * derxy(l,i,2) !   line of epsilon
          viscous(l,1,2,3*i-1) = 0.5 * derxy(l,i,1) ! 3rd index:
          viscous(l,1,2,3*i)   = 0.0                !   column of epsilon
          viscous(l,1,3,3*i-2) = 0.5 * derxy(l,i,3) ! 4th index:
          viscous(l,1,3,3*i-1) = 0.0                !   elemental vel dof
          viscous(l,1,3,3*i)   = 0.5 * derxy(l,i,1)
          viscous(l,2,1,3*i-2) = 0.5 * derxy(l,i,2)
          viscous(l,2,1,3*i-1) = 0.5 * derxy(l,i,1)
          viscous(l,2,1,3*i)   = 0.0
          viscous(l,2,2,3*i-2) = 0.0
          viscous(l,2,2,3*i-1) = derxy(l,i,2)
          viscous(l,2,2,3*i)   = 0.0
          viscous(l,2,3,3*i-2) = 0.0
          viscous(l,2,3,3*i-1) = 0.5 * derxy(l,i,3)
          viscous(l,2,3,3*i)   = 0.5 * derxy(l,i,2)
          viscous(l,3,1,3*i-2) = 0.5 * derxy(l,i,3)
          viscous(l,3,1,3*i-1) = 0.0
          viscous(l,3,1,3*i)   = 0.5 * derxy(l,i,1)
          viscous(l,3,2,3*i-2) = 0.0
          viscous(l,3,2,3*i-1) = 0.5 * derxy(l,i,3)
          viscous(l,3,2,3*i)   = 0.5 * derxy(l,i,2)
          viscous(l,3,3,3*i-2) = 0.0
          viscous(l,3,3,3*i-1) = 0.0
          viscous(l,3,3,3*i)   = derxy(l,i,3)

C================= build stabilisation operatior ===================
C         viscous term  - grad * epsilon(u):
C             /                                                \
C             |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
C           1 |                                                |
C         - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
C           2 |                                                |
C             |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
C             \                                                /
C            with N_x .. x-line of N
C                 N_y .. y-line of N
          viscs2(l,1,3*i-2) = - derxy2(l,i,1)
     &                        - 0.5 * ( derxy2(l,i,2) + derxy2(l,i,3) )
          viscs2(l,1,3*i-1) = - 0.5 *  derxy2(l,i,4)
          viscs2(l,1,3*i)   = - 0.5 *  derxy2(l,i,5)
          viscs2(l,2,3*i-2) = - 0.5 *  derxy2(l,i,4)
          viscs2(l,2,3*i-1) = - 0.5 * ( derxy2(l,i,1) + derxy2(l,i,3) )
     &                        - derxy2(l,i,2)
          viscs2(l,2,3*i)   = - 0.5 *  derxy2(l,i,6)
          viscs2(l,3,3*i-2) = - 0.5 *  derxy2(l,i,5)
          viscs2(l,3,3*i-1) = - 0.5 *  derxy2(l,i,6)
          viscs2(l,3,3*i)   = - 0.5 * ( derxy2(l,i,1) + derxy2(l,i,2) )
     &                        - derxy2(l,i,3)

        enddo ! the loop

C       convective part u_old * grad (funct)
        do l=1,sizevec(5) ! the loop
          if(flagvec(2).ne.0) then
            conv_c(l,i) = derxy(l,i,1) * aleconvint(l,1)
     &                  + derxy(l,i,2) * aleconvint(l,2)
     &                  + derxy(l,i,3) * aleconvint(l,3)
          else
C         u_old_x * N,x  +  u_old_y * N,y  +  u_old_z * N,z
C           with  N .. form function matrix
            conv_c(l,i) = derxy(l,i,1) * velint(l,1)
     &                  + derxy(l,i,2) * velint(l,2)
     &                  + derxy(l,i,3) * velint(l,3)
          endif
        enddo ! the loop
      enddo   ! nodes i
C 
C     now build single residual terms
      do ri=1,sizevec(2)     ! row index
        do l=1,sizevec(5)  ! the loop
          rix = ri*4-3
          riy = ri*4-2
          riz = ri*4-1
C         ************ integrate element residuum vector ***************
C         simple parts, which are not partially integrated
          eforce(l,rix) = eforce(l,rix) 
     &                  + funct(ri) * rhsint(l,1) * invtime * fac(l)
          eforce(l,riy) = eforce(l,riy) 
     &                  + funct(ri) * rhsint(l,2) * invtime * fac(l)
          eforce(l,riz) = eforce(l,riz) 
     &                  + funct(ri) * rhsint(l,3) * invtime * fac(l)
C         viscous forces integrated by parts
          eforce(l,rix) = eforce(l,rix) - twovisc * fac(l) *
     &                    ( viscous(l,1,1,ri*3-2) * eps_u(l,1,1)
     &                     +viscous(l,1,2,ri*3-2) * eps_u(l,1,2)
     &                     +viscous(l,1,3,ri*3-2) * eps_u(l,1,3)
     &                     +viscous(l,2,1,ri*3-2) * eps_u(l,2,1)
     &                     +viscous(l,2,2,ri*3-2) * eps_u(l,2,2)
     &                     +viscous(l,2,3,ri*3-2) * eps_u(l,2,3)
     &                     +viscous(l,3,1,ri*3-2) * eps_u(l,3,1)
     &                     +viscous(l,3,2,ri*3-2) * eps_u(l,3,2)
     &                     +viscous(l,3,3,ri*3-2) * eps_u(l,3,3) )
          eforce(l,riy) = eforce(l,riy) - twovisc * fac(l) *
     &                    ( viscous(l,1,1,ri*3-1) * eps_u(l,1,1)
     &                     +viscous(l,1,2,ri*3-1) * eps_u(l,1,2)
     &                     +viscous(l,1,3,ri*3-1) * eps_u(l,1,3)
     &                     +viscous(l,2,1,ri*3-1) * eps_u(l,2,1)
     &                     +viscous(l,2,2,ri*3-1) * eps_u(l,2,2)
     &                     +viscous(l,2,3,ri*3-1) * eps_u(l,2,3)
     &                     +viscous(l,3,1,ri*3-1) * eps_u(l,3,1)
     &                     +viscous(l,3,2,ri*3-1) * eps_u(l,3,2)
     &                     +viscous(l,3,3,ri*3-1) * eps_u(l,3,3) )
          eforce(l,riz) = eforce(l,riz) - twovisc * fac(l) *
     &                    ( viscous(l,1,1,ri*3) * eps_u(l,1,1)
     &                     +viscous(l,1,2,ri*3) * eps_u(l,1,2)
     &                     +viscous(l,1,3,ri*3) * eps_u(l,1,3)
     &                     +viscous(l,2,1,ri*3) * eps_u(l,2,1)
     &                     +viscous(l,2,2,ri*3) * eps_u(l,2,2)
     &                     +viscous(l,2,3,ri*3) * eps_u(l,2,3)
     &                     +viscous(l,3,1,ri*3) * eps_u(l,3,1)
     &                     +viscous(l,3,2,ri*3) * eps_u(l,3,2)
     &                     +viscous(l,3,3,ri*3) * eps_u(l,3,3) )
C         pressure forces integrated by parts
          eforce(l,rix) = eforce(l,rix) 
     &                  + presint(l) * derxy(l,ri,1) * fac(l)
          eforce(l,riy) = eforce(l,riy)
     &                  + presint(l) * derxy(l,ri,2) * fac(l)
          eforce(l,riz) = eforce(l,riz)
     &                  + presint(l) * derxy(l,ri,3) * fac(l)
C         stabilisation part - impulse stabilisation
          eforce(l,rix) = eforce(l,rix)
     &                  + tau_M(l) * conv_c(l,ri) * resid(l,1)
          eforce(l,riy) = eforce(l,riy)
     &                  + tau_M(l) * conv_c(l,ri) * resid(l,2)
          eforce(l,riz) = eforce(l,riz)
     &                  + tau_M(l) * conv_c(l,ri) * resid(l,3)
C 
          eforce(l,rix) = eforce(l,rix) + tau_Mp(l) * twovisc *
     &                  ( viscs2(l,1,ri*3-2) * resid(l,1)
     &                   +viscs2(l,2,ri*3-2) * resid(l,2)
     &                   +viscs2(l,3,ri*3-2) * resid(l,3) )
          eforce(l,riy) = eforce(l,riy) + tau_Mp(l) * twovisc *
     &                  ( viscs2(l,1,ri*3-1) * resid(l,1)
     &                   +viscs2(l,2,ri*3-1) * resid(l,2)
     &                   +viscs2(l,3,ri*3-1) * resid(l,3) ) 
          eforce(l,riz) = eforce(l,riz) + tau_Mp(l) * twovisc *
     &                  ( viscs2(l,1,ri*3) * resid(l,1)
     &                   +viscs2(l,2,ri*3) * resid(l,2)
     &                   +viscs2(l,3,ri*3) * resid(l,3) )
C         stabilisation part - continuity stabilistation
          eforce(l,rix) = eforce(l,rix) - tau_C(l) * timefac *
     &                   ( vderxy(l,1,1)
     &                   + vderxy(l,2,2)
     &                   + vderxy(l,3,3) ) * derxy(l,ri,1)
          eforce(l,riy) = eforce(l,riy) - tau_C(l) * timefac *
     &                   ( vderxy(l,1,1)
     &                   + vderxy(l,2,2)
     &                   + vderxy(l,3,3) ) * derxy(l,ri,2)
          eforce(l,riz) = eforce(l,riz) - tau_C(l) * timefac *
     &                   ( vderxy(l,1,1)
     &                   + vderxy(l,2,2)
     &                   + vderxy(l,3,3) ) * derxy(l,ri,3)
C } /* end loop over rows */
        enddo ! the loop
      enddo   ! loop over rows
      return
      end subroutine ! f3fcalresvec
