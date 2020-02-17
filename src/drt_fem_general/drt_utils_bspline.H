/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of B-splines

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/


#ifndef BSPLINE_H
#define BSPLINE_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "../linalg/linalg_serialdensevector.H"


namespace DRT
{
  namespace NURBS
  {
    namespace UTILS
    {
      /*!
            \class BsplinePolynomial
      \brief A class which allows to evaluate all bspline
             polynomials which are nonzero in the 'central'
             knotvector interval
      */
      class BsplinePolynomial
      {
       public:
        //--------------------------------------------------
        /*! \brief constructor

            \param degree           int (i)
                                    degree of polynomial, for
                                    consistency checks

            \param local_knotvector Epetra_SerialDenseVector (i)
                                    knot range that contains
                                    the compact support of
                                    all bsplines which are
                                    not zero in x


            The knot support looks like this:

            \verbatim

            |<----degree----->|     |<----degree----->|
            |                 |     |                 |
            |                 |     |                 |

            +-----+-----+-----+-----+-----+-----+-----+

                              |     |
                              |     |
                              |<--->|

                        element containing x

            \endverbatim

        */
        //--------------------------------------------------
        BsplinePolynomial(const int degree, const Epetra_SerialDenseVector local_knotvector);
        //--------------------------------------------------
        //! \brief destructor
        //--------------------------------------------------
        virtual ~BsplinePolynomial();
        //--------------------------------------------------
        //! \brief copy constructor
        //!
        //! \param old const BsplinePolynomial (i)
        //!            bspline polynomial to copy
        //--------------------------------------------------
        BsplinePolynomial(const BsplinePolynomial& old);

        //--------------------------------------------------
        /*! \brief compute a bspline value at position x

         \param bspline_value double (o)
                              the bspline value computed

         \param  x            double (i)
                              position to evaluate bspline

         \param  ldofid       int    (i)
                              which bspline to evaluate,
                              ldofid is the local dof id,
                              i.e. the number of the
                              bspline we want to evaluate.
                              0 is the leftmost nonzero
                              bspline, degree the right
                              most. All of them have
                              degree+1 nonzero intervals


         example for the use of ldofid: (degree is 3)

         \verbatim
                                ^
                     ****       ^        +-----------+
                    *    *      ^        | ldofid==0 |
                   *      *     ^        +-----------+
                 **        **   ^
              ***            ***^
          +***---+-----+-----+--***+-----+-----+-----+
                                ^
                                ^
                           **** ^        +-----------+
                          *    *^        | ldofid==1 |
                         *      *        +-----------+
                       **       ^**
                    ***         ^  ***
          +-----+***---+-----+-----+--***+-----+-----+
                                ^
                                ^
          +-----------+         ^****
          | ldofid==2 |         *    *
          +-----------+        *^     *
                             ** ^      **
                          ***   ^        ***
          +-----+-----+***---+-----+-----+--***+-----+
                                ^
                                ^
          +-----------+         ^      ****
          | ldofid==3 |         ^     *    *
          +-----------+         ^    *      *
                                ^  **        **
                                ***            ***
          +-----+-----+-----+***---+-----+-----+--***+
                                ^
                                ^
                                x

         \endverbatim


         Recursion formula used:

         \verbatim

                  x - x                x     - x
          p            i     p-1        i+p+1         p-1
         N (x) = -------- * N   (x) + ------------ * N   (x)
          i      x   - x     i        x     - x       i+1
                  i+p   i              i+p+1   i+1

         \endverbatim

         The recursion is started from the following functions
         (depending on ldofid).

         \verbatim

                   +-
                   |
            0      |  1   for x contained in interval i
           N (x) = |
            i      |  0   otherwise
                   |
                   +-

         \endverbatim

         In the example above (degree 3), we have:

         \verbatim


          ldofid = 0:
                            +-----+
                            |     |
                            |     |
                            |     |
          +-----+-----+-----+--x--+



          ldofid = 1:
                            +-----+
                            |     |
                            |     |
                            |     |
                +-----+-----+--x--+-----+


          ldofid = 2:
                            +-----+
                            |     |
                            |     |
                            |     |
                      +-----+--x--+-----+-----+


          ldofid = 3:
                            +-----+
                            |     |
                            |     |
                            |     |
                            +--x--+-----+-----+-----+

                            |     |
                            |     |
                            |<--->|

                      element containing x

         \endverbatim



         the recursion is done in Aitken-Neville-style:

         \verbatim
                |        |        |        |
                | rr==0  | rr==1  | rr==2  |
                |        |        |        |

              N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
                |       /|       /|       /
                |      / |      / |      /
                |     /  |     /  |     /
                |    /   |    /   |    /
                |   /    |   /    |   /
                |  /     |  /     |  /
                | /      | /      | /
              N(0,1)  N(1,1)   N(2,1)             p == 1
                |       /|       /
                |      / |      /
                |     /  |     /
                |    /   |    /
                |   /    |   /
                |  /     |  /
                | /      | /
              N(0,2)  N(1,2)                      p == 2
                |       /
                |      /
                |     /
                |    /
                |   /
                |  /
                | /
              N(0,3)                              p == 3

         \endverbatim

         memory is reused, i.e. in the end, N(0,3) is
         contained in bspline[0]

        */
        template <typename scalar_type>
        void EvaluateBspline(scalar_type& bspline_value, const scalar_type x, const int ldofid)
        {
          //                        ^
          //             ****       ^        +-----------+
          //            *    *      ^        | ldofid==0 |
          //           *      *     ^        +-----------+
          //         **        **   ^
          //      ***            ***^
          //  +***---+-----+-----+--***+-----+-----+-----+
          //                        ^
          //                        ^
          //                   **** ^        +-----------+
          //                  *    *^        | ldofid==1 |
          //                 *      *        +-----------+
          //               **       ^**
          //            ***         ^  ***
          //  +-----+***---+-----+-----+--***+-----+-----+
          //                        ^
          //                        ^
          //  +-----------+         ^****
          //  | ldofid==2 |         *    *
          //  +-----------+        *^     *
          //                     ** ^      **
          //                  ***   ^        ***
          //  +-----+-----+***---+-----+-----+--***+-----+
          //                        ^
          //                        ^
          //  +-----------+         ^      ****
          //  | ldofid==3 |         ^     *    *
          //  +-----------+         ^    *      *
          //                        ^  **        **
          //                        ***            ***
          //  +-----+-----+-----+***---+-----+-----+--***+
          //                        ^
          //                        ^
          //                        x

          /*****************************************************************
           *  WARNING: here was a consistency check which was removed       *
           *           in revision 19453. It guaranteed that only           *
           *           evaluations of points located within the element     *
           *           are valid. But, the contact algorithm requires eval. *
           *           outside the element domain due to GP-projections !!! *
           ******************************************************************/

          // define the vector of values at x of all initial polynomials
          std::vector<scalar_type> bspline(degree_ + 1);

          // The nonzero initial bspline polynomial and the intervals
          // that define the compact support of the bspline number lid
          // of given degree:
          //
          //
          //  ldofid = 0:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //  +-----+-----+-----+--x--+
          //
          //  all other initial bsplines are zero at x
          //  (empty intervals indicate the compact support of
          //   the other bspline polynomials which contribute
          //   to the bspline value associated with the given
          //   dof lid. The union of all intervals defines the
          //   support of the computed bspline of degree 3)
          //
          //
          //  ldofid = 1:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //        +-----+-----+--x--+-----+
          //
          //
          //  ldofid = 2:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //              +-----+--x--+-----+-----+
          //
          //
          //  ldofid = 3:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //                    +--x--+-----+-----+-----+
          //
          //                    |     |
          //                    |     |
          //                    |<--->|
          //
          //              element containing x
          //

          for (int rr = 0; rr < degree_ + 1; rr++)
          {
            bspline[rr] = 0;
          }
          bspline[degree_ - ldofid] = 1;

          //        |        |        |        |
          //        | rr==0  | rr==1  | rr==2  |
          //        |        |        |        |
          //
          //      N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
          //        |       /|       /|       /
          //        |      / |      / |      /
          //        |     /  |     /  |     /
          //        |    /   |    /   |    /
          //        |   /    |   /    |   /
          //        |  /     |  /     |  /
          //        | /      | /      | /
          //      N(0,1)  N(1,1)   N(2,1)             p == 1
          //        |       /|       /
          //        |      / |      /
          //        |     /  |     /
          //        |    /   |    /
          //        |   /    |   /
          //        |  /     |  /
          //        | /      | /
          //      N(0,2)  N(1,2)                      p == 2
          //        |       /
          //        |      /
          //        |     /
          //        |    /
          //        |   /
          //        |  /
          //        | /
          //      N(0,3)                              p == 3
          //
          //
          //
          // memory is reused, i.e. in the end, N(0,3) is contained
          // in bspline[0]
          //

          // loop all rows in the upper table
          for (int p = 0; p < degree_; ++p)
          {
// do computation of bspline values of specified degree,
// corresponding to one row in the scheme above
#pragma clang loop vectorize(disable)  // prevent clang from vectorizing
            for (int rr = 0; rr < degree_ - p; ++rr)
            {
              // id of first bspline function of this combination
              int first = ldofid + rr;

              scalar_type fact1;
              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              if (fabs(myknotvector_(first + p + 1) - myknotvector_(first)) < 10e-9)
              {
                fact1 = 0;
              }
              else
              {
                fact1 = (x - myknotvector_(first));
                fact1 /= (myknotvector_(first + p + 1) - myknotvector_(first));
              }

              scalar_type fact2;
              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              if (fabs(myknotvector_(first + p + 2) - myknotvector_(first + 1)) < 10e-9)
              {
                fact2 = 0;
              }
              else
              {
                fact2 = (myknotvector_(first + p + 2) - x);
                fact2 /= (myknotvector_(first + p + 2) - myknotvector_(first + 1));
              }
              // do the actual bspline recursion --- memory is reused!
              bspline[rr] = fact1 * bspline[rr] + fact2 * bspline[rr + 1];
            }
          }

          // set the output
          bspline_value = bspline[0];

          return;
        }

        //--------------------------------------------------
        /*! \brief Compute ldofid's Bspline value at point x
                   In addiditon, compute its first derivative

         \param bsplineval   double   (o)
                             bspline value

         \param bsplineder   double   (o)
                             first derivative of bspline

         \param x            double   (i)
                             position to evaluate bspline

         \param ldofid       int      (i)
                             which bspline to evaluate,
                             ldofid is the local dof id,
                             i.e. the number of the
                             bspline we want to evaluate.
                             0 is the leftmost nonzero
                             bspline, degree the right
                             most. All of them have
                             degree+1 nonzero intervals

         Recursion formulas used:

         \verbatim

                  x - x                x     - x
          p            i     p-1        i+p+1         p-1
         N (x) = -------- * N   (x) + ------------ * N   (x)
          i      x   - x     i        x     - x       i+1
                  i+p   i              i+p+1   i+1

           p          p        p-1           p         p-1
         N'  (x) = -------- * N   (x) - ----------- * N   (x)
           i       x   - x     i        x     - x      i+1
                    i+p   i              i+p+1   i+1

         \endverbatim

         mind that the recursion scheme is changed
         (we have to descent into a branch for the
          first derivatives):

         \verbatim

                |        |        |        |
                | rr==0  | rr==1  | rr==2  |
                |        |        |        |

              N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
                |       /|       /|       /
                |      / |      / |      /
                |     /  |     /  |     /
                |    /   |    /   |    /
                |   /    |   /    |   /
                |  /     |  /     |  /
                | /      | /      | /
              N(0,1)  N(1,1)   N(2,1)             p == 1
                |       /|       /
                |      / |      /
                |     /  |     /
                |    /   |    /
                |   /    |   /
                |  /     |  /
                | /      | /
              N(0,2)  N(1,2)                      p == 2



         ----------------------------------------------------



              N(0,2)  N(1,2)                      p == 2   +--
                |\      /|                                 |
                | \    / |                                 | branch
                |  \  /  |                                 |
                |   \/   |                                 | for first
                |   /\   |                                 |
                |  /  \  |                                 | derivatives
                | /    \ |                                 |
             val(0,3) der(0,3)                    p == 3   +--



         memory is reused on the first level. For the last
         step, we have additional memory to be able to access
         N(0,3) twice


         \endverbatim

        */
        template <typename scalar_type>
        void EvaluateBsplineAndDeriv(
            scalar_type& bsplineval, scalar_type& bsplineder, const scalar_type x, const int ldofid)
        {
          //                        ^
          //             ****       ^        +-----------+
          //            *    *      ^        | ldofid==0 |
          //           *      *     ^        +-----------+
          //         **        **   ^
          //      ***            ***^
          //  +***---+-----+-----+--***+-----+-----+-----+
          //                        ^
          //                        ^
          //                   **** ^        +-----------+
          //                  *    *^        | ldofid==1 |
          //                 *      *        +-----------+
          //               **       ^**
          //            ***         ^  ***
          //  +-----+***---+-----+-----+--***+-----+-----+
          //                        ^
          //                        ^
          //  +-----------+         ^****
          //  | ldofid==2 |         *    *
          //  +-----------+        *^     *
          //                     ** ^      **
          //                  ***   ^        ***
          //  +-----+-----+***---+-----+-----+--***+-----+
          //                        ^
          //                        ^
          //  +-----------+         ^      ****
          //  | ldofid==3 |         ^     *    *
          //  +-----------+         ^    *      *
          //                        ^  **        **
          //                        ***            ***
          //  +-----+-----+-----+***---+-----+-----+--***+
          //                        ^
          //                        ^
          //                        x

          /*****************************************************************
           *  WARNING: here was a consistency check which was removed       *
           *           in revision 19453. It guaranteed that only           *
           *           evaluations of points located within the element     *
           *           are valid. But, the contact algorithm requires eval. *
           *           outside the element domain due to GP-projections !!! *
           ******************************************************************/

          // define the vector of values at x of all initial polynomials
          std::vector<scalar_type> bspline(degree_ + 1);

          // The nonzero initial bspline polynomial and the intervals
          // that define the compact support of the bspline number lid
          // of given degree:
          //
          //
          //  ldofid = 0:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //  +-----+-----+-----+--x--+
          //
          //  all other initial bsplines are zero at x
          //  (empty intervals indicate the compact support of
          //   the other bspline polynomials which contribute
          //   to the bspline value associated with the given
          //   dof lid. The union of all intervals defines the
          //   support of the computed bspline of degree 3)
          //
          //
          //  ldofid = 1:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //        +-----+-----+--x--+-----+
          //
          //
          //  ldofid = 2:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //              +-----+--x--+-----+-----+
          //
          //
          //  ldofid = 3:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //                    +--x--+-----+-----+-----+
          //
          //                    |     |
          //                    |     |
          //                    |<--->|
          //
          //              element containing x
          //

          // initial values for recursion
          //
          //           +-
          //           |
          //    0      |  1   for x contained in interval i
          //   N (x) = |
          //    i      |  0   otherwise
          //           |
          //           +-
          //

          for (int rr = 0; rr < degree_ + 1; rr++)
          {
            bspline[rr] = 0;
          }
          bspline[degree_ - ldofid] = 1;

          //        |        |        |        |
          //        | rr==0  | rr==1  | rr==2  |
          //        |        |        |        |
          //
          //      N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
          //        |       /|       /|       /
          //        |      / |      / |      /
          //        |     /  |     /  |     /
          //        |    /   |    /   |    /
          //        |   /    |   /    |   /
          //        |  /     |  /     |  /
          //        | /      | /      | /
          //      N(0,1)  N(1,1)   N(2,1)             p == 1
          //        |       /|       /
          //        |      / |      /
          //        |     /  |     /
          //        |    /   |    /
          //        |   /    |   /
          //        |  /     |  /
          //        | /      | /
          //      N(0,2)  N(1,2)                      p == 2
          //
          //
          //
          // ----------------------------------------------------
          //
          //
          //
          //      N(0,2)  N(1,2)                      p == 2   +--
          //        |\      /|                                 |
          //        | \    / |                                 | branch
          //        |  \  /  |                                 |
          //        |   \/   |                                 | for first
          //        |   /\   |                                 |
          //        |  /  \  |                                 | derivatives
          //        | /    \ |                                 |
          //     val(0,3) der(0,3)                    p == 3   +--
          //
          //
          //
          // memory is reused on the first level. For the last
          // step, we have additional memory tobe able to access
          // N(0,3) twice
          //

          // loop all rows in the upper table up to the last
          // but one. Both arguments are still required to compute
          // the derivatives, so do not throw them away
          // (or overwrite)
          for (int p = 0; p < degree_ - 1; ++p)
          {
// do computation of bspline values of specified degree,
// corresponding to one row in the scheme above
#pragma clang loop vectorize(disable)  // prevent clang from vectorizing
            // clang-3.5 runs into an FPE e.g. for the test case contact2D_nurbs9_dual_consistent-p1
            // if vectorization is enabled due to "0.0/0.0" (whose result does not propagate due
            // to binary operations on vectorized data)
            for (int rr = 0; rr < degree_ - p; ++rr)
            {
              // id of first bspline function of this combination
              int i = ldofid + rr;

              // recursion for the computation of the basis
              // function
              //
              //          x - x                x     - x
              //  p            i     p-1        i+p+1         p-1
              // N (x) = -------- * N   (x) + ------------ * N   (x)
              //  i      x   - x     i        x     - x       i+1
              //          i+p   i              i+p+1   i+1
              //
              //        |        |           |            |
              //        +--------+           +------------+
              //           fact1                  fact2
              //

              scalar_type fact1;
              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              if (fabs(myknotvector_(i + p + 1) - myknotvector_(i)) < 10e-9)
              {
                fact1 = 0;
              }
              else
              {
                fact1 = (x - myknotvector_(i));
                fact1 /= (myknotvector_(i + p + 1) - myknotvector_(i));
              }

              scalar_type fact2;
              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              if (fabs(myknotvector_(i + p + 2) - myknotvector_(i + 1)) < 10e-9)
              {
                fact2 = 0;
              }
              else
              {
                fact2 = (myknotvector_(i + p + 2) - x);
                fact2 /= (myknotvector_(i + p + 2) - myknotvector_(i + 1));
              }
              // do the actual bspline recursion --- memory is reused!
              bspline[rr] = fact1 * bspline[rr] + fact2 * bspline[rr + 1];
            }
          }

          //---------------------------------------------------
          // do computation of bspline value in the last level
          // corresponding to one row in the scheme above

          scalar_type fact1;
          // the first part of the if statement allows to
          // enforce interpolation using multiple nodes
          // the second part is a part of the bspline recursion
          if (fabs(myknotvector_(ldofid + degree_) - myknotvector_(ldofid)) < 10e-9)
          {
            fact1 = 0;
          }
          else
          {
            fact1 = (x - myknotvector_(ldofid));
            fact1 /= (myknotvector_(ldofid + degree_) - myknotvector_(ldofid));
          }

          scalar_type fact2;
          // the first part of the if statement allows to
          // enforce interpolation using multiple nodes
          // the second part is a part of the bspline recursion
          if (fabs(myknotvector_(ldofid + degree_ + 1) - myknotvector_(ldofid + 1)) < 10e-9)
          {
            fact2 = 0;
          }
          else
          {
            fact2 = (myknotvector_(ldofid + degree_ + 1) - x);
            fact2 /= (myknotvector_(ldofid + degree_ + 1) - myknotvector_(ldofid + 1));
          }
          // do the actual bspline recursion --- memory is reused!
          bsplineval = fact1 * bspline[0] + fact2 * bspline[1];

          //---------------------------------------------------
          // do computation of bspline derivatives from the
          // last level corresponding to one row in the scheme
          // above

          if (degree_ > 0)
          {
            //
            //
            //   p          p        p-1           p         p-1
            // N'  (x) = -------- * N   (x) - ----------- * N   (x)
            //   i       x   - x     i        x     - x      i+1
            //            i+p   i              i+p+1   i+1
            //
            //          |        |           |            |
            //          +--------+           +------------+
            //             fact1                  fact2

            // the first part of the if statement allows to
            // enforce interpolation using multiple nodes
            // the second part computes fact1 in the equation
            // above
            if (fabs(myknotvector_(ldofid + degree_) - myknotvector_(ldofid)) < 10e-9)
            {
              fact1 = 0;
            }
            else
            {
              fact1 = degree_;
              fact1 /= (myknotvector_(ldofid + degree_) - myknotvector_(ldofid));
            }

            // the first part of the if statement allows to
            // enforce interpolation using multiple nodes
            // the second part is a part of the bspline recursion
            // see above
            if (fabs(myknotvector_(ldofid + degree_ + 1) - myknotvector_(ldofid + 1)) < 10e-9)
            {
              fact2 = 0;
            }
            else
            {
              fact2 = degree_;
              fact2 /= (myknotvector_(ldofid + degree_ + 1) - myknotvector_(ldofid + 1));
            }

            // compute the actual bspline derivative formula
            bsplineder = fact1 * bspline[0] - fact2 * bspline[1];
          }
          else
          {
            // piecewise constants get 0 derivatives
            bsplineder = 0;
          }

          return;
        }

        //--------------------------------------------------
        /*! \brief Compute ldofid's Bspline value at point x
                   In addiditon, compute its first and
                   second derivative

        \param bsplineval  double  (o)
                           bspline value

        \param bsplineder  double  (o)
                           first derivative of bspline

        \param bsplineder2 double  (o)
                           second derivative of bspline

        \param x           double  (i)
                           position to evaluate bspline

        \param ldofid      int     (i)
                           which bspline to evaluate,
                           ldofid is the local dof id,
                           i.e. the number of the
                           bspline we want to evaluate.
                           0 is the leftmost nonzero
                           bspline, degree the right
                           most. All of them have
                           degree+1 nonzero intervals

         (Recursion) formulas:

         \verbatim

         Recursion for Basis functions:

                  x - x                x     - x
          p            i     p-1        i+p+1         p-1
         N (x) = -------- * N   (x) + ------------ * N   (x)
          i      x   - x     i        x     - x       i+1
                  i+p   i              i+p+1   i+1

         Computation of first derivatives from basis functions
         of the previous level:

           p          p        p-1           p         p-1
         N'  (x) = -------- * N   (x) - ----------- * N   (x)
           i       x   - x     i        x     - x      i+1
                    i+p   i              i+p+1   i+1

         Computation of second derivatives:

        -------------------------------------------------------

           p          p        p-1           p         p-1
         N'' (x) = -------- * N'  (x) - ----------- * N'  (x)
           i       x   - x     i        x     - x      i+1
                    i+p   i              i+p+1   i+1

         with

           p-1         p-1       p-2          p-1        p-2
         N'  (x) = ---------- * N   (x) - ----------- * N   (x)
           i       x     - x     i        x   - x        i+1
                    i+p-1   i              i+p   i+1

         \endverbatim


         \verbatim


                |        |        |        |
                | rr==0  | rr==1  | rr==2  |
                |        |        |        |

              N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
                |       /|       /|       /
                |      / |      / |      /
                |     /  |     /  |     /
                |    /   |    /   |    /
                |   /    |   /    |   /
                |  /     |  /     |  /
                | /      | /      | /
              N(0,1)  N(1,1)   N(2,1)             p == 1



         ====================================================



              N(0,1)  N(1,1)   N(2,1)             p == 1   +--
                |       /|       /                         |
                |      / |      /                          | branch
                |     /  |     /                           |
                |    /   |    /                            | for second
                |   /    |   /                             |
                |  /     |  /                              | derivatives
                | /      | /                               |
              N(0,2)  N(1,2)                      p == 2   +--



         ====================================================



              N(0,2)  N(1,2)                      p == 2   +--
                |\      /|                                 |
                | \    / |                                 | branch
                |  \  /  |                                 |
                |   \/   |                                 | for first
                |   /\   |                                 |
                |  /  \  |                                 | derivatives
                | /    \ |                                 |
             val(0,3) der(0,3)                    p == 3   +--

         \endverbatim

        */
        template <typename scalar_type>
        void EvaluateBsplineFirstAndSecondDeriv(scalar_type& bsplineval, scalar_type& bsplineder,
            scalar_type& bsplineder2, const scalar_type x, const int ldofid)
        {
          /*

                    x - x                x     - x
            p            i     p-1        i+p+1         p-1
           N (x) = -------- * N   (x) + ------------ * N   (x)
            i      x   - x     i        x     - x       i+1
                    i+p   i              i+p+1   i+1

             p          p        p-1           p         p-1
           N'  (x) = -------- * N   (x) - ----------- * N   (x)
             i       x   - x     i        x     - x      i+1
                      i+p   i              i+p+1   i+1

          -------------------------------------------------------

             p          p        p-1           p         p-1
           N'' (x) = -------- * N'  (x) - ----------- * N'  (x)
             i       x   - x     i        x     - x      i+1
                      i+p   i              i+p+1   i+1

             p-1         p-1       p-2          p-1        p-2
           N'  (x) = ---------- * N   (x) - ----------- * N   (x)
             i       x     - x     i        x   - x        i+1
                      i+p-1   i              i+p   i+1


          */

          //                        ^
          //             ****       ^        +-----------+
          //            *    *      ^        | ldofid==0 |
          //           *      *     ^        +-----------+
          //         **        **   ^
          //      ***            ***^
          //  +***---+-----+-----+--***+-----+-----+-----+
          //                        ^
          //                        ^
          //                   **** ^        +-----------+
          //                  *    *^        | ldofid==1 |
          //                 *      *        +-----------+
          //               **       ^**
          //            ***         ^  ***
          //  +-----+***---+-----+-----+--***+-----+-----+
          //                        ^
          //                        ^
          //  +-----------+         ^****
          //  | ldofid==2 |         *    *
          //  +-----------+        *^     *
          //                     ** ^      **
          //                  ***   ^        ***
          //  +-----+-----+***---+-----+-----+--***+-----+
          //                        ^
          //                        ^
          //  +-----------+         ^      ****
          //  | ldofid==3 |         ^     *    *
          //  +-----------+         ^    *      *
          //                        ^  **        **
          //                        ***            ***
          //  +-----+-----+-----+***---+-----+-----+--***+
          //                        ^
          //                        ^
          //                        x

          /*****************************************************************
           *  WARNING: here was a consistency check which was removed       *
           *           in revision 19453. It guaranteed that only           *
           *           evaluations of points located within the element     *
           *           are valid. But, the contact algorithm requires eval. *
           *           outside the element domain due to GP-projections !!! *
           ******************************************************************/

          // The nonzero initial bspline polynomial and the intervals
          // that define the compact support of the bspline number lid
          // of given degree:
          //
          //
          //  ldofid = 0:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //  +-----+-----+-----+--x--+
          //
          //  all other initial bsplines are zero at x
          //  (empty intervals indicate the compact support of
          //   the other bspline polynomials which contribute
          //   to the bspline value associated with the given
          //   dof lid. The union of all intervals defines the
          //   support of the computed bspline of degree 3)
          //
          //
          //  ldofid = 1:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //        +-----+-----+--x--+-----+
          //
          //
          //  ldofid = 2:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //              +-----+--x--+-----+-----+
          //
          //
          //  ldofid = 3:
          //                    +-----+
          //                    |     |
          //                    |     |
          //                    |     |
          //                    +--x--+-----+-----+-----+
          //
          //                    |     |
          //                    |     |
          //                    |<--->|
          //
          //              element containing x
          //

          // initial values for recursion
          //
          //           +-
          //           |
          //    0      |  1   for x contained in interval i
          //   N (x) = |
          //    i      |  0   otherwise
          //           |
          //           +-
          //

          for (int rr = 0; rr < degree_plus_one_; rr++)
          {
            bspline_[rr] = 0;
          }
          bspline_[degree_ - ldofid] = 1;

          //        |        |        |        |
          //        | rr==0  | rr==1  | rr==2  |
          //        |        |        |        |
          //
          //      N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
          //        |       /|       /|       /
          //        |      / |      / |      /
          //        |     /  |     /  |     /
          //        |    /   |    /   |    /
          //        |   /    |   /    |   /
          //        |  /     |  /     |  /
          //        | /      | /      | /
          //      N(0,1)  N(1,1)   N(2,1)             p == 1
          //
          //
          //
          // ====================================================
          //
          //
          //
          //      N(0,1)  N(1,1)   N(2,1)             p == 1   +--
          //        |       /|       /                         |
          //        |      / |      /                          | branch
          //        |     /  |     /                           |
          //        |    /   |    /                            | for second
          //        |   /    |   /                             |
          //        |  /     |  /                              | derivatives
          //        | /      | /                               |
          //      N(0,2)  N(1,2)                      p == 2   +--
          //
          //
          //
          // ====================================================
          //
          //
          //
          //      N(0,2)  N(1,2)                      p == 2   +--
          //        |\      /|                                 |
          //        | \    / |                                 | branch
          //        |  \  /  |                                 |
          //        |   \/   |                                 | for first
          //        |   /\   |                                 |
          //        |  /  \  |                                 | derivatives
          //        | /    \ |                                 |
          //     val(0,3) der(0,3)                    p == 3   +--
          //
          //
          //
          // memory is reused on the first level. For the last
          // step, we have additional memory to be able to access
          // N(0,3) twice
          //


          // loop all rows in the upper table up to the last
          // but one. Both arguments are still required to compute
          // the derivatives, so do not throw them away
          // (or overwrite)
          for (int p = 0; p < degree_ - 2; ++p)
          {
// do computation of bspline values of specified degree,
// corresponding to one row in the scheme above
#pragma clang loop vectorize(disable)  // prevent clang from vectorizing
            for (int rr = 0; rr < degree_ - p; ++rr)
            {
              // id of first bspline function of this combination
              const int i = ldofid + rr;

              // recursion for the computation of the basis
              // function
              //
              //          x - x                x     - x
              //  p            i     p-1        i+p+1         p-1
              // N (x) = -------- * N   (x) + ------------ * N   (x)
              //  i      x   - x     i        x     - x       i+1
              //          i+p   i              i+p+1   i+1
              //
              //        |        |           |            |
              //        +--------+           +------------+
              //         fact_[0]               fact_[1]
              //

              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              double dx = myknotvector_(i + p + 1) - myknotvector_(i);

              if (fabs(dx) < 10e-9)
              {
                fact_[0] = 0;
              }
              else
              {
                fact_[0] = (x - myknotvector_(i)) / dx;
              }

              dx = myknotvector_(i + p + 2) - myknotvector_(i + 1);
              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              if (fabs(dx) < 10e-9)
              {
                fact_[1] = 0;
              }
              else
              {
                fact_[1] = (myknotvector_(i + p + 2) - x) / dx;
              }
              // do the actual bspline recursion --- memory is reused!
              bspline_[rr] = fact_[0] * bspline_[rr] + fact_[1] * bspline_[rr + 1];
            }
          }

          //
          // ====================================================
          //

          //---------------------------------------------------
          // do computation of both bspline derivatives
          // from the p-1 level

          if (degree_ > 1)
          {
            for (int rr = 0; rr < 2; ++rr)
            {
              // id of first bspline function of this combination
              const int i = ldofid + rr;
              const int i_plus_degree = i + degree_;

              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part computes fact_[0] in the equation
              // above

              const double dxi = myknotvector_(i_plus_degree - 1) - myknotvector_(i);

              if (fabs(dxi) < 10e-9)
              {
                fact_[2] = 0;
                fact_[0] = 0;
              }
              else
              {
                fact_[2] = (degree_ - 1) / dxi;
                fact_[0] = (x - myknotvector_(i)) / dxi;
              }

              const double dxip = myknotvector_(i_plus_degree) - myknotvector_(i + 1);

              // the first part of the if statement allows to
              // enforce interpolation using multiple nodes
              // the second part is a part of the bspline recursion
              // see above
              if (fabs(dxip) < 10e-9)
              {
                fact_[3] = 0;
                fact_[1] = 0;
              }
              else
              {
                fact_[3] = (degree_ - 1) / dxip;
                fact_[1] = (myknotvector_(i_plus_degree) - x) / dxip;
              }

              //
              //
              //
              //     p-1         p-1       p-2          p-1        p-2
              //   N'  (x) = ---------- * N   (x) - ----------- * N   (x)
              //     i       x     - x     i        x   - x        i+1
              //              i+p-1   i              i+p   i+1
              //             |        |           |            |
              //             +--------+           +------------+
              //              fact_[2]               fact_[3]

              // compute the actual bspline derivative formula
              pmo_deriv_[rr] = fact_[2] * bspline_[rr] - fact_[3] * bspline_[rr + 1];
              //----------------------------------------------
              // bspline[rr] for this level will be destroyed
              // NOW and is replaced by the next levels value!
              //----------------------------------------------

              // recursion for the computation of the basis
              // function
              //
              //            x - x                  x   - x
              //  p-1            i       p-2        i+p           p-2
              // N   (x) = ---------- * N   (x) + ------------ * N   (x)
              //  i        x     - x     i         x   - x       i+1
              //            i+p-1   i               i+p   i+1
              //
              //           |        |            |            |
              //           +--------+            +------------+
              //            fact_[0]                fact_[1]
              //

              // do the actual bspline recursion --- memory is reused!
              bspline_[rr] = fact_[0] * bspline_[rr] + fact_[1] * bspline_[rr + 1];
            }
          }
          else
          {
            // piecewise constants get 0 derivatives
            pmo_deriv_[0] = 0;
            pmo_deriv_[1] = 0;
          }

          //
          // ====================================================
          //

          //---------------------------------------------------
          // do computation of bspline value in the last level

          // the first part of the if statement allows to
          // enforce interpolation using multiple nodes
          // the second part is a part of the bspline recursion
          const double dxilast = myknotvector_(ldofid + degree_) - myknotvector_(ldofid);

          if (fabs(dxilast) < 10e-9)
          {
            fact_[0] = 0;
            fact_[2] = 0;
          }
          else
          {
            fact_[0] = (x - myknotvector_(ldofid)) / dxilast;
            fact_[2] = degree_ / dxilast;
          }

          // the first part of the if statement allows to
          // enforce interpolation using multiple nodes
          // the second part is a part of the bspline recursion
          const double dxiplast =
              myknotvector_(ldofid + degree_plus_one_) - myknotvector_(ldofid + 1);

          if (fabs(dxiplast) < 10e-9)
          {
            fact_[1] = 0;
            fact_[3] = 0;
          }
          else
          {
            fact_[1] = (myknotvector_(ldofid + degree_plus_one_) - x) / dxiplast;
            fact_[3] = degree_ / dxiplast;
          }
          // do the actual bspline recursion --- memory is reused!
          bsplineval = fact_[0] * bspline_[0] + fact_[1] * bspline_[1];

          //---------------------------------------------------
          // evaluate the second derivatives

          //   p          p        p-1           p         p-1
          // N'' (x) = -------- * N'  (x) - ----------- * N'  (x)
          //   i       x   - x     i        x     - x      i+1
          //            i+p   i              i+p+1   i+1
          //
          //          |        |           |            |
          //          +--------+           +------------+
          //           fact_[2]               fact_[3]

          bsplineder2 = fact_[2] * pmo_deriv_[0] - fact_[3] * pmo_deriv_[1];

          //---------------------------------------------------
          // do computation of bspline derivatives from the
          // last level corresponding to one row in the scheme
          // above

          if (degree_ > 0)
          {
            //
            //
            //   p          p        p-1           p         p-1
            // N'  (x) = -------- * N   (x) - ----------- * N   (x)
            //   i       x   - x     i        x     - x      i+1
            //            i+p   i              i+p+1   i+1
            //
            //          |        |           |            |
            //          +--------+           +------------+
            //           fact_[2]               fact_[3]

            // compute the actual bspline derivative formula
            bsplineder = fact_[2] * bspline_[0] - fact_[3] * bspline_[1];
          }
          else
          {
            // piecewise constants get 0 derivatives
            bsplineder = 0;
          }

          return;
        }

        //--------------------------------------------------
        //! \brief Print some information on the bspline
        //--------------------------------------------------
        void PrintBspline()
        {
          std::cout << "------------------\n";
          std::cout << "Bspline knotvector\n";
          std::cout << myknotvector_;
          std::cout << std::endl;
          std::cout << "Bspline degree: " << degree_ << "\n";
          std::cout << std::endl;
          std::cout << "Allows to evaluate bsplines at positions ";
          std::cout << "in interval [";
          std::cout << myknotvector_(degree_);
          std::cout << ",";
          std::cout << myknotvector_(degree_ + 1);
          std::cout << "]\n";

          return;
        };

       private:
        //! the part of the knotvector where bspline
        //! polynomials of the center-interval are non-zero
        Epetra_SerialDenseVector myknotvector_;

        //! a working array for the construction of bsplines
        std::vector<double> bspline_;

        //! the degree of the bspline-polynomial under
        //! consideration
        int degree_;
        //! the degree of the bspline-polynomial under
        //! consideration plus 1
        int degree_plus_one_;

        //! temporary doubles
        double fact_[4];

        //! both bspline derivatives from the p-1 level
        double pmo_deriv_[2];

        //--------------------------------------------------
        //! \brief throw an errormessage
        //--------------------------------------------------
        void Throwerror(const std::string errormessage);
      };

    }  // namespace UTILS
  }    // namespace NURBS
}  // namespace DRT

#endif
