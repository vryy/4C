/*!----------------------------------------------------------------------
\file drt_utils_bspline.cpp

\brief Specification of B-splines

<pre>
\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_utils_bspline.H"
#include "../drt_lib/drt_dserror.H"

//--------------------------------------------------
// Constructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(
    const int degree, const Epetra_SerialDenseVector local_knotvector)
    : myknotvector_(local_knotvector),
      bspline_(degree + 1),
      degree_(degree),
      degree_plus_one_(degree + 1)
{
  return;
}

//--------------------------------------------------
// Copy constructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(const BsplinePolynomial& old)
    : degree_(old.degree_)
{
  myknotvector_ = old.myknotvector_;
  return;
}

//--------------------------------------------------
// Destructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::~BsplinePolynomial() { return; }

//--------------------------------------------------
// Compute ldofid's Bspline value at point x
//--------------------------------------------------
void DRT::NURBS::UTILS::BsplinePolynomial::EvaluateBspline(
    double& bspline_value, const double x, const int ldofid)
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
  std::vector<double> bspline(degree_ + 1);

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

      double fact1;
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

      double fact2;
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
// Compute ldofid's Bspline value at point x
// In addiditon, compute its first derivative
//--------------------------------------------------
void DRT::NURBS::UTILS::BsplinePolynomial::EvaluateBsplineAndDeriv(
    double& bsplineval, double& bsplineder, const double x, const int ldofid)
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
  std::vector<double> bspline(degree_ + 1);

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

      double fact1;
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

      double fact2;
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

  double fact1;
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

  double fact2;
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
// Compute ldofid's Bspline value at point x
// In addiditon, compute its first and second
// derivative
//--------------------------------------------------
void DRT::NURBS::UTILS::BsplinePolynomial::EvaluateBsplineFirstAndSecondDeriv(
    double& bsplineval, double& bsplineder, double& bsplineder2, const double x, const int ldofid)
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
  const double dxiplast = myknotvector_(ldofid + degree_plus_one_) - myknotvector_(ldofid + 1);

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

void DRT::NURBS::UTILS::BsplinePolynomial::Throwerror(const std::string errormessage)
{
  // give some information on bspline
  PrintBspline();

  // and the throw the error and exit with a
  // sigsegv
  dserror(errormessage);

  return;
}
