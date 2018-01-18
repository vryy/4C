/*---------------------------------------------------------------------*/
/*!
\file least_squares.H

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include<vector>
#include "Epetra_SerialDenseSolver.h"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
namespace GEO
{
  namespace CUT
  {
  /*!
  \brief This class solves the system of equations using linear least squares method
  */
    class LeastSquares
    {
      public:
        LeastSquares(std::vector<std::vector<double> >matri,Epetra_SerialDenseVector sourc)
                :matri_(matri),
                 sourc_(sourc)
        {
        }
        ~LeastSquares(){};

        /*!
        \brief Performs the linear least squares. The resulting equations are solved using non-iterative method
        because the least squares produces full matrix
        */
        Epetra_SerialDenseVector linear_least_square();

      private:

        /*!
        \brief Construct the square matrix by multiplying the transpose of the original system matrix
        */
        Epetra_SerialDenseMatrix get_square_matrix(Epetra_SerialDenseVector &rhs);
        /*std::vector<double> ConjugateGradient(std::vector<std::vector<double> >coeff, std::vector<double> rhs);
        std::vector<double> multiply(std::vector<std::vector<double> > mat, std::vector<double> ve);
        double multi_vec(std::vector<double> mm1, std::vector<double> mm2);
        double maxAbsolute(std::vector<double>a);*/

        std::vector<std::vector<double> > matri_;
        Epetra_SerialDenseVector sourc_;
        Epetra_SerialDenseVector unknown_;
    };
  }
}

#endif
