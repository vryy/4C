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
class LeastSquares
{
        public:
                LeastSquares(std::vector<std::vector<double> >matri,Epetra_SerialDenseVector sourc)
                        :matri_(matri),
                         sourc_(sourc)
                {
                }
                ~LeastSquares(){};
                Epetra_SerialDenseVector linear_least_square();

        private:
                std::vector<std::vector<double> > matri_;
                Epetra_SerialDenseVector sourc_;
                Epetra_SerialDenseVector unknown_;
                Epetra_SerialDenseMatrix get_square_matrix(Epetra_SerialDenseVector &rhs);
                std::vector<double> ConjugateGradient(std::vector<std::vector<double> >coeff, std::vector<double> rhs);
                std::vector<double> multiply(std::vector<std::vector<double> > mat, std::vector<double> ve);
                double multi_vec(std::vector<double> mm1, std::vector<double> mm2);
		double maxAbsolute(std::vector<double>a);
};
}
}

#endif
