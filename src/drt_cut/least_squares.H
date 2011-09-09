#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include<vector>
namespace GEO
{
        namespace CUT
        {
class LeastSquares
{
        public:
                LeastSquares(std::vector<std::vector<double> >matri,std::vector<double> sourc)
                        :matri_(matri),
                         sourc_(sourc)
                {
                }
                ~LeastSquares(){};
                std::vector<double> linear_least_square();

        private:
                std::vector<std::vector<double> > matri_;
                std::vector<double> sourc_;
                std::vector<double> unknown_;
                std::vector<std::vector<double> > get_square_matrix(std::vector<double> &rhs);
                std::vector<double> ConjugateGradient(std::vector<std::vector<double> >coeff, std::vector<double> rhs);
                std::vector<double> multiply(std::vector<std::vector<double> > mat, std::vector<double> ve);
                double multi_vec(std::vector<double> mm1, std::vector<double> mm2);
};
}
}

#endif
