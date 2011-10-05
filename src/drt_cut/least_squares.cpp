#include "least_squares.H"
#include <iostream>
#include <cmath>

//premultiplying the matrix with its transpose to get the square matrix
//the source terms also get multiplied
std::vector<std::vector<double> > GEO::CUT::LeastSquares::get_square_matrix(std::vector<double> &rhs)
{
        std::vector<std::vector<double> > sqr;
        sqr.resize(matri_[0].size());
        rhs.resize(matri_[0].size());
        for(unsigned i=0;i<matri_[0].size();i++)
                sqr[i].resize(matri_[0].size());

        for(unsigned i=0;i<matri_[0].size();i++)
        {
                for(unsigned j=0;j<matri_[0].size();j++)
                {
                        sqr[i][j] = 0.0;
                        for(unsigned k=0;k<matri_.size();k++)
                                sqr[i][j] += matri_[k][i]*matri_[k][j]; 
                }
        }

/*      for(unsigned i=0;i<matri_[0].size();i++)
        {
                for(unsigned j=0;j<matri_[0].size();j++)
                        std::cout<<sqr[i][j]<<"\t";
                std::cout<<"\n";
        }*/

        for(unsigned i=0;i<matri_[0].size();i++)
        {
                rhs[i] = 0.0;
                for(unsigned j=0;j<matri_.size();j++)
                        rhs[i] += matri_[j][i]*sourc_[j];
        }

/*      for(unsigned i=0;i<matri_[0].size();i++)
                std::cout<<rhs[i]<<std::endl;*/
        return sqr;
}

//solves the system of equations using conjugate gradient method
//this can be used only when the matrix is symmetric
//since in our case we multiplied the A with its transpose the resulting system is always symmetric
std::vector<double> GEO::CUT::LeastSquares::ConjugateGradient(std::vector<std::vector<double> >coeff, std::vector<double> rhs)
{
        std::vector<double> resi = rhs;
        std::vector<double> resi_new = rhs;
        std::vector<double> p = rhs;
        std::vector<double> soln,tempo;
        soln.resize(rhs.size());
        tempo.resize(rhs.size());

/*      for(unsigned i=0;i<rhs.size();i++)
                std::cout<<soln[i]<<std::endl;*/
	int ittno = 0;
	double epsi = 1e-15;
        while(1)
        {
		ittno++;
                tempo = multiply(coeff,p);
                double alpha = multi_vec(resi,resi)/multi_vec(p,tempo); 
//double alpha = multi_vec(resi,resi);
                
                for(unsigned i=0;i<rhs.size();i++)
                {
                        soln[i] = soln[i]+alpha*p[i];
                        resi_new[i] = resi[i]-alpha*tempo[i];
                }
                double conv_check = multi_vec(resi_new,resi_new);
          //    std::cout<<sqrt(conv_check)<<std::endl;
                if(sqrt(conv_check)<epsi || ittno>4500)
                        break;
		if(ittno>150)
			epsi = 1e-12;
		if(ittno>200)
			epsi = 1e-10;
		if(ittno>350)
			epsi = 1e-08;
                double beta = multi_vec(resi_new,resi_new)/multi_vec(resi,resi);
                for(unsigned i=0;i<rhs.size();i++)
                {
                        p[i] = resi_new[i]+beta*p[i];
                        resi[i] = resi_new[i];
                }
        }
	std::cout<<ittno<<std::endl;//blockkk or remove
        return soln;
}

//multiply a n*n matrix and a n vector
std::vector<double> GEO::CUT::LeastSquares::multiply(std::vector<std::vector<double> > mat, std::vector<double> ve)
{
        std::vector<double> resu(ve.size());
        for(unsigned i=0;i<ve.size();i++)
        {
                resu[i] = 0.0;
                for(unsigned j=0;j<ve.size();j++)
                        resu[i] += mat[i][j]*ve[j];
        }
        return resu;
}

//multiply two vectors
double GEO::CUT::LeastSquares::multi_vec(std::vector<double> mm1, std::vector<double> mm2)
{
        double resu = 0.0;
        for(unsigned i=0;i<mm1.size();i++)
                resu += mm1[i]*mm2[i];
        return resu;
}

//solve the rectangular system with linear least squares
std::vector<double> GEO::CUT::LeastSquares::linear_least_square()
{
        std::vector<std::vector<double> > sqr;
        std::vector<double> rhs;
        sqr = get_square_matrix(rhs);

        unknown_ = ConjugateGradient(sqr, rhs);
        return unknown_;
}
