/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/

#include "least_squares.H"
#include <iostream>
#include <cmath>

#if 0
//solves the system of equations using conjugate gradient method
//this can be used only when the matrix is symmetric
//since in our case we multiplied the A with its transpose the resulting system is always symmetric
/*std::vector<double> GEO::CUT::LeastSquares::ConjugateGradient(std::vector<std::vector<double> >coeff, std::vector<double> rhs)
{
        std::vector<double> resi = rhs;
        std::vector<double> resi_new = rhs;
        std::vector<double> p = rhs;
        std::vector<double> soln,tempo;
        soln.resize(rhs.size());
        tempo.resize(rhs.size());

//      for(unsigned i=0;i<rhs.size();i++)
//                std::cout<<soln[i]<<std::endl;
  int ittno = 0;
  double epsi = 1e-12;
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
                double conv_check =  maxAbsolute(resi_new);
          //    std::cout<<sqrt(conv_check)<<std::endl;
                if(conv_check<epsi || ittno>900)
                        break;
    if(ittno>250)
      epsi = 1e-08;
    if(ittno>300)
      epsi = 1e-06;
    if(ittno>450)
      epsi = 1e-05;
                double beta = multi_vec(resi_new,resi_new)/multi_vec(resi,resi);
                for(unsigned i=0;i<rhs.size();i++)
                {
                        p[i] = resi_new[i]+beta*p[i];
                        resi[i] = resi_new[i];
                }
    std::cout<<"conv_check"<<conv_check<<"\n";
        }
  std::cout<<"least_iter = "<<ittno<<"\t"<<"epsi = "<<epsi<<"\n";
        return soln;
}*/

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



double GEO::CUT::LeastSquares::maxAbsolute(std::vector<double>a)
{
  double maxx=0.0;
  for(unsigned i=0;i<a.size();i++)
  {
    if(fabs(a[i])>maxx)
      maxx = fabs(a[i]);
  }
  return maxx;
}
#endif

// solve the rectangular system with linear least squares
Epetra_SerialDenseVector GEO::CUT::LeastSquares::linear_least_square()
{
  Epetra_SerialDenseMatrix sqr(matri_[0].size(), matri_[0].size());
  Epetra_SerialDenseVector rhs(matri_[0].size());
  sqr = get_square_matrix(rhs);
  unknown_.Size(matri_[0].size());

  Epetra_SerialDenseSolver solve_for_GPweights;
  solve_for_GPweights.SetMatrix(sqr);
  solve_for_GPweights.SetVectors(unknown_, rhs);
  solve_for_GPweights.FactorWithEquilibration(true);
  int err2 = solve_for_GPweights.Factor();
  int err = solve_for_GPweights.Solve();
  if ((err != 0) && (err2 != 0))
    dserror(
        "Computation of Gauss weights failed, Ill"
        "conditioned matrix in least square");


  /*  Epetra_SerialDenseMatrix matt(sqr.size(),sqr.size());
    Epetra_SerialDenseVector unn(sqr.size());
    Epetra_SerialDenseVector rrr(sqr.size());
    for(unsigned i=0;i<sqr.size();i++)
    {
      for(unsigned j=0;j<sqr.size();j++)
        matt(j,i) = sqr[i][j];
      unn(i) = 0.0;
      rrr(i) = rhs[i];
    }
    unknown_ = ConjugateGradient(sqr, rhs);*/

  return unknown_;
}

// premultiplying the matrix with its transpose to get the square matrix
// the source terms also get multiplied
Epetra_SerialDenseMatrix GEO::CUT::LeastSquares::get_square_matrix(Epetra_SerialDenseVector &rhs)
{
  Epetra_SerialDenseMatrix sqr(matri_[0].size(), matri_[0].size());

  for (unsigned i = 0; i < matri_[0].size(); i++)
  {
    for (unsigned j = 0; j < matri_[0].size(); j++)
    {
      sqr(j, i) = 0.0;

      // it is sqr(j,i) because the Epetra elements are column ordered first
      for (unsigned k = 0; k < matri_.size(); k++) sqr(j, i) += matri_[k][i] * matri_[k][j];
    }
  }

  for (unsigned i = 0; i < matri_[0].size(); i++)
  {
    rhs(i) = 0.0;
    for (unsigned j = 0; j < matri_.size(); j++) rhs(i) += matri_[j][i] * sourc_(j);
  }

  return sqr;
}
