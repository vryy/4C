/*----------------------------------------------------------------------*/
/*!
 * \file inv_analysis.H
 *
 *\brief inv analysis
\level 2

<pre>
\maintainer Sebastian Kehl
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INV_ANALYSIS_H
#define INV_ANALYSIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

// forward declarations

namespace STR
{
  class TimInt;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class Solver;
}

namespace DRT
{
  class Condition;
  class Discretization;
}

namespace MAT
{
  namespace PAR
  {
    class Material;
  }
}

namespace STR
{

/// Set material parameters
void SetMaterialParameters(int prob,
    Epetra_SerialDenseVector& p_cur,
    std::set<int>& mymatset,
    std::set<int>& myehmatset);


/// Inverse Analysis of Structures
///
/// my long and detailed description
/// including references to papers
///
/// \author rausch
class InvAnalysis {

private:

  /// my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  /// vector layout from the discretization to construct matching vectors and matrices
  //Epetra_Map* dofrowmap_;

  /// linear algebraic solver
  Teuchos::RCP<LINALG::Solver> solver_;

  /// context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  ///  time integrator for structural dynamics
  Teuchos::RCP<STR::TimInt> sti_;

  ///  number of measurment points the
  double tstep_;                          // timestep

  int problem_type_;        // Neumann or Dirichlet Problem
  int np_;                                // number of fitted parameters
  int nmp_;                                // number of timesteps

  Epetra_SerialDenseVector mcurve_;// measured displacment of the experiments (target value)

  double mu_;
  double tol_mu_;

  double error_;
  double error_o_;
  double tol_;

  std::string filename_;

  double kappa_multi_;

  int numb_run_;

  Epetra_SerialDenseVector p_;      // parameter vector

  /// @name Storage
  //@{
  Epetra_SerialDenseMatrix p_s_;
  Epetra_SerialDenseMatrix delta_p_s_;
  Epetra_SerialDenseMatrix ccurve_s_;
  Epetra_SerialDenseMatrix par_s_;
  Epetra_SerialDenseVector mu_s_;
  Epetra_SerialDenseVector error_s_;
  //@}

  std::vector<DRT::Condition*> surfneum_;
  std::vector<DRT::Condition*> surfdir_;

  int reset_out_count_;

  std::vector<std::set<int> > matset_;
  std::vector<std::set<int> > eh_matset_;


public:

  /// standard constructor
  InvAnalysis(
    Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<LINALG::Solver> solver,  ///< algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  /// Read in routine for the material parameters
  void ReadInParameters();

  /// Set material parameters
  void SetParameters(Epetra_SerialDenseVector p_cur);

  /// get current state of calculated curve
  Epetra_SerialDenseVector GetCalculatedCurve();

  /// coarse fitting of parameters
  void PrecondParameters();

  /// fitting the material parameters using the Levenberg-Marquart algorithm
  void CalcNewParameters(Epetra_SerialDenseMatrix cmatrix,  std::vector<double> inc);

  Epetra_SerialDenseVector CalcCvector(bool outputtofile);

  ///
  void Integrate();

  ///
  void PrintStorage(Epetra_SerialDenseMatrix cmatrix,  Epetra_SerialDenseVector delta_p);
  void PrintFile();

  /// Re-initialized micro-scale
  void MultiInvAnaInit();


};  // class InvAnalysis

}  // namespace STR

#endif /*INV_ANALYSIS_H_*/
