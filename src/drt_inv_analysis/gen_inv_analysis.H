/*----------------------------------------------------------------------*/
/*!
 * \file gen_inv_analysis.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gee
            089 - 289-15239
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef GEN_INV_ANALYSIS_H
#define GEN_INV_ANALYSIS_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_SerialDenseVector.h>
#include <math.h>
#include "../drt_adapter/adapter_structure.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_structure/strugenalpha.H"
#include "../drt_structure/strtimint.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_matelast/elast_summand.H"


namespace STR
{

/// Inverse Analysis of Structures
///
/// my long and detailed description
/// including references to papers
///
/// \author rausch
class GenInvAnalysis {

private:

  /// my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  /// linear algebraic solver
  Teuchos::RCP<LINALG::Solver> solver_;

  /// context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  ///  time integrator for structural dynamics
  Teuchos::RCP<STR::TimInt> sti_;

  ///  number of measurment points the
  double tstep_;                    // timestep

  int problem_type_;		    // Neumann or Dirichlet Problem
  int np_;                          // number of fitted parameters (default 3, 4 with kappa)
  int nmp_;                         // number of fitted parameters (default 3, 4 with kappa)

  Epetra_SerialDenseVector mcurve_; // measured displacment of the experiments (target value)

  double mu_;

  double error_;
  double error_o_;
  double tol_;

  string filename_;

  double kappa_multi_;

  int numb_run_;

  Epetra_SerialDenseVector p_;	   // parameter vector

  /// @name Storage
  //@{
  Epetra_SerialDenseMatrix p_s_;
  Epetra_SerialDenseMatrix delta_p_s_;
  Epetra_SerialDenseMatrix par_s_;
  Epetra_SerialDenseVector mu_s_;
  Epetra_SerialDenseVector error_s_;
  //@}

  int reset_out_count_;

  // Michael's inv_analysis stuff
  int                  nnodes_;        // number of nodes observed
  //int                  ndirs_;         // number of directions per node observed
  int                  ndofs_;         // total number of dofs monitored
  int                  nsteps_;        // number of time steps observed
  vector<int>          nodes_;         // list of node gids observed
  vector<vector<int> > dofs_;          // list of dofs on each node that are monitored

public:

  /// standard constructor
  GenInvAnalysis(
    Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<LINALG::Solver> solver,  ///< algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  /// Read in routine for the material parameters
  void ReadInParameters();

  /// Set material parameters
  void SetParameters(Epetra_SerialDenseVector p_cur);

  /// get current state of calculated curve
  Epetra_SerialDenseVector GetCalculatedCurve(Epetra_Vector& disp);

  /// coarse fitting of parameters
  void PrecondParameters();

  /// fitting the material parameters using the Levenberg-Marquart algorithm
  void CalcNewParameters(Epetra_SerialDenseMatrix& cmatrix, vector<double>& perturb);


  ///
  Epetra_SerialDenseVector CalcCvector(bool outputtofile);

  ///
  void Integrate();

  ///
  void PrintStorage(Epetra_SerialDenseVector delta_p);

};  // class GenInvAnalysis

}  // namespace STR

#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
