/*----------------------------------------------------------------------*/
/*!
 * \file inv_analysis.H

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/rausch
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INV_ANALYSIS_H
#define INV_ANALYSIS_H

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
class InvAnalysis {

private:

  /// my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  /// vector layout from the discretization to construct matching vectors and matrices
  Epetra_Map* dofrowmap_;

  /// linear algebraic solver
  Teuchos::RCP<LINALG::Solver> solver_;

  /// context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  ///  time integrator for structural dynamics
  Teuchos::RCP<STR::TimInt> sti_;

  ///  number of measurment points the
  int mp_;
  double tstep_;                          // timestep

  vector<double>   disp_;		        // displacment vector in all three direction
  vector<int>      bc_nodes_;             // nodes with bc
  vector<int>      bc_dof_;               // nodes with bc
  vector<double>   load_;		        // load vector in all three direction
  int problem_type_;		        // Neumann or Dirichlet Problem
  int pp_;                                // number of fitted parameters (default 3, 4 with kappa)

  Epetra_SerialDenseVector ccurve_;       // calculated displacments
  Epetra_SerialDenseVector ccurve_o_;	// calculated displacment of the previous loop
  Epetra_SerialDenseVector mcurve_;	// measured displacment of the experiments (target value)
  Epetra_SerialDenseVector rcurve_;	// final_disp_ - measured_disp

  double mu_;
  double mup_;
  double mum_;

  int precond_;

  double error_;
  double error_o_;
  double tol_;

  int numb_run_;
  int direction_;                       // gives the direction of the stretch

  Epetra_SerialDenseVector p_;	      // parameter vector
  Epetra_SerialDenseVector p_neg_;	      // parameter vector
  Epetra_SerialDenseVector p_o_;        // old parameter vector
  Epetra_SerialDenseVector delta_p_;    // change of parameter

  /// @name Storage
  //@{
  Epetra_SerialDenseMatrix p_s_;
  Epetra_SerialDenseMatrix delta_p_s_;
  Epetra_SerialDenseMatrix ccurve_s_;
  Epetra_SerialDenseMatrix par_s_;
  Epetra_SerialDenseVector mu_s_;
  Epetra_SerialDenseVector error_s_;
  //@}

  vector<DRT::Condition*> surfneum_;
  vector<DRT::Condition*> surfdir_;

public:

  /// standard constructor
  InvAnalysis(
    Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<LINALG::Solver> solver,  ///< algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  /// reset parameters for new fitting loop
  void ResetParameters();

  /// get current state of calculated curve
  void GetCalculatedCurve(
    const Teuchos::RCP<Epetra_Vector> value,  ///< either
                                              ///< displacements on
                                              ///< NBC or
                                              ///< reaction forces on DBC

    int numb  ///< time step
    );

  /// coarse fitting of parameters
  void PrecondParameters();

  /// fitting the material parameters using the Levenberg-Marquart algorithm
  void CalcNewParameters();

  ///
  void Evaluate();

  ///
  void Integrate();

  ///
  void PrintStorage();

};  // class InvAnalysis

}  // namespace STR

#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
