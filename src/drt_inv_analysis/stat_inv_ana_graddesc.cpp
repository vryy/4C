/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_graddesc.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "stat_inv_ana_graddesc.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_inpar/drt_validparameters.H"

// needed to deal with materials
#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_material.H"

#include "objective_funct.H"
#include "timint_adjoint.H"
#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::INVANA::StatInvAnaGradDesc::StatInvAnaGradDesc(Teuchos::RCP<DRT::Discretization> dis):
  StatInvAnalysis(dis),
stepsize_(0.0),
maxiter_(0)
{
  const Teuchos::ParameterList& invap = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // max number of iterations
  maxiter_ = invap.get<int>("MAXITER");

  //set stepsize for gradient scheme
  stepsize_ = invap.get<double>("STEPSIZE");

}

/*----------------------------------------------------------------------*/
/* decide for an optimization algorithm*/
void STR::INVANA::StatInvAnaGradDesc::Optimize()
{

  objgrad_o_->Scale(0.0);

  SolveForwardProblem();

  SolveAdjointProblem();

  EvaluateGradient();



  Epetra_SerialDenseVector vnorm(matman_->GetParams()->NumVectors());
  for (int i=0; i<matman_->GetParams()->NumVectors(); i++)
    vnorm(i) = 1.0e16;


  int i=0;

  while (vnorm.Norm1() > 1.0e-10 && i<maxiter_)
  {
    UpdateOptStep();

    // summary
    objgrad_->NormInf(vnorm.Values());

    cout << "this steps optimality evaluation: ERROR: " << objval_ << "   GRADNORM: " << vnorm.Norm1() << " -> next stepsize: " << stepsize_ << endl;

    cout << "************************** END OF STEP " << i << " *******************************" << endl;

    i=i+1;
  }

  //cout << "************************** FINAL SET OF MATERIAL PARAMETERS **************************" << endl;
  //cout << *(matman_->GetParams()) << endl;
  return;

}

/*----------------------------------------------------------------------*/
/* do the update of the parameter vector */
void STR::INVANA::StatInvAnaGradDesc::UpdateOptStep()
{
  objgrad_o_->Scale(1.0, *objgrad_);

  //scale the gradient
  double normgrad;
  objgrad_->Norm1(&normgrad);

  objgrad_->Scale(-1.0/normgrad*stepsize_);

  matman_->UpdateParams(objgrad_);

  SolveForwardProblem();

  //get actual value of the objective
  objval_ = objfunct_->Evaluate(dis_,matman_);

  //this is an extremly simple timeadaptive gradient descent scheme
  if (objval_ < objval_o_)
  {
    if (stepsize_ < 1.0e4) stepsize_ = stepsize_*1.2;

    objval_o_ = objval_;

    SolveAdjointProblem();
    EvaluateGradient();
  }
  else
  {
    matman_->ResetParams();
    if (stepsize_ > 1.0e-4) stepsize_ = stepsize_*0.5;
    objgrad_->Scale(1.0,*objgrad_o_);

  }

  return;

}
