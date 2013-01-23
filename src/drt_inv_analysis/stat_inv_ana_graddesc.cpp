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
  //const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& invap = DRT::Problem::Instance()->StatInverseAnalysisParams();
  maxiter_ = invap.get<int>("MAXITER");

  //set stepsize for gradient scheme
  stepsize_ = invap.get<double>("STEPSIZE");

}

/*----------------------------------------------------------------------*/
/* decide for an optimization algorithm*/
void STR::INVANA::StatInvAnaGradDesc::Optimize()
{
  for (int i=0; i<maxiter_; i++)
  {
    Epetra_MultiVector objgrad(*params_);

    matman_->UpdateParams(params_);

    SolveForwardProblem();

    SolveAdjointProblem();

    EvaluateGradient(&objgrad);

    UpdateOptStep(&objgrad,i);
  }

  return;
}

/*----------------------------------------------------------------------*/
/* do the update of the parameter vector */
void STR::INVANA::StatInvAnaGradDesc::UpdateOptStep(Epetra_MultiVector* objgrad, int nstep)
{

  objval_o_ = objval_;
  //get actual value of the objective
  objfunct_->Evaluate(dis_,&objval_);

  if (nstep == 0)
    return;

  //this is an extremly simple timeadaptive gradient descent scheme
  if (objval_ <= objval_o_)
  {
    //scale the gradient
    Epetra_SerialDenseVector normvec(objgrad->NumVectors());
    objgrad->NormInf(normvec.Values());
    double normgrad = normvec.NormInf();
    cout << "gradnorm: " << normgrad << endl;
    objgrad->Scale(-1.0/normgrad);

    params_->Update(stepsize_,*objgrad,1.0);
    stepsize_ = stepsize_*1.2;
  }
  else
  {
    stepsize_ = stepsize_*0.5;
  }

  cout << "this steps optimality evaluation: " << objval_ << " step: " << stepsize_ << endl;
  cout << *params_ << endl;

}
