/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_analysis.cpp

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "stat_inv_analysis.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
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
STR::INVANA::StatInvAnalysis::StatInvAnalysis(Teuchos::RCP<DRT::Discretization> dis):
discret_(dis),
dofrowmap_(NULL)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();


  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");
  else
    dofrowmap_ = discret_->DofRowMap();

  // this is supposed to be the number of simulation steps in the primal AND the dual problem
  msteps_ = sdyn.get<int>("NUMSTEP");
  double timestep = sdyn.get<double>("TIMESTEP");

  // initialize "state" vectors
  dis_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  disdual_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));

  // initialize the vector of time steps according to the structural dynamic params
  time_ = Teuchos::rcp(new std::vector<double>(msteps_,0.0));
  for (int i=0; i<=msteps_-1; i++)
  {
    (*time_)[i] = (i+1)*timestep;
  }

  // set up an objective function
  objfunct_ = Teuchos::rcp(new STR::INVANA::ObjectiveFunct(discret_, msteps_, time_));


  // set up the material parameter handler:
  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvMatParametrization>(statinvp,"PARAMETRIZATION"))
  {
    case INPAR::STR::stat_inv_mp_smoothkernel:
    {
      dserror("no parametrization based on gaussian kernels yet!");
      //matman_ = Teuchos::rcp(new STR::INVANA::MatParManKernelPara(discret_));
    }
    break;
    case INPAR::STR::stat_inv_mp_elementwise:
    {
      matman_ = Teuchos::rcp(new STR::INVANA::MatParManager(discret_));
    }
    break;
    default:
      dserror("choose a valid method of parametrizing the material parameter field");
    break;
  }

  //set these infeasibly high:
  objval_ = 1.0e17;
  objval_o_ = 1.0e16;

  objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->GetParams())));
  objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*(matman_->GetParams())));

}

/*----------------------------------------------------------------------*/
/* MStep EpetraVector to EpetraMultiVector */
void STR::INVANA::StatInvAnalysis::MStepEpetraToEpetraMulti(Teuchos::RCP<TimIntMStep<Epetra_Vector> > mstepvec,
                                                            Teuchos::RCP<Epetra_MultiVector> multivec)
{
  for (int i=0; i<msteps_; i++)
    (*multivec)(i)->Update(1.0,*(*mstepvec)(-(msteps_-i-1)),0.0);
}

/*----------------------------------------------------------------------*/
/* Mstep double to std::vector<double> */
void STR::INVANA::StatInvAnalysis::MStepDToStdVecD(Teuchos::RCP<TimIntMStep<double> > mstepvec,
                                                   Teuchos::RCP<std::vector<double> > stdvec)
{
  for (int i=0; i<msteps_; i++)
    (*stdvec)[i] = *(*mstepvec)(-(msteps_-i-1));

}

/*----------------------------------------------------------------------*/
/* solve a forward problem */
void STR::INVANA::StatInvAnalysis::SolveForwardProblem()
{

  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
    case INPAR::STR::dyna_statics:
    {
      ADAPTER::StructureBaseAlgorithm adapterbase(sdyn,const_cast<Teuchos::ParameterList&>(sdyn), discret_);
      ADAPTER::Structure& structadaptor = const_cast<ADAPTER::Structure&>(adapterbase.StructureField());

      // do restart
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        structadaptor.ReadRestart(restart);
      }
      structadaptor.Integrate();
      //output_ = structadaptor.DiscWriter();

      // get displacement and time
      MStepEpetraToEpetraMulti(structadaptor.DispMStep(), dis_);
      MStepDToStdVecD(structadaptor.TimeMStep(), time_);

      break;
    }
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_expleuler:
    case INPAR::STR::dyna_centrdiff:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_euma:
    case INPAR::STR::dyna_euimsto:
      dserror("return of multistep-variables only for static analysis (so far)");
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* solve a dual problem */
void STR::INVANA::StatInvAnalysis::SolveAdjointProblem()
{
  //Setup RHS for the adjoints
  Teuchos::RCP<Epetra_MultiVector> objgrad;
  objgrad = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));
  objfunct_->EvaluateGradient(dis_,objgrad);

  //initilize adjoint time integration with RHS as input
  STR::TimIntAdjoint timintadj = TimIntAdjoint(discret_);
  timintadj.SetupAdjoint(objgrad, dis_);

  // adjoint time integration
  timintadj.Integrate();

  // get the solution
  disdual_ = timintadj.ExtractSolution();

//  //append dual solution to the existing forward problem output:
//  for (int i=0; i<msteps_; i++)
//  {
//    output_->NewStep(msteps_+i+1,0.1);
//    output_->WriteVector("displacement", Teuchos::rcp((*disdual_)(i),false));
//  }
//
//  std::cout << "done" << std::endl;

}

/*----------------------------------------------------------------------*/
/* evaluate the objective functions gradient using adjoints */
void STR::INVANA::StatInvAnalysis::EvaluateGradient()
{
  //zero out gradient vector initially
  objgrad_->Scale(0.0);

  //loop the time steps
  for (int j=0; j<msteps_; j++)
  {
    discret_->SetState("displacement", Teuchos::rcp((*dis_)(j),false));
    discret_->SetState("residual displacement",  Teuchos::rcp((*dis_)(j),false));

    matman_->Evaluate(objgrad_, Teuchos::rcp((*disdual_)(j),false));
  }

  //matman_->AddRegularizationGradient(Teuchos::rcp(&objgrad_,false));

}
