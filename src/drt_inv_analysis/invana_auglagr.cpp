/*----------------------------------------------------------------------*/
/*!
 * \file invana_auglagr.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "invana_auglagr.H"
#include "invana_utils.H"

#include "matpar_manager.H"
#include "objective_funct.H"
#include "regularization_base.H"

#include "../drt_adapter/ad_str_structure.H"
#include "timint_adjoint.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_timintmstep.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_inpar/drt_validparameters.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*/
/* standard constructor                                      keh 10/13  */
/*----------------------------------------------------------------------*/
STR::INVANA::InvanaAugLagr::InvanaAugLagr():
  InvanaBase(),
dis_(Teuchos::null),
disdual_(Teuchos::null),
time_(Teuchos::null),
msteps_(0),
fprestart_(0)
{
  const Teuchos::ParameterList& sdyn  = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // this is supposed to be the number of simulation steps in the primal AND the dual problem
  msteps_ = sdyn.get<int>("NUMSTEP");
  double timestep = sdyn.get<double>("TIMESTEP");

  // initialize the vector of time steps according to the structural dynamic params
  time_ = Teuchos::rcp(new std::vector<double>(msteps_,0.0));
  for (int i=0; i<=msteps_-1; i++)
  {
    (*time_)[i] = (i+1)*timestep;
  }

  fprestart_ = invp.get<int>("FPRESTART"); //this is supposed to be the forward problem restart

}

void STR::INVANA::InvanaAugLagr::Setup()
{
  if (not Discret()->Filled() || not Discret()->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");

  // initialize "state" vectors
  dis_ = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()),msteps_,true));
  disdual_ = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()),msteps_,true));

  //output for the forward problem
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::string prefix = DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
  size_t pos = filename.rfind('/');
  size_t pos2 = prefix.rfind('-');
  std::string filenameout = filename.substr(0,pos+1) + prefix.substr(0,pos2) + "_forward" + filename.substr(pos+1+prefix.length());

  Teuchos::RCP<IO::OutputControl> controlfile =
    Teuchos::rcp(new IO::OutputControl(
      Discret()->Comm(),
      DRT::Problem::Instance()->ProblemName(),
      DRT::Problem::Instance()->SpatialApproximation(),
      DRT::Problem::Instance()->OutputControlFile()->InputFileName(),
      filenameout,
      DRT::Problem::Instance()->NDim(),
      fprestart_,
      DRT::Problem::Instance()->OutputControlFile()->FileSteps(),
      DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_BIN")
    )
  );
  // give the discretization another controlfile for output
  Discret()->Writer()->SetOutput(controlfile);

  return;
}

/*----------------------------------------------------------------------*/
/* MStep EpetraVector to EpetraMultiVector                   keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::InvanaAugLagr::MStepEpetraToEpetraMulti(Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > mstepvec,
                                                            Teuchos::RCP<Epetra_MultiVector> multivec)
{
  for (int i=0; i<msteps_; i++)
    (*multivec)(i)->Update(1.0,*(*mstepvec)(-(msteps_-i-1)),0.0);
}

/*----------------------------------------------------------------------*/
/* Mstep double to std::vector<double>                       keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::InvanaAugLagr::MStepDToStdVecD(Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > mstepvec,
                                                   Teuchos::RCP<std::vector<double> > stdvec)
{
  for (int i=0; i<msteps_; i++)
    (*stdvec)[i] = *(*mstepvec)(-(msteps_-i-1));

}

/*----------------------------------------------------------------------*/
/* solve primal problem                                      keh 10/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::InvanaAugLagr::SolveForwardProblem()
{
  // use the same control file for every run since usually the last one is of interest
  Discret()->Writer()->OverwriteResultFile();

  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
  {
    case INPAR::STR::dyna_statics:
    {
      ADAPTER::StructureBaseAlgorithm adapterbase(sdyn,const_cast<Teuchos::ParameterList&>(sdyn), Discret());
      ADAPTER::Structure& structadaptor = const_cast<ADAPTER::Structure&>(adapterbase.StructureField());

      // do restart but the one which is explicitly given in the INVERSE ANALYSIS section
      if (fprestart_)
      {
        dserror("Restarting from within a timestep of the forward problem needs some tweaking first!");
        structadaptor.ReadRestart(fprestart_);
      }
      structadaptor.Integrate();

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
/* solve dual problem                                       keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::InvanaAugLagr::SolveAdjointProblem()
{
  //Setup RHS for the adjoints
  Teuchos::RCP<Epetra_MultiVector> objgrad;
  objgrad = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()),msteps_,true));
  ObjectiveFunct()->EvaluateGradient(dis_,objgrad);

  //initialize adjoint time integration with RHS as input
  STR::TimIntAdjoint timintadj = TimIntAdjoint(Discret(), *time_);
  timintadj.SetupAdjoint(objgrad, dis_);

  // adjoint time integration
  timintadj.Integrate();

  // get the solution
  disdual_ = timintadj.ExtractSolution();

}

void STR::INVANA::InvanaAugLagr::Evaluate(const Epetra_MultiVector& sol, double* val, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  Matman()->ReplaceParams(sol);
  ResetDiscretization();

  if ( gradient != Teuchos::null or val!=NULL )
  {
    SolveForwardProblem();
    EvaluateError(sol,val);

    if (gradient != Teuchos::null)
    {
      SolveAdjointProblem();
      EvaluateGradient(sol,gradient);
    }
  }
}

/*----------------------------------------------------------------------*/
/* evaluate gradient of the objective function                          */
/* using the adjoint equations                              keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::InvanaAugLagr::EvaluateGradient(const Epetra_MultiVector& sol, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  //zero out gradient vector initially
  gradient->Scale(0.0);

  Teuchos::RCP<Epetra_Vector> zeros;
  zeros = LINALG::CreateVector(*(Discret()->DofRowMap()), true);
  zeros->Scale(0.0);

  //loop the time steps
  bool sumaccrosprocs = false;
  for (int j=0; j<msteps_; j++)
  {
    Discret()->SetState(0, "displacement", Teuchos::rcp((*dis_)(j),false));
    Discret()->SetState(0, "residual displacement", zeros);
    Discret()->SetState(0, "dual displacement", Teuchos::rcp((*disdual_)(j),false));

    if (j==msteps_-1) sumaccrosprocs = true;

    Matman()->Evaluate((*time_)[j], gradient, sumaccrosprocs);
  }

  if (Regman() != Teuchos::null)
    Regman()->EvaluateGradient(sol,gradient);

}

void STR::INVANA::InvanaAugLagr::ResetDiscretization()
{
  Teuchos::ParameterList p;
  p.set("action","calc_struct_reset_all");
  Discret()->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
}

/*----------------------------------------------------------------------*/
/* Evaluate the objective function                          keh 10/13   */
/*----------------------------------------------------------------------*/
void STR::INVANA::InvanaAugLagr::EvaluateError(const Epetra_MultiVector& sol, double* val)
{
  ObjectiveFunct()->Evaluate(dis_,*val);

  Regman()->Evaluate(sol,val);
}
