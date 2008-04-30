/*!----------------------------------------------------------------------
\file micromaterialgp_static.cpp

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "micromaterialgp_static.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_validparameters.H"

#include "../drt_stru_multi/microstatic.H"

#include "../drt_io/io_control.H"

using namespace std;
using namespace Teuchos;
using namespace IO;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

extern struct _MATERIAL    *mat;


RefCountPtr<MicroStatic> MAT::MicroMaterialGP::microstatic_;
int MAT::MicroMaterialGP::microstaticcounter_;


/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterialGP::MicroMaterialGP(const int gp, const int ele_ID)
  : gp_(gp),
    ele_ID_(ele_ID)
{
  RefCountPtr<DRT::Problem> microproblem = DRT::Problem::Instance(1);
  RefCountPtr<DRT::Discretization> microdis = microproblem->Dis(0, 0);
  dism_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  dis_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  lastalpha_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldalpha_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldfeas_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldKaainv_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldKda_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  microstaticcounter_ += 1;

  // if class handling microscale simulations is not yet initialized
  // -> set up

  if (microstatic_ == null)
  {
    MAT::MicroMaterialGP::SetUpMicroStatic();
  }

  // create and initialize "empty" EAS history map
  EasInit();

  RefCountPtr<DRT::Discretization> actdis = DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  // Check for surface stress conditions due to interfacial phenomena
  vector<DRT::Condition*> surfstresscond(0);
  actdis->GetCondition("SurfaceStress",surfstresscond);
  if (surfstresscond.size())
  {
    surf_stress_man_=rcp(new DRT::SurfStressManager(*actdis));
  }


  // set up micro output
  //
  // Get the macro output prefix and insert element and gauss point
  // identifier. We use the original name here and rely on our (micro)
  // OutputControl object below to act just like the macro (default)
  // OutputControl. In particular we assume that there are always micro and
  // macro control files on restart.
  RCP<IO::OutputControl> macrocontrol = DRT::Problem::Instance(0)->OutputControlFile();
  std::string microprefix = macrocontrol->RestartName();

  unsigned pos = microprefix.rfind('-');
  if (pos!=std::string::npos)
  {
    std::string number = microprefix.substr(pos+1);
    std::string prefix = microprefix.substr(0,pos);

    ostringstream s;
    s << prefix << "_el" << ele_ID_ << "_gp" << gp_ << "-" << number;
    microprefix = s.str();
  }
  else
  {
    ostringstream s;
    s << microprefix << "_el" << ele_ID_ << "_gp" << gp_;
    microprefix = s.str();
  }

  RCP<OutputControl> microcontrol =
    rcp(new OutputControl(actdis->Comm(),
                          DRT::Problem::Instance(1)->ProblemType(),
                          "micro-input-file-not-known",
                          microprefix,
                          genprob.ndim,
                          genprob.restart,
                          macrocontrol->FileSteps()));
  micro_output_ = rcp(new DiscretizationWriter(actdis,microcontrol));

  // do restart if demanded from input file
  istep_ = 0;
  if (genprob.restart)
  {
    istep_ = genprob.restart;
    microstatic_->ReadRestart(istep_, dis_, lastalpha_, microprefix);
    // both dis_ and dism_ are the same
    dism_->Update(1.0, *dis_, 0.0);
    // both lastalpha and oldalpha are the same
    *oldalpha_ = *lastalpha_;
    microstatic_->SetOldState(dis_, dism_, surf_stress_man_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);
  }
  else
  {
    micro_output_->WriteMesh(istep_, microstatic_->GetTime());
  }

  // we are using the same structural dynamic parameters as on the
  // macroscale, so to avoid checking the equivalence of the reader
  // GiD sections we simply ask the macroproblem for its parameters.
  // note that the microscale solver is currently always UMFPACK
  // (hard coded)!
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

//   timen_ = microstatic_->GetTime();
  // we set the microscale time _ALWAYS_ to 0 (also in restart case)
  // in order to handle the implicit update query!!!
  timen_ = 0.;

  dt_    = sdyn.get<double>("TIMESTEP");
}

/// destructor

MAT::MicroMaterialGP::~MicroMaterialGP()
{
  microstaticcounter_ -= 1;
  if (microstaticcounter_==0)
    microstatic_ = Teuchos::null;
}


/// Set up microscale generalized alpha

void MAT::MicroMaterialGP::SetUpMicroStatic()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // currently taking the parameters of the macroscale problem here!!!
  // -> it is generally no problem to take the ones of the microscale,
  // but then one should check that the two inputfiles are in sync at
  // least for the dynamic parameters so that e.g. dt is the same for
  // both problems. output options could/should be different, but the
  // output interval (= output every nstep) is part of StructuralDynamicsParams

  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // always choose UMFPACK as microstructural solver
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  solveparams->set("solver","umfpack");
  solveparams->set("symmetric",false);
  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a static "time integrator"
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> params = rcp(new ParameterList());
  MicroStatic::SetDefaults(*params);

  params->set<double>("beta",sdyn.get<double>("BETA"));
  params->set<double>("gamma",sdyn.get<double>("GAMMA"));
  params->set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
  params->set<double>("alpha f",sdyn.get<double>("ALPHA_F"));
  params->set<string>("convcheck", sdyn.get<string>("CONV_CHECK"));
  params->set<double>("total time",0.0);
  params->set<double>("delta time",sdyn.get<double>("TIMESTEP"));
  params->set<int>   ("step",0);
  params->set<int>   ("nstep",sdyn.get<int>("NUMSTEP"));
  params->set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
  params->set<int>   ("num iterations",-1);

  params->set<double>("tolerance residual",sdyn.get<double>("TOLRES"));
  params->set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
  params->set<bool>  ("print to screen",true);

  params->set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
  params->set<int>   ("io disp every nstep",sdyn.get<int>("RESEVRYDISP"));

  switch (Teuchos::getIntegralValue<STRUCT_STRESS_TYP>(ioflags,"STRUCT_STRESS"))
  {
  case struct_stress_none:
    params->set<string>("io structural stress", "none");
    break;
  case struct_stress_cauchy:
    params->set<string>("io structural stress", "cauchy");
    break;
  case struct_stress_pk:
    params->set<string>("io structural stress", "2PK");
    break;
  default:
    params->set<string>("io structural stress", "none");
    break;
  }

  params->set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

  switch (Teuchos::getIntegralValue<STRUCT_STRAIN_TYP>(ioflags,"STRUCT_STRAIN"))
  {
  case struct_strain_none:
    params->set<string>("io structural strain", "none");
    break;
  case struct_strain_ea:
    params->set<string>("io structural strain", "euler_almansi");
    break;
  case struct_strain_gl:
    params->set<string>("io structural strain", "green_lagrange");
    break;
  default:
    params->set<string>("io structural strain", "none");
    break;
  }

  params->set<int>   ("restart",probtype.get<int>("RESTART"));
  params->set<int>   ("write restart every",sdyn.get<int>("RESTARTEVRY"));

  microstatic_ = rcp(new MicroStatic(params,actdis,solver));
}

void MAT::MicroMaterialGP::EasInit()
{
  // create the parameters for the discretization
  ParameterList p;
  // action for elements
  p.set("action","eas_init_multi");
  p.set("lastalpha", lastalpha_);
  p.set("oldalpha", oldalpha_);
  p.set("oldfeas", oldfeas_);
  p.set("oldKaainv", oldKaainv_);
  p.set("oldKda", oldKda_);

  RefCountPtr<DRT::Discretization> actdis = DRT::Problem::Instance(1)->Dis(genprob.numsf,0);
  actdis->Evaluate(p,null,null,null,null,null);

  return;
}


/// perform microscale simulation

void MAT::MicroMaterialGP::PerformMicroSimulation(const Epetra_SerialDenseMatrix* defgrd,
                                                  Epetra_SerialDenseVector* stress,
                                                  Epetra_SerialDenseMatrix* cmat,
                                                  double* density,
                                                  const double time,
                                                  string action)
{
  // this is a comparison of two doubles, but since timen_ is always a
  // copy of time as long as we are within a time step, any variation
  // must be due to an update of macroscale time!
  // -> no update for a new time step if timen_=0
  // note that timen_ is set in the following if statement!
  if (time != timen_ && timen_ != 0.)
  {
    microstatic_->UpdateNewTimeStep(dis_, dism_, oldalpha_, lastalpha_);
  }

  // set displacements and EAS data of last step
  microstatic_->SetOldState(dis_, dism_, surf_stress_man_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);

  // check if we have to update absolute time and step number
  // in case of restart, timen_ is set to the current total time in
  // the constructor, whereas "time" handed on to this function is
  // still 0. since the macroscopic update of total time is done after
  // the initializing phase in which this function is called the first
  // time
  if (time != timen_ && time != 0.)
  {
    // Microscale data should be output when macroscale is entering a
    // new timestep, not in every macroscopic iteration! Therefore
    // output is written in the beginning of a microscopic step if
    // necessary at all. Problem: we don't get any output for the very
    // last time step since the macro-program finishes and the
    // micro-program is not called again to write output.

    // We don't want to write results after just having constructed
    // the StruGenAlpha class which corresponds to a total time of 0.
    // (in the calculation phase, total time is instantly set to the first
    // time step)
    if (timen_ != 0.)
      microstatic_->Output(micro_output_, timen_, istep_, dt_);
    timen_ = time;
    istep_++;
  }

  // set current absolute time, step number
  microstatic_->SetTime(timen_, istep_);

  microstatic_->Predictor(defgrd);
  microstatic_->FullNewton();
  microstatic_->StaticHomogenization(stress, cmat, density, defgrd);

  // note that it is not necessary to save displacements and EAS data
  // explicitly since we dealt with RCP's -> any update in class
  // microstatic and the elements, respectively, inherently updates the
  // micromaterialgp_static data!

  // save calculated displacements
  //dism_ = microstatic_->ReturnNewDism();

  // clear displacements in MicroStruGenAlpha for next usage
  microstatic_->ClearState();
}

#endif
