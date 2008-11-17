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

extern struct _MATERIAL    *mat;


std::map<int, RefCountPtr<STRUMULTI::MicroStatic> > MAT::MicroMaterialGP::microstaticmap_;
std::map<int, int> MAT::MicroMaterialGP::microstaticcounter_;


/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterialGP::MicroMaterialGP(const int gp, const int ele_ID, const bool eleowner,
                                      const double time, const int microdisnum)
  : gp_(gp),
    ele_ID_(ele_ID),
    microdisnum_(microdisnum)
{
  RefCountPtr<DRT::Problem> microproblem = DRT::Problem::Instance(microdisnum_);
  RefCountPtr<DRT::Discretization> microdis = microproblem->Dis(0, 0);
  dism_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  disn_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  dis_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  lastalpha_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldalpha_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldfeas_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldKaainv_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  oldKda_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);

  // if class handling microscale simulations is not yet initialized
  // -> set up

  if (microstaticmap_.find(microdisnum_) == microstaticmap_.end() or microstaticmap_[microdisnum_] == Teuchos::null)
  {
    // create "time integration" class for this microstructure
    MAT::MicroMaterialGP::SetUpMicroStatic();
    // create a counter of macroscale GP associated with this "time integration" class
    // note that the counter is immediately updated afterwards!
    microstaticcounter_[microdisnum_] = 0;
  }

  microstaticcounter_[microdisnum] += 1;

  // create and initialize "empty" EAS history map
  EasInit();

  // Check for surface stress conditions due to interfacial phenomena
  vector<DRT::Condition*> surfstresscond(0);
  microdis->GetCondition("SurfaceStress",surfstresscond);
  if (surfstresscond.size())
  {
    surf_stress_man_=rcp(new UTILS::SurfStressManager(*microdis));
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
  std::string micronewprefix = macrocontrol->NewOutputFileName();

  // figure out how the file we restart from is called on the microscale
  unsigned pos = microprefix.rfind('-');
  if (pos!=std::string::npos)
  {
    std::string number = microprefix.substr(pos+1);
    std::string prefix = microprefix.substr(0,pos);
    ostringstream s;
    s << prefix << "_el" << ele_ID_ << "_gp" << gp_;
    microprefix = s.str();
    s << "-" << number;
    restartname_ = s.str();
  }
  else
  {
    ostringstream s;
    s << microprefix << "_el" << ele_ID_ << "_gp" << gp_;
    restartname_ = s.str();
  }

  if (eleowner)
  {
    // figure out how the new output file is called on the microscale
    // note: the trailing number must be the same as on the macroscale
    std::string newfilename;
    unsigned pos = micronewprefix.rfind('-');
    if (pos!=std::string::npos)
    {
      std::string number = micronewprefix.substr(pos+1);
      std::string prefix = micronewprefix.substr(0,pos);
      ostringstream s;
      s << prefix << "_el" << ele_ID_ << "_gp" << gp_ << "-" << number;
      newfilename = s.str();
    }
    else
    {
      ostringstream s;
      s << micronewprefix << "_el" << ele_ID_ << "_gp" << gp_;
      newfilename = s.str();
    }

    RCP<OutputControl> microcontrol =
      rcp(new OutputControl(microdis->Comm(),
                            DRT::Problem::Instance(microdisnum_)->ProblemType(),
                            microproblem->SpatialApproximation(),
                            "micro-input-file-not-known",
                            restartname_,
                            newfilename,
                            genprob.ndim,
                            genprob.restart,
                            macrocontrol->FileSteps()));

    micro_output_ = rcp(new DiscretizationWriter(microdis,microcontrol));
  }

  istep_ = 0;

  if (eleowner) micro_output_->WriteMesh(istep_, microstaticmap_[microdisnum_]->GetTime());


  // we are using the same structural dynamic parameters as on the
  // macroscale, so to avoid checking the equivalence of the reader
  // GiD sections we simply ask the macroproblem for its parameters.
  // note that the microscale solver is currently always UMFPACK
  // (hard coded)!
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // we set the microscale time _ALWAYS_ to 0 (also in restart case)
  // in order to handle the implicit update query!!!
  timen_ = 0.;

  dt_    = sdyn.get<double>("TIMESTEP");

  // check whether we are using modified Newton as a nonlinear solver
  // on the macroscale or not
  if (Teuchos::getIntegralValue<int>(sdyn,"NLNSOL")==STRUCT_DYNAMIC::modnewton)
    mod_newton_ = true;
  else
    mod_newton_ = false;

  build_stiff_ = true;
}

/// destructor

MAT::MicroMaterialGP::~MicroMaterialGP()
{
  microstaticcounter_[microdisnum_] -= 1;
  if (microstaticcounter_[microdisnum_]==0)
    microstaticmap_[microdisnum_] = Teuchos::null;
}


/// Read restart

void MAT::MicroMaterialGP::ReadRestart()
{
  istep_ = genprob.restart;
  microstaticmap_[microdisnum_]->ReadRestart(istep_, dis_, lastalpha_, surf_stress_man_, restartname_);
  // both dism_ and disn_ equal dis_
  dism_->Update(1.0, *dis_, 0.0);
  disn_->Update(1.0, *dis_, 0.0);
  // both lastalpha and oldalpha are the same
  *oldalpha_ = *lastalpha_;
}

/// Set up microscale generalized alpha

void MAT::MicroMaterialGP::SetUpMicroStatic()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance(microdisnum_)->Dis(genprob.numsf,0);

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
//   RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
//   solveparams->set("solver","umfpack");
//   solveparams->set("symmetric",false);
//   RefCountPtr<LINALG::Solver> solver =
//     rcp(new LINALG::Solver(solveparams,actdis->Comm(),
//                            DRT::Problem::Instance()->ErrorFile()->Handle()));
//   actdis->ComputeNullSpaceIfNecessary(*solveparams);
  RefCountPtr<LINALG::Solver> solver = rcp (new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                                                               actdis->Comm(),
                                                               DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // create a static "time integrator"
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> params = rcp(new ParameterList());
  STRUMULTI::MicroStatic::SetDefaults(*params);

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

  params->set<bool>  ("ADAPTCONV",getIntegralValue<int>(sdyn,"ADAPTCONV")==1);
  params->set<double>("ADAPTCONV_BETTER",sdyn.get<double>("ADAPTCONV_BETTER"));

  microstaticmap_[microdisnum_] = rcp(new STRUMULTI::MicroStatic(params,actdis,solver));
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

  RefCountPtr<DRT::Discretization> actdis = DRT::Problem::Instance(microdisnum_)->Dis(genprob.numsf,0);
  actdis->Evaluate(p,null,null,null,null,null);

  return;
}


/// perform microscale simulation

void MAT::MicroMaterialGP::PerformMicroSimulation(LINALG::Matrix<3,3>* defgrd,
                                                  LINALG::Matrix<6,1>* stress,
                                                  LINALG::Matrix<6,6>* cmat,
                                                  double* density,
                                                  const double time,
                                                  const bool eleowner)
{
  // select corresponding "time integration class" for this microstructure
  RefCountPtr<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  // this is a comparison of two doubles, but since timen_ is always a
  // copy of time as long as we are within a time step, any variation
  // must be due to an update of macroscale time!
  // -> no update for a new time step if timen_=0
  // note that timen_ is set in the following if statement!
  if (time != timen_ && timen_ != 0.)
  {
    microstatic->UpdateNewTimeStep(dis_, dism_, disn_, oldalpha_, lastalpha_, surf_stress_man_);
    microstatic->SetNewStep(true);           // this is needed for possible surface stresses
                                             // (comparison of A_old and A_new, which are
                                             // nearly the same in the beginning of a new
                                             // time step due to our choice of predictors)
  }

  if (time != timen_)
    build_stiff_ = true;

  // set displacements and EAS data of last step
  microstatic->SetOldState(dis_, dism_, disn_, surf_stress_man_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);

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
    // Note: only the process owning the corresponding macro-element
    // is writing output here!

    // We don't want to write results after just having constructed
    // the StruGenAlpha class which corresponds to a total time of 0.
    // (in the calculation phase, total time is instantly set to the first
    // time step)
    if (timen_ != 0. && eleowner)
      microstatic->Output(micro_output_, timen_, istep_, dt_);
    timen_ = time;
    istep_++;
  }



  // set current absolute time, step number
  microstatic->SetTime(timen_, istep_);

  microstatic->Predictor(defgrd);
  microstatic->SetNewStep(false);
  microstatic->FullNewton();
  microstatic->StaticHomogenization(stress, cmat, density, defgrd, mod_newton_, build_stiff_);

  // note that it is not necessary to save displacements and EAS data
  // explicitly since we dealt with RCP's -> any update in class
  // microstatic and the elements, respectively, inherently updates the
  // micromaterialgp_static data!

  // save calculated displacements
  //dism_ = microstatic->ReturnNewDism();

  // clear displacements in MicroStruGenAlpha for next usage
  microstatic->ClearState();
}

#endif
