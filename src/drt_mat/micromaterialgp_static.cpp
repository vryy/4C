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
#include "../drt_inpar/inpar_structure.H"

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


std::map<int, RefCountPtr<STRUMULTI::MicroStatic> > MAT::MicroMaterialGP::microstaticmap_;
std::map<int, int> MAT::MicroMaterialGP::microstaticcounter_;


/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterialGP::MicroMaterialGP(const int gp, const int ele_ID, const bool eleowner,
                                      const double time, const int microdisnum, const double V0)
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
    microstaticmap_[microdisnum_] = rcp(new STRUMULTI::MicroStatic(microdisnum_, V0));
    // create a counter of macroscale GP associated with this "time integration" class
    // note that the counter is immediately updated afterwards!
    microstaticcounter_[microdisnum_] = 0;
  }

  microstaticcounter_[microdisnum] += 1;

  // create and initialize "empty" EAS history map
  EasInit();

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
  size_t pos = microprefix.rfind('-');
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

  // figure out how the new output file is called on the microscale
  // note: the trailing number must be the same as on the macroscale
  std::string newfilename;
  size_t posn = micronewprefix.rfind('-');
  if (posn!=std::string::npos)
  {
    std::string number = micronewprefix.substr(posn+1);
    std::string prefix = micronewprefix.substr(0,posn);
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

  if (eleowner)
  {
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

  // we are using some parameters from the macroscale input file
  // (e.g. time step size, alphaf etc. which need to be consistent in
  // both micro- and macroscale input file) whereas individual
  // parameters for the microscale can be used e.g. wrt output
  // options, kind of predictor etc.
  const Teuchos::ParameterList& sdyn_macro = DRT::Problem::Instance()->StructuralDynamicParams();

  // Initialize SurfStressManager for handling surface stress conditions due to interfacial phenomena
  // Note that this has to be done _after_ finding the output file name!
  // Note also that we are using the macroscale parameterlist here
  // because the SurfStressManager needs to know alphaf here
  surf_stress_man_=rcp(new UTILS::SurfStressManager(microdis, sdyn_macro, newfilename));

  istep_ = 0;

  if (eleowner) micro_output_->WriteMesh(istep_, microstaticmap_[microdisnum_]->GetTime());


  // we set the microscale time _ALWAYS_ to 0 (also in restart case)
  // in order to handle the implicit update query!!!
  timen_ = 0.;
  dt_    = sdyn_macro.get<double>("TIMESTEP");

  // check whether we are using modified Newton as a nonlinear solver
  // on the macroscale or not
  if (Teuchos::getIntegralValue<INPAR::STR::NonlinSolTech>(sdyn_macro,"NLNSOL")==INPAR::STR::soltech_newtonmod)
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
                                                  const double dt,
                                                  const bool eleowner)
{
  // select corresponding "time integration class" for this microstructure
  RefCountPtr<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  // check if we have to update absolute time and step number
  // in case of restart, timen_ is set to the current total time in
  // the constructor, whereas "time" handed on to this function is
  // still 0. since the macroscopic update of total time is done after
  // the initializing phase in which this function is called the first
  // time.
  // note that actually time and timen_ need to be _exactly_ the same
  // while we are in the same macro time step, since timen_ is always
  // a copy of time.
  if (fabs(time - timen_) > 0.5*dt)  // this is an arbitrarily chosen bound
  {
    build_stiff_ = true;

    if (timen_ > 0.5*dt)  // confer above comment (actually if (timen_!=0.))
    {
      microstatic->UpdateNewTimeStep(dis_, dism_, disn_, oldalpha_, lastalpha_, surf_stress_man_);

      // set displacements and EAS data of last step
      microstatic->SetOldState(dis_, dism_, disn_, surf_stress_man_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);

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
      if (eleowner)
      {
        microstatic->Output(micro_output_, timen_, istep_, dt_);
      }
      timen_ = time;
      dt_ = dt;
      istep_++;
    }
    else
    {
      // set displacements and EAS data of last step
      microstatic->SetOldState(dis_, dism_, disn_, surf_stress_man_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);
      timen_ = time;
      dt_ = dt;
      istep_++;
    }
  }
  else
  {
    // set displacements and EAS data of last step
    microstatic->SetOldState(dis_, dism_, disn_, surf_stress_man_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);
  }

  // set current absolute time, time step size and step number
  microstatic->SetTime(timen_, dt_, istep_);

  microstatic->Predictor(defgrd);
  microstatic->FullNewton();
  microstatic->StaticHomogenization(stress, cmat, density, defgrd, mod_newton_, build_stiff_);

  // note that it is not necessary to save displacements and EAS data
  // explicitly since we dealt with RCP's -> any update in class
  // microstatic and the elements, respectively, inherently updates the
  // micromaterialgp_static data!

  // clear displacements in MicroStruGenAlpha for next usage
  microstatic->ClearState();
}

#endif
