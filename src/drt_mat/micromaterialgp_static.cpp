/*!----------------------------------------------------------------------
\file micromaterialgp_static.cpp

<pre>
Maintainer: Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef D_SOLID3

#include "micromaterialgp_static.H"
#include "../drt_stru_multi/microstatic.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_shw6.H"
#include "../drt_inpar/inpar_structure.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


std::map<int, Teuchos::RCP<STRUMULTI::MicroStatic> > MAT::MicroMaterialGP::microstaticmap_;
std::map<int, int> MAT::MicroMaterialGP::microstaticcounter_;


/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterialGP::MicroMaterialGP(const int gp, const int ele_ID, const bool eleowner,
                                      const int microdisnum, const double V0)
  : gp_(gp),
    ele_ID_(ele_ID),
    microdisnum_(microdisnum)
{
  Teuchos::RCP<DRT::Problem> microproblem = DRT::Problem::Instance(microdisnum_);
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->Dis(0, 0);
  dism_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  disn_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  dis_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  lastalpha_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  oldalpha_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  oldfeas_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  oldKaainv_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  oldKda_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);

  // we are using some parameters from the macroscale input file
  // (e.g. time step size, alphaf etc. which need to be consistent in
  // both micro- and macroscale input file) whereas individual
  // parameters for the microscale can be used e.g. wrt output
  // options, kind of predictor etc.
  const Teuchos::ParameterList& sdyn_macro = DRT::Problem::Instance()->StructuralDynamicParams();

  dt_    = sdyn_macro.get<double>("TIMESTEP");
  step_  = 0;
  stepn_ = step_ + 1;
  time_  = 0.;
  timen_ = time_ + dt_;

  // if class handling microscale simulations is not yet initialized
  // -> set up

  if (microstaticmap_.find(microdisnum_) == microstaticmap_.end() or microstaticmap_[microdisnum_] == Teuchos::null)
  {
    // create "time integration" class for this microstructure
    microstaticmap_[microdisnum_] = Teuchos::rcp(new STRUMULTI::MicroStatic(microdisnum_, V0));
    // create a counter of macroscale GP associated with this "time integration" class
    // note that the counter is immediately updated afterwards!
    microstaticcounter_[microdisnum_] = 0;
  }

  microstaticcounter_[microdisnum] += 1;

  // create and initialize "empty" EAS history map (if necessary)
  EasInit();

  std::string newfilename;
  NewResultFile(eleowner, newfilename);

  // Initialize SurfStressManager for handling surface stress conditions due to interfacial phenomena
  // Note that this has to be done _after_ finding the output file name!
  // Note also that we are using the macroscale parameterlist here
  // because the SurfStressManager needs to know alphaf here
  surf_stress_man_ = Teuchos::rcp(new UTILS::SurfStressManager(microdis, sdyn_macro, newfilename));

  // check whether we are using modified Newton as a nonlinear solver
  // on the macroscale or not
  if (DRT::INPUT::IntegralValue<INPAR::STR::NonlinSolTech>(sdyn_macro,"NLNSOL")==INPAR::STR::soltech_newtonmod)
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
  step_ = genprob.restart;
  microstaticmap_[microdisnum_]->ReadRestart(step_, dis_, lastalpha_, surf_stress_man_, restartname_);
  // both dism_ and disn_ equal dis_
  dism_->Update(1.0, *dis_, 0.0);
  disn_->Update(1.0, *dis_, 0.0);
  // both lastalpha and oldalpha are the same
  *oldalpha_ = *lastalpha_;
}


/// New resultfile

void MAT::MicroMaterialGP::NewResultFile(bool eleowner, std::string& newfilename)
{
  // set up micro output
  //
  // Get the macro output prefix and insert element and gauss point
  // identifier. We use the original name here and rely on our (micro)
  // OutputControl object below to act just like the macro (default)
  // OutputControl. In particular we assume that there are always micro and
  // macro control files on restart.
  Teuchos::RCP<IO::OutputControl> macrocontrol = DRT::Problem::Instance(0)->OutputControlFile();
  std::string microprefix = macrocontrol->RestartName();
  std::string micronewprefix = macrocontrol->NewOutputFileName();

  Teuchos::RCP<DRT::Problem> microproblem = DRT::Problem::Instance(microdisnum_);
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->Dis(0, 0);

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
    Teuchos::RCP<IO::OutputControl> microcontrol =
      Teuchos::rcp(new IO::OutputControl(microdis->Comm(),
                            "structure",
                            microproblem->SpatialApproximation(),
                            "micro-input-file-not-known",
                            restartname_,
                            newfilename,
                            genprob.ndim,
                            genprob.restart,
                            macrocontrol->FileSteps()));

    micro_output_ = Teuchos::rcp(new IO::DiscretizationWriter(microdis,microcontrol));

    micro_output_->WriteMesh(step_, time_);
  }

  return;
}


void MAT::MicroMaterialGP::EasInit()
{
  Teuchos::RCP<DRT::Discretization> discret = (DRT::Problem::Instance(microdisnum_))->Dis(0, 0);

  for (int lid=0; lid<discret->ElementRowMap()->NumMyElements(); ++lid)
  {
    DRT::Element* actele = discret->lRowElement(lid);

    if (actele->ElementType()==DRT::ELEMENTS::So_hex8Type::Instance() or
        actele->ElementType()==DRT::ELEMENTS::So_shw6Type::Instance())
    {
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","multi_eas_init");
      p.set("lastalpha", lastalpha_);
      p.set("oldalpha", oldalpha_);
      p.set("oldfeas", oldfeas_);
      p.set("oldKaainv", oldKaainv_);
      p.set("oldKda", oldKda_);

      Epetra_SerialDenseMatrix elematrix1;
      Epetra_SerialDenseMatrix elematrix2;
      Epetra_SerialDenseVector elevector1;
      Epetra_SerialDenseVector elevector2;
      Epetra_SerialDenseVector elevector3;
      vector<int> lm;

      actele->Evaluate(p,*discret,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
    }
  }

  return;
}



void MAT::MicroMaterialGP::ResetTimeAndStep()
{
  time_  = 0.0;
  timen_ = time_ + dt_;
  step_  = 0;
  stepn_ = step_ + 1;
}


/// perform microscale simulation

void MAT::MicroMaterialGP::PerformMicroSimulation(LINALG::Matrix<3,3>* defgrd,
                                                  LINALG::Matrix<6,1>* stress,
                                                  LINALG::Matrix<6,6>* cmat,
                                                  double* density)
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  // set displacements and EAS data of last step
  microstatic->SetState(dis_, dism_, disn_, surf_stress_man_, stress_, strain_, plstrain_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);

  // set current time, time step size and step number
  microstatic->SetTime(time_, timen_, dt_, step_, stepn_);

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


void MAT::MicroMaterialGP::Update()
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  microstatic->UpdateNewTimeStep(dis_, dism_, disn_, oldalpha_, lastalpha_, surf_stress_man_);

  time_  = timen_;
  timen_ += dt_;
  step_  = stepn_;
  stepn_++;

  // in case of modified Newton, the stiffness matrix needs to be rebuilt at
  // the beginning of the new time step
  build_stiff_ = true;
}


void MAT::MicroMaterialGP::PrepareOutput()
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  stress_ = Teuchos::rcp(new std::vector<char>());
  strain_ = Teuchos::rcp(new std::vector<char>());
  plstrain_ = Teuchos::rcp(new std::vector<char>());

  microstatic->SetState(dis_, dism_, disn_, surf_stress_man_, stress_, strain_, plstrain_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);
  microstatic->SetTime(time_, timen_, dt_, step_, stepn_);
  microstatic->PrepareOutput();
}


void MAT::MicroMaterialGP::Output()
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  // set displacements and EAS data of last step
  microstatic->SetState(dis_, dism_, disn_, surf_stress_man_, stress_, strain_, plstrain_, lastalpha_, oldalpha_, oldfeas_, oldKaainv_, oldKda_);
  microstatic->Output(micro_output_, time_, step_, dt_);

  // we don't need these containers anymore
  stress_ = Teuchos::null;
  strain_ = Teuchos::null;
  plstrain_ = Teuchos::null;
}

#endif //#ifdef D_SOLID3
#endif
