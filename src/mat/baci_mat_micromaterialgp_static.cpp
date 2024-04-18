/*----------------------------------------------------------------------*/
/*! \file
\brief
class for handling of micro-macro transitions

\level 3


*----------------------------------------------------------------------*/


#include "baci_mat_micromaterialgp_static.hpp"

#include "baci_global_data.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_io.hpp"
#include "baci_io_control.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_so3_hex8.hpp"
#include "baci_so3_shw6.hpp"
#include "baci_stru_multi_microstatic.hpp"

#include <filesystem>

FOUR_C_NAMESPACE_OPEN


std::map<int, Teuchos::RCP<STRUMULTI::MicroStatic>> MAT::MicroMaterialGP::microstaticmap_;
std::map<int, int> MAT::MicroMaterialGP::microstaticcounter_;


/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterialGP::MicroMaterialGP(
    const int gp, const int ele_ID, const bool eleowner, const int microdisnum, const double V0)
    : gp_(gp), ele_ID_(ele_ID), microdisnum_(microdisnum)
{
  GLOBAL::Problem* microproblem = GLOBAL::Problem::Instance(microdisnum_);
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->GetDis("structure");
  dis_ = CORE::LINALG::CreateVector(*microdis->DofRowMap(), true);
  disn_ = CORE::LINALG::CreateVector(*microdis->DofRowMap(), true);
  lastalpha_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>);
  oldalpha_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>);
  oldfeas_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>);
  oldKaainv_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>);
  oldKda_ = Teuchos::rcp(new std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>);

  // data must be consistent between micro and macro input file
  const Teuchos::ParameterList& sdyn_macro = GLOBAL::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& sdyn_micro = microproblem->StructuralDynamicParams();

  dt_ = sdyn_macro.get<double>("TIMESTEP");
  microdis->Comm().Broadcast(&dt_, 1, 0);
  step_ = 0;
  stepn_ = step_ + 1;
  time_ = 0.;
  timen_ = time_ + dt_;

  // if class handling microscale simulations is not yet initialized
  // -> set up

  if (microstaticmap_.find(microdisnum_) == microstaticmap_.end() or
      microstaticmap_[microdisnum_] == Teuchos::null)
  {
    // create "time integration" class for this microstructure
    microstaticmap_[microdisnum_] = Teuchos::rcp(new STRUMULTI::MicroStatic(microdisnum_, V0));
    // create a counter of macroscale GP associated with this "time integration" class
    // note that the counter is immediately updated afterwards!
    microstaticcounter_[microdisnum_] = 0;
  }

  microstaticcounter_[microdisnum] += 1;
  density_ = (microstaticmap_[microdisnum_])->Density();

  // create and initialize "empty" EAS history map (if necessary)
  EasInit();

  std::string newfilename;
  NewResultFile(eleowner, newfilename);

  // check whether we are using modified Newton as a nonlinear solver
  // on the macroscale or not
  if (CORE::UTILS::IntegralValue<INPAR::STR::NonlinSolTech>(sdyn_micro, "NLNSOL") ==
      INPAR::STR::soltech_newtonmod)
    mod_newton_ = true;
  else
    mod_newton_ = false;

  build_stiff_ = true;
}

/// destructor

MAT::MicroMaterialGP::~MicroMaterialGP()
{
  microstaticcounter_[microdisnum_] -= 1;
  if (microstaticcounter_[microdisnum_] == 0) microstaticmap_[microdisnum_] = Teuchos::null;
}


/// Read restart

void MAT::MicroMaterialGP::ReadRestart()
{
  step_ = GLOBAL::Problem::Instance()->Restart();
  microstaticmap_[microdisnum_]->ReadRestart(step_, dis_, lastalpha_, restartname_);

  *oldalpha_ = *lastalpha_;

  disn_->Update(1.0, *dis_, 0.0);
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
  Teuchos::RCP<IO::OutputControl> macrocontrol = GLOBAL::Problem::Instance(0)->OutputControlFile();
  std::string microprefix = macrocontrol->RestartName();
  std::string micronewprefix = macrocontrol->NewOutputFileName();

  GLOBAL::Problem* microproblem = GLOBAL::Problem::Instance(microdisnum_);
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->GetDis("structure");

  if (microdis->Comm().MyPID() == 0)
  {
    // figure out prefix of micro-scale restart files
    restartname_ = NewResultFilePath(microprefix);

    // figure out new prefix for micro-scale output files
    newfilename = NewResultFilePath(micronewprefix);
  }

  // restart file name and new output file name are sent to supporting procs
  if (microdis->Comm().NumProc() > 1)
  {
    {
      // broadcast restartname_ for micro scale
      int length = restartname_.length();
      std::vector<int> name(restartname_.begin(), restartname_.end());
      int err = microdis->Comm().Broadcast(&length, 1, 0);
      if (err) dserror("communication error");
      name.resize(length);
      err = microdis->Comm().Broadcast(name.data(), length, 0);
      if (err) dserror("communication error");
      restartname_.assign(name.begin(), name.end());
    }

    {
      // broadcast newfilename for micro scale
      int length = newfilename.length();
      std::vector<int> name(newfilename.begin(), newfilename.end());
      int err = microdis->Comm().Broadcast(&length, 1, 0);
      if (err) dserror("communication error");
      name.resize(length);
      err = microdis->Comm().Broadcast(name.data(), length, 0);
      if (err) dserror("communication error");
      newfilename.assign(name.begin(), name.end());
    }
  }

  if (eleowner)
  {
    const int ndim = GLOBAL::Problem::Instance()->NDim();
    const int restart = GLOBAL::Problem::Instance()->Restart();
    bool adaptname = true;
    // in case of restart, the new output file name is already adapted
    if (restart) adaptname = false;

    Teuchos::RCP<IO::OutputControl> microcontrol =
        Teuchos::rcp(new IO::OutputControl(microdis->Comm(), "Structure",
            microproblem->SpatialApproximationType(), "micro-input-file-not-known", restartname_,
            newfilename, ndim, restart, macrocontrol->FileSteps(),
            CORE::UTILS::IntegralValue<bool>(microproblem->IOParams(), "OUTPUT_BIN"), adaptname));

    micro_output_ = Teuchos::rcp(new IO::DiscretizationWriter(
        microdis, microcontrol, microproblem->SpatialApproximationType()));
    micro_output_->SetOutput(microcontrol);

    micro_output_->WriteMesh(step_, time_);
  }
}

std::string MAT::MicroMaterialGP::NewResultFilePath(const std::string& newprefix)
{
  std::string newfilename;

  // create path from string to extract only filename prefix
  const std::filesystem::path path(newprefix);
  const std::string newfileprefix = path.filename().string();

  const size_t posn = newfileprefix.rfind('-');
  if (posn != std::string::npos)
  {
    std::string number = newfileprefix.substr(posn + 1);
    std::string prefix = newfileprefix.substr(0, posn);

    // recombine path and file
    const std::filesystem::path parent_path(path.parent_path());
    const std::filesystem::path filen_name(prefix);
    const std::filesystem::path recombined_path = parent_path / filen_name;

    std::ostringstream s;
    s << recombined_path.string() << "_el" << ele_ID_ << "_gp" << gp_ << "-" << number;
    newfilename = s.str();
  }
  else
  {
    std::ostringstream s;
    s << newprefix << "_el" << ele_ID_ << "_gp" << gp_;
    newfilename = s.str();
  }
  return newfilename;
}


void MAT::MicroMaterialGP::EasInit()
{
  Teuchos::RCP<DRT::Discretization> discret =
      (GLOBAL::Problem::Instance(microdisnum_))->GetDis("structure");

  for (int lid = 0; lid < discret->ElementRowMap()->NumMyElements(); ++lid)
  {
    DRT::Element* actele = discret->lRowElement(lid);

    if (actele->ElementType() == DRT::ELEMENTS::SoHex8Type::Instance() or
        actele->ElementType() == DRT::ELEMENTS::SoShw6Type::Instance())
    {
      // create the parameters for the discretization
      Teuchos::ParameterList p;
      // action for elements
      p.set("action", "multi_eas_init");
      p.set("lastalpha", lastalpha_);
      p.set("oldalpha", oldalpha_);
      p.set("oldfeas", oldfeas_);
      p.set("oldKaainv", oldKaainv_);
      p.set("oldKda", oldKda_);

      CORE::LINALG::SerialDenseMatrix elematrix1;
      CORE::LINALG::SerialDenseMatrix elematrix2;
      CORE::LINALG::SerialDenseVector elevector1;
      CORE::LINALG::SerialDenseVector elevector2;
      CORE::LINALG::SerialDenseVector elevector3;
      std::vector<int> lm;

      actele->Evaluate(p, *discret, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
    }
  }

  return;
}



void MAT::MicroMaterialGP::ResetTimeAndStep()
{
  time_ = 0.0;
  timen_ = time_ + dt_;
  step_ = 0;
  stepn_ = step_ + 1;
}


/// perform microscale simulation

void MAT::MicroMaterialGP::PerformMicroSimulation(CORE::LINALG::Matrix<3, 3>* defgrd,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat)
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  // set displacements and EAS data of last step
  microstatic->SetState(dis_, disn_, stress_, strain_, plstrain_, lastalpha_, oldalpha_, oldfeas_,
      oldKaainv_, oldKda_);

  // set current time, time step size and step number
  microstatic->SetTime(time_, timen_, dt_, step_, stepn_);

  microstatic->Predictor(defgrd);
  microstatic->FullNewton();
  microstatic->StaticHomogenization(stress, cmat, defgrd, mod_newton_, build_stiff_);

  // note that it is not necessary to save displacements and EAS data
  // explicitly since we dealt with Teuchos::RCP's -> any update in class
  // microstatic and the elements, respectively, inherently updates the
  // micromaterialgp_static data!

  // clear displacements in MicroStruGenAlpha for next usage
  microstatic->ClearState();
}


void MAT::MicroMaterialGP::Update()
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  time_ = timen_;
  timen_ += dt_;
  step_ = stepn_;
  stepn_++;

  dis_->Update(1.0, *disn_, 0.0);

  GLOBAL::Problem* microproblem = GLOBAL::Problem::Instance(microdisnum_);
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->GetDis("structure");
  const Epetra_Map* elemap = microdis->ElementRowMap();

  for (int i = 0; i < elemap->NumMyElements(); ++i) (*lastalpha_)[i] = (*oldalpha_)[i];

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

  microstatic->SetState(dis_, disn_, stress_, strain_, plstrain_, lastalpha_, oldalpha_, oldfeas_,
      oldKaainv_, oldKda_);
  microstatic->SetTime(time_, timen_, dt_, step_, stepn_);
  microstatic->PrepareOutput();
}


void MAT::MicroMaterialGP::Output()
{
  // select corresponding "time integration class" for this microstructure
  Teuchos::RCP<STRUMULTI::MicroStatic> microstatic = microstaticmap_[microdisnum_];

  // set displacements and EAS data of last step
  microstatic->SetState(dis_, disn_, stress_, strain_, plstrain_, lastalpha_, oldalpha_, oldfeas_,
      oldKaainv_, oldKda_);
  microstatic->Output(micro_output_, time_, step_, dt_);

  // we don't need these containers anymore
  stress_ = Teuchos::null;
  strain_ = Teuchos::null;
  plstrain_ = Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
