/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all structural time integration strategies.


\level 3

*/
/*-----------------------------------------------------------*/


#include "baci_structure_new_timint_base.H"

#include "baci_beaminteraction_str_model_evaluator.H"
#include "baci_comm_utils.H"
#include "baci_inpar_contact.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_gmsh.H"
#include "baci_io_pstream.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_blocksparsematrix.H"
#include "baci_structure_new_dbc.H"
#include "baci_structure_new_enum_lists.H"
#include "baci_structure_new_factory.H"
#include "baci_structure_new_integrator.H"
#include "baci_structure_new_model_evaluator_data.H"
#include "baci_structure_new_model_evaluator_factory.H"
#include "baci_structure_new_resulttest.H"
#include "baci_structure_new_timint_basedataio_monitor_dbc.H"
#include "baci_structure_new_timint_basedataio_runtime_vtk_output.H"
#include "baci_structure_new_timint_basedataio_runtime_vtp_output.H"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Base::Base()
    : StructureNew(),
      isinit_(false),
      issetup_(false),
      isrestarting_(false),
      state_is_insync_with_noxgroup_(true),
      dataio_(Teuchos::null),
      datasdyn_(Teuchos::null),
      dataglobalstate_(Teuchos::null),
      int_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null)
{
  Epetra_Object::SetTracebackMode(1);
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Init(const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate)
{
  // ---------------------------------------------------------------------------
  // We need to call Setup() after Init()
  // ---------------------------------------------------------------------------
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // initialize the data container ptrs
  // ---------------------------------------------------------------------------
  dataio_ = dataio;
  datasdyn_ = datasdyn;
  dataglobalstate_ = dataglobalstate;

  // ---------------------------------------------------------------------------
  // set isInit flag
  // ---------------------------------------------------------------------------
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Setup()
{
  CheckInit();

  // ---------------------------------------------------------------------------
  // Create the Dirichlet Boundary Condition handler
  // ---------------------------------------------------------------------------
  dbc_ptr_ = STR::BuildDbc(DataSDyn());
  /* FixMe It would be sufficient to use a constant discretization,
   * unfortunately this wasn't considered during the implementation of the
   * discretization routines. Therefore many methods need a slight modification
   * (most times adding a "const" should fix the problem).          hiermeier */
  Teuchos::RCP<DRT::Discretization> discret_ptr = DataGlobalState().GetDiscret();
  dbc_ptr_->Init(discret_ptr, DataGlobalState().GetFreactNp(), Teuchos::rcp(this, false));
  dbc_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // Create the explicit/implicit integrator
  // ---------------------------------------------------------------------------
  int_ptr_ = STR::BuildIntegrator(DataSDyn());
  int_ptr_->Init(
      DataSDynPtr(), DataGlobalStatePtr(), DataIOPtr(), dbc_ptr_, Teuchos::rcp(this, false));
  int_ptr_->Setup();
  int_ptr_->PostSetup();
  // Initialize and Setup the input/output writer for every Newton iteration
  dataio_->InitSetupEveryIterationWriter(this, DataSDyn().GetNoxParams());

  // Initialize the output of system energy
  if (dataio_->GetWriteEnergyEveryNStep())
  {
    SelectEnergyTypesToBeWritten();

    if (dataglobalstate_->GetMyRank() == 0) InitializeEnergyFileStreamAndWriteHeaders();
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Reset()
{
  dserror(
      "Reset of all class variables is not yet implemented for "
      "the modelevaluator!");
  // ModelEvaluator().Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResetStep()
{
  CheckInitSetup();

  int_ptr_->ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::Base::NotFinished() const
{
  CheckInitSetup();

  // check for early stopping
  const bool early_stop = int_ptr_->EarlyStopping();
  if (early_stop)
  {
    // Simulation is finished regardless of the simulation time or the time step.
    return false;
  }

  // check the current time
  const double& timenp = dataglobalstate_->GetTimeNp();
  const double& timemax = datasdyn_->GetTimeMax();
  const double& dt = (*dataglobalstate_->GetDeltaTime())[0];
  // check the step counter
  const int& stepnp = dataglobalstate_->GetStepNp();
  const int& stepmax = datasdyn_->GetStepMax();

  return (timenp <= timemax + 1.0e-8 * dt and stepnp <= stepmax);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SetRestart(int stepn, double timen, Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
    Teuchos::RCP<std::vector<char>> elementdata, Teuchos::RCP<std::vector<char>> nodedata)
{
  CheckInitSetup();

  dserror("SetRestartState() is deprecated, use the ReadRestart() routine instead!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& STR::TIMINT::Base::GetMassDomainMap() const
{
  CheckInitSetup();
  return dataglobalstate_->GetMassMatrix()->DomainMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MapExtractor> STR::TIMINT::Base::GetDBCMapExtractor()
{
  CheckInitSetup();
  return dbc_ptr_->GetDBCMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MapExtractor> STR::TIMINT::Base::GetDBCMapExtractor() const
{
  CheckInitSetup();
  return dbc_ptr_->GetDBCMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::LocsysManager> STR::TIMINT::Base::LocsysManager()
{
  CheckInitSetup();
  return dbc_ptr_->LocSysManagerPtr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Generic& STR::TIMINT::Base::ModelEvaluator(
    INPAR::STR::ModelType mtype) const
{
  return Integrator().ModelEval().Evaluator(mtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::TIMINT::Base::ModelEvaluator(INPAR::STR::ModelType mtype)
{
  return Integrator().ModelEval().Evaluator(mtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::Base::TimIntParam() const
{
  CheckInitSetup();
  return int_ptr_->GetIntParam();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResizeMStepTimAda()
{
  CheckInitSetup();
  // resize time and stepsize fields
  const double& timen = dataglobalstate_->GetTimeN();
  dataglobalstate_->GetMultiTime()->Resize(-1, 0, timen);
  const double& dtn = (*dataglobalstate_->GetDeltaTime())[0];
  dataglobalstate_->GetDeltaTime()->Resize(-1, 0, dtn);

  // resize state vectors, AB2 is a 2-step method, thus we need two
  // past steps at t_{n} and t_{n-1}
  const Epetra_Map* dofrowmap_ptr = dataglobalstate_->DofRowMapView();
  dataglobalstate_->GetMultiDis()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->GetMultiVel()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->GetMultiAcc()->Resize(-1, 0, dofrowmap_ptr, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Update()
{
  CheckInitSetup();
  int_ptr_->PreUpdate();
  int_ptr_->UpdateStructuralEnergy();
  int_ptr_->UpdateStepState();
  UpdateStepTime();
  SetNumberOfNonlinearIterations();
  int_ptr_->UpdateStepElement();
  int_ptr_->PostUpdate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::UpdateStepTime()
{
  CheckInitSetup();
  double& timenp = dataglobalstate_->GetTimeNp();
  int& stepnp = dataglobalstate_->GetStepNp();
  int& stepn = dataglobalstate_->GetStepN();

  // --------------------------------------------------------------------------
  // update old time and step variables
  // --------------------------------------------------------------------------
  dataglobalstate_->GetMultiTime()->UpdateSteps(timenp);
  stepn = stepnp;

  // --------------------------------------------------------------------------
  // update the new time and step variables
  // --------------------------------------------------------------------------
  // get current time step size
  const double& dtn = (*dataglobalstate_->GetDeltaTime())[0];
  dataglobalstate_->GetDeltaTime()->UpdateSteps(dtn);
  timenp += dtn;
  stepnp += 1;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SetNumberOfNonlinearIterations()
{
  int nlniter = 0;

  if (DataSDyn().GetNoxParams().isSublist("Output"))
  {
    const Teuchos::ParameterList& nox_output = DataSDyn().GetNoxParams().sublist("Output");
    if (nox_output.isParameter("Nonlinear Iterations"))
      nlniter = nox_output.get<int>("Nonlinear Iterations");
  }

  dataglobalstate_->SetNlnIterationNumber(nlniter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SelectEnergyTypesToBeWritten()
{
  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();

  // decide which types of energy contributions shall be written separately
  const std::set<enum INPAR::STR::ModelType>& mtypes = datasdyn_->GetModelTypes();

  std::set<enum INPAR::STR::ModelType>::const_iterator model_iter;
  for (model_iter = mtypes.begin(); model_iter != mtypes.end(); ++model_iter)
  {
    switch (*model_iter)
    {
      case INPAR::STR::model_structure:
      {
        evaldata.InsertEnergyTypeToBeConsidered(STR::internal_energy);
        evaldata.InsertEnergyTypeToBeConsidered(STR::kinetic_energy);
        break;
      }
      case INPAR::STR::model_beaminteraction:
      {
        STR::MODELEVALUATOR::BeamInteraction const beaminteraction_evaluator =
            dynamic_cast<STR::MODELEVALUATOR::BeamInteraction const&>(
                int_ptr_->ModelEvalPtr()->Evaluator(INPAR::STR::model_beaminteraction));

        if (beaminteraction_evaluator.HaveSubModelType(
                INPAR::BEAMINTERACTION::submodel_beamcontact))
        {
          evaldata.InsertEnergyTypeToBeConsidered(STR::beam_contact_penalty_potential);
        }
        if (beaminteraction_evaluator.HaveSubModelType(INPAR::BEAMINTERACTION::submodel_potential))
        {
          evaldata.InsertEnergyTypeToBeConsidered(STR::beam_interaction_potential);
        }
        if (beaminteraction_evaluator.HaveSubModelType(
                INPAR::BEAMINTERACTION::submodel_crosslinking))
        {
          evaldata.InsertEnergyTypeToBeConsidered(STR::beam_to_beam_link_internal_energy);
          evaldata.InsertEnergyTypeToBeConsidered(STR::beam_to_beam_link_kinetic_energy);
        }
        if (beaminteraction_evaluator.HaveSubModelType(
                INPAR::BEAMINTERACTION::submodel_spherebeamlink))
        {
          evaldata.InsertEnergyTypeToBeConsidered(STR::beam_to_sphere_link_internal_energy);
          evaldata.InsertEnergyTypeToBeConsidered(STR::beam_to_sphere_link_kinetic_energy);
        }
        break;
      }
      default:
      {
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::InitializeEnergyFileStreamAndWriteHeaders()
{
  auto& evaldata = int_ptr_->EvalData();

  dataio_->SetupEnergyOutputFile();

  // write column headers to file
  dataio_->GetEnergyOutputStream() << std::setw(12) << "#timestep," << std::setw(24) << "time,";

  for (const auto& energy_data : evaldata.GetEnergyData())
  {
    dataio_->GetEnergyOutputStream()
        << std::setw(36) << STR::EnergyType2String(energy_data.first) + ",";
  }

  dataio_->GetEnergyOutputStream() << std::setw(24) << "total_energy" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> STR::TIMINT::Base::CreateFieldTest()
{
  CheckInitSetup();
  Teuchos::RCP<STR::ResultTest> resulttest = Teuchos::rcp(new STR::ResultTest());
  resulttest->Init(GetDataGlobalState(), Integrator().EvalData());
  resulttest->Setup();

  return resulttest;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
    Teuchos::RCP<Epetra_Vector> disnp, Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> accnp, Teuchos::RCP<std::vector<char>> elementdata,
    Teuchos::RCP<std::vector<char>> nodedata)
{
  CheckInitSetup();
  // at some point we have to create a copy
  *step = dataglobalstate_->GetStepN();
  *time = dataglobalstate_->GetTimeN();
  Teuchos::RCP<const DRT::Discretization> discret_ptr =
      Teuchos::rcp_dynamic_cast<const DRT::Discretization>(dataglobalstate_->GetDiscret());
  *elementdata = *(discret_ptr->PackMyElements());
  *nodedata = *(discret_ptr->PackMyNodes());

  // get restart data is only for simple structure problems
  // hence if the model set is larger than one, we throw an error
  if (datasdyn_->GetModelTypes().size() > 1)
    dserror("The GetRestartData routine supports the structural model case ONLY!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::PrepareOutput(bool force_prepare_timestep)
{
  CheckInitSetup();
  // --- stress, strain and optional quantity calculation ---------------------
  if ((dataio_->IsWriteResultsEnabled() && force_prepare_timestep) ||
      dataio_->WriteResultsForThisStep(dataglobalstate_->GetStepNp()))
  {
    int_ptr_->DetermineStressStrain();
    int_ptr_->DetermineOptionalQuantity();

    if (dataio_->IsWriteCurrentEleVolume())
    {
      Teuchos::RCP<Epetra_Vector> elevolumes = Teuchos::null;
      Teuchos::RCP<const Epetra_Vector> disnp = dataglobalstate_->GetDisNp();

      int_ptr_->DetermineElementVolumes(*disnp, elevolumes);
      int_ptr_->EvalData().SetElementVolumeData(elevolumes);
    }
  }
  if ((dataio_->IsRuntimeVtkOutputEnabled() && force_prepare_timestep) ||
      dataio_->WriteRuntimeVtkResultsForThisStep(dataglobalstate_->GetStepNp()) ||
      dataio_->WriteRuntimeVtpResultsForThisStep(dataglobalstate_->GetStepNp()))
  {
    int_ptr_->RuntimePreOutputStepState();
  }
  // --- energy calculation ---------------------------------------------------
  if ((dataio_->GetWriteEnergyEveryNStep() and
          (force_prepare_timestep ||
              dataglobalstate_->GetStepNp() % dataio_->GetWriteEnergyEveryNStep() == 0)))
  {
    STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();
    evaldata.ClearValuesForAllEnergyTypes();

    int_ptr_->DetermineEnergy();

    // sum processor-local values of all separate contributions into global value
    double energy_local = 0.0;
    double energy_global = 0.0;

    for (const auto& energy_data : evaldata.GetEnergyData())
    {
      energy_local = energy_data.second;

      dataglobalstate_->GetComm().SumAll(&energy_local, &energy_global, 1);

      evaldata.SetValueForEnergyType(energy_global, energy_data.first);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Output(bool forced_writerestart)
{
  CheckInitSetup();
  OutputStep(forced_writerestart);
  // write Gmsh output
  WriteGmshStrucOutputStep();
  int_ptr_->PostOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputStep(bool forced_writerestart)
{
  CheckInitSetup();
  // special treatment is necessary when restart is forced
  if (forced_writerestart)
  {
    // reset possible history data on element level
    ResetStep();
    // restart has already been written or simulation has just started
    if (dataio_->ShouldWriteRestartForStep(dataglobalstate_->GetStepN()) or
        dataglobalstate_->GetStepN() == DRT::Problem::Instance()->Restart())
      return;
    // if state already exists, add restart information
    if (dataio_->WriteResultsForThisStep(dataglobalstate_->GetStepN()))
    {
      AddRestartToOutputState();
      return;
    }
  }

  /* This flag indicates whether some form of output has already been written in the current time
   * step. It is passed along subroutines and prevents repeated initialization of output writer,
   * printing of state vectors, or similar.
   */
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if (forced_writerestart || dataio_->ShouldWriteRestartForStep(dataglobalstate_->GetStepN()))
  {
    OutputRestart(datawritten);
    dataio_->SetLastWrittenResults(dataglobalstate_->GetStepN());
  }

  // output results (not necessary if restart in same step)
  if (dataio_->IsWriteState() and dataio_->WriteResultsForThisStep(dataglobalstate_->GetStepN()) and
      (not datawritten))
  {
    NewIOStep(datawritten);
    OutputState();
    dataio_->SetLastWrittenResults(dataglobalstate_->GetStepN());
  }

  // output results during runtime ( not used for restart so far )
  if (dataio_->WriteRuntimeVtkResultsForThisStep(dataglobalstate_->GetStepN()) or
      dataio_->WriteRuntimeVtpResultsForThisStep(dataglobalstate_->GetStepN()))
  {
    RuntimeOutputState();
  }

  // write reaction forces
  if (dataio_->ShouldWriteReactionForcesForThisStep(dataglobalstate_->GetStepN()))
  {
    OutputReactionForces();
  }

  // output stress, strain and optional quantity
  if (dataio_->ShouldWriteStressStrainForThisStep(dataglobalstate_->GetStepN()))
  {
    NewIOStep(datawritten);
    OutputStressStrain();
    OutputOptionalQuantity();
  }

  if (dataio_->WriteResultsForThisStep(dataglobalstate_->GetStepN()) and
      dataio_->IsWriteCurrentEleVolume())
  {
    NewIOStep(datawritten);
    IO::DiscretizationWriter& iowriter = *(dataio_->GetOutputPtr());
    OutputElementVolume(iowriter);
  }

  // output energy
  if (dataio_->ShouldWriteEnergyForThisStep(dataglobalstate_->GetStepN()))
  {
    OutputEnergy();
  }

  // print error norms
  OutputErrorNorms();

  //  OutputVolumeMass();

  // ToDo output of nodal positions in current configuration
  //  OutputNodalPositions();

  // ToDo write output on micro-scale (multi-scale analysis)
  //  if (HaveMicroMat())
  //    dserror("OutputMicro() is not yet implemented!"); // OutputMicro();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::NewIOStep(bool& datawritten)
{
  if (not datawritten)
  {
    // Make new step
    dataio_->GetOutputPtr()->NewStep(dataglobalstate_->GetStepN(), dataglobalstate_->GetTimeN());

    datawritten = true;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputState()
{
  CheckInitSetup();
  IO::DiscretizationWriter& iowriter = *(dataio_->GetOutputPtr());

  OutputState(iowriter, dataio_->IsFirstOutputOfRun());

  dataio_->SetFirstOutputOfRun(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputDebugState(IO::DiscretizationWriter& iowriter, bool write_owner) const
{
  OutputState(iowriter, write_owner);

  // write element volumes as additional debugging information, if activated
  if (dataio_->IsWriteCurrentEleVolume()) OutputElementVolume(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputState(IO::DiscretizationWriter& iowriter, bool write_owner) const
{
  // owner of elements is just written once because it does not change during
  // simulation (so far)
  iowriter.WriteElementData(write_owner);
  iowriter.WriteNodeData(write_owner);

  int_ptr_->OutputStepState(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::RuntimeOutputState()
{
  CheckInitSetup();
  int_ptr_->RuntimeOutputStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputReactionForces()
{
  CheckInitSetup();
  IO::DiscretizationWriter& iowriter = *(dataio_->GetOutputPtr());
  int_ptr_->MonitorDbc(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputElementVolume(IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();

  iowriter.WriteVector("current_ele_volumes",
      Teuchos::rcpFromRef(evaldata.CurrentElementVolumeData()), IO::elementvector);

  evaldata.SetElementVolumeData(Teuchos::null);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputStressStrain()
{
  CheckInitSetup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();

  // ---------------------------------------------------------------------------
  // write stress output
  // ---------------------------------------------------------------------------
  std::string text = "";
  if (dataio_->GetStressOutputType() != INPAR::STR::stress_none)
  {
    switch (dataio_->GetStressOutputType())
    {
      case INPAR::STR::stress_cauchy:
        text = "gauss_cauchy_stresses_xyz";
        break;
      case INPAR::STR::stress_2pk:
        text = "gauss_2PK_stresses_xyz";
        break;
      default:
        dserror("Requested stress type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.StressData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.StressDataPtr() = Teuchos::null;

  // ---------------------------------------------------------------------------
  // write coupling stress output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->GetCouplingStressOutputType() != INPAR::STR::stress_none)
  {
    switch (dataio_->GetCouplingStressOutputType())
    {
      case INPAR::STR::stress_cauchy:
        text = "gauss_cauchy_coupling_stresses_xyz";
        break;
      case INPAR::STR::stress_2pk:
        text = "gauss_2PK_coupling_stresses_xyz";
        break;
      default:
        dserror("Requested coupling stress type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.CouplingStressData(), *(Discretization()->ElementRowMap()));
  }

  // ---------------------------------------------------------------------------
  // write strain output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->GetStrainOutputType() != INPAR::STR::strain_none)
  {
    switch (dataio_->GetStrainOutputType())
    {
      case INPAR::STR::strain_none:
        break;
      case INPAR::STR::strain_ea:
        text = "gauss_EA_strains_xyz";
        break;
      case INPAR::STR::strain_gl:
        text = "gauss_GL_strains_xyz";
        break;
      case INPAR::STR::strain_log:
        text = "gauss_LOG_strains_xyz";
        break;
      default:
        dserror("Requested strain type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.StrainData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.StrainDataPtr() = Teuchos::null;

  // ---------------------------------------------------------------------------
  // write plastic strain output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->GetPlasticStrainOutputType() != INPAR::STR::strain_none)
  {
    switch (dataio_->GetPlasticStrainOutputType())
    {
      case INPAR::STR::strain_ea:
        text = "gauss_pl_EA_strains_xyz";
        break;
      case INPAR::STR::strain_gl:
        text = "gauss_pl_GL_strains_xyz";
        break;
      default:
        dserror("Requested plastic strain type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.PlasticStrainData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.PlasticStrainDataPtr() = Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputEnergy() const
{
  CheckInitSetup();

  if (dataglobalstate_->GetMyRank() == 0)
  {
    std::ostream& energy_output_stream = dataio_->GetEnergyOutputStream();

    energy_output_stream << std::setw(11) << dataglobalstate_->GetStepN() << std::setw(1) << ","
                         << std::scientific << std::setprecision(14) << std::setw(23)
                         << dataglobalstate_->GetTimeN() << std::setw(1) << ",";

    STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();

    double total_energy = 0.0;

    for (const auto& energy_data : evaldata.GetEnergyData())
    {
      energy_output_stream << std::setw(35) << energy_data.second << std::setw(1) << ",";
      total_energy += energy_data.second;
    }

    energy_output_stream << std::setw(24) << total_energy << std::endl;

    IO::cout(IO::verbose) << "\n\nOutput for energy written to file!" << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputOptionalQuantity()
{
  CheckInitSetup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();

  // ---------------------------------------------------------------------------
  // write optional quantity output
  // ---------------------------------------------------------------------------
  std::string text = "";
  if (dataio_->GetOptQuantityOutputType() != INPAR::STR::optquantity_none)
  {
    switch (dataio_->GetOptQuantityOutputType())
    {
      case INPAR::STR::optquantity_membranethickness:
        text = "gauss_membrane_thickness";
        break;
      default:
        dserror("Requested optional quantity type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.OptQuantityData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.OptQuantityDataPtr() = Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputRestart(bool& datawritten)
{
  CheckInitSetup();

  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();
  // write restart output, please
  if (dataglobalstate_->GetStepN() != 0)
    output_ptr->WriteMesh(dataglobalstate_->GetStepN(), dataglobalstate_->GetTimeN());
  NewIOStep(datawritten);

  output_ptr->WriteElementData(dataio_->IsFirstOutputOfRun());
  output_ptr->WriteNodeData(dataio_->IsFirstOutputOfRun());
  dataio_->SetFirstOutputOfRun(false);

  // add velocity and acceleration if necessary
  output_ptr->WriteVector("velocity", dataglobalstate_->GetVelN());
  output_ptr->WriteVector("acceleration", dataglobalstate_->GetAccN());

  /* Add the restart information of the different time integrators and model
   * evaluators. */
  int_ptr_->WriteRestart(*output_ptr);

  // info dedicated to user's eyes staring at standard out
  if ((dataglobalstate_->GetMyRank() == 0) and (dataio_->GetPrint2ScreenEveryNStep() > 0) and
      (StepOld() % dataio_->GetPrint2ScreenEveryNStep() == 0))
  {
    IO::cout << "====== Restart for field 'Structure' written in step "
             << dataglobalstate_->GetStepN() << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::AddRestartToOutputState()
{
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();

  // force output of velocity and acceleration in case it is not written previously by the model
  // evaluators
  if (!dataio_->IsWriteVelAcc())
  {
    output_ptr->WriteVector("velocity", dataglobalstate_->GetVelN());
    output_ptr->WriteVector("acceleration", dataglobalstate_->GetAccN());
  }

  /* Add the restart information of the different time integrators and model
   * evaluators. */
  int_ptr_->WriteRestart(*output_ptr, true);

  // finally add the missing mesh information, order is important here
  output_ptr->WriteMesh(DataGlobalState().GetStepN(), DataGlobalState().GetTimeN());

  // info dedicated to user's eyes staring at standard out
  if ((dataglobalstate_->GetMyRank() == 0) and (dataio_->GetPrint2ScreenEveryNStep() > 0) and
      (StepOld() % dataio_->GetPrint2ScreenEveryNStep() == 0))
  {
    IO::cout << "====== Restart for field 'Structure' written in step "
             << dataglobalstate_->GetStepN() << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::WriteGmshStrucOutputStep()
{
  CheckInitSetup();
  if (!dataio_->IsGmsh()) return;

  const std::string filename = IO::GMSH::GetFileName(
      "struct", dataglobalstate_->GetStepNp(), false, dataglobalstate_->GetMyRank());
  std::ofstream gmshfilecontent(filename.c_str());

  // add 'View' to Gmsh postprocessing file
  gmshfilecontent << "View \" "
                  << "struct displacement \" {" << std::endl;
  // draw vector field 'struct displacement' for every element
  IO::GMSH::VectorFieldDofBasedToGmsh(Discretization(), Dispn(), gmshfilecontent, 0, true);
  gmshfilecontent << "};" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ReadRestart(const int stepn)
{
  CheckInit();
  // that restarting flag
  isrestarting_ = true;

  // create an input/output reader
  IO::DiscretizationReader ioreader(Discretization(), stepn);
  dataglobalstate_->GetStepN() = stepn;
  dataglobalstate_->GetStepNp() = stepn + 1;
  dataglobalstate_->GetMultiTime() =
      Teuchos::rcp(new BACI::TIMINT::TimIntMStep<double>(0, 0, ioreader.ReadDouble("time")));
  const double& timen = dataglobalstate_->GetTimeN();
  const double& dt = (*dataglobalstate_->GetDeltaTime())[0];
  dataglobalstate_->GetTimeNp() = timen + dt;

  // ---------------------------------------------------------------------------
  // The order is important at this point!
  // (0) read element and node data --> new discretization
  // (1) Setup() the model evaluator and time integrator
  // (2) read and possibly overwrite the general dynamic state
  // (3) read specific time integrator and model evaluator data
  // ---------------------------------------------------------------------------
  // (0) read element and node data
  ioreader.ReadHistoryData(stepn);

  // (1) Setup() the model evaluator and time integrator
  /* Since we call a redistribution on the structural discretization, we have to
   * setup the structural time integration strategy at this point and not as
   * usually during the adapter call.                         hiermeier 05/16 */
  Setup();

  // (2) read (or overwrite) the general dynamic state
  Teuchos::RCP<Epetra_Vector>& velnp = dataglobalstate_->GetVelNp();
  ioreader.ReadVector(velnp, "velocity");
  dataglobalstate_->GetMultiVel()->UpdateSteps(*velnp);
  Teuchos::RCP<Epetra_Vector>& accnp = dataglobalstate_->GetAccNp();
  ioreader.ReadVector(accnp, "acceleration");
  dataglobalstate_->GetMultiAcc()->UpdateSteps(*accnp);

  // (3) read specific time integrator (forces, etc.) and model evaluator data
  int_ptr_->ReadRestart(ioreader);
  int_ptr_->PostSetup();  // compute here the equilibrium system to account for initial
                          // displacement/velocity.

  // short screen output
  if (dataglobalstate_->GetMyRank() == 0)
    IO::cout << "====== Restart of the structural simulation from step " << stepn << IO::endl;

  // end of restarting
  isrestarting_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> STR::TIMINT::Base::Discretization()
{
  CheckInit();
  return dataglobalstate_->GetDiscret();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SetActionType(const DRT::ELEMENTS::ActionType& action)
{
  CheckInitSetup();
  int_ptr_->EvalData().SetActionType(action);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Base::GroupId() const
{
  Teuchos::RCP<CORE::COMM::Communicators> group = DRT::Problem::Instance()->GetCommunicators();
  return group->GroupId();
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputErrorNorms()
{
  // get parameters
  const Teuchos::ParameterList& params = DRT::Problem::Instance()->StructuralDynamicParams();

  // get error calculation info
  const auto calcerr = DRT::INPUT::IntegralValue<INPAR::STR::CalcError>(params, "CALCERROR");

  switch (calcerr)
  {
    case INPAR::STR::no_error_calculation:
    {
      return;
    }
    case INPAR::STR::byfunct:
    {
      // initialize variables
      Teuchos::RCP<CORE::LINALG::SerialDenseVector> norms =
          Teuchos::rcp(new CORE::LINALG::SerialDenseVector(3));
      norms->putScalar(0.0);

      // call discretization to evaluate error norms
      Teuchos::ParameterList p;
      p.set("action", "calc_struct_errornorms");
      Discretization()->ClearState();
      Discretization()->SetState("displacement", Dispnp());
      Discretization()->EvaluateScalars(p, norms);
      Discretization()->ClearState();

      // print error
      if (dataglobalstate_->GetMyRank() == 0)
      {
        {
          std::cout.precision(8);
          std::cout << std::endl;
          std::cout << "---- Error norm for analytical solution -------------------" << std::endl;
          std::cout << "| absolute L_2 displacement error norm:   " << sqrt((*norms)(0))
                    << std::endl;
          std::cout << "-----------------------------------------------------------" << std::endl;
          std::cout << std::endl;
        }

        // print last error in a seperate file

        // append error of the last time step to the error file
        if ((dataglobalstate_->GetStepN() == datasdyn_->GetStepMax()) or
            (dataglobalstate_->GetTimeN() == datasdyn_->GetTimeMax()))  // write results to file
        {
          std::ostringstream temp;
          const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
          const std::string fname = simulation + ".abserror";

          std::ofstream f;
          f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
          f << "#| " << simulation << "\n";
          f << "#| Step | Time | abs. L2-error displacement |\n";
          f << dataglobalstate_->GetStepN() << " " << dataglobalstate_->GetTimeN() << " "
            << sqrt((*norms)(0)) << "\n";
          f.flush();
          f.close();
        }

        std::ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation + "_time.abserror";

        if (dataglobalstate_->GetStepN() == 1)
        {
          std::ofstream f;
          f.open(fname.c_str());
          f << "#| Step | Time | abs. L2-error displacement |\n";
          f << std::setprecision(10) << dataglobalstate_->GetStepN() << " " << std::setw(1)
            << std::setprecision(5) << dataglobalstate_->GetTimeN() << std::setw(1)
            << std::setprecision(6) << " " << sqrt((*norms)(0)) << "\n";
          f.flush();
          f.close();
        }
        else
        {
          std::ofstream f;
          f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
          f << std::setprecision(10) << dataglobalstate_->GetStepN() << " " << std::setw(3)
            << std::setprecision(5) << dataglobalstate_->GetTimeN() << std::setw(1)
            << std::setprecision(6) << " " << sqrt((*norms)(0)) << "\n";
          f.flush();
          f.close();
        }
      }
    }
    break;
    default:
      dserror("Cannot calculate error. Unknown type of analytical test problem");
      break;
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::PostUpdate() { int_ptr_->PostUpdate(); }

void STR::TIMINT::Base::PostTimeLoop() { int_ptr_->PostTimeLoop(); }

bool STR::TIMINT::Base::HasFinalStateBeenWritten() const
{
  return dataio_->GetLastWrittenResults() == dataglobalstate_->GetStepN();
}

BACI_NAMESPACE_CLOSE
