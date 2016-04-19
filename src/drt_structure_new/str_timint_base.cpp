/*-----------------------------------------------------------*/
/*!
\file str_timint_base.cpp

\maintainer Michael Hiermeier

\date Aug 12, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_base.H"
#include "str_factory.H"
#include "str_model_evaluator_factory.H"
#include "str_dbc.H"
#include "str_integrator.H"
#include "str_resulttest.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io.H"

#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Base::Base()
    : StructureNew(),
      isinit_(false),
      issetup_(false),
      dataio_(Teuchos::null),
      datasdyn_(Teuchos::null),
      dataglobalstate_(Teuchos::null),
      int_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate
    )
{
  // ---------------------------------------------------------------------------
  // We need to call Setup() after Init()
  // ---------------------------------------------------------------------------
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // initilize the data container ptrs
  // ---------------------------------------------------------------------------
  dataio_ = dataio;
  datasdyn_ = datasdyn;
  dataglobalstate_ = dataglobalstate;

  // ---------------------------------------------------------------------------
  // set isInit flag
  // ---------------------------------------------------------------------------
  isinit_ = true;

  // good bye
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Setup()
{
  CheckInit();

  // ---------------------------------------------------------------------------
  // Create the Dirichlet Boundary Condition handler
  // ---------------------------------------------------------------------------
  dbc_ptr_ = Teuchos::rcp(new STR::Dbc());
  /* FixMe It would be sufficient to use a constant discretization,
   * unfortunately this wasn't considered during the implementation of the
   * discretization routines. Therefore many methods need a slight modification
   * (most times adding a "const" should fix the problem).          hiermeier */
  Teuchos::RCP<DRT::Discretization> discret_ptr =
      DataGlobalState().GetMutableDiscret();
  dbc_ptr_->Init(discret_ptr,DataGlobalState().GetMutableFreactNp(),
      Teuchos::rcp(this,false));
  dbc_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // Create the explicit/implicit integrator
  // ---------------------------------------------------------------------------
  int_ptr_ = STR::BuildIntegrator(DataSDyn());
  int_ptr_->Init(DataSDynPtr(),DataGlobalStatePtr(),DataIOPtr(),dbc_ptr_,
      Teuchos::rcp(this,false));
  int_ptr_->Setup();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Reset()
{
  dserror("Reset of all class variables is not yet implemented for "
      "the modelevaluator!");
  // ModelEvaluator().Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResetStep()
{
  dserror("ResetStep of all \"n+1\" class variables is not yet implemented for "
      "the modelevaluator!");
  // ModelEvaluator().ResetStep();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::Base::NotFinished() const
{
  CheckInitSetup();
  // check the current time
  const double& timenp = dataglobalstate_->GetTimeNp();
  const double& timemax = datasdyn_->GetTimeMax();
  const double& dt = (*dataglobalstate_->GetDeltaTime())[0];
  // check the step counter
  const int& stepnp  = dataglobalstate_->GetStepNp();
  const int& stepmax = datasdyn_->GetStepMax();

  return (timenp <= timemax + 1.0e-8*dt and stepnp <= stepmax);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ReadRestart(const int stepn)
{
  CheckInitSetup();
  Teuchos::RCP<DRT::Discretization> actdis = DataGlobalState().GetMutableDiscret();
  IO::DiscretizationReader reader(actdis, stepn);
  if (stepn != reader.ReadInt("step"))
    dserror("Time step on file not equal to given step");

//  step_ = step;
//  stepn_ = step_ + 1;
//  time_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
//  timen_ = (*time_)[0] + (*dt_)[0];
  // set step and time variables
  dataglobalstate_->GetMutableStepN()  = stepn;
  dataglobalstate_->GetMutableStepNp() = stepn + 1;
  dataglobalstate_->GetMutableMultiTime() =
      Teuchos::rcp(new ::TIMINT::TimIntMStep<double>(0, 0, reader.ReadDouble("time")));
  dataglobalstate_->GetMutableTimeNp() =
      dataglobalstate_->GetTimeN() + (*dataglobalstate_->GetDeltaTime())[0];

  // TODO: restart for model--evaluators...
  dserror("FixMe: Restart model evaluators is still missing!");

//  ModelEvaluator().ReadRestart();

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SetRestart(
    int stepn,
    double timen,
    Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln,
    Teuchos::RCP<Epetra_Vector> accn,
    Teuchos::RCP<std::vector<char> > elementdata,
    Teuchos::RCP<std::vector<char> > nodedata)
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
const Teuchos::RCP<const LINALG::MapExtractor> STR::TIMINT::Base::GetDBCMapExtractor()
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
double STR::TIMINT::Base::TimIntParam() const
{
  CheckInitSetup();
  return int_ptr_->GetIntParam();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::Base::HaveModelType(const INPAR::STR::ModelType& modeltype)
const
{
  CheckInitSetup();
  return (datasdyn_->GetModelTypes().find(modeltype)!=
      datasdyn_->GetModelTypes().end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::Base::HaveEleTech(const INPAR::STR::EleTech& eletech)
const
{
  CheckInitSetup();
  return (datasdyn_->GetElementTechnologies().find(eletech) !=
      datasdyn_->GetElementTechnologies().end());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResizeMStepTimAda()
{
  CheckInitSetup();
  // resize time and stepsize fields
  const double& timen = dataglobalstate_->GetTimeN();
  dataglobalstate_->GetMutableMultiTime()->Resize(-1,0,timen);
  const double& dtn = (*dataglobalstate_->GetMutableDeltaTime())[0];
  dataglobalstate_->GetMutableDeltaTime()->Resize(-1,0,dtn);

  // resize state vectors, AB2 is a 2-step method, thus we need two
  // past steps at t_{n} and t_{n-1}
  const Epetra_Map* dofrowmap_ptr = dataglobalstate_->DofRowMapView();
  dataglobalstate_->GetMutableMultiDis()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->GetMutableMultiVel()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->GetMutableMultiAcc()->Resize(-1, 0, dofrowmap_ptr, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Update()
{
  CheckInitSetup();
  int_ptr_->PreUpdate();
  int_ptr_->UpdateStepState();
  UpdateStepTime();
  int_ptr_->UpdateStepElement();
  int_ptr_->PostUpdate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::UpdateStepTime()
{
  CheckInitSetup();
  double& timenp = dataglobalstate_->GetMutableTimeNp();
  int& stepnp = dataglobalstate_->GetMutableStepNp();
  int& stepn  = dataglobalstate_->GetMutableStepN();

  // --------------------------------------------------------------------------
  // update old time and step variables
  // --------------------------------------------------------------------------
  dataglobalstate_->GetMutableMultiTime()->UpdateSteps(timenp);
  stepn = stepnp;

  // --------------------------------------------------------------------------
  // update the new time and step variables
  // --------------------------------------------------------------------------
  // get current time step size
  const double& dtn = (*dataglobalstate_->GetDeltaTime())[0];
  timenp += dtn;
  stepnp += 1;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> STR::TIMINT::Base::CreateFieldTest()
{
  CheckInitSetup();
  Teuchos::RCP<STR::ResultTest> resulttest = Teuchos::rcp(new STR::ResultTest());
  resulttest->Init(GetDataGlobalState());
  resulttest->Setup();

  return resulttest;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::GetRestartData(
    Teuchos::RCP<int> step,
    Teuchos::RCP<double> time,
    Teuchos::RCP<Epetra_Vector> disnp,
    Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> accnp,
    Teuchos::RCP<std::vector<char> > elementdata,
    Teuchos::RCP<std::vector<char> > nodedata)
{
  CheckInitSetup();
  // at some point we have to create a copy
  *step = dataglobalstate_->GetStepN();
  *time = dataglobalstate_->GetTimeN();
  disnp = Teuchos::rcp(new Epetra_Vector(*dataglobalstate_->GetDisNp()));
  velnp = Teuchos::rcp(new Epetra_Vector(*dataglobalstate_->GetVelNp()));
  accnp = Teuchos::rcp(new Epetra_Vector(*dataglobalstate_->GetAccNp()));
  *elementdata = *(dataglobalstate_->GetDiscret()->PackMyElements());
  *nodedata = *(dataglobalstate_->GetDiscret()->PackMyNodes());

  // get restart data is only for simple structure problems
  // hence if the model set is large than one, we throw an error
  if (datasdyn_->GetModelTypes().size()>1)
    dserror("The GetRestartData routine supports the structural model case ONLY!");

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::PrepareOutput()
{
  CheckInitSetup();
  // --- stress and strain calculation -----------------------------------------
  if (dataio_->GetWriteResultsEveryNStep()
      and ((dataio_->GetStressOutputType() != INPAR::STR::stress_none)
      or   (dataio_->GetCouplingStressOutputType() != INPAR::STR::stress_none)
      or   (dataio_->GetStrainOutputType() != INPAR::STR::strain_none)
      or   (dataio_->GetPlasticStrainOutputType() != INPAR::STR::strain_none))
      and (dataglobalstate_->GetStepN()%dataio_->GetWriteResultsEveryNStep() == 0))
  {
    int_ptr_->DetermineStressStrain();
  }
 // --- energy calculation ------------------------------------------------------
  if (     dataio_->GetWriteEnergyEveryNStep()
      and (dataglobalstate_->GetStepN()%dataio_->GetWriteEnergyEveryNStep() == 0))
  {
    int_ptr_->DetermineEnergy();
  }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Output(bool forced_writerestart)
{
  CheckInitSetup();
  PreOutput();
  OutputStep(forced_writerestart);
  // write Gmsh output
  writeGmshStrucOutputStep();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputStep(bool forced_writerestart)
{
  CheckInitSetup();
  // print iterations instead of steps
  if (dataio_->IsOutputEveryIter())
  {
    dserror("OutputEveryIter is not yet implemented!");
    // OutputEveryIter();
    return;
  }

  // special treatment is necessary when restart is forced
  if(forced_writerestart)
  {
    // reset possible history data on element level
    ResetStep();
    // restart has already been written or simulation has just started
    if((    dataio_->GetWriteRestartEveryNStep()
        and (dataglobalstate_->GetStepN()%dataio_->GetWriteRestartEveryNStep() == 0))
        or   dataglobalstate_->GetStepN() == DRT::Problem::Instance()->Restart())
      return;
    // if state already exists, add restart information
    if(    dataio_->GetWriteResultsEveryNStep()
       and (dataglobalstate_->GetStepN()%dataio_->GetWriteResultsEveryNStep() == 0))
    {
      dserror("AddRestartToOutputState() is not yet implemented!");
//      AddRestartToOutputState();
      return;
    }
  }

  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if (    (dataio_->GetWriteRestartEveryNStep()
      and (dataglobalstate_->GetStepN()%dataio_->GetWriteRestartEveryNStep() == 0)
      and  dataglobalstate_->GetStepN() != 0)
      or   forced_writerestart )
  {
    dserror("OutputRestart() is not yet implemented!");
    NewIOStep(datawritten);
//    OutputRestart();
  }

  // output results (not necessary if restart in same step)
  if ( dataio_->IsWriteState()
       and dataio_->GetWriteResultsEveryNStep()
       and (dataglobalstate_->GetStepN()%dataio_->GetWriteResultsEveryNStep() == 0)
       and (not datawritten) )
  {
    NewIOStep(datawritten);
    OutputState();
  }

  // output stress & strain
  if ( dataio_->GetWriteResultsEveryNStep()
       and ( (dataio_->GetStressOutputType() != INPAR::STR::stress_none)
             or (dataio_->GetCouplingStressOutputType() != INPAR::STR::stress_none)
             or (dataio_->GetStrainOutputType() != INPAR::STR::strain_none)
             or (dataio_->GetPlasticStrainOutputType() != INPAR::STR::strain_none) )
       and (dataglobalstate_->GetStepN()%dataio_->GetWriteResultsEveryNStep() == 0) )
  {
//    dserror("OutputStressStrain() is not yet implemented!");
//    // has to be done by the corresponding model evaluators
//    NewIOStep(datawritten);
//    OutputStressStrain();
  }

  // output energy
  if (     dataio_->GetWriteEnergyEveryNStep()
      and (dataglobalstate_->GetStepN()%dataio_->GetWriteEnergyEveryNStep() == 0) )
  {
    dserror("OutputEnergy() is not yet implemented!");
//    OutputEnergy();
  }

  // ToDo print error norms
//  OutputErrorNorms();

//  OutputVolumeMass();

  // ToDo output of nodal positions in current configuration
//  OutputNodalPositions();

  // ToDo write output on micro-scale (multi-scale analysis)
//  if (HaveMicroMat())
//    dserror("OutputMicro() is not yet implemented!"); // OutputMicro();

  // write patient specific output
  if (     dataio_->GetWriteResultsEveryNStep()
      and (dataglobalstate_->GetStepN()%dataio_->GetWriteResultsEveryNStep() == 0))
  {
    // ToDo OutputPatspec()
    // ToDo OutputCell()
  }

  // what's next?
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::NewIOStep(bool& datawritten)
{
  if (datawritten)
    return;

  // Make new step
  if (not datawritten)
    dataio_->GetMutableOutputPtr()->NewStep(dataglobalstate_->GetStepN(),
        dataglobalstate_->GetTimeN());

  datawritten = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputState()
{
  CheckInitSetup();
  int_ptr_->OutputStepState();
  // owner of elements is just written once because it does not change during simulation (so far)
  IO::DiscretizationWriter& owriter = *(dataio_->GetMutableOutputPtr());
  owriter.WriteElementData(dataio_->IsFirstOutputOfRun());
  owriter.WriteNodeData(dataio_->IsFirstOutputOfRun());
  dataio_->SetFirstOutputOfRun(false);

}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::writeGmshStrucOutputStep()
{
  CheckInitSetup();
  if (!dataio_->IsGmsh())
    return;

  const std::string filename = IO::GMSH::GetFileName("struct",
      dataglobalstate_->GetStepNp(), false, dataglobalstate_->GetMyRank());
  std::ofstream gmshfilecontent(filename.c_str());

  // add 'View' to Gmsh postprocessing file
  gmshfilecontent << "View \" " << "struct displacement \" {" << std::endl;
  // draw vector field 'struct displacement' for every element
  IO::GMSH::VectorFieldDofBasedToGmsh(dataglobalstate_->GetMutableDiscret(),
      Dispn(),gmshfilecontent,0,true);
  gmshfilecontent << "};" << std::endl;
}
