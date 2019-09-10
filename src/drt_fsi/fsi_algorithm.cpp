/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all FSI algorithms

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------*/


#include "fsi_algorithm.H"
#include "fsi_str_model_evaluator_partitioned.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these entries
// define the order in which the filters handle the Discretizations, which in
// turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*/
FSI::Algorithm::Algorithm(const Epetra_Comm& comm)
    : AlgorithmBase(comm, DRT::Problem::Instance()->FSIDynamicParams()),
      adapterbase_ptr_(Teuchos::null),
      use_old_structure_(false)
{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Algorithm::~Algorithm() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Setup()
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // access the fsi dynamic params
  Teuchos::ParameterList& fsidyn =
      const_cast<Teuchos::ParameterList&>(DRT::Problem::Instance()->FSIDynamicParams());

  // build and register fsi model evaluator
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> fsi_model_ptr =
      Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedFSI());

  // todo FIX THIS !!!!
  // Decide whether to use old structural time integration or new structural time integration.
  // This should be removed as soon as possible! We need to clean up crack fsi first!
  // Also all structural elements need to be adapted first!
  // Then, we can switch the 3 remaining fsi tests using the old time integration to the new one,
  // i.e.: crfsi_inclDomain_suppFluid_interfaceBuild
  //       fsi_dc3D_part_ait_ga_ost_xwall
  //       fsi_ow3D_mtr_drt
  // build structure
  if (sdyn.get<std::string>("INT_STRATEGY") == "Standard")
  {
    adapterbase_ptr_ = ADAPTER::STR::BuildStructureAlgorithm(sdyn);
    adapterbase_ptr_->Init(fsidyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
    adapterbase_ptr_->RegisterModelEvaluator("Partitioned Coupling Model", fsi_model_ptr);
    adapterbase_ptr_->Setup();
    structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::FSIStructureWrapper>(
        adapterbase_ptr_->StructureField());

    // set pointer in FSIStructureWrapper
    structure_->SetModelEvaluatorPtr(
        Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedFSI>(fsi_model_ptr));

    if (structure_ == Teuchos::null)
      dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");
  }
  else if (sdyn.get<std::string>("INT_STRATEGY") ==
           "Old")  // todo this is the part that should be removed !
  {
    if (Comm().MyPID() == 0)
      std::cout << "\n"
                << " USING OLD STRUCTURAL TIME INEGRATION! FIX THIS! THIS IS ONLY SUPPOSED TO BE "
                   "TEMPORARY!"
                   "\n"
                << std::endl;

    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(
        new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams(),
            const_cast<Teuchos::ParameterList&>(sdyn), structdis));
    structure_ =
        Teuchos::rcp_dynamic_cast<::ADAPTER::FSIStructureWrapper>(structure->StructureField());
    structure_->Setup();

    if (structure_ == Teuchos::null)
      dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

    use_old_structure_ = true;
  }
  else
    dserror(
        "Unknown time integration requested!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCUTRAL DYNAMIC section!\n"
        "If you want to use yet unsupported elements or you want to do crack simulation,\n"
        "set INT_STRATEGY to Old in ---STRUCUTRAL DYNAMIC section!");

  Teuchos::RCP<::ADAPTER::FluidMovingBoundaryBaseAlgorithm> MBFluidbase =
      Teuchos::rcp(new ADAPTER::FluidMovingBoundaryBaseAlgorithm(
          DRT::Problem::Instance()->FSIDynamicParams(), "FSICoupling"));
  fluid_ = MBFluidbase->MBFluidField();

  coupsf_ = Teuchos::rcp(new ADAPTER::Coupling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  double time = MBFluidField()->ReadRestart(step);
  SetTimeStep(time, step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  StructureField()->PrepareTimeStep();
  MBFluidField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Update()
{
  StructureField()->Update();
  MBFluidField()->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::PrepareOutput() { StructureField()->PrepareOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField()->Output();
  MBFluidField()->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::Coupling& FSI::Algorithm::StructureFluidCoupling() { return *coupsf_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const ADAPTER::Coupling& FSI::Algorithm::StructureFluidCoupling() const { return *coupsf_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::StructToFluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::FluidToStruct(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}
