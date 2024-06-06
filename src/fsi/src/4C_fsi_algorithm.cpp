/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all FSI algorithms


\level 1
*/
/*----------------------------------------------------------------------*/


#include "4C_fsi_algorithm.hpp"

#include "4C_adapter_str_factory.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fsi_str_model_evaluator_partitioned.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these entries
// define the order in which the filters handle the Discretizations, which in
// turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*/
FSI::Algorithm::Algorithm(const Epetra_Comm& comm)
    : AlgorithmBase(comm, Global::Problem::Instance()->FSIDynamicParams()),
      adapterbase_ptr_(Teuchos::null),
      use_old_structure_(false)
{
  // empty constructor
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Setup()
{
  // access the structural discretization
  Teuchos::RCP<Discret::Discretization> structdis =
      Global::Problem::Instance()->GetDis("structure");

  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();

  // access the fsi dynamic params
  Teuchos::ParameterList& fsidyn =
      const_cast<Teuchos::ParameterList&>(Global::Problem::Instance()->FSIDynamicParams());

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
    adapterbase_ptr_ = Adapter::build_structure_algorithm(sdyn);
    adapterbase_ptr_->Init(fsidyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
    adapterbase_ptr_->register_model_evaluator("Partitioned Coupling Model", fsi_model_ptr);
    adapterbase_ptr_->Setup();
    structure_ = Teuchos::rcp_dynamic_cast<Adapter::FSIStructureWrapper>(
        adapterbase_ptr_->structure_field());

    // set pointer in FSIStructureWrapper
    structure_->set_model_evaluator_ptr(
        Teuchos::rcp_dynamic_cast<STR::MODELEVALUATOR::PartitionedFSI>(fsi_model_ptr));

    if (structure_ == Teuchos::null)
      FOUR_C_THROW("cast from Adapter::Structure to Adapter::FSIStructureWrapper failed");
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

    Teuchos::RCP<Adapter::StructureBaseAlgorithm> structure = Teuchos::rcp(
        new Adapter::StructureBaseAlgorithm(Global::Problem::Instance()->FSIDynamicParams(),
            const_cast<Teuchos::ParameterList&>(sdyn), structdis));
    structure_ =
        Teuchos::rcp_dynamic_cast<Adapter::FSIStructureWrapper>(structure->structure_field());
    structure_->Setup();

    if (structure_ == Teuchos::null)
      FOUR_C_THROW("cast from Adapter::Structure to Adapter::FSIStructureWrapper failed");

    use_old_structure_ = true;
  }
  else
    FOUR_C_THROW(
        "Unknown time integration requested!\n"
        "Set parameter INT_STRATEGY to Standard in ---STRUCUTRAL DYNAMIC section!\n"
        "If you want to use yet unsupported elements or you want to do crack simulation,\n"
        "set INT_STRATEGY to Old in ---STRUCUTRAL DYNAMIC section!");

  Teuchos::RCP<Adapter::FluidMovingBoundaryBaseAlgorithm> MBFluidbase =
      Teuchos::rcp(new Adapter::FluidMovingBoundaryBaseAlgorithm(
          Global::Problem::Instance()->FSIDynamicParams(), "FSICoupling"));
  fluid_ = MBFluidbase->MBFluidField();

  coupsf_ = Teuchos::rcp(new Core::Adapter::Coupling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::read_restart(int step)
{
  structure_field()->read_restart(step);
  double time = MBFluidField()->read_restart(step);
  SetTimeStep(time, step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::prepare_time_step()
{
  increment_time_and_step();

  print_header();

  structure_field()->prepare_time_step();
  MBFluidField()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::update()
{
  structure_field()->Update();
  MBFluidField()->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::prepare_output(bool force_prepare)
{
  structure_field()->prepare_output(force_prepare);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  structure_field()->Output();
  MBFluidField()->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::struct_to_fluid(Teuchos::RCP<Epetra_Vector> iv)
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::fluid_to_struct(Teuchos::RCP<Epetra_Vector> iv)
{
  return coupsf_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::Adapter::Coupling& FSI::Algorithm::structure_fluid_coupling() { return *coupsf_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Core::Adapter::Coupling& FSI::Algorithm::structure_fluid_coupling() const { return *coupsf_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::struct_to_fluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Algorithm::fluid_to_struct(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}

FOUR_C_NAMESPACE_CLOSE
