/*---------------------------------------------------------------------*/
/*! \file
\brief global algorithm control class for all immersed algorithms

\level 2

*/
/*---------------------------------------------------------------------*/
#include "4C_immersed_problem_dyn.hpp"

#include "4C_adapter_str_fsiwrapper_immersed.hpp"
#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_discretization_fem_general_utils_createdis.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_immersed_problem_immersed_base.hpp"
#include "4C_immersed_problem_immersed_partitioned_fsi_dirichletneumann.hpp"
#include "4C_inpar_immersed.hpp"
#include "4C_rebalance_binning_based.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN



void immersed_problem_drt()
{
  ///////////////////////////////////////////////////////////////////////
  // General declarations and variables
  ///////////////////////////////////////////////////////////////////////
  // declare general ParameterList that can be handed into Algorithm
  Teuchos::ParameterList params;
  // get pointer to global problem
  Global::Problem* problem = Global::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  ///////////////////////////////////////////////////////////////////////
  // Get Parameters
  ///////////////////////////////////////////////////////////////////////
  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->immersed_method_params();
  // choose algorithm
  int coupling = Core::UTILS::IntegralValue<int>(immersedmethodparams, "COUPALGO");
  int scheme = Core::UTILS::IntegralValue<int>(immersedmethodparams, "SCHEME");


  ///////////////////////////////////////////////////////////////////////
  // Query algorithms
  ///////////////////////////////////////////////////////////////////////
  switch (coupling)
  {
    case Inpar::Immersed::partitioned:
    {
      switch (Global::Problem::Instance()->GetProblemType())
      {
        case Core::ProblemType::immersed_fsi:
        {
          // fill discretizations
          problem->GetDis("structure")->fill_complete(false, false, false);
          problem->GetDis("fluid")->fill_complete(false, false, false);

          // SAFETY FIRST
          {
            // check if INODE is defined in input file
            int gid = problem->GetDis("fluid")->ElementRowMap()->GID(0);
            Core::Nodes::ImmersedNode* inode = dynamic_cast<Core::Nodes::ImmersedNode*>(
                (problem->GetDis("fluid")->gElement(gid)->Nodes()[0]));

            if (inode == nullptr)
              FOUR_C_THROW(
                  "dynamic cast from Node to ImmersedNode failed.\n"
                  "Make sure you defined INODE instead of NODE in your input file.");
          }

          {
            // check if structural predictor ConstDisVelAcc is chosen in input file
            if (problem->structural_dynamic_params().get<std::string>("PREDICT") !=
                "ConstDisVelAcc")
              FOUR_C_THROW(
                  "Invalid structural predictor for immersed fsi!\n"
                  "Choose ConstDisVelAcc as predictor in ---STRUCTURAL DYNAMIC section.\n"
                  "Structural state projected onto fluid in new time step should be the same as in "
                  "previous time step.");
          }

          Teuchos::RCP<Immersed::ImmersedPartitionedFSIDirichletNeumann> algo = Teuchos::null;
          if (scheme == Inpar::Immersed::dirichletneumann)
            algo = Teuchos::rcp(new Immersed::ImmersedPartitionedFSIDirichletNeumann(comm));
          else
          {
            algo = Teuchos::null;
            FOUR_C_THROW("unknown coupling scheme");
          }

          // init algo
          algo->Init(params);

          // ghost structure redundantly on all procs
          Core::Rebalance::GhostDiscretizationOnAllProcs(problem->GetDis("structure"));

          // setup algo
          algo->Setup();

          // PARTITIONED FSI ALGORITHM

          // read restart step
          const int restart = Global::Problem::Instance()->Restart();
          if (restart)
            algo->read_restart(restart);
          else
            // additional setup for structural search tree, etc.
            algo->setup_structural_discretization();

          algo->Timeloop(algo);

          if (immersedmethodparams.get<std::string>("TIMESTATS") == "endofsim")
          {
            Teuchos::TimeMonitor::summarize();
            Teuchos::TimeMonitor::zeroOutTimers();
          }

          // create result tests for single fields
          Global::Problem::Instance()->AddFieldTest(algo->MBFluidField()->CreateFieldTest());
          Global::Problem::Instance()->AddFieldTest(algo->structure_field()->CreateFieldTest());

          // do the actual testing
          Global::Problem::Instance()->TestAll(comm);

          break;
        }  // case Core::ProblemType::immersed_fsi

        default:
        {
          FOUR_C_THROW("no valid problem type specified");
          break;
        }  // default
        break;
      }  // switch problemtype
      break;
    }  // case partitioned (default)
    case Inpar::Immersed::monolithic:
    {
      FOUR_C_THROW(
          "Monolithic solution scheme not implemented for immersed problems, yet.\n "
          "Make sure that the parameter COUPALGO is set to 'partitioned'");
      break;
    }  // case monolithic
  }    // end switch(coupling)

}  // immersed_problem_drt()

FOUR_C_NAMESPACE_CLOSE
