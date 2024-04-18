/*---------------------------------------------------------------------*/
/*! \file
\brief global algorithm control class for all immersed algorithms

\level 2

*/
/*---------------------------------------------------------------------*/
#include "baci_immersed_problem_dyn.hpp"

#include "baci_adapter_str_fsiwrapper_immersed.hpp"
#include "baci_ale_utils_clonestrategy.hpp"
#include "baci_fsi_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_immersed_problem_immersed_base.hpp"
#include "baci_immersed_problem_immersed_partitioned_fsi_dirichletneumann.hpp"
#include "baci_inpar_immersed.hpp"
#include "baci_lib_utils_createdis.hpp"
#include "baci_rebalance_binning_based.hpp"

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
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  ///////////////////////////////////////////////////////////////////////
  // Get Parameters
  ///////////////////////////////////////////////////////////////////////
  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();
  // choose algorithm
  int coupling = CORE::UTILS::IntegralValue<int>(immersedmethodparams, "COUPALGO");
  int scheme = CORE::UTILS::IntegralValue<int>(immersedmethodparams, "SCHEME");


  ///////////////////////////////////////////////////////////////////////
  // Query algorithms
  ///////////////////////////////////////////////////////////////////////
  switch (coupling)
  {
    case INPAR::IMMERSED::partitioned:
    {
      switch (GLOBAL::Problem::Instance()->GetProblemType())
      {
        case GLOBAL::ProblemType::immersed_fsi:
        {
          // fill discretizations
          problem->GetDis("structure")->FillComplete(false, false, false);
          problem->GetDis("fluid")->FillComplete(false, false, false);

          // SAFETY FIRST
          {
            // check if INODE is defined in input file
            int gid = problem->GetDis("fluid")->ElementRowMap()->GID(0);
            DRT::ImmersedNode* inode = dynamic_cast<DRT::ImmersedNode*>(
                (problem->GetDis("fluid")->gElement(gid)->Nodes()[0]));

            if (inode == nullptr)
              dserror(
                  "dynamic cast from Node to ImmersedNode failed.\n"
                  "Make sure you defined INODE instead of NODE in your input file.");
          }

          {
            // check if structural predictor ConstDisVelAcc is chosen in input file
            if (problem->StructuralDynamicParams().get<std::string>("PREDICT") != "ConstDisVelAcc")
              dserror(
                  "Invalid structural predictor for immersed fsi!\n"
                  "Choose ConstDisVelAcc as predictor in ---STRUCTURAL DYNAMIC section.\n"
                  "Structural state projected onto fluid in new time step should be the same as in "
                  "previous time step.");
          }

          Teuchos::RCP<IMMERSED::ImmersedPartitionedFSIDirichletNeumann> algo = Teuchos::null;
          if (scheme == INPAR::IMMERSED::dirichletneumann)
            algo = Teuchos::rcp(new IMMERSED::ImmersedPartitionedFSIDirichletNeumann(comm));
          else
          {
            algo = Teuchos::null;
            dserror("unknown coupling scheme");
          }

          // init algo
          algo->Init(params);

          // ghost structure redundantly on all procs
          CORE::REBALANCE::GhostDiscretizationOnAllProcs(problem->GetDis("structure"));

          // setup algo
          algo->Setup();

          // PARTITIONED FSI ALGORITHM

          // read restart step
          const int restart = GLOBAL::Problem::Instance()->Restart();
          if (restart)
            algo->ReadRestart(restart);
          else
            // additional setup for structural search tree, etc.
            algo->SetupStructuralDiscretization();

          algo->Timeloop(algo);

          if (immersedmethodparams.get<std::string>("TIMESTATS") == "endofsim")
          {
            Teuchos::TimeMonitor::summarize();
            Teuchos::TimeMonitor::zeroOutTimers();
          }

          // create result tests for single fields
          GLOBAL::Problem::Instance()->AddFieldTest(algo->MBFluidField()->CreateFieldTest());
          GLOBAL::Problem::Instance()->AddFieldTest(algo->StructureField()->CreateFieldTest());

          // do the actual testing
          GLOBAL::Problem::Instance()->TestAll(comm);

          break;
        }  // case GLOBAL::ProblemType::immersed_fsi

        default:
        {
          dserror("no valid problem type specified");
          break;
        }  // default
        break;
      }  // switch problemtype
      break;
    }  // case partitioned (default)
    case INPAR::IMMERSED::monolithic:
    {
      dserror(
          "Monolithic solution scheme not implemented for immersed problems, yet.\n "
          "Make sure that the parameter COUPALGO is set to 'partitioned'");
      break;
    }  // case monolithic
  }    // end switch(coupling)

}  // immersed_problem_drt()

FOUR_C_NAMESPACE_CLOSE
