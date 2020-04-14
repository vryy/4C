/*---------------------------------------------------------------------*/
/*! \file
\brief global algorithm control class for all immersed algorithms

\level 2

\maintainer Jonas Eichinger
*/
/*---------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>
#include "immersed_problem_dyn.H"
#include "immersed_base.H"
#include "immersed_partitioned_fsi_dirichletneumann.H"

#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_immersed.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_ale/ale_utils_clonestrategy.H"



void immersed_problem_drt()
{
  ///////////////////////////////////////////////////////////////////////
  // General declarations and variables
  ///////////////////////////////////////////////////////////////////////
  // declare general ParameterList that can be handed into Algorithm
  Teuchos::ParameterList params;
  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  ///////////////////////////////////////////////////////////////////////
  // Get Parameters
  ///////////////////////////////////////////////////////////////////////
  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();
  // choose algorithm
  int coupling = DRT::INPUT::IntegralValue<int>(immersedmethodparams, "COUPALGO");
  int scheme = DRT::INPUT::IntegralValue<int>(immersedmethodparams, "SCHEME");


  ///////////////////////////////////////////////////////////////////////
  // Query algorithms
  ///////////////////////////////////////////////////////////////////////
  switch (coupling)
  {
    case INPAR::IMMERSED::partitioned:
    {
      switch (DRT::Problem::Instance()->GetProblemType())
      {
        case prb_immersed_fsi:
        {
          // fill discretizations
          problem->GetDis("structure")->FillComplete(false, false, false);
          problem->GetDis("fluid")->FillComplete(false, false, false);

          // SAFETY FIRST
          {
            // check if INODE is defined in input file
            int gid = problem->GetDis("fluid")->ElementRowMap()->GID(0);
            IMMERSED::ImmersedNode* inode = dynamic_cast<IMMERSED::ImmersedNode*>(
                (problem->GetDis("fluid")->gElement(gid)->Nodes()[0]));

            if (inode == NULL)
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
          DRT::UTILS::GhostDiscretizationOnAllProcs(problem->GetDis("structure"));

          // setup algo
          algo->Setup();

          // PARTITIONED FSI ALGORITHM

          // read restart step
          const int restart = DRT::Problem::Instance()->Restart();
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
          DRT::Problem::Instance()->AddFieldTest(algo->MBFluidField()->CreateFieldTest());
          DRT::Problem::Instance()->AddFieldTest(algo->StructureField()->CreateFieldTest());

          // do the actual testing
          DRT::Problem::Instance()->TestAll(comm);

          break;
        }  // case prb_immersed_fsi

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
