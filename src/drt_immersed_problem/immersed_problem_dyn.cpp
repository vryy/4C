/*!----------------------------------------------------------------------
\file immersed_problem_dyn.cpp

\brief global algorithm control class for all immersed algorithms

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_problem_dyn.H"
#include "immersed_base.H"
#include "immersed_partitioned.H"
#include "immersed_partitioned_fsi.H"
#include "immersed_partitioned_cellmigration.H"
#include "immersed_partitioned_fsi_dirichletneumann.H"
#include "immersed_partitioned_fsi_dirichletneumann_ale.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_immersed.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_poroelast/poroelast_monolithic.H"
#include "../drt_poroelast/poro_scatra_part_2wc.H"
#include "../drt_poroelast/poro_utils_clonestrategy.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include <Teuchos_TimeMonitor.hpp>


void immersed_problem_drt()
{
  // get pointer to global problem
  DRT::Problem* problem = DRT::Problem::Instance();
  // get communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // get parameterlist for immersed method
  const Teuchos::ParameterList& immersedmethodparams = problem->ImmersedMethodParams();
  // choose algorithm
  int coupling = DRT::INPUT::IntegralValue<int>(immersedmethodparams,"COUPALGO");
  int scheme   = DRT::INPUT::IntegralValue<int>(immersedmethodparams,"SCHEME");
  switch (coupling)
  {
  case INPAR::IMMERSED::partitioned:
  {
    switch(DRT::Problem::Instance()->ProblemType())
    {
    case prb_immersed_fsi:
    {
      // fill discretizations
      problem->GetDis("structure")->FillComplete();
      problem->GetDis("fluid")    ->FillComplete();

      Teuchos::RCP<IMMERSED::ImmersedPartitionedFSIDirichletNeumann> algo = Teuchos::null;
      if(scheme == INPAR::IMMERSED::dirichletneumann)
        algo = Teuchos::rcp(new IMMERSED::ImmersedPartitionedFSIDirichletNeumann(comm));
      else
      {
        algo = Teuchos::null;
        dserror("unknown coupling scheme");
      }

      // PARTITIONED FSI ALGORITHM
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        algo->ReadRestart(restart);
      }

      algo->SetupStructuralDiscretization();
      algo->Timeloop(algo);

      if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
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
    }// case prb_immersed_fsi
    case prb_immersed_ale_fsi:
    {
      // fill discretizations
      problem->GetDis("structure")->FillComplete();
      problem->GetDis("fluid")    ->FillComplete();

      Teuchos::RCP<IMMERSED::ImmersedPartitionedFSIDirichletNeumann> algo = Teuchos::null;
      if(scheme == INPAR::IMMERSED::dirichletneumann)
        algo = Teuchos::rcp(new IMMERSED::ImmersedPartitionedFSIDirichletNeumannALE(comm));
      else
      {
        algo = Teuchos::null;
        dserror("unknown coupling scheme");
      }

      // PARTITIONED FSI ALGORITHM
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        algo->ReadRestart(restart);
      }

      algo->SetupStructuralDiscretization();
      algo->Timeloop(algo);

      if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
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
    }// case prb_immersed_fsi
    case prb_immersed_cell:
    {
      // fill discretizations
      Teuchos::RCP<DRT::Discretization> structdis    = problem->GetDis("structure");
      Teuchos::RCP<DRT::Discretization> porofluiddis = problem->GetDis("porofluid");
      Teuchos::RCP<DRT::Discretization> scatradis    = problem->GetDis("scatra");

      structdis   ->FillComplete(true,true,true);
      porofluiddis->FillComplete(true,true,true);

      if (porofluiddis->NumGlobalNodes() == 0)
      {
        // Create the fluid discretization within the porous medium
        DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastImmersedCloneStrategy>(structdis,porofluiddis);

        //set material pointers
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,porofluiddis);
      }
      else
      {
        dserror("Porofluiddis already exists! That's quite impossible!");
      }

      scatradis->FillComplete(true,true,true);

      if (scatradis->NumGlobalNodes()==0)
      {
        // fill scatra discretization by cloning structure discretization
        DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroScatraCloneStrategy>(structdis,scatradis);

        // set implementation type
        for(int i=0; i<scatradis->NumMyColElements(); ++i)
        {
          DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if(element == NULL)
            dserror("Invalid element type!");
          else
            element->SetImplType(DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(problem->PoroScatraControlParams(),"SCATRATYPE"));
        }

        // assign materials. Order is important here!
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(structdis,scatradis);
        POROELAST::UTILS::SetMaterialPointersMatchingGrid(porofluiddis,scatradis);
      }
      else
      dserror("Structure AND ScaTra discretization present. This is not supported.");

      problem->GetDis("cell")->FillComplete(true,true,true);

      Teuchos::RCP<IMMERSED::ImmersedPartitionedCellMigration> algo =
          Teuchos::rcp(new IMMERSED::ImmersedPartitionedCellMigration(comm));

      algo->SetupImmersedDiscretization();
      algo->SetupBackgroundDiscretization();

      algo->Timeloop(algo);

      if(immersedmethodparams.get<std::string>("TIMESTATS")=="endofsim")
      {
        Teuchos::TimeMonitor::summarize();
        Teuchos::TimeMonitor::zeroOutTimers();
      }

      // create result tests for single fields
      DRT::Problem::Instance()->AddFieldTest(algo->CellField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->PoroField()->FluidField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->PoroField()->StructureField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->PoroField()->ScaTraFieldBase()->CreateScaTraFieldTest());

      // do the actual testing
      DRT::Problem::Instance()->TestAll(comm);

      break;
    }// case prb_cell_migration
    default:
    {
      dserror("no valid problem type specified");
      break;
    }// default
    break;
    }//switch problemtype
    break;
  } // case partitioned (default)
  case INPAR::IMMERSED::monolithic:
  {
    dserror("Monolithic solution scheme not implemented for immersed problems, yet.\n "
            "Make sure that the parameter COUPALGO is set to 'partitioned'");
    break;
  }// case monolithic
  }// end switch(coupling)

}// immersed_problem_drt()
