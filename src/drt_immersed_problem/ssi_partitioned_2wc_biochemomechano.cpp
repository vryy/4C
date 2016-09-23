/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_biochemomechano.cpp

\brief specialization of ssi2wc for biochemo-mechano coupled active cell material

\level 3

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#include "ssi_partitioned_2wc_biochemomechano.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_poroelast/poro_scatra_base.H"

#include "../drt_immersed_problem/immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | constructor                                              rauch 01/16 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_BIOCHEMOMECHANO::SSI_Part2WC_BIOCHEMOMECHANO(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams)
  : SSI_Part2WC(comm, globaltimeparams),
    myrank_(comm.MyPID()),
    exchange_manager_(NULL)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call he setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
}


/*----------------------------------------------------------------------*
 | Setup this object                                        rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_BIOCHEMOMECHANO::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
{
  int returnvar=0;

  // call setup of base class
  returnvar =
      SSI::SSI_Part2WC::Init(comm,globaltimeparams,scatraparams,structparams,struct_disname,scatra_disname);

  return returnvar;
}


/*----------------------------------------------------------------------*
 | Setup this specific object                               rauch 08/16 |
 *----------------------------------------------------------------------*/
bool SSI::SSI_Part2WC_BIOCHEMOMECHANO::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
{
  int returnvar=0;

  // call standard setup
  returnvar =
      Init(comm,globaltimeparams,scatraparams,structparams,struct_disname,scatra_disname);

  // get pointer poroelast-scatra interaction subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase> >("RCPToPoroScatra");

  return returnvar;
}



/*----------------------------------------------------------------------*
 | Initialize this class                                    rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::Setup()
{
  specialized_structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(StructureField());
  if(specialized_structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  Teuchos::RCP<Epetra_MultiVector> phinp = scatra_->ScaTraField()->Phinp();
  const Epetra_Map* colmap = scatra_->ScaTraField()->Discretization()->DofColMap(0);
  phi_ = LINALG::CreateVector(*colmap,true);
  LINALG::Export(*phinp,*phi_);

  // Construct vectors for exchange between Structure --> Scatra
  const Epetra_Map* elementcolmap = scatra_->ScaTraField()->Discretization()->ElementColMap();

  ratesactin_ = LINALG::CreateVector(*elementcolmap,true);
  ratesactin_.reset(new Epetra_MultiVector(*elementcolmap,8));

  rates_ = LINALG::CreateVector(*elementcolmap,true);
  rates_.reset(new Epetra_MultiVector(*elementcolmap,8));

  // get pointer to the ImmersedFieldExchangeManager
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();

  // Set pointers to multivectors for the rates
  exchange_manager_->SetPointerToRates(rates_);
  exchange_manager_->SetPointerToRatesActin(ratesactin_);

  // Set pointer to the concentrations
  exchange_manager_->SetPointerToPhinps(phi_);

  return;
}


/*----------------------------------------------------------------------*
| update time step and print to screen                      rauch 01/16 |
------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::UpdateAndOutput()
{
  SSI::SSI_Part2WC::UpdateAndOutput();

  // also for dummy poroscatra problem
  poroscatra_subproblem_->PrepareTimeStep(false);
  poroscatra_subproblem_->PrepareOutput();
  poroscatra_subproblem_->Update();
  poroscatra_subproblem_->Output();
}


/*----------------------------------------------------------------------*
 | convergence check                                        rauch 01/16 |
 *----------------------------------------------------------------------*/
bool SSI::SSI_Part2WC_BIOCHEMOMECHANO::ConvergenceCheck(int itnum)
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //     | scalar+1 |_2
  //
  // AND
  //
  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //             dt * n

  // variables to save different L2 - Norms
  // define L2-norm of incremental scalar and scalar
  double scaincnorm_L2(0.0);
  double scanorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  scaincnp_->Update(1.0,*(scatra_->ScaTraField()->Phinp()),-1.0);
  dispincnp_->Update(1.0,*(structure_->Dispnp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  scatra_->ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n";
    std::cout<<"***********************************************************************************\n";
    std::cout<<"    OUTER BIOCHEMO-MECHANO STRESS FIBER MODEL ITERATION STEP    \n";
    std::cout<<"***********************************************************************************\n";
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]  -|-  scalar-inc  -|-  disp-inc      -|-  scalar-rel-inc  -|-  disp-rel-inc  -|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E      |  %10.3E        |  %10.3E      |",
         itnum,itmax_,ittol_,scaincnorm_L2/Dt()/sqrt(scaincnp_->GlobalLength()),dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()),scaincnorm_L2/scanorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("\n");
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
  }

  // converged
  if (  ( ((scaincnorm_L2/scanorm_L2) <= ittol_) and ((dispincnorm_L2/dispnorm_L2) <= ittol_) ) or ( ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))<=ittol_) and ((scaincnorm_L2/Dt()/sqrt(scaincnp_->GlobalLength()))<=ittol_) ) )
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 )
    {
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                                                      |\n", itnum,itmax_);
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ( (itnum==itmax_) and (  ( ((scaincnorm_L2/scanorm_L2) > ittol_) and ((dispincnorm_L2/dispnorm_L2) > ittol_) ) or ( ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))>ittol_) and ((scaincnorm_L2/Dt()/sqrt(scaincnp_->GlobalLength()))>ittol_) ) ))
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in itemax steps!                                                                      |\n");
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned SSI solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}


/*----------------------------------------------------------------------*
 | Update values of Phinp                                   rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::UpdateScalars()
{

  const Epetra_Map* colmap = scatra_->ScaTraField()->Discretization()->DofColMap(0);
  // Get pointers to Phi
  Teuchos::RCP<Epetra_MultiVector> pointstoPhinp = exchange_manager_->GetPointerToPhinps();
  // Get current value of Phi
  Teuchos::RCP<Epetra_MultiVector> phinp = scatra_->ScaTraField()->Phinp();
  Teuchos::RCP<Epetra_MultiVector> tmp = LINALG::CreateVector(*colmap,true);
  LINALG::Export(*phinp,*tmp);

  int err = pointstoPhinp->Update(1.0, *tmp, 0.0);

  if(err!=0)
    dserror(" Epetra_Vector update threw error code %i ",err);

}


/*----------------------------------------------------------------------*
 | Solve structure filed                                    rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n CELL SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_-> Solve();
}


/*----------------------------------------------------------------------*
 | Solve Scatra field                                       rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n*******************************\n  BIOCHEMO TRANSPORT SOLVER \n*******************************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();
  UpdateScalars();

}


/*----------------------------------------------------------------------*/
//prepare time step                                         rauch 01/16 |
/*----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  if(printheader)
    PrintHeader();

  SetScatraSolution(scatra_->ScaTraField()->Phin());
  structure_-> PrepareTimeStep();

  SetStructSolution(structure_->Dispn(),structure_->Veln());
  scatra_->ScaTraField()->PrepareTimeStep();

  UpdateScalars();
}


