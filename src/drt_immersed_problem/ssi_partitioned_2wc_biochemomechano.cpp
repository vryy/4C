/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_biochemomechano.cpp

\brief specialization of ssi2wc for biochemo-mechano coupled active cell material

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "ssi_partitioned_2wc_biochemomechano.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_immersed_problem/immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | constructor                                               rauch 01/16 |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_BIOCHEMOMECHANO::SSI_Part2WC_BIOCHEMOMECHANO(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
  : SSI_Part2WC(comm, globaltimeparams, scatraparams, structparams,struct_disname,scatra_disname)
{

  // set communicator
  myrank_ = comm.MyPID();

  specialized_structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(StructureField());
  if(specialized_structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // Initialize pointers for exchange between ScaTra and Structure
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();

  Teuchos::RCP<Epetra_MultiVector> phinp = scatra_->ScaTraField()->Phinp();
  const Epetra_Map* colmap = scatra_->ScaTraField()->Discretization()->DofColMap(0);
  Teuchos::RCP<Epetra_MultiVector> tmp = LINALG::CreateVector(*colmap,true);
  LINALG::Export(*phinp,*tmp);

  // Set pointer for the first time
  exchange_manager_->SetPointerToPhinps(tmp);


  // Initialize Pointers for exchange between Structure --> ScaTra
  const Epetra_Map* elementcolmap = scatra_->ScaTraField()->Discretization()->ElementColMap();

  Teuchos::RCP<Epetra_MultiVector> ratesactin = LINALG::CreateVector(*elementcolmap,true);
  ratesactin.reset(new Epetra_MultiVector(*elementcolmap,8));

  Teuchos::RCP<Epetra_MultiVector> rates = LINALG::CreateVector(*elementcolmap,true);
  rates.reset(new Epetra_MultiVector(*elementcolmap,8));
  // Set pointer for the first time
  exchange_manager_->SetPointerToRates(rates);
  exchange_manager_->SetPointerToRatesActin(ratesactin);


  return;
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
 | Solve Scatra field                                       rauch 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_BIOCHEMOMECHANO::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
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
