/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.cpp

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "../drt_adapter/adapter_coupling.H"

#include "../drt_fsi/fsi_matrixtransform.H"
//
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
//
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_mapextractor.H"

#include <Teuchos_TimeMonitor.hpp>

#include "scatra_timint_meshtying_strategy_s2i_elch.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch* elchtimint
    ) :
MeshtyingStrategyS2I(elchtimint)
{
  return;
} // SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch


/*--------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions (electrochemistry)   fang 10/14 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying() const
{
  // time measurement: evaluate condition 'S2ICoupling'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'S2ICoupling'");

  switch(s2imortartype_)
  {
  case INPAR::SCATRA::s2i_mortar_none:
  {
    // check matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix = scatratimint_->BlockSystemMatrix();
    if(blocksystemmatrix == Teuchos::null)
      dserror("System matrix is not a block matrix!");

    // create parameter list
    Teuchos::ParameterList condparams;

    // action for elements
    condparams.set<int>("action",SCATRA::bd_calc_s2icoupling);
    condparams.set<int>("scatratype",scatratimint_->ScaTraType());
    condparams.set<bool>("isale",scatratimint_->IsALE());
    condparams.set<double>("frt",ElchTimInt()->FRT()); // factor F/RT

    // set global state vectors according to time-integration scheme
    scatratimint_->Discretization()->ClearState();
    scatratimint_->AddTimeIntegrationSpecificVectors();

    // fill interface state vector imasterphinp_ with transformed master dof values and add to discretization
    imasterphinp_->Update(1.,*(icoup_->MasterToSlave(maps_->ExtractVector(*(scatratimint_->Phiafnp()),2))),0.);
    scatratimint_->Discretization()->SetState("imasterphinp",imasterphinp_);

    // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
    islavematrix_->Zero();
    imastermatrix_->Zero();
    islaveresidual_->PutScalar(0.);
    scatratimint_->Discretization()->EvaluateCondition(condparams,islavematrix_,imastermatrix_,islaveresidual_,Teuchos::null,Teuchos::null,"S2ICouplingSlave");
    scatratimint_->Discretization()->ClearState();

    // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
    islavematrix_->Complete();
    blocksystemmatrix->Matrix(1,1).Add(*islavematrix_,false,1.,1.);

    // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global system matrix
    imastermatrix_->Complete();
    (*islavetomastercoltransform_)(imastermatrix_->RowMap(),imastermatrix_->ColMap(),*imastermatrix_,1.,
        ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(1,2));

    // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global system matrix
    (*islavetomasterrowtransform_)(*islavematrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,1));

    // derive linearizations of master fluxes w.r.t. master dofs and assemble into global system matrix
    (*islavetomasterrowcoltransform_)(*imastermatrix_,-1.,ADAPTER::CouplingSlaveConverter(*icoup_),ADAPTER::CouplingSlaveConverter(*icoup_),blocksystemmatrix->Matrix(2,2),true,true);

    // assemble slave residuals into global residual vector
    maps_->AddVector(islaveresidual_,1,scatratimint_->Residual());

    // transform master residuals and assemble into global residual vector
    maps_->AddVector(icoup_->SlaveToMaster(islaveresidual_),2,scatratimint_->Residual(),-1.);

    break;
  }

  default:
  {
    dserror("Not yet implemented!");
    break;
  }
  }

  return;
} // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying
