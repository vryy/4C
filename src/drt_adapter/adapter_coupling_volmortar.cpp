/*!----------------------------------------------------------------------
\file adapter_coupling_volmortar.cpp

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  farah 10/13|
 *----------------------------------------------------------------------*/
#include "adapter_coupling_volmortar.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_volmortar/volmortar_coupling.H"
#include "../drt_volmortar/volmortar_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_particle/binning_strategy.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_io/io_gmsh.H"
#include"../drt_inpar/inpar_volmortar.H"

/*----------------------------------------------------------------------*
 |  ctor                                                     farah 10/13|
 *----------------------------------------------------------------------*/
ADAPTER::MortarVolCoupl::MortarVolCoupl()
{
  //empty...
}

/*----------------------------------------------------------------------*
 |  setup                                                    farah 10/13|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::Setup(Teuchos::RCP<DRT::Discretization> dis1, // on Omega_1
                                    Teuchos::RCP<DRT::Discretization> dis2, // on Omega_2
                                    std::vector<int>* coupleddof12,
                                    std::vector<int>* coupleddof21,
                                    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy)
{
  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  const int dim = DRT::Problem::Instance()->NDim();

  const Teuchos::ParameterList& params = DRT::Problem::Instance()->VolmortarParams();

  // create vector of discr.
  std::vector<Teuchos::RCP<DRT::Discretization> > dis;
  dis.push_back(dis1);
  dis.push_back(dis2);

  //binning strategy for parallel redistribution
  Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy = Teuchos::null;

  std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
  std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

  // redistribute discr. with help of binning strategy
  if(dis1->Comm().NumProc()>1)
  {
    /// binning strategy is created and parallel redistribution is performed
    binningstrategy = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis,stdelecolmap,stdnodecolmap));
  }

  if(materialstrategy==Teuchos::null)
    materialstrategy= Teuchos::rcp(new VOLMORTAR::UTILS::DefaultMaterialStrategy() );
  // create coupling instance
  Teuchos::RCP<VOLMORTAR::VolMortarCoupl> coupdis =
      Teuchos::rcp(new VOLMORTAR::VolMortarCoupl(dim,dis1,dis2,coupleddof12,coupleddof21,materialstrategy));

  //-----------------------
  // Evaluate volmortar coupling:
  if(DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::CouplingType>(params,"COUPLINGTYPE") ==
      INPAR::VOLMORTAR::couplingtype_volmortar)
    coupdis->EvaluateVolmortar();
  //-----------------------
  // consistent interpolation (NO VOLMORTAR)
  else if (DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::CouplingType>(params,"COUPLINGTYPE")==
      INPAR::VOLMORTAR::couplingtype_coninter)
    coupdis->EvaluateConsistentInterpolation();
  //-----------------------
  else
    dserror("ERROR: Chosen coupling not implemented!!!");

  // get the P operators
  P12_ = coupdis->GetPMatrix12();
  P21_ = coupdis->GetPMatrix21();

  if(dis1->Comm().NumProc()>1)
  {
    /// revert extended ghosting
    if (not DRT::INPUT::IntegralValue<int>(params, "KEEP_EXTENDEDGHOSTING"))
      binningstrategy->RevertExtendedGhosting(dis,stdelecolmap,stdnodecolmap);
  }

  /***********************************************************
   * Assign materials                                        *
   ***********************************************************/
  //assign materials from one discretization to the other
  coupdis->AssignMaterials();

  return;
}

/*----------------------------------------------------------------------*
 |  AssignMaterials                                          vuong 09/14|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::AssignMaterials(
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2,
    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy)
{
  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  const int dim = DRT::Problem::Instance()->NDim();

  if(materialstrategy==Teuchos::null)
    materialstrategy= Teuchos::rcp(new VOLMORTAR::UTILS::DefaultMaterialStrategy() );
  // create coupling instance
  Teuchos::RCP<VOLMORTAR::VolMortarCoupl> coupdis =
      Teuchos::rcp(new VOLMORTAR::VolMortarCoupl(dim,dis1,dis2,NULL,NULL,materialstrategy));

  //assign materials from one discretization to the other
  coupdis->AssignMaterials();

  return;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_2 --> Omega_1                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::MortarVolCoupl::ApplyVectorMapping12(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  Teuchos::RCP<Epetra_Vector> mapvec = LINALG::CreateVector(P12_->RowMap(),true);
  P12_->Multiply(false,*vec,*mapvec);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_1 --> Omega_2                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::MortarVolCoupl::ApplyVectorMapping21(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  Teuchos::RCP<Epetra_Vector> mapvec = LINALG::CreateVector(P21_->RowMap(),true);
  P21_->Multiply(false,*vec,*mapvec);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_2 --> Omega_1                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::MortarVolCoupl::ApplyMatrixMapping12(
    Teuchos::RCP<const LINALG::SparseMatrix> mat) const
{
  return LINALG::MLMultiply(*mat,false,*P12_,false,false,false,true);
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_1 --> Omega_2                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::MortarVolCoupl::ApplyMatrixMapping21(
    Teuchos::RCP<const LINALG::SparseMatrix> mat) const
{
  return LINALG::MLMultiply(*mat,false,*P21_,false,false,false,true);
}
