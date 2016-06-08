/*!----------------------------------------------------------------------
\file adapter_coupling_volmortar.cpp

\brief adapter for the volmortar framework

\level 2

<pre>
\maintainer Philipp Farah
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
ADAPTER::MortarVolCoupl::MortarVolCoupl() :
  issetup_(false),
  P12_(Teuchos::null),
  P21_(Teuchos::null)
{
  //empty...
}


/*----------------------------------------------------------------------*
 |  setup                                                    farah 10/13|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::Setup(Teuchos::RCP<DRT::Discretization> dis1, // masterdis - on Omega_1
                                    Teuchos::RCP<DRT::Discretization> dis2, // slavedis  - on Omega_2
                                    std::vector<int>* coupleddof12,
                                    std::vector<int>* coupleddof21,
                                    std::pair<int,int>* dofsets12,
                                    std::pair<int,int>* dofsets21,
                                    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy,
                                    bool redistribute,
                                    bool createauxdofs)
{
  // get problem dimension (2D or 3D)
  const int dim = DRT::Problem::Instance()->NDim();

  const Teuchos::ParameterList& params =
      DRT::Problem::Instance()->VolmortarParams();

  if ((dis1->NumDofSets() == 1) and (dis2->NumDofSets() == 1) and createauxdofs)
  {
    if(coupleddof12==NULL or coupleddof21==NULL)
      dserror("ERROR: No coupling dofs for volmortar algorithm specified!");

    CreateAuxDofsets(dis1, dis2, coupleddof12, coupleddof21);
  }

  // create vector of discr.
  std::vector<Teuchos::RCP<DRT::Discretization> > dis;
  dis.push_back(dis1);
  dis.push_back(dis2);

  //binning strategy for parallel redistribution
  Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy = Teuchos::null;

  std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
  std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

  // redistribute discr. with help of binning strategy
  if(dis1->Comm().NumProc()>1 and redistribute)
  {
    /// binning strategy is created and parallel redistribution is performed
    binningstrategy = Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis,stdelecolmap,stdnodecolmap));
  }

  // create material strategy
  if(materialstrategy==Teuchos::null)
    materialstrategy= Teuchos::rcp(new VOLMORTAR::UTILS::DefaultMaterialStrategy() );

  // create coupling instance
  Teuchos::RCP<VOLMORTAR::VolMortarCoupl> coupdis =
      Teuchos::rcp(new VOLMORTAR::VolMortarCoupl(
          dim,
          dis1,
          dis2,
          coupleddof12,
          coupleddof21,
          dofsets12,
          dofsets21,
          materialstrategy));

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

  if(dis1->Comm().NumProc()>1 and redistribute)
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

  // set flag
  issetup_=true;

  // bye
  return;
}


/*----------------------------------------------------------------------*
 |  Create Auxiliary dofsets for multiphysics                farah 06/15|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::CreateAuxDofsets(
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2,
    std::vector<int>* coupleddof12,
    std::vector<int>* coupleddof21)
{
  //first call FillComplete for single discretizations.
  //This way the physical dofs are numbered successively
  dis1->FillComplete();
  dis2->FillComplete();

  //build auxiliary dofsets, i.e. pseudo dofs on each discretization
  if (dis2->BuildDofSetAuxProxy(coupleddof21->size(), 0, 0, true) != 1)
    dserror("unexpected dof sets in fluid field");
  if (dis1->BuildDofSetAuxProxy(coupleddof12->size(), 0, 0, true ) != 1)
    dserror("unexpected dof sets in structure field");

  //call AssignDegreesOfFreedom also for auxiliary dofsets
  //note: the order of FillComplete() calls determines the gid numbering!
  // 1. dofs 1
  // 2. dofs 2
  // 3. auxiliary dofs 1
  // 4. auxiliary dofs 2
  dis1->FillComplete(true, false,false);
  dis2->FillComplete(true, false,false);

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
      Teuchos::rcp(new VOLMORTAR::VolMortarCoupl(dim,dis1,dis2,NULL,NULL,NULL,NULL,materialstrategy));

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
  // safety check
  CheckSetup();

  Teuchos::RCP<Epetra_Vector> mapvec = LINALG::CreateVector(P12_->RowMap(),true);
  int err = P12_->Multiply(false,*vec,*mapvec);
  if(err!=0)
    dserror("ERROR: Matrix multiply returned error code %i", err);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_1 --> Omega_2                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::MortarVolCoupl::ApplyVectorMapping21(
    Teuchos::RCP<const Epetra_Vector> vec) const
{
  // safety check
  CheckSetup();

  Teuchos::RCP<Epetra_Vector> mapvec = LINALG::CreateVector(P21_->RowMap(),true);
  int err = P21_->Multiply(false,*vec,*mapvec);
  if(err!=0)
    dserror("ERROR: Matrix multiply returned error code %i", err);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_2 --> Omega_1                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::MortarVolCoupl::ApplyMatrixMapping12(
    Teuchos::RCP<const LINALG::SparseMatrix> mat) const
{
  // safety check
  CheckSetup();

  return LINALG::MLMultiply(*mat,false,*P12_,false,false,false,true);
}

/*----------------------------------------------------------------------*
 |  ApplyMapping from Omega_1 --> Omega_2                    farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::MortarVolCoupl::ApplyMatrixMapping21(
    Teuchos::RCP<const LINALG::SparseMatrix> mat) const
{
  // safety check
  CheckSetup();

  return LINALG::MLMultiply(*mat,false,*P21_,false,false,false,true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_Vector> mv) const
{
  // safety check
  CheckSetup();

  //create vector
  Teuchos::RCP<Epetra_Vector> sv = LINALG::CreateVector(P21_->RowMap(),true);
  //project
  MasterToSlave(mv,sv);

  return sv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_MultiVector> mv,
    Teuchos::RCP<Epetra_MultiVector> sv) const
{
#ifdef DEBUG
  if (not mv->Map().PointSameAs(P21_->RowMap()))
    dserror("master dof map vector expected");
  if (not sv->Map().PointSameAs(P21_->RowMap()))
    dserror("slave dof map vector expected");
  if (sv->NumVectors()!=mv->NumVectors())
    dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

  // safety check
  CheckSetup();

  //slave vector with auxiliary dofmap
  Epetra_MultiVector sv_aux(P21_->RowMap(),sv->NumVectors());

  //project
  int err = P21_->Multiply(false,*mv,sv_aux);
  if(err!=0)
    dserror("ERROR: Matrix multiply returned error code %i", err);

  //copy from auxiliary to physical map (needed for coupling in fluid ale algorithm)
  std::copy(sv_aux.Values(), sv_aux.Values()+(sv_aux.MyLength()*sv_aux.NumVectors()), sv->Values());

  // in contrast to the ADAPTER::Coupling class we do not need to export here, as
  // the binning has (or should have) guaranteed the same distribution of master and slave dis
  // on all procs
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ADAPTER::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_MultiVector> mv) const
{
  // safety check
  CheckSetup();

  //create vector
  Teuchos::RCP<Epetra_MultiVector> sv =
    Teuchos::rcp(new Epetra_MultiVector(P21_->RowMap(),mv->NumVectors()));
  //project
  MasterToSlave(mv,sv);

  return sv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::MortarVolCoupl::SlaveToMaster(
    Teuchos::RCP<const Epetra_Vector> sv) const
{
  // safety check
  CheckSetup();

  //create vector
  Teuchos::RCP<Epetra_Vector> mv = LINALG::CreateVector(P12_->RowMap(),true);
  //project
  SlaveToMaster(sv,mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ADAPTER::MortarVolCoupl::SlaveToMaster(
    Teuchos::RCP<const Epetra_MultiVector> sv) const
{
  // safety check
  CheckSetup();

  //create vector
  Teuchos::RCP<Epetra_MultiVector> mv =
    Teuchos::rcp(new Epetra_MultiVector(P12_->RowMap(),sv->NumVectors()));
  //project
  SlaveToMaster(sv,mv);

  return mv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::SlaveToMaster(
    Teuchos::RCP<const Epetra_MultiVector> sv,
    Teuchos::RCP<Epetra_MultiVector> mv) const
{
#ifdef DEBUG
  if (not mv->Map().PointSameAs(P12_->RowMap()))
    dserror("master dof map vector expected");
  if (not sv->Map().PointSameAs(P12_->RowMap()))
    dserror("slave dof map vector expected");
  if (sv->NumVectors()!=mv->NumVectors())
    dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

  // safety check
  CheckSetup();

  //master vector with auxiliary dofmap
  Epetra_MultiVector mv_aux(P12_->RowMap(),mv->NumVectors());

  //project
  int err = P12_->Multiply(false,*sv,mv_aux);
  if(err!=0)
    dserror("ERROR: Matrix multiply returned error code %i", err);

  //copy from auxiliary to physical map (needed for coupling in fluid ale algorithm)
  std::copy(mv_aux.Values(), mv_aux.Values()+(mv_aux.MyLength()*mv_aux.NumVectors()), mv->Values());

  // in contrast to the ADAPTER::Coupling class we do not need to export here, as
  // the binning has (or should have) guaranteed the same distribution of master and slave dis
  // on all procs
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>  ADAPTER::MortarVolCoupl::MasterDofMap() const
{
  // safety check
  CheckSetup();

  return Teuchos::rcpFromRef(P12_->RowMap());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::MortarVolCoupl::SlaveDofMap() const
{
  // safety check
  CheckSetup();

  return Teuchos::rcpFromRef(P21_->RowMap());
}
