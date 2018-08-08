/*!----------------------------------------------------------------------
\file adapter_coupling_volmortar.cpp

\brief adapter for the volmortar framework

\level 2

<pre>
\maintainer Alexander Popp
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  farah 10/13|
 *----------------------------------------------------------------------*/
#include "adapter_coupling_volmortar.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_volmortar/volmortar_coupling.H"
#include "../drt_volmortar/volmortar_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_io/io_gmsh.H"
#include"../drt_inpar/inpar_volmortar.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

/*----------------------------------------------------------------------*
 |  ctor                                                     farah 10/13|
 *----------------------------------------------------------------------*/
ADAPTER::MortarVolCoupl::MortarVolCoupl() :
  issetup_(false),
  isinit_(false),
  P12_(Teuchos::null),
  P21_(Teuchos::null),
  masterdis_(Teuchos::null),
  slavedis_(Teuchos::null),
  coupleddof12_(NULL),
  coupleddof21_(NULL),
  dofsets12_(NULL),
  dofsets21_(NULL),
  materialstrategy_(Teuchos::null)
{
  //empty...
}


/*----------------------------------------------------------------------*
 |  init                                                     farah 10/13|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::Init(Teuchos::RCP<DRT::Discretization> dis1,  // masterdis - on Omega_1
                                    Teuchos::RCP<DRT::Discretization> dis2, // slavedis  - on Omega_2
                                    std::vector<int>* coupleddof12,
                                    std::vector<int>* coupleddof21,
                                    std::pair<int,int>* dofsets12,
                                    std::pair<int,int>* dofsets21,
                                    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy,
                                    bool createauxdofs)
{
  // Note : We need to make sure that the parallel distribution of discretizations
  //        is the same externally! The best thing is if you do this in your *_dyn.cpp,
  //        i.e., your global control algorithm.

  // reset the setup flag
  issetup_=false;

  // set pointers to discretizations
  masterdis_=dis1;
  slavedis_=dis2;

  // set various pointers
  coupleddof12_ = coupleddof12;
  coupleddof21_ = coupleddof21;
  dofsets12_    = dofsets12;
  dofsets21_    = dofsets21;
  materialstrategy_ = materialstrategy;

  if ((dis1->NumDofSets() == 1) and (dis2->NumDofSets() == 1) and createauxdofs)
  {
    if(coupleddof12==NULL or coupleddof21==NULL)
      dserror("ERROR: No coupling dofs for volmortar algorithm specified!");

    CreateAuxDofsets(dis1, dis2, coupleddof12, coupleddof21);
  }

  // set flag
  isinit_=true;

  // bye
  return;
}


/*----------------------------------------------------------------------*
 |  setup                                                    rauch 08/16|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::Setup()
{
  CheckInit();

  // get problem dimension (2D or 3D)
  const int dim = DRT::Problem::Instance()->NDim();

  // get volmortar params
  const Teuchos::ParameterList& params =
      DRT::Problem::Instance()->VolmortarParams();

  // create material strategy
  if(materialstrategy_.is_null())
    materialstrategy_= Teuchos::rcp(new VOLMORTAR::UTILS::DefaultMaterialStrategy() );

  // create coupling instance
  Teuchos::RCP<VOLMORTAR::VolMortarCoupl> coupdis =
      Teuchos::rcp(new VOLMORTAR::VolMortarCoupl(
          dim,
          masterdis_,
          slavedis_,
          coupleddof12_,
          coupleddof21_,
          dofsets12_,
          dofsets21_,
          materialstrategy_));

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

  /***********************************************************
   * Assign materials                                        *
   ***********************************************************/
  //assign materials from one discretization to the other
  coupdis->AssignMaterials();

  // validate flag issetup_
  issetup_=true;

  return;
}


/*----------------------------------------------------------------------*
 |  redistribute                                             rauch 08/16|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::Redistribute()
{
  CheckInit();

  // create vector of discr.
  std::vector<Teuchos::RCP<DRT::Discretization> > dis;
  dis.push_back(masterdis_);
  dis.push_back(slavedis_);

  DRT::UTILS::RedistributeDiscretizationsByBinning(dis,false);

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
  // add proxy of velocity related degrees of freedom to scatra discretization
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(coupleddof21->size(), 0, 0, true));
  if (dis2->AddDofSet(dofsetaux) != 1)
    dserror("unexpected dof sets in fluid field");
  dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(coupleddof12->size(), 0, 0, true));
  if (dis1->AddDofSet(dofsetaux) != 1)
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
  CheckInit();

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
  CheckInit();

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
  CheckInit();

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
  CheckInit();

  return LINALG::MLMultiply(*mat,false,*P21_,false,false,false,true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::MortarVolCoupl::MasterToSlave(
    Teuchos::RCP<const Epetra_Vector> mv) const
{
  // safety check
  CheckSetup();
  CheckInit();

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
  if (not mv->Map().PointSameAs(P21_->DomainMap()))
    dserror("master dof map vector expected");
  if (not sv->Map().PointSameAs(P21_->RowMap()))
    dserror("slave dof map vector expected");
  if (sv->NumVectors()!=mv->NumVectors())
    dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

  // safety check
  CheckSetup();
  CheckInit();

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
  CheckInit();

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
  CheckInit();

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
  CheckInit();

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
  if (not sv->Map().PointSameAs(P21_->RowMap()))
    dserror("slave dof map vector expected");
  if (sv->NumVectors()!=mv->NumVectors())
    dserror("column number mismatch %d!=%d",sv->NumVectors(),mv->NumVectors());
#endif

  // safety check
  CheckSetup();
  CheckInit();

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
  CheckInit();

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
