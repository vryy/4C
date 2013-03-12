/*----------------------------------------------------------------------*/
/*!
\file adapter_coupling_mortar.cpp

\brief A class providing coupling capabilities based on mortar methods

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#include "adapter_coupling_mortar.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/inpar_fluid.H"

#define ALLDOF

ADAPTER::CouplingMortar::CouplingMortar()
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Setup(DRT::Discretization& masterdis,
    DRT::Discretization& slavedis, DRT::Discretization& aledis,
    const Epetra_Comm& comm, bool structslave)
{
  // initialize maps for row nodes
  std::map<int, DRT::Node*> masternodes;
  std::map<int, DRT::Node*> slavenodes;

  // initialize maps for column nodes
  std::map<int, DRT::Node*> mastergnodes;
  std::map<int, DRT::Node*> slavegnodes;

  //initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element> > masterelements;
  std::map<int, Teuchos::RCP<DRT::Element> > slaveelements;

  // Fill maps based on condition for master side
  DRT::UTILS::FindConditionObjects(masterdis, masternodes, mastergnodes, masterelements,
      "FSICoupling");

  // Fill maps based on condition for slave side
  DRT::UTILS::FindConditionObjects(slavedis, slavenodes, slavegnodes, slaveelements,
      "FSICoupling");

  // get mortar coupling parameters
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MortarCouplingParams();

  // check for invalid parameter values
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_dual)
    dserror("Mortar coupling adapter only works for dual shape functions");
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(input,"INTTYPE") != INPAR::MORTAR::inttype_segments)
    dserror("Mortar coupling adapter only works for segment-based integration");

  // check for parallel redistribution (only if more than 1 proc)
  bool parredist = false;
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    if (comm.NumProc()>1) parredist = true;

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  // IMPORTANT: We assume that all nodes have 'dim' DoF, that have to be considered for coupling.
  //            Possible pressure DoF are not transferred to MortarInterface.
  const int dim = DRT::Problem::Instance()->NDim();

  // create an empty mortar interface
  // (To be on the safe side we still store all interface nodes and elements
  // fully redundant here in the mortar ADAPTER. This makes applications such
  // as SlidingALE much easier, whereas it would not be needed for others.)
  // FIXME: We should use the input flag (popp, 11/2012)
  INPAR::MORTAR::RedundantStorage redundant = INPAR::MORTAR::redundant_all;
  if (redundant != INPAR::MORTAR::redundant_all)
    dserror("Mortar coupling adapter only works for redundant slave and master storage");
  Teuchos::RCP<MORTAR::MortarInterface> interface = Teuchos::rcp(new MORTAR::MortarInterface(0, comm, dim, input, redundant));

  // feeding master nodes to the interface including ghosted nodes
  // only consider the first 'dim' dofs
  std::map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    std::vector<int> dofids(dim);
    for (int k=0;k<dim;++k) dofids[k] = masterdis.Dof(node)[k];
    Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                    dim, dofids, false));

    interface->AddMortarNode(mrtrnode);
  }

  // feeding slave nodes to the interface including ghosted nodes
  for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    std::vector<int> dofids(dim);
    for (int k=0;k<dim;++k) dofids[k] = slavedis.Dof(node)[k];
    Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                    dim, dofids, true));

    interface->AddMortarNode(mrtrnode);
  }

  // max master element ID needed for unique eleIDs in interface discretization
  // will be used as offset for slave elements
  int EleOffset = masterdis.ElementRowMap()->MaxAllGID()+1;

  // feeding master elements to the interface
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator elemiter;
  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                new MORTAR::MortarElement(ele->Id(), ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), false));

    interface->AddMortarElement(mrtrele);
  }

  // feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                new MORTAR::MortarElement(ele->Id() + EleOffset, ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), true));

    interface->AddMortarElement(mrtrele);
  }

  // finalize the contact interface construction
  interface->FillComplete();

  // store old row maps (before parallel redistribution)
  slavedofrowmap_  = Teuchos::rcp(new Epetra_Map(*interface->SlaveRowDofs()));
  masterdofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->MasterRowDofs()));

  // print parallel distribution
  interface->PrintParallelDistribution(1);

  //**********************************************************************
  // PARALLEL REDISTRIBUTION OF INTERFACE
  //**********************************************************************
  if (parredist && comm.NumProc()>1)
  {
    // redistribute optimally among all procs
    interface->Redistribute();

    // call fill complete again
    interface->FillComplete();

    // print parallel distribution again
    interface->PrintParallelDistribution(1);
  }
  //**********************************************************************

  // create binary search tree
  interface->CreateSearchTree();

  // all the following stuff has to be done once in setup
  // in order to get initial D_ and M_

  // interface displacement (=0) has to be merged from slave and master discretization
  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(masterdofrowmap_,slavedofrowmap_, false);
  Teuchos::RCP<Epetra_Vector> dispn = LINALG::CreateVector(*dofrowmap, true);

  // set displacement state in mortar interface
  interface->SetState("displacement", dispn);

  // print message
  if(comm.MyPID()==0)
  {
    std::cout << "\nPerforming mortar coupling...............";
    fflush(stdout);
  }

  //in the following two steps MORTAR does all the work
  interface->Initialize();
  interface->Evaluate();

  // print message
  if(comm.MyPID()==0) std::cout << "done!" << std::endl;

  // preparation for AssembleDM
  // (Note that redistslave and redistmaster are the slave and master row maps
  // after parallel redistribution. If no redistribution was performed, they
  // are of course identical to slavedofrowmap_/masterdofrowmap_!)
  Teuchos::RCP<Epetra_Map> redistslave  = interface->SlaveRowDofs();
  Teuchos::RCP<Epetra_Map> redistmaster = interface->MasterRowDofs();
  Teuchos::RCP<LINALG::SparseMatrix> dmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 10));
  Teuchos::RCP<LINALG::SparseMatrix> mmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 100));
  interface->AssembleDM(*dmatrix, *mmatrix);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*redistmaster, *redistslave);
  D_ = dmatrix;
  M_ = mmatrix;

  // Build Dinv
  Dinv_ = Teuchos::rcp(new LINALG::SparseMatrix(*D_));

  // extract diagonal of invd into diag
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*redistslave,true);
  Dinv_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  diag->Reciprocal(*diag);
  Dinv_->ReplaceDiagonalValues(*diag);
  Dinv_->Complete( D_->RangeMap(), D_->DomainMap() );
  DinvM_ = MLMultiply(*Dinv_,*M_,false,false,true);

  // store interface
  interface_ = interface;

  // mesh initialization (for rotational invariance)
  if (!structslave)
    MeshInit(masterdis,slavedis,aledis,redistmaster,redistslave,comm,structslave);

  // only for parallel redistribution case
  if (parredist)
  {
    // transform everything back to old distribution
    D_     = MORTAR::MatrixRowColTransform(D_,slavedofrowmap_,slavedofrowmap_);
    M_     = MORTAR::MatrixRowColTransform(M_,slavedofrowmap_,masterdofrowmap_);
    Dinv_  = MORTAR::MatrixRowColTransform(Dinv_,slavedofrowmap_,slavedofrowmap_);
    DinvM_ = MORTAR::MatrixRowColTransform(DinvM_,slavedofrowmap_,masterdofrowmap_);
  }

  // check for overlap of slave and Dirichlet boundaries
  // (this is not allowed in order to avoid over-constraint)
  bool overlap = false;
  Teuchos::ParameterList p;
  p.set("total time", 0.0);
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps = Teuchos::rcp(new LINALG::MapExtractor());
  Teuchos::RCP<Epetra_Vector > temp = LINALG::CreateVector(*(slavedis.DofRowMap()), true);
  slavedis.EvaluateDirichlet(p,temp,Teuchos::null,Teuchos::null,Teuchos::null,dbcmaps);

  // loop over all slave row nodes of the interface
  for (int j=0;j<interface_->SlaveRowNodes()->NumMyElements();++j)
  {
    int gid = interface_->SlaveRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // check if this node's dofs are in dbcmap
    for (int k=0;k<mtnode->NumDof();++k)
    {
      int currdof = mtnode->Dofs()[k];
      int lid = (dbcmaps->CondMap())->LID(currdof);

      // found slave node with dbc
      if (lid>=0)
      {
        overlap = true;
        break;
      }
    }
  }

  // print warning message to screen
  if (overlap && comm.MyPID()==0)
  {
    std::cout << RED << "\nWARNING: Slave boundary and Dirichlet boundary conditions overlap!" << std::endl;
    std::cout << "This leads to over-constraint, so you might encounter some problems!" << END_COLOR << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Setup
(
    DRT::Discretization& dis
)
{
  // initialize maps for row nodes
  std::map<int, DRT::Node*> masternodes;
  std::map<int, DRT::Node*> slavenodes;

  // initialize maps for column nodes
  std::map<int, DRT::Node*> mastergnodes;
  std::map<int, DRT::Node*> slavegnodes;

  //initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element> > masterelements;
  std::map<int, Teuchos::RCP<DRT::Element> > slaveelements;

  // Fill maps based on condition for master side
  DRT::UTILS::FindConditionObjects(dis, masternodes, mastergnodes, masterelements,
      "FSICoupling");

  // Fill maps based on condition for slave side
  DRT::UTILS::FindConditionObjects(dis, slavenodes, slavegnodes, slaveelements,
      "FSICoupling");

  // get mortar coupling parameters
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MortarCouplingParams();

  // check for invalid parameter values
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_dual)
    dserror("Mortar coupling adapter only works for dual shape functions");
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(input,"INTTYPE") != INPAR::MORTAR::inttype_segments)
    dserror("Mortar coupling adapter only works for segment-based integration");

  // check for parallel redistribution (only if more than 1 proc)
  bool parredist = false;
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    if (dis.Comm().NumProc()>1) parredist = true;

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  // IMPORTANT: We assume that all nodes have 'dim' DoF, that have to be considered for coupling.
  //            Possible pressure DoF are not transferred to MortarInterface.
  const int dim = DRT::Problem::Instance()->NDim();

  // create an empty mortar interface
  // (To be on the safe side we still store all interface nodes and elements
  // fully redundant here in the mortar ADAPTER. This makes applications such
  // as SlidingALE much easier, whereas it would not be needed for others.)
  // FIXME: We should use the input flag (popp, 11/2012)
  INPAR::MORTAR::RedundantStorage redundant = INPAR::MORTAR::redundant_all;
  if (redundant != INPAR::MORTAR::redundant_all)
    dserror("Mortar coupling adapter only works for redundant slave and master storage");
  Teuchos::RCP<MORTAR::MortarInterface> interface = Teuchos::rcp(new MORTAR::MortarInterface(0, dis.Comm(), dim, input, redundant));

  int NodeOffset = dis.NodeRowMap()->MaxAllGID()+1;
  int DofOffset = dis.DofRowMap()->MaxAllGID()+1;

  // feeding master nodes to the interface including ghosted nodes
  // only consider the first 'dim' dofs
  std::map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    std::vector<int> dofids(dim);
    for (int k=0;k<dim;++k) dofids[k] = dis.Dof(node)[k];
    Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                    dim, dofids, false));

    interface->AddMortarNode(mrtrnode);
  }

  // feeding slave nodes to the interface including ghosted nodes
  for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    std::vector<int> dofids(dim);
    for (int k=0;k<dim;++k) dofids[k] = dis.Dof(node)[k]+DofOffset;
    Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                new MORTAR::MortarNode(node->Id()+NodeOffset, node->X(), node->Owner(),
                    dim, dofids, true));

    interface->AddMortarNode(mrtrnode);
  }

  // max master element ID needed for unique eleIDs in interface discretization
  // will be used as offset for slave elements
  int EleOffset = dis.ElementRowMap()->MaxAllGID()+1;

  // feeding master elements to the interface
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator elemiter;
  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                new MORTAR::MortarElement(ele->Id(), ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), false));

    interface->AddMortarElement(mrtrele);
  }

  // feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    std::vector<int> nidsoff;
    for(int i=0; i<ele->NumNode(); i++)
    {
      nidsoff.push_back(ele->NodeIds()[ele->NumNode()-1-i]+NodeOffset);
    }

    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                new MORTAR::MortarElement(ele->Id() + EleOffset, ele->Owner(), ele->Shape(),
                    ele->NumNode(), &(nidsoff[0]), true));

    interface->AddMortarElement(mrtrele);
  }

  // finalize the contact interface construction
  interface->FillComplete();
  // store old row maps (before parallel redistribution)
  slavedofrowmap_  = Teuchos::rcp(new Epetra_Map(*interface->SlaveRowDofs()));
  masterdofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->MasterRowDofs()));

  // print parallel distribution
  interface->PrintParallelDistribution(1);

//  //**********************************************************************
//  // PARALLEL REDISTRIBUTION OF INTERFACE
//  //**********************************************************************
//  if (parredist && dis.Comm().NumProc()>1)
//  {
//    // redistribute optimally among all procs
//    interface->Redistribute();
//
//    // call fill complete again
//    interface->FillComplete();
//
//    // print parallel distribution again
//    interface->PrintParallelDistribution(1);
//  }
//  //**********************************************************************

  // create binary search tree
  interface->CreateSearchTree();

  // all the following stuff has to be done once in setup
  // in order to get initial D_ and M_

  // interface displacement (=0) has to be merged from slave and master discretization
  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(masterdofrowmap_,slavedofrowmap_, false);
  Teuchos::RCP<Epetra_Vector> dispn = LINALG::CreateVector(*dofrowmap, true);

  // set displacement state in mortar interface
  interface->SetState("displacement", dispn);

  // print message
  if(dis.Comm().MyPID()==0)
  {
    std::cout << "\nPerforming mortar coupling...............";
    fflush(stdout);
  }

  //in the following two steps MORTAR does all the work
  interface->Initialize();
  interface->Evaluate();

  // print message
  if(dis.Comm().MyPID()==0) std::cout << "done!" << std::endl;

  // preparation for AssembleDM
  // (Note that redistslave and redistmaster are the slave and master row maps
  // after parallel redistribution. If no redistribution was performed, they
  // are of course identical to slavedofrowmap_/masterdofrowmap_!)
  Teuchos::RCP<Epetra_Map> redistslave  = interface->SlaveRowDofs();
  Teuchos::RCP<Epetra_Map> redistmaster = interface->MasterRowDofs();
  Teuchos::RCP<LINALG::SparseMatrix> dmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 10));
  Teuchos::RCP<LINALG::SparseMatrix> mmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 100));
  interface->AssembleDM(*dmatrix, *mmatrix);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*redistmaster, *redistslave);
  D_ = dmatrix;
  M_ = mmatrix;

  // Build Dinv
  Dinv_ = Teuchos::rcp(new LINALG::SparseMatrix(*D_));

  // extract diagonal of invd into diag
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*redistslave,true);
  Dinv_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  diag->Reciprocal(*diag);
  Dinv_->ReplaceDiagonalValues(*diag);
  Dinv_->Complete( D_->RangeMap(), D_->DomainMap() );
  DinvM_ = MLMultiply(*Dinv_,*M_,false,false,true);

  // store interface
  interface_ = interface;

  // only for parallel redistribution case
  if (parredist)
  {
    // transform everything back to old distribution
    D_     = MORTAR::MatrixRowColTransform(D_,slavedofrowmap_,slavedofrowmap_);
    M_     = MORTAR::MatrixRowColTransform(M_,slavedofrowmap_,masterdofrowmap_);
    Dinv_  = MORTAR::MatrixRowColTransform(Dinv_,slavedofrowmap_,slavedofrowmap_);
    DinvM_ = MORTAR::MatrixRowColTransform(DinvM_,slavedofrowmap_,masterdofrowmap_);
  }



  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool ADAPTER::CouplingMortar::Setup(DRT::Discretization& dis,
    const Epetra_Comm& comm, int meshtyingoption, bool structslave)
{
  bool pcoupled = true;
  int myrank  = dis.Comm().MyPID();

  // initialize maps for row nodes
  std::map<int, DRT::Node*> masternodes;
  std::map<int, DRT::Node*> slavenodes;

  // initialize maps for column nodes
  std::map<int, DRT::Node*> mastergnodes;
  std::map<int, DRT::Node*> slavegnodes;

  //initialize maps for elements
  std::map<int, Teuchos::RCP<DRT::Element> > masterelements;
  std::map<int, Teuchos::RCP<DRT::Element> > slaveelements;

  std::vector<DRT::Condition*> conds;
  std::vector<DRT::Condition*> conds_master(0);
  std::vector<DRT::Condition*> conds_slave(0);
  const string& condname = "Mortar";
  dis.GetCondition(condname, conds);

  for (unsigned i=0; i<conds.size(); i++)
  {
    const string* side = conds[i]->Get<string>("Side");

    if (*side == "Master")
      conds_master.push_back(conds[i]);

    if(*side == "Slave")
      conds_slave.push_back(conds[i]);
  }

  // Fill maps based on condition for master side
  DRT::UTILS::FindConditionObjects(dis, masternodes, mastergnodes, masterelements,conds_master);

  // Fill maps based on condition for slave side
  DRT::UTILS::FindConditionObjects(dis, slavenodes, slavegnodes, slaveelements,conds_slave);

  // get mortar coupling parameters
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MortarCouplingParams();

  // check for invalid parameter values
  if(meshtyingoption != INPAR::FLUID::sps_coupled and meshtyingoption != INPAR::FLUID::sps_pc)
    if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_dual)
      dserror("Condensation works only for dual shape functions");
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(input,"INTTYPE") != INPAR::MORTAR::inttype_segments)
    dserror("Mortar coupling adapter only works for segment-based integration");

  if (myrank == 0)
  {
    if (DRT::INPUT::IntegralValue<int>(input,"SHAPEFCN")==INPAR::MORTAR::shape_standard)
    {
      if(comm.MyPID()==0)
        std::cout << "Shape functions:  standard" << std::endl;
    }
    else if (DRT::INPUT::IntegralValue<int>(input,"SHAPEFCN")==INPAR::MORTAR::shape_dual)
    {
      if(comm.MyPID()==0)
        std::cout << "Shape functions:  dual" << std::endl;
    }
    else
    {
      dserror("ERROR: Only standard or dual shape functions in this mortar adapter");
    }
  }

  // check for parallel redistribution (only if more than 1 proc)
  //bool parredist = false;
  //if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
  //  if (comm.NumProc()>1) parredist = true;

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  // IMPORTANT: We assume that all nodes have 'dim' DoF, that have to be considered for coupling.

  const int dim = DRT::Problem::Instance()->NDim();

  // create an empty mortar interface
  // (To be on the safe side we still store all interface nodes and elements
  // fully redundant here in the mortar ADAPTER. This makes applications such
  // as SlidingALE much easier, whereas it would not be needed for others.)
  // FIXME: We should use the input flag (popp, 11/2012)
  INPAR::MORTAR::RedundantStorage redundant = INPAR::MORTAR::redundant_all;
  if (redundant != INPAR::MORTAR::redundant_all)
    dserror("Mortar coupling adapter only works for redundant slave and master storage");
  Teuchos::RCP<MORTAR::MortarInterface> interface = Teuchos::rcp(new MORTAR::MortarInterface(0, comm, dim, input, redundant));

  //  Pressure DoF are also transferred to MortarInterface
#ifdef ALLDOF
  int dof = dis.NumDof(dis.lRowNode(0));
  pcoupled = true;

  if(myrank==0)
  {
    std::cout << "All dof's are coupled!! " << std::endl << std::endl;
#else
  //  Pressure DoF are not transferred to MortarInterface
  int dof = dis.NumDof(dis.lRowNode(0))-1;
  pcoupled = false;

  if(myrank==0)
  {
    std::cout << "Warning: NOT all dof's are coupled!! " << std::endl << std::endl;
#endif
  }

  if(meshtyingoption != INPAR::FLUID::coupling_iontransport_laplace)
  {
    // feeding master nodes to the interface including ghosted nodes
    std::map<int, DRT::Node*>::const_iterator nodeiter;
    for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
    {
      DRT::Node* node = nodeiter->second;
      std::vector<int> dofids(dof);
      for (int k=0;k<dof;++k) dofids[k] = dis.Dof(node)[k];
      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                  new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                      dof, dofids, false));

      interface->AddMortarNode(mrtrnode);
    }

    // feeding slave nodes to the interface including ghosted nodes
    for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
    {
      DRT::Node* node = nodeiter->second;
      std::vector<int> dofids(dof);
      for (int k=0;k<dof;++k) dofids[k] = dis.Dof(node)[k];
      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                  new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                      dof, dofids, true));

      interface->AddMortarNode(mrtrnode);
    }
  }
  else
  {
    // only the potential is coupled (at position dof-1)
    // feeding master nodes to the interface including ghosted nodes
    std::map<int, DRT::Node*>::const_iterator nodeiter;
    for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
    {
      DRT::Node* node = nodeiter->second;
      std::vector<int> dofids(1);
      for (int k=0;k<1;++k) dofids[k] = dis.Dof(node)[dis.NumDof(dis.lRowNode(0))-1];
      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                  new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                      1, dofids, false));

      interface->AddMortarNode(mrtrnode);
    }

    // feeding slave nodes to the interface including ghosted nodes
    for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
    {
      DRT::Node* node = nodeiter->second;
      std::vector<int> dofids(1);
      for (int k=0;k<1;++k) dofids[k] = dis.Dof(node)[dis.NumDof(dis.lRowNode(0))-1];
      Teuchos::RCP<MORTAR::MortarNode> mrtrnode = Teuchos::rcp(
                  new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                      1, dofids, true));

      interface->AddMortarNode(mrtrnode);
    }
  }

  // max master element ID needed for unique eleIDs in interface discretization
  // will be used as offset for slave elements
  int EleOffset = dis.ElementRowMap()->MaxAllGID()+1;

  // feeding master elements to the interface
  std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator elemiter;
  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                new MORTAR::MortarElement(ele->Id(), ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), false));

    interface->AddMortarElement(mrtrele);
  }

  // feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    Teuchos::RCP<DRT::Element> ele = elemiter->second;
    Teuchos::RCP<MORTAR::MortarElement> mrtrele = Teuchos::rcp(
                new MORTAR::MortarElement(ele->Id() + EleOffset, ele->Owner(), ele->Shape(),
                    ele->NumNode(), ele->NodeIds(), true));

    interface->AddMortarElement(mrtrele);
  }

  // TODO: Difference between condensed and saddlePointproblem
  //       parallel distribution
  // finalize the contact interface construction
  if(meshtyingoption != INPAR::FLUID::sps_coupled or meshtyingoption != INPAR::FLUID::sps_pc) // Saddle point problem
  {
    interface->FillComplete(dis.DofRowMap()->MaxAllGID());
  }
  else
    interface->FillComplete();

  // store old row maps (before parallel redistribution)
  slavedofrowmap_  = Teuchos::rcp(new Epetra_Map(*interface->SlaveRowDofs()));
  masterdofrowmap_ = Teuchos::rcp(new Epetra_Map(*interface->MasterRowDofs()));


  // print parallel distribution
  interface->PrintParallelDistribution(1);

  //**********************************************************************
  // PARALLEL REDISTRIBUTION OF INTERFACE
  //**********************************************************************
  /*
  if (parredist && comm.NumProc()>1)
  {
    // redistribute optimally among all procs
    interface->Redistribute();

    // call fill complete again
    interface->FillComplete();

    // print parallel distribution again
    interface->PrintParallelDistribution(1);
  }
  */
  //**********************************************************************

  // create binary search tree
  interface->CreateSearchTree();

  // all the following stuff has to be done once in setup
  // in order to get initial D_ and M_

  // interface displacement (=0) has to be merged from slave and master discretization
  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(masterdofrowmap_,slavedofrowmap_, false);
  Teuchos::RCP<Epetra_Vector> dispn = LINALG::CreateVector(*dofrowmap, true);

  // set displacement state in mortar interface
  interface->SetState("displacement", dispn);

  // print message
  if(comm.MyPID()==0)
  {
    std::cout << "\nPerforming mortar coupling...............";
    fflush(stdout);
  }

  //in the following two steps MORTAR does all the work
  interface->Initialize();
  interface->Evaluate();

  // print message
  if(comm.MyPID()==0) std::cout << "done!" << std::endl << std::endl;

  // preparation for AssembleDM
  // (Note that redistslave and redistmaster are the slave and master row maps
  // after parallel redistribution. If no redistribution was performed, they
  // are of course identical to slavedofrowmap_/masterdofrowmap_!)
  Teuchos::RCP<Epetra_Map> redistslave  = interface->SlaveRowDofs();
  Teuchos::RCP<Epetra_Map> redistmaster = interface->MasterRowDofs();
  Teuchos::RCP<LINALG::SparseMatrix> dmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 10));
  Teuchos::RCP<LINALG::SparseMatrix> mmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 100));
  interface->AssembleDM(*dmatrix, *mmatrix);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*redistmaster, *redistslave);
  D_ = dmatrix;
  M_ = mmatrix;

  // Build Dinv
  Dinv_ = Teuchos::rcp(new LINALG::SparseMatrix(*D_));

  // extract diagonal of invd into diag
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*redistslave,true);
  Dinv_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  diag->Reciprocal(*diag);
  Dinv_->ReplaceDiagonalValues(*diag);
  Dinv_->Complete( D_->RangeMap(), D_->DomainMap() );
  DinvM_ = MLMultiply(*Dinv_,*M_,false,false,true);

  if(meshtyingoption != INPAR::FLUID::sps_coupled or meshtyingoption != INPAR::FLUID::sps_pc)  // Saddle point problem
  {
    // make numbering of LM dofs consecutive and unique across N interfaces
    int offset_if = 0;

    // build Lagrange multiplier dof map
    interface->UpdateLagMultSets(offset_if);

    glmdofrowmap_ = Teuchos::null;
    glmdofrowmap_ = LINALG::MergeMap(glmdofrowmap_, interface->LagMultDofs());

    offset_if = glmdofrowmap_->NumGlobalElements();
    if (offset_if < 0) offset_if = 0;

    // first setup
    Teuchos::RCP<LINALG::SparseMatrix> constrmt = Teuchos::rcp(new LINALG::SparseMatrix(*(dis.DofRowMap()),100,false,true));
    constrmt->Add(*D_,true,1.0,1.0);
    constrmt->Add(*M_,true,-1.0,1.0);
    constrmt->Complete(*slavedofrowmap_,*(dis.DofRowMap()));

    conmatrix_ = MORTAR::MatrixColTransformGIDs(constrmt,glmdofrowmap_);

  }
  // store interface
  interface_ = interface;

  // mesh initialization (for rotational invariance)
  //MeshInit(masterdis,slavedis,aledis,redistmaster,redistslave,comm,structslave);

  // only for parallel redistribution case
  /*
  if (parredist)
  {
    // transform everything back to old distribution
    D_     = MORTAR::MatrixRowColTransform(D_,slavedofrowmap_,slavedofrowmap_);
    M_     = MORTAR::MatrixRowColTransform(M_,slavedofrowmap_,masterdofrowmap_);
    Dinv_  = MORTAR::MatrixRowColTransform(Dinv_,slavedofrowmap_,slavedofrowmap_);
    DinvM_ = MORTAR::MatrixRowColTransform(DinvM_,slavedofrowmap_,masterdofrowmap_);
  }
  */

  // check for overlap of slave and Dirichlet boundaries
  // (this is not allowed in order to avoid over-constraint)

  bool overlap = false;
  Teuchos::ParameterList p;
  p.set("total time", 0.0);
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps = Teuchos::rcp(new LINALG::MapExtractor());
  Teuchos::RCP<Epetra_Vector > temp = LINALG::CreateVector(*(dis.DofRowMap()), true);
  dis.EvaluateDirichlet(p,temp,Teuchos::null,Teuchos::null,Teuchos::null,dbcmaps);

  // loop over all slave row nodes of the interface
  for (int j=0;j<interface_->SlaveRowNodes()->NumMyElements();++j)
  {
    int gid = interface_->SlaveRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // check if this node's dofs are in dbcmap
    for (int k=0;k<mtnode->NumDof();++k)
    {
      int currdof = mtnode->Dofs()[k];
      int lid = (dbcmaps->CondMap())->LID(currdof);

      // found slave node with dbc
      if (lid>=0)
      {
        overlap = true;
        break;
      }
    }
  }

  // print warning message to screen
  if (overlap && comm.MyPID()==0)
  {
    dserror("Slave boundary and Dirichlet boundary conditions overlap!\n"
        "This leads to over-constraint");
  }

  return pcoupled;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::MeshInit(DRT::Discretization& masterdis,
    DRT::Discretization& slavedis, DRT::Discretization& aledis,
    Teuchos::RCP<Epetra_Map> masterdofrowmap, Teuchos::RCP<Epetra_Map> slavedofrowmap,
    const Epetra_Comm& comm, bool structslave)
{
  // problem dimension
  const int dim = DRT::Problem::Instance()->NDim();

  //**********************************************************************
  // (0) check constraints in reference configuration
  //**********************************************************************
  // build global vectors of slave and master coordinates
  Teuchos::RCP<Epetra_Vector> xs = LINALG::CreateVector(*slavedofrowmap,true);
  Teuchos::RCP<Epetra_Vector> xm = LINALG::CreateVector(*masterdofrowmap,true);

  // loop over all slave row nodes
  for (int j=0; j<interface_->SlaveRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->SlaveRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // prepare assembly
    Epetra_SerialDenseVector val(dim);
    std::vector<int> lm(dim);
    std::vector<int> lmowner(dim);

    for (int k=0;k<dim;++k)
    {
      val[k] = mtnode->X()[k];
      lm[k] = mtnode->Dofs()[k];
      lmowner[k] = mtnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*xs,val,lm,lmowner);
  }

  // loop over all master row nodes
  for (int j=0; j<interface_->MasterRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->MasterRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // prepare assembly
    Epetra_SerialDenseVector val(dim);
    std::vector<int> lm(dim);
    std::vector<int> lmowner(dim);

    for (int k=0;k<dim;++k)
    {
      val[k] = mtnode->X()[k];
      lm[k] = mtnode->Dofs()[k];
      lmowner[k] = mtnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*xm,val,lm,lmowner);
  }

  // compute g-vector at global level
  Teuchos::RCP<Epetra_Vector> Dxs = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap));
  D_->Multiply(false,*xs,*Dxs);
  Teuchos::RCP<Epetra_Vector> Mxm = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap));
  M_->Multiply(false,*xm,*Mxm);
  Teuchos::RCP<Epetra_Vector > gold = LINALG::CreateVector(*slavedofrowmap, true);
  gold->Update(1.0,*Dxs,1.0);
  gold->Update(-1.0,*Mxm,1.0);
  double gnorm = 0.0;
  gold->Norm2(&gnorm);

  // no need to do mesh initialization if g already very small
  if (gnorm < 1.0e-12) return;

  // print message
  if(comm.MyPID()==0)
  {
    std::cout << "Performing mesh initialization...........";
    fflush(stdout);
  }

  //**********************************************************************
  // (1) get master positions on global level
  //**********************************************************************
  // fill Xmaster first
  Teuchos::RCP<Epetra_Vector> Xmaster = LINALG::CreateVector(*masterdofrowmap, true);

  // loop over all master row nodes on the current interface
  for (int j=0; j<interface_->MasterRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->MasterRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // do assembly (overwrite duplicate nodes)
    for (int k=0;k<dim;++k)
    {
      int dof = mtnode->Dofs()[k];
      (*Xmaster)[(Xmaster->Map()).LID(dof)] = mtnode->X()[k];
    }
  }

  //**********************************************************************
  // (2) solve for modified slave positions on global level
  //**********************************************************************
  // initialize modified slave positions
  Teuchos::RCP<Epetra_Vector> Xslavemod = LINALG::CreateVector(*slavedofrowmap,true);

  // this is trivial for dual Lagrange multipliers
  DinvM_->Multiply(false,*Xmaster,*Xslavemod);


  //**********************************************************************
  // (3) perform mesh initialization node by node
  //**********************************************************************
  // export Xslavemod to fully overlapping column map for current interface
  Teuchos::RCP<Epetra_Map> fullsdofs  = LINALG::AllreduceEMap(*(interface_->SlaveRowDofs()));
  Teuchos::RCP<Epetra_Map> fullsnodes = LINALG::AllreduceEMap(*(interface_->SlaveRowNodes()));
  Epetra_Vector Xslavemodcol(*fullsdofs,false);
  LINALG::Export(*Xslavemod,Xslavemodcol);

  // loop over all slave nodes on the current interface
  for (int j=0; j<fullsnodes->NumMyElements(); ++j)
  {
    // get global ID of current node
    int gid = fullsnodes->GID(j);

    // be careful to modify BOTH mtnode in interface discret ...
    // (check if the node is available on this processor)
    bool isininterfacecolmap = false;
    int ilid = interface_->SlaveColNodes()->LID(gid);
    if (ilid>=0) isininterfacecolmap = true;
    DRT::Node* node = NULL;
    MORTAR::MortarNode* mtnode = NULL;
    if (isininterfacecolmap)
    {
      node = interface_->Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      mtnode = static_cast<MORTAR::MortarNode*>(node);
    }

    // ... AND standard node in underlying slave discret
    // (check if the node is available on this processor)
    bool isinproblemcolmap = false;
    int lid = slavedis.NodeColMap()->LID(gid);
    if (lid>=0) isinproblemcolmap = true;
    DRT::Node* pnode = NULL;
    if (isinproblemcolmap)
    {
      pnode = slavedis.gNode(gid);
      if (!pnode) dserror("ERROR: Cannot find node with gid %",gid);
    }

    // ... AND standard node in ALE discret if fluid=slave
    // (check if the node is available on this processor)
    bool isinproblemcolmap2 = false;
    int lid2 = aledis.NodeColMap()->LID(gid);
    if (lid2>=0) isinproblemcolmap2 = true;
    DRT::Node* alenode = NULL;
    if (isinproblemcolmap2)
    {
      alenode = aledis.gNode(gid);
      if (!structslave && !alenode) dserror("ERROR: Cannot find node with gid %",gid);
    }

    // new nodal position and problem dimension
    double Xnew[3] = {0.0, 0.0, 0.0};
    double Xnewglobal[3] = {0.0, 0.0, 0.0};

    //******************************************************************
    // compute new nodal position
    //******************************************************************
    // first sort out procs that do not know of mtnode
    if (isininterfacecolmap)
    {
      // owner processor of this node will do computation
      if (comm.MyPID()==mtnode->Owner())
      {
        // get corresponding entries from Xslavemod
        int numdof = mtnode->NumDof();
        if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

        // find DOFs of current node in Xslavemod and extract this node's position
        std::vector<int> locindex(numdof);

        for (int dof=0;dof<numdof;++dof)
        {
          locindex[dof] = (Xslavemodcol.Map()).LID(mtnode->Dofs()[dof]);
          if (locindex[dof]<0) dserror("ERROR: Did not find dof in map");
          Xnew[dof] = Xslavemodcol[locindex[dof]];
        }

        // check is mesh distortion is still OK
        // (throw a dserror if length of relocation is larger than 80%
        // of an adjacent element edge -> see Puso, IJNME, 2004)
        double limit = 0.8;
        double relocation = 0.0;
        if (dim==2)
        {
          relocation = sqrt((Xnew[0]-mtnode->X()[0])*(Xnew[0]-mtnode->X()[0])
                           +(Xnew[1]-mtnode->X()[1])*(Xnew[1]-mtnode->X()[1]));
        }
        else if (dim==3)
        {
          relocation = sqrt((Xnew[0]-mtnode->X()[0])*(Xnew[0]-mtnode->X()[0])
                           +(Xnew[1]-mtnode->X()[1])*(Xnew[1]-mtnode->X()[1])
                           +(Xnew[2]-mtnode->X()[2])*(Xnew[2]-mtnode->X()[2]));
        }
        else dserror("ERROR: Problem dimension must be either 2 or 3!");
        bool isok = mtnode->CheckMeshDistortion(relocation,limit);
        if (!isok) dserror("ERROR: Mesh distortion generated by relocation is too large!");
      }
    }

    // communicate new position Xnew to all procs
    // (we can use SumAll here, as Xnew will be zero on all processors
    // except for the owner processor of the current node)
    comm.SumAll(&Xnew[0],&Xnewglobal[0],3);

    // const_cast to force modifed X() into mtnode
    // const_cast to force modifed xspatial() into mtnode
    // const_cast to force modifed X() into pnode
    // const_cast to force modifed X() into alenode if fluid=slave
    // (remark: this is REALLY BAD coding)
    for (int k=0;k<dim;++k)
    {
      // modification in interface discretization
      if (isininterfacecolmap)
      {
        const_cast<double&>(mtnode->X()[k])        = Xnewglobal[k];
        const_cast<double&>(mtnode->xspatial()[k]) = Xnewglobal[k];
      }

      // modification in problem discretization
      if (isinproblemcolmap)
        const_cast<double&>(pnode->X()[k])       = Xnewglobal[k];

      // modification in ALE discretization
      if (isinproblemcolmap2 && !structslave)
        const_cast<double&>(alenode->X()[k])     = Xnewglobal[k];
    }
  }

  //**********************************************************************
  // (4) re-evaluate constraints in reference configuration
  //**********************************************************************
  // build global vectors of slave and master coordinates
  xs = LINALG::CreateVector(*slavedofrowmap,true);
  xm = LINALG::CreateVector(*masterdofrowmap,true);

  // loop over all slave row nodes
  for (int j=0; j<interface_->SlaveRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->SlaveRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // prepare assembly
    Epetra_SerialDenseVector val(dim);
    std::vector<int> lm(dim);
    std::vector<int> lmowner(dim);

    for (int k=0;k<dim;++k)
    {
      val[k] = mtnode->X()[k];
      lm[k] = mtnode->Dofs()[k];
      lmowner[k] = mtnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*xs,val,lm,lmowner);
  }

  // loop over all master row nodes
  for (int j=0; j<interface_->MasterRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->MasterRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

    // prepare assembly
    Epetra_SerialDenseVector val(dim);
    std::vector<int> lm(dim);
    std::vector<int> lmowner(dim);

    for (int k=0;k<dim;++k)
    {
      val[k] = mtnode->X()[k];
      lm[k] = mtnode->Dofs()[k];
      lmowner[k] = mtnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*xm,val,lm,lmowner);
  }

  // compute g-vector at global level
  Dxs = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap));
  D_->Multiply(false,*xs,*Dxs);
  Mxm = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap));
  M_->Multiply(false,*xm,*Mxm);
  Teuchos::RCP<Epetra_Vector > gnew = LINALG::CreateVector(*slavedofrowmap, true);
  gnew->Update(1.0,*Dxs,1.0);
  gnew->Update(-1.0,*Mxm,1.0);
  gnew->Norm2(&gnorm);

  // error if g is still non-zero
  if (gnorm > 1.0e-12) dserror("ERROR: Mesh initialization was not successful!");

  //**********************************************************************
  // (5) re-initialize finite elements (if slave=structure)
  //**********************************************************************
  // if slave=fluid, we are lucky because fluid elements do not
  // need any re-initialization (unlike structural elements)
  if (structslave) DRT::ParObjectFactory::Instance().InitializeElements(slavedis);

  // print message
  if (comm.MyPID()==0) std::cout << "done!" << std::endl;

  return;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Evaluate(Teuchos::RCP<Epetra_Vector> idisp)
{
  // set new displacement state in mortar interface
  interface_->SetState("displacement", idisp);
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Evaluate(Teuchos::RCP<Epetra_Vector> idispma, Teuchos::RCP<Epetra_Vector> idispsl)
{

  const Epetra_BlockMap stdmap = idispsl->Map();
  idispsl->ReplaceMap(*slavedofrowmap_);

  Teuchos::RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*masterdofrowmap_,*slavedofrowmap_, true);
  Teuchos::RCP<Epetra_Import> msimpo = Teuchos::rcp(new Epetra_Import(*dofrowmap,*masterdofrowmap_));
  Teuchos::RCP<Epetra_Import> slimpo = Teuchos::rcp(new Epetra_Import(*dofrowmap,*slavedofrowmap_));

  Teuchos::RCP<Epetra_Vector> idispms = LINALG::CreateVector(*dofrowmap,true);

  idispms -> Import(*idispma,*msimpo,Add);
  idispms -> Import(*idispsl,*slimpo,Add);

  // set new displacement state in mortar interface
  interface_->SetState("displacement", idispms);

  Evaluate();

  idispsl->ReplaceMap(stdmap);

  return;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Evaluate()
{
  // check for parallel redistribution
  bool parredist = false;
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MortarCouplingParams();
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(input,"PARALLEL_REDIST") != INPAR::MORTAR::parredist_none)
    parredist = true;

  // in the following two steps MORTAR does all the work for new interface displacements
  interface_->Initialize();
  interface_->Evaluate();

  // preparation for AssembleDM
  // (Note that redistslave and redistmaster are the slave and master row maps
  // after parallel redistribution. If no redistribution was performed, they
  // are of course identical to slavedofrowmap_/masterdofrowmap_!)
  Teuchos::RCP<Epetra_Map> redistslave  = interface_->SlaveRowDofs();
  Teuchos::RCP<Epetra_Map> redistmaster = interface_->MasterRowDofs();
  Teuchos::RCP<LINALG::SparseMatrix> dmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 10));
  Teuchos::RCP<LINALG::SparseMatrix> mmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*redistslave, 100));
  interface_->AssembleDM(*dmatrix, *mmatrix);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*redistmaster, *redistslave);
  D_ = dmatrix;
  M_ = mmatrix;

  // Build Dinv
  Dinv_ = Teuchos::rcp(new LINALG::SparseMatrix(*D_));

  // extract diagonal of invd into diag
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*redistslave,true);
  Dinv_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  diag->Reciprocal(*diag);
  Dinv_->ReplaceDiagonalValues(*diag);
  Dinv_->Complete( D_->RangeMap(), D_->DomainMap() );
  DinvM_ = MLMultiply(*Dinv_,*M_,false,false,true);

  // only for parallel redistribution case
  if (parredist)
  {
    // transform everything back to old distribution
    D_     = MORTAR::MatrixRowColTransform(D_,slavedofrowmap_,slavedofrowmap_);
    M_     = MORTAR::MatrixRowColTransform(M_,slavedofrowmap_,masterdofrowmap_);
    Dinv_  = MORTAR::MatrixRowColTransform(Dinv_,slavedofrowmap_,slavedofrowmap_);
    DinvM_ = MORTAR::MatrixRowColTransform(DinvM_,slavedofrowmap_,masterdofrowmap_);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::CouplingMortar::MasterToSlave
(
  Teuchos::RCP<Epetra_Vector> mv
) const
{
  dsassert( masterdofrowmap_->SameAs( mv->Map() ),
      "Vector with master dof map expected" );

  Epetra_Vector tmp = Epetra_Vector(M_->RowMap());

  if (M_->Multiply(false, *mv, tmp))
    dserror( "M*mv multiplication failed" );

  Teuchos::RCP<Epetra_Vector> sv = Teuchos::rcp( new Epetra_Vector( *slavedofrowmap_ ) );

  if ( Dinv_->Multiply( false, tmp, *sv ) )
    dserror( "D^{-1}*v multiplication failed" );

  return sv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::CouplingMortar::SlaveToMaster
(
  Teuchos::RCP<Epetra_Vector> sv
) const
{
  Epetra_Vector tmp = Epetra_Vector(M_->RangeMap());
  std::copy(sv->Values(), sv->Values() + sv->MyLength(), tmp.Values());

  Teuchos::RCP<Epetra_Vector> mv = Teuchos::rcp(new Epetra_Vector(*masterdofrowmap_));
  if (M_->Multiply(true, tmp, *mv))
    dserror( "M^{T}*sv multiplication failed" );

  return mv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::CouplingMortar::MasterToSlave
(
  Teuchos::RCP<const Epetra_Vector> mv
) const
{
  dsassert( masterdofrowmap_->SameAs( mv->Map() ),
      "Vector with master dof map expected" );

  Epetra_Vector tmp = Epetra_Vector(M_->RowMap());

  if (M_->Multiply(false, *mv, tmp))
    dserror( "M*mv multiplication failed" );

  Teuchos::RCP<Epetra_Vector> sv = Teuchos::rcp( new Epetra_Vector( *slavedofrowmap_ ) );

  if ( Dinv_->Multiply( false, tmp, *sv ) )
    dserror( "D^{-1}*v multiplication failed" );

  return sv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::CouplingMortar::SlaveToMaster
(
  Teuchos::RCP<const Epetra_Vector> sv
) const
{
  Epetra_Vector tmp = Epetra_Vector(M_->RangeMap());
  std::copy(sv->Values(), sv->Values() + sv->MyLength(), tmp.Values());

  Teuchos::RCP<Epetra_Vector> mv = Teuchos::rcp(new Epetra_Vector(*masterdofrowmap_));
  if (M_->Multiply(true, tmp, *mv))
    dserror( "M^{T}*sv multiplication failed" );

  return mv;
}

