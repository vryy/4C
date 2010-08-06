/*----------------------------------------------------------------------*/
/*!
 \file adapter_coupling_mortar.cpp

 \brief

 <pre>
 Maintainer: Thomas Kloeppel
 kloeppel@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15257
 </pre>
 */
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_coupling_mortar.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

using namespace std;
extern struct _GENPROB genprob;

ADAPTER::CouplingMortar::CouplingMortar()
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Setup(DRT::Discretization& masterdis,
    DRT::Discretization& slavedis, DRT::Discretization& aledis,
    Epetra_Comm& comm, bool structslave)
{
  // initialize maps for row nodes
  map<int, DRT::Node*> masternodes;
  map<int, DRT::Node*> slavenodes;

  // initialize maps for column nodes
  map<int, DRT::Node*> mastergnodes;
  map<int, DRT::Node*> slavegnodes;

  //initialize maps for elements
  map<int, RefCountPtr<DRT::Element> > masterelements;
  map<int, RefCountPtr<DRT::Element> > slaveelements;

  // Fill maps based on condition for master side
  DRT::UTILS::FindConditionObjects(masterdis, masternodes, mastergnodes, masterelements,
      "FSICoupling");

  // Fill maps based on condition for slave side
  DRT::UTILS::FindConditionObjects(slavedis, slavenodes, slavegnodes, slaveelements,
      "FSICoupling");

  // get structural dynamics parameter
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->MeshtyingAndContactParams();

  // check for invalid paramater values
  if (Teuchos::getIntegralValue<INPAR::MORTAR::ShapeFcn>(input,"SHAPEFCN") != INPAR::MORTAR::shape_dual)
    dserror("Mortar coupling adapter only works for dual shape functions");

  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  // IMPORTANT: We assume that all nodes have 'dim' DoF, that have to be considered for coupling.
  //            Possible pressure DoF are not transferred to MortarInterface.
  const int dim = genprob.ndim;
  RCP<MORTAR::MortarInterface> interface = rcp(new MORTAR::MortarInterface(0, comm, dim, input));

  // feeding master nodes to the interface including ghosted nodes
  vector<int> masterdofs;
  map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    RCP<MORTAR::MortarNode> mrtrnode = rcp(
                new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                    dim, masterdis.Dof(node), false));

    interface->AddMortarNode(mrtrnode);
  }

  // build master dof row map
  for (nodeiter = masternodes.begin(); nodeiter != masternodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    vector<int> dofs = masterdis.Dof(node);

    dofs.resize(dim);
    masterdofs.insert(masterdofs.end(), dofs.begin(), dofs.end());

  }

  //build Epetra_Map of master row dofs
  masterdofrowmap_ = rcp(new Epetra_Map(-1, masterdofs.size(), &masterdofs[0], 0, comm));

  //feeding slave nodes to the interface including ghosted nodes
  vector<int> slaverowdofs;
  vector<int> slavecoldofs;

  for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;

    RCP<MORTAR::MortarNode> mrtrnode = rcp(
                new MORTAR::MortarNode(node->Id(), node->X(), node->Owner(),
                    dim, slavedis.Dof(node), true));

    interface->AddMortarNode(mrtrnode);

    vector<int> dofs = slavedis.Dof(node);
    dofs.resize(dim);
    slavecoldofs.insert(slavecoldofs.end(), dofs.begin(), dofs.end());
  }

  // build slave dof row map
  for (nodeiter = slavenodes.begin(); nodeiter != slavenodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    vector<int> dofs = slavedis.Dof(node);
    dofs.resize(dim);
    slaverowdofs.insert(slaverowdofs.end(), dofs.begin(), dofs.end());

  }
  slavenodes.clear();

  //build map of slave dofs and slave mortar dofs without pressure
  slavedofrowmap_ = rcp(
      new Epetra_Map(-1, slaverowdofs.size(), &slaverowdofs[0], 0, comm));
  slavedofcolmap_ = rcp(
      new Epetra_Map(-1, slavecoldofs.size(), &slavecoldofs[0], 0, comm));

  //feeding master elements to the interface
  map<int, RefCountPtr<DRT::Element> >::const_iterator elemiter;

  // max master element ID needed for unique eleIDs in interface discretization
  // will be used as offset for slave elements
  int EleOffset = masterdis.ElementRowMap()->MaxAllGID()+1;

  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    RefCountPtr<DRT::Element> ele = elemiter->second;

    RCP<MORTAR::MortarElement>
        mrtrele =
            rcp(
                new MORTAR::MortarElement(ele->Id(), ele->Owner(), ele->Shape(), ele->NumNode(), ele->NodeIds(), false));

    interface->AddMortarElement(mrtrele);

  }

  //feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    RefCountPtr<DRT::Element> ele = elemiter->second;
    RCP<MORTAR::MortarElement>
        mrtrele =
            rcp(
                new MORTAR::MortarElement(ele->Id() + EleOffset, ele->Owner(), ele->Shape(), ele->NumNode(), ele->NodeIds(), true));

    interface->AddMortarElement(mrtrele);
  }

  // finalize the contact interface construction
  interface->FillComplete();

  // create binary search tree
  interface->CreateSearchTree();

  // all the following stuff has to be done once in setup
  // in order to get initial D_ and M_

  // interface displacement (=0) has to be merged from slave and master discretization
  RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(masterdofrowmap_,slavedofrowmap_, true);
  RCP<Epetra_Vector> dispn = LINALG::CreateVector(*dofrowmap, true);

  interface->SetState("displacement", dispn);
  //in the following two steps MORTAR does all the work
  interface->Initialize();
  interface->Evaluate();

  // preparation for AssembleDM
  RCP<Epetra_Map> slavedofrowmap = interface->SlaveRowDofs();
  RCP<Epetra_Map> masterdofrowmap = interface->MasterRowDofs();

  RCP<LINALG::SparseMatrix> dmatrix = rcp(
      new LINALG::SparseMatrix(*slavedofrowmap, 10));

  RCP<LINALG::SparseMatrix> mmatrix = rcp(
      new LINALG::SparseMatrix(*slavedofrowmap, 100));

  interface->AssembleDM(*dmatrix, *mmatrix);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*masterdofrowmap, *slavedofrowmap);

  D_ = dmatrix;

  M_ = mmatrix;

  //Build Dinv
  Dinv_ = rcp(new LINALG::SparseMatrix(*D_));

  RCP<Epetra_Vector> diag = LINALG::CreateVector(*slavedofrowmap,true);

  // extract diagonal of invd into diag
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
  MeshInit(masterdis,slavedis,aledis,masterdofrowmap,slavedofrowmap,comm,structslave);

  // check for overlap of slave and Dirichlet boundaries
  // (this is not allowed in order to avoid over-constraint)
  bool overlap = false;
	Teuchos::ParameterList p;
	p.set("total time", 0.0);
	RCP<LINALG::MapExtractor> dbcmaps = rcp(new LINALG::MapExtractor());
	RCP<Epetra_Vector > temp = LINALG::CreateVector(*(slavedis.DofRowMap()), true);
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
  	cout << RED << "WARNING: Slave boundary and Dirichlet boundary conditions overlap!" << endl;
  	cout << "This leads to over-constraint, so you might encounter some problems!" << END_COLOR << endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::MeshInit(DRT::Discretization& masterdis,
    DRT::Discretization& slavedis, DRT::Discretization& aledis,
    RCP<Epetra_Map> masterdofrowmap, RCP<Epetra_Map> slavedofrowmap,
    Epetra_Comm& comm, bool structslave)
{
	// problem dimension
	int dim = genprob.ndim;

	//**********************************************************************
	// (0) check constraints in reference configuration
	//**********************************************************************
	// intialize and assemble g-vector
	RCP<Epetra_Vector> gold = LINALG::CreateVector(*slavedofrowmap, true);
	interface_->AssembleG(*gold);
	double gnorm = 0.0;
	gold->Norm2(&gnorm);

	// no need to do mesh initialization if g already very small
	if (gnorm < 1.0e-12) return;

  // print message
  if(comm.MyPID()==0)
  {
    cout << "Performing mesh initialization...........";
    fflush(stdout);
  }

	//**********************************************************************
	// (1) get master positions on global level
	//**********************************************************************
	// fill Xmaster first
	RCP<Epetra_Vector> Xmaster = LINALG::CreateVector(*masterdofrowmap, true);

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
	RCP<Epetra_Vector> Xslavemod = LINALG::CreateVector(*slavedofrowmap,true);

	// this is trivial for dual Lagrange multipliers
	DinvM_->Multiply(false,*Xmaster,*Xslavemod);


	//**********************************************************************
	// (3) perform mesh initialization node by node
	//**********************************************************************
	// export Xslavemod to standard column map overlap for current interface
	Epetra_Vector Xslavemodcol(*(interface_->SlaveColDofs()),false);
	LINALG::Export(*Xslavemod,Xslavemodcol);

	// loop over all slave column nodes on the current interface
	for (int j=0; j<interface_->SlaveColNodes()->NumMyElements(); ++j)
	{
		int gid = interface_->SlaveColNodes()->GID(j);

		// be careful to modify BOTH mtnode in interface discret ...
		DRT::Node* node = interface_->Discret().gNode(gid);
		if (!node) dserror("ERROR: Cannot find node with gid %",gid);
		MORTAR::MortarNode* mtnode = static_cast<MORTAR::MortarNode*>(node);

		// ... AND standard node in underlying slave discret
		DRT::Node* pnode = slavedis.gNode(gid);
		if (!pnode) dserror("ERROR: Cannot find node with gid %",gid);

		//... AND standard node in ALE discret if fluid=slave
		DRT::Node* alenode = aledis.gNode(gid);
		if (!structslave && !alenode) dserror("ERROR: Cannot find node with gid %",gid);

		// new nodal position
		double Xnew[3] = {0.0, 0.0, 0.0};

		// get corresponding entries from Xslavemod
		int numdof = mtnode->NumDof();
		if (dim!=numdof) dserror("ERROR: Inconsisteny Dim <-> NumDof");

		// find DOFs of current node in Xslavemod and extract this node's position
		vector<int> locindex(numdof);

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

		// const_cast to force modifed X() into mtnode
		// const_cast to force modifed xspatial() into mtnode
		// const_cast to force modifed X() into pnode
		// const_cast to force modifed X() into alenode if fluid=slave
		// (remark: this is REALLY BAD coding)
		for (int k=0;k<dim;++k)
		{
			const_cast<double&>(mtnode->X()[k])                    = Xnew[k];
			const_cast<double&>(mtnode->xspatial()[k])             = Xnew[k];
			const_cast<double&>(pnode->X()[k])                     = Xnew[k];
			if (!structslave) const_cast<double&>(alenode->X()[k]) = Xnew[k];
		}
	}

	//**********************************************************************
	// (4) re-evaluate constraints in reference configuration
	//**********************************************************************
	// intialize and assemble g-vector
	RCP<Epetra_Vector > gnew = LINALG::CreateVector(*slavedofrowmap, true);
	interface_->AssembleG(*gnew);
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
  if (comm.MyPID()==0) cout << "done!\n" << endl;

  return;

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Evaluate(RCP<Epetra_Vector> idisp)
{
  interface_->SetState("displacement", idisp);
  //in the following two steps MORTAR does all the work for new interface displacements
  interface_->Initialize();
  interface_->Evaluate();

  //Preparation for AssembleDMG
  RCP<Epetra_Map> slavedofrowmap = interface_->SlaveRowDofs();
  RCP<Epetra_Map> masterdofrowmap = interface_->MasterRowDofs();

  RCP<LINALG::SparseMatrix> dmatrix = rcp(
      new LINALG::SparseMatrix(*slavedofrowmap, 10));

  RCP<LINALG::SparseMatrix> mmatrix = rcp(
      new LINALG::SparseMatrix(*slavedofrowmap, 100));

  interface_->AssembleDM(*dmatrix, *mmatrix);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*masterdofrowmap, *slavedofrowmap);

  D_ = dmatrix;

  M_ = mmatrix;

  //Build Dinv
  Dinv_ = rcp(new LINALG::SparseMatrix(*D_));

  RCP<Epetra_Vector> diag = LINALG::CreateVector(*slavedofrowmap,true);

  // extract diagonal of invd into diag
  Dinv_->ExtractDiagonalCopy(*diag);

  // set zero diagonal values to dummy 1.0
  for (int i=0;i<diag->MyLength();++i)
    if ((*diag)[i]==0.0) (*diag)[i]=1.0;

  // scalar inversion of diagonal values
  diag->Reciprocal(*diag);
  Dinv_->ReplaceDiagonalValues(*diag);

  Dinv_->Complete( D_->RangeMap(), D_->DomainMap() );

  DinvM_ = MLMultiply(*Dinv_,*M_,false,false,true);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> ADAPTER::CouplingMortar::MasterToSlave
(
  RefCountPtr<Epetra_Vector> mv
) const
{
  dsassert( masterdofrowmap_->SameAs( mv->Map() ),
      "Vector with master dof map expected" );

  Epetra_Vector tmp = Epetra_Vector(M_->RowMap());

  if (M_->Multiply(false, *mv, tmp))
    dserror( "M*mv multiplication failed" );

  RefCountPtr<Epetra_Vector> sv = rcp( new Epetra_Vector( *slavedofrowmap_ ) );

  if ( Dinv_->Multiply( false, tmp, *sv ) )
    dserror( "D^{-1}*v multiplication failed" );

  return sv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> ADAPTER::CouplingMortar::SlaveToMaster
(
  RefCountPtr<Epetra_Vector> sv
) const
{
  Epetra_Vector tmp = Epetra_Vector(M_->RangeMap());
  copy(sv->Values(), sv->Values() + sv->MyLength(), tmp.Values());

  RefCountPtr<Epetra_Vector> mv = rcp(new Epetra_Vector(*masterdofrowmap_));
  if (M_->Multiply(true, tmp, *mv))
    dserror( "M^{T}*sv multiplication failed" );

  return mv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> ADAPTER::CouplingMortar::MasterToSlave
(
  RefCountPtr<const Epetra_Vector> mv
) const
{
  dsassert( masterdofrowmap_->SameAs( mv->Map() ),
      "Vector with master dof map expected" );

  Epetra_Vector tmp = Epetra_Vector(M_->RowMap());

  if (M_->Multiply(false, *mv, tmp))
    dserror( "M*mv multiplication failed" );

  RefCountPtr<Epetra_Vector> sv = rcp( new Epetra_Vector( *slavedofrowmap_ ) );

  if ( Dinv_->Multiply( false, tmp, *sv ) )
    dserror( "D^{-1}*v multiplication failed" );

  return sv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> ADAPTER::CouplingMortar::SlaveToMaster
(
  RefCountPtr<const Epetra_Vector> sv
) const
{
  Epetra_Vector tmp = Epetra_Vector(M_->RangeMap());
  copy(sv->Values(), sv->Values() + sv->MyLength(), tmp.Values());

  RefCountPtr<Epetra_Vector> mv = rcp(new Epetra_Vector(*masterdofrowmap_));
  if (M_->Multiply(true, tmp, *mv))
    dserror( "M^{T}*sv multiplication failed" );

  return mv;
}

RCP<LINALG::SparseMatrix> ADAPTER::CouplingMortar::GetMortarTrafo()
{
  return DinvM_;
}

#endif // CCADISCRET
