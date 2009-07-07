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
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_contact/drt_contact_manager.H"
#include "../drt_contact/drt_contact_interface.H"
#include "../drt_contact/contactdefines.H"
//#include "../drt_lib/linalg_ana.H"
#include "../drt_io/io.H"
#include "../drt_lib/linalg_utils.H"

using namespace std;
extern struct _GENPROB genprob;

ADAPTER::CouplingMortar::CouplingMortar()
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Setup(const DRT::Discretization& masterdis,
    const DRT::Discretization& slavedis, Epetra_Comm& comm)
{
  // initialize maps for nodes and elements
  map<int, DRT::Node*> masternodes;
  map<int, DRT::Node*> slavenodes;
  
  map<int, DRT::Node*> mastergnodes;
  map<int, DRT::Node*> slavegnodes;

  map<int, RefCountPtr<DRT::Element> > masterelements;
  map<int, RefCountPtr<DRT::Element> > slaveelements;


  DRT::UTILS::FindInterfaceObjects(masterdis, masternodes, mastergnodes, masterelements,
      "FSICoupling");

  DRT::UTILS::FindInterfaceObjects(slavedis, slavenodes, slavegnodes, slaveelements,
      "FSICoupling");
  //	parameter list for contact definition 
  // const Teuchos::ParameterList& input = DRT::Problem::Instance()->StructuralContactParams();
  
  Teuchos::ParameterList tmpinput;
  tmpinput.set<string> ("search algorithm", "binarytree");
  tmpinput.set<double> ("search parameter", 0.3);
  // get problem dimension (2D or 3D) and initialize (CONTACT::) interface
  const int dim = genprob.ndim;
  RCP<CONTACT::Interface> interface = rcp(
      new CONTACT::Interface(0, comm, dim, tmpinput));

  
  //feeding master nodes to the interface
  vector<int> masterdofs;
  map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = mastergnodes.begin(); nodeiter != mastergnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;

    RCP<CONTACT::CNode> cnode = rcp(
                new CONTACT::CNode(node->Id(), node->X(), node->Owner(), 
                    masterdis.NumDof(node), masterdis.Dof(node), false, true));

    interface->AddCNode(cnode);
  }
  
  for (nodeiter = masternodes.begin(); nodeiter != masternodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    vector<int> dofs = masterdis.Dof(node);
    masterdofs.insert(masterdofs.end(), dofs.begin(), dofs.end());

  }

  //build map of master dofs
  masterdofmap_ = rcp(
      new Epetra_Map(-1, masterdofs.size(), &masterdofs[0], 0, comm));

  //feeding slave nodes to the interface       
  vector<int> slavedofs;
  vector<int> slavemortardofs;
  int dofoffset = masterdofmap_->MaxAllGID() + 1;

  
  for (nodeiter = slavegnodes.begin(); nodeiter != slavegnodes.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;

    RCP<CONTACT::CNode> cnode = rcp(
                new CONTACT::CNode(node->Id(), node->X(), node->Owner(), 
                    slavedis.NumDof(node)-1, slavedis.Dof(node), true, true));

    interface->AddCNode(cnode);
  }
  
  for (nodeiter = slavenodes.begin(); nodeiter != slavenodes.end(); ++nodeiter)
  {
    // We expect to have a pressure dof at each node. Mortar
    // couples just the displacements, so remove the pressure dof.

    DRT::Node* node = nodeiter->second;
    vector<int> dofs = slavedis.Dof(node);
    dofs.resize(dofs.size() - 1);
    slavedofs.insert(slavedofs.end(), dofs.begin(), dofs.end());

  }
  slavenodes.clear();

  //build map of slave dofs and slave mortar dofs without pressure
  slavedofmap_ = rcp(
      new Epetra_Map(-1, slavedofs.size(), &slavedofs[0], 0, comm));

  //feeding master elements to the interface
  map<int, RefCountPtr<DRT::Element> >::const_iterator elemiter;
  
  // max master element ID needed for unique eleIDs in interface discretization
  // will be used as offset for slave elements
  int EleOffset = masterdis.ElementRowMap()->MaxAllGID()+1; 

  for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
  {
    RefCountPtr<DRT::Element> ele = elemiter->second;

    RCP<CONTACT::CElement>
        cele =
            rcp(
                new CONTACT::CElement(ele->Id(), DRT::Element::element_contact, ele->Owner(), ele->Shape(), ele->NumNode(), ele->NodeIds(), false));

    interface->AddCElement(cele);

  }

  //feeding slave elements to the interface
  for (elemiter = slaveelements.begin(); elemiter != slaveelements.end(); ++elemiter)
  {
    RefCountPtr<DRT::Element> ele = elemiter->second;
    RCP<CONTACT::CElement>
        cele =
            rcp(
                new CONTACT::CElement(ele->Id() + EleOffset, DRT::Element::element_contact, ele->Owner(), ele->Shape(), ele->NumNode(), ele->NodeIds(), true));

    interface->AddCElement(cele);
  }

  //finalize the contact interface construction	
  interface->FillComplete();
  
  // all the following stuff has to be done once in setup 
  // in order to get initial D_ and M_ 

  // interface displacement of the first time step (=0) has to be merged from
  // slave and master discretization
  RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(masterdofmap_,slavedofmap_, true);

  RCP<Epetra_Vector> dispn = LINALG::CreateVector(*dofrowmap, true);

  interface->SetState("displacement", dispn);

  //in the following two steps do the hard work; thanks to CONTACT
  interface->Initialize();
  interface->Evaluate();

  //Preparation for AssembleDMG
  RCP<Epetra_Map> slavedofrowmap = interface->SlaveRowDofs();
  RCP<Epetra_Map> masterdofrowmap = interface->MasterRowDofs();

  RCP<LINALG::SparseMatrix> dmatrix = rcp(
      new LINALG::SparseMatrix(*slavedofrowmap, 10));

  RCP<LINALG::SparseMatrix> mmatrix = rcp(
      new LINALG::SparseMatrix(*slavedofrowmap, 100));

  RCP<Epetra_Map> gwighting = interface->SlaveRowNodes();

  RCP<Epetra_Vector> g = rcp(new Epetra_Vector(*gwighting, true));

  interface->AssembleDMG(*dmatrix, *mmatrix, *g);

  // Complete() global Mortar matrices
  dmatrix->Complete();
  mmatrix->Complete(*masterdofrowmap, *slavedofrowmap);

  D_ = dmatrix->EpetraMatrix();

  M_ = mmatrix->EpetraMatrix();

  //Build Dinv 		
  Dinv_ = rcp(
      new Epetra_CrsMatrix(Copy, D_->DomainMap(), D_->RangeMap(), 1, true));

  const Epetra_Map& Dinvmap = Dinv_->RowMap();
  const Epetra_Map& Dmap = D_->RowMap();
  for (int row = 0; row < Dinvmap.NumMyElements(); ++row)
  {
    int rowgid = Dinvmap.GID(row);
    int colgid = Dmap.GID(row);

    int numentries;
    double *values;
    int *indices;

    // do some validation
    
    if (D_->ExtractMyRowView(row, numentries, values, indices))
      dserror("ExtractMyRowView failed");
    double value = 1.0/values[0];
    if ( Dinv_->InsertGlobalValues(rowgid, 1, &value, &colgid) )
        dserror( "InsertGlobalValues failed" );

  }

  if ( Dinv_->FillComplete( D_->RangeMap(), D_->DomainMap() ) )
    dserror( "Filling failed" );

  //store interface	
  interface_ = interface;  
  

  
  return;
}

//the next function is not used up to now
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingMortar::Evaluate()
{

  dserror("Evaluate method for mortar coupling not implemented yet!");

}

/*----------------------------------------------------------------------*
 * Compute sv = D^{-1}M(mv)
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> ADAPTER::CouplingMortar::MasterToSlave(RefCountPtr<
    Epetra_Vector> mv)
{
  
  dsassert( masterdofmap_->SameAs( mv->Map() ),
      "Vector with master dof map expected" );

  Epetra_Vector tmp = Epetra_Vector(M_->RowMap());

  if (M_->Multiply(false, *mv, tmp))
    dserror( "M*mv multiplication failed" );

    RefCountPtr<Epetra_Vector> sv = rcp( new Epetra_Vector( *slavedofmap_ ) );

    if ( Dinv_->Multiply( false, tmp, *sv ) )
    dserror( "D^{-1}*v multiplication failed" );

    return sv;
  }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> ADAPTER::CouplingMortar::SlaveToMaster(RefCountPtr<
    Epetra_Vector> sv)
{
  Epetra_Vector tmp = Epetra_Vector(M_->RangeMap());
  copy(sv->Values(), sv->Values() + sv->MyLength(), tmp.Values());

  RefCountPtr<Epetra_Vector> mv = rcp(new Epetra_Vector(*masterdofmap_));
  if (M_->Multiply(true, tmp, *mv))
    dserror( "M^{T}*sv multiplication failed" );

    return mv;
  }

#endif // CCADISCRET
