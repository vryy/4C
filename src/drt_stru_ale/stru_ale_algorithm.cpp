/*----------------------------------------------------------------------*/
/*!
\file stru_ale_algorithm.cpp

\brief  Basis of all structure approaches with ale
        (Lagrangian step followed by Eulerian step )
<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               mgit 04/11 |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                   mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "stru_ale_algorithm.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_w1/wall1.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_ale/ale.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*
 | constructor (public)                                      mgit 05/11 |
 *----------------------------------------------------------------------*/
STRU_ALE::Algorithm::Algorithm(const Epetra_Comm& comm)
 : AlgorithmBase(comm,DRT::Problem::Instance()->StructuralDynamicParams()),
   StructureBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams(), DRT::Problem::Instance()->GetDis("structure")),
   AleBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams())
{
  // contact/meshtying manager
  cmtman_= StructureField().ContactManager();
}

/*----------------------------------------------------------------------*
 | destructor (public)                                       mgit 05/11 |
 *----------------------------------------------------------------------*/
STRU_ALE::Algorithm::~Algorithm()
{
}

/*----------------------------------------------------------------------*
 | time loop                                                 mgit 05/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::TimeLoop()
{

  // time loop
  while (NotFinished())
  {
    // counter and print header
    IncrementTimeAndStep();
    PrintHeader();

    /********************************************************************/
    /* LAGRANGE STEP                                                    */  
    /* structural lagrange step with contact                            */
    /********************************************************************/
    // predict and solve structural system
    StructureField().PrepareTimeStep();
    StructureField().Solve();

    // calculate stresses, strains, energies
    StructureField().PrepareOutput();

    // write output to screen and files
    StructureField().Update();

    // write output to screen and files
    StructureField().Output();
    
    /********************************************************************/
    /* EULERIAN STEP                                                    */  
    /* 1. mesh displacements due to wear from ALE system                */
    /* 2. mapping of results from "old" to "new" mesh                   */
    /********************************************************************/
    
    // 1.----------------------------------------------------------------- 
    // wear as interface displacements in ale dofs 
    RCP<Epetra_Vector> idisale;
    InterfaceDisp(idisale);
    
    // system of equation
    AleField().BuildSystemMatrix();

    // application of interface displacements as dirichlet conditions
    AleField().ApplyInterfaceDisplacements(idisale);
    
    // solution
    AleField().Solve();

    // 2.-----------------------------------------------------------------     
    // mesh displacement from solution of ALE field in structural dofs
    // FIXGIT: Has to be done with transformation of vector
    RCP<Epetra_Vector> idis = Teuchos::rcp(new Epetra_Vector(*(StructureField().Discretization()->DofRowMap()),true));
    for (int i=0; i<idis->MyLength(); ++i)
      (*idis)[i]=(*(AleField().ExtractDisplacement()))[i];
    
    // application of mesh displacements to structural field, 
    // mapping of results
    ApplyMeshDisplacement(idis);

  }  // time loop
}  // STRUE_ALE::Algorithm::TimeLoop()

/*----------------------------------------------------------------------*
 | Vector of interface displacements in ALE dofs              mgit 07/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::InterfaceDisp(Teuchos::RCP<Epetra_Vector>& disinterface)
{
  // FIXGIT: From global slave vector
  // FIXGIT: Perhaps master nodes
  
  // get vector of unweighted wear
  RCP<Epetra_Vector> realwear = cmtman_->GetStrategy().ContactWear();
  
  // map of slave dofs
  RCP<Epetra_Map> slavedofs;
   
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // dimension of the problem
  int dim = strategy.Dim();

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // slave nodes
    const RCP<Epetra_Map> slavenodes = interface[m]->SlaveRowNodes();

    // define local variables
    int slavecountnodes = 0;
    vector<int> myslavealedofs((slavenodes->NumMyElements())*dim);

    // loop over all slave nodes of the interface
    for (int i=0;i<slavenodes->NumMyElements();++i)
    {
      int gid = slavenodes->GID(i);
      DRT::Node* node = StructureField().Discretization()->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: ConvertMaps: Node ownership inconsistency!");

      // ale dofs
      if (dim == 2)
      {  
        myslavealedofs[slavecountnodes*dim] = (AleField().Discretization()->Dof(0,node))[0];
        myslavealedofs[slavecountnodes*dim+1] = (AleField().Discretization()->Dof(0,node))[1];
      }
      else
        dserror("ERROR: 3D not yet implemented.");
      
      ++slavecountnodes;
    }  

    // resize the temporary vectors
    myslavealedofs.resize(slavecountnodes*dim);

    // communicate slavecuntnodes
    int gslavecountnodes;
    Comm().SumAll(&slavecountnodes,&gslavecountnodes,1);

    // create slave node map and active dof map
    slavedofs = Teuchos::rcp(new Epetra_Map(gslavecountnodes*dim,slavecountnodes*dim,&myslavealedofs[0],0,Comm()));
  }
  
  // additional spatial displacements
  // FIXGIT: check has to be done with nodal normal
  // FIXGIT: transformation of results
  disinterface = Teuchos::rcp(new Epetra_Vector(*slavedofs,true));
  for (int i=0; i<disinterface->MyLength(); ++i)
  {
    if (i%2 > 0 and (*cstrategy.ContactWear())[i] > 0.0)
      (*disinterface)[i] = 0;
    else
      (*disinterface)[i] = -(*cstrategy.ContactWear())[i]; 
   }
  
  return;
}  // STRU_ALE::Algorithm::InterfaceDisp()

/*----------------------------------------------------------------------*
 | Application of mesh displacement                           mgit 07/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector>& disale)
{
  // vector of current spatial displacements
  RCP<Epetra_Vector> dispn = StructureField().ExtractDispn();

  // additional spatial displacements
  RCP<Epetra_Vector> disadditional = Teuchos::rcp(new Epetra_Vector(dispn->Map()),true);

  // material displacements
  RCP<Epetra_Vector> dismat = Teuchos::rcp(new Epetra_Vector(dispn->Map()),true);
  
  // set state
  (StructureField().Discretization())->SetState(0,"displacement",dispn);

  // set state
  (StructureField().Discretization())->SetState(0,"material displacement",StructureField().DispMat());
 
  // loop over all row nodes to fill graph
    for (int k=0;k<StructureField().Discretization()->NumMyRowNodes();++k)
    {
      int gid = StructureField().Discretization()->NodeRowMap()->GID(k);
      
      DRT::Node* node = StructureField().Discretization()->gNode(gid);
      DRT::Element** ElementPtr = node->Elements();
        
      int numelement = node->NumElement();
        
      double XMat[2];
      double XMesh[2];
      
      XMat[0] = node->X()[0];
      XMat[1] = node->X()[1];
        
      int locid = (dispn->Map()).LID(2*gid);
  
      XMesh[0]=node->X()[0]+(*dispn)[locid]+(*disale)[locid];
      XMesh[1]=node->X()[1]+(*dispn)[locid+1]+(*disale)[locid+1];
      
      AdvectionMap(&XMat[0],&XMat[1],&XMesh[0],&XMesh[1],ElementPtr,numelement);

      (*disadditional)[locid]   = (*disale)[locid];
      (*disadditional)[locid+1] = (*disale)[locid+1];

      (*dismat)[locid] = XMat[0]-node->X()[0];
      (*dismat)[locid+1] = XMat[1]-node->X()[1];
    }

    // update spatial displacements
    dispn->Update(1.0,*disadditional,1.0);

    // apply material displacements to structural field
    StructureField().ApplyDisMat(dismat);
    
    return;

}  // STRU_ALE::Algorithm::SubtractWear()

/*----------------------------------------------------------------------*
 | material coordinates evaluated from spatial ones          mgit 05/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::AdvectionMap(double* XMat1,
                                       double* XMat2,
                                       double* XMesh1,
                                       double* XMesh2,
                                       DRT::Element** ElementPtr,
                                       int numelements)
{
  // found element the spatial coordinate lies in
  bool found = false;

  // loop over adjacent elements
  for(int jele = 0; jele < numelements; jele++)
  {
    // get element
    DRT::Element* actele = ElementPtr[jele];

    // get element location vector, dirichlet flags and ownerships
    // get element location vector
    // FIXGIT: Why this "1"
    DRT::Element::LocationArray la(1);
    actele->LocationVector(*(StructureField().Discretization()),la,false);

    // get state
    RCP<const Epetra_Vector> disp = (StructureField().Discretization())->GetState("displacement");
    RCP<const Epetra_Vector> dispmat = (StructureField().Discretization())->GetState("material displacement");

//#ifdef D_WALL1
    // structure with ale only for wall element so far
    if (actele->ElementType() != DRT::ELEMENTS::Wall1Type::Instance())
      dserror ("Construction of advection map only for wall element so far.");

    // cast element to wall element
    DRT::ELEMENTS::Wall1* w1ele = static_cast<DRT::ELEMENTS::Wall1*>(actele);

    // checks if the spatial coordinate lies within this element
    // if yes, returns the material displacements
    w1ele->AdvectionMapElement(XMat1,XMat2,XMesh1,XMesh2,disp,dispmat,la,found);
//#endif
    // leave when element is found
    if (found == true)
      return;
  }

  // error if element is not found
  if (found == false)
    cout << "STRU_ALE::Algorithm::AdvectionMap: Particle tracking not successful." << endl;

  // bye
  return;

}  // STUE_ALE::Algorithm::AdvectionMap()

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     mgit 05/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::ReadRestart(int step)
{
  StructureField().ReadRestart(step);
  SetTimeStep(StructureField().GetTime(),step);

  return;
}
/*----------------------------------------------------------------------*/
