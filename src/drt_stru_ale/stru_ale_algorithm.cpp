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

#include "Epetra_SerialComm.h"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "../drt_adapter/adapter_coupling.H"
/*----------------------------------------------------------------------*
 | constructor (public)                                     farah 05/13 |
 *----------------------------------------------------------------------*/
STRU_ALE::Algorithm::Algorithm(const Epetra_Comm& comm)
 : AlgorithmBase(comm,DRT::Problem::Instance()->StructuralDynamicParams())
{
  /*--------------------------------------------------------------------*
   | first create structure then ale --> important for discretization   |
   | numbering and therefore for the post_drt_ensight.cpp               |
   *--------------------------------------------------------------------*/

  // create structure
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams(), const_cast<Teuchos::ParameterList&>(DRT::Problem::Instance()->StructuralDynamicParams()), DRT::Problem::Instance()->GetDis("structure")));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  // create ale
  Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams()));
  ale_ = ale->AleFieldrcp();

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // contact/meshtying manager
  cmtman_= StructureField()->ContactManager();

  const int ndim = DRT::Problem::Instance()->NDim();

  // create ale-struct coupling
  const Epetra_Map* structdofmap = StructureField()->Discretization()->NodeRowMap();
  const Epetra_Map* aledofmap   = AleField().Discretization()->NodeRowMap();

  // if there are two identical nodes (i.e. for initial contact) the nodes matching creates an error !!!
  coupalestru_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupalestru_->SetupCoupling(*AleField().Discretization(),
                       *StructureField()->Discretization(),
                       *aledofmap,
                       *structdofmap,
                        ndim);

  //create interface coupling
  coupstrualei_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupstrualei_->SetupConditionCoupling(*StructureField()->Discretization(),
                                 StructureField()->Interface()->AleWearCondMap(),
                                 *AleField().Discretization(),
                                 AleField().Interface()->Map(AleField().Interface()->cond_ale_wear),
                                 "AleWear",
                                 ndim);
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
    StructureField()->PrepareTimeStep();
    StructureField()->Solve();

    // calculate stresses, strains, energies
    StructureField()->PrepareOutput();

    // update at time step
    StructureField()->Update();

    // write output to screen and files
    StructureField()->Output();
    
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

    // couplig of struct/mortar and ale dofs
    DispCoupling(idisale);

    // application of interface displacements as dirichlet conditions
    AleField().ApplyInterfaceDisplacements(idisale);
    
    // solution
    AleField().Solve();

    // output
    //FIX: output creates error by writing in resultfile
    //AleField().Output();

    // 2.-----------------------------------------------------------------     
    // application of mesh displacements to structural field, 
    ApplyMeshDisplacement();

  }  // time loop
}  // STRUE_ALE::Algorithm::TimeLoop()


/*----------------------------------------------------------------------*
 | Perform Coupling from struct/mortar to ale dofs          farah 05/13 |
 | This is necessary due to the parallel redistribution                 |
 | of the contact interface                                             |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::DispCoupling(Teuchos::RCP<Epetra_Vector>& disinterface)
{
  //Teuchos::RCP<Epetra_Vector> aledofs = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->Map(AleField().Interface()->cond_ale_wear)),true);
  Teuchos::RCP<Epetra_Vector> strudofs = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->AleWearCondMap()),true);

  // change the parallel distribution from mortar interface to structure
  LINALG::Export(*disinterface,*strudofs);

  // perform coupling
  disinterface.reset();
  disinterface = coupstrualei_->MasterToSlave(strudofs);

  return;
}

/*----------------------------------------------------------------------*
 | Vector of interface displacements in ALE dofs             farah 05/13|
 | Currently just for 1 interface                                       |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::InterfaceDisp(Teuchos::RCP<Epetra_Vector>& disinterface)
{
  // FIXGIT: From global slave vector
  // FIXGIT: Perhaps master nodes
  
  // get vector of unweighted wear
  RCP<Epetra_Vector> realwear = cmtman_->GetStrategy().ContactWear();
  
  // map of slave dofs
  //RCP<Epetra_Map> slavedofs;

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  // this currently works only for one interface yet
  if (interface.size()>1)
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");

  // get slave row dofs as map
  Teuchos::RCP<Epetra_Map> slavedofs = interface[0]->SlaveRowDofs();

  // dimension of the problem
  //int dim = strategy.Dim();

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // do nothing
  }
  
  // additional spatial displacements
  // FIXGIT: check has to be done with nodal normal
  // FIXGIT: transformation of results
  // FIX: is here the wear vector multiplied with the nodal normal or should it be done?
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
void STRU_ALE::Algorithm::ApplyMeshDisplacement()
{
  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  RCP<Epetra_Vector> aledisp = AleField().ExtractDispnp();
  RCP<Epetra_Vector> disale = coupalestru_->MasterToSlave(aledisp);

  // vector of current spatial displacements
  RCP<Epetra_Vector> dispn = StructureField()->ExtractDispn();

  // additional spatial displacements
  //RCP<Epetra_Vector> disadditional = Teuchos::rcp(new Epetra_Vector(dispn->Map()),true);

  // material displacements
  RCP<Epetra_Vector> dismat = Teuchos::rcp(new Epetra_Vector(dispn->Map()),true);
  
  // set state
  (StructureField()->Discretization())->SetState(0,"displacement",dispn);

  // set state
  (StructureField()->Discretization())->SetState(0,"material displacement",StructureField()->DispMat());
 
  // loop over all row nodes to fill graph
  for (int k=0;k<StructureField()->Discretization()->NumMyRowNodes();++k)
  {
    int gid = StructureField()->Discretization()->NodeRowMap()->GID(k);

    DRT::Node* node = StructureField()->Discretization()->gNode(gid);
    DRT::Element** ElementPtr = node->Elements();
      
    int numelement = node->NumElement();
      
    double XMat[2];
    double XMesh[2];

    XMat[0] = node->X()[0];
    XMat[1] = node->X()[1];
      
    int locid = (dispn->Map()).LID(2*gid);

    // reference node position + displacement t_n + delta displacement t_n+1
    XMesh[0]=node->X()[0]+(*dispn)[locid]+(*disale)[locid];
    XMesh[1]=node->X()[1]+(*dispn)[locid+1]+(*disale)[locid+1];

    // create updated  XMat --> via nonlinear interpolation between nodes (like gp projection)
    AdvectionMap(&XMat[0],&XMat[1],&XMesh[0],&XMesh[1],ElementPtr,numelement);

    // copy disale? so we copy just the wear? i think this is nonsense
    //(*disadditional)[locid]   = (*disale)[locid];
    //(*disadditional)[locid+1] = (*disale)[locid+1];

    // create delta displacement in material configuration
    (*dismat)[locid] = XMat[0]-node->X()[0];
    (*dismat)[locid+1] = XMat[1]-node->X()[1];
  }

  // update spatial displacements
  //dispn->Update(1.0,*disadditional,1.0);
  dispn->Update(1.0,*disale,1.0);

  //StructureField()->Dispnp() = dispn;

  // apply material displacements to structural field
  StructureField()->ApplyDisMat(dismat);
  //cout << "dismat= " << *dismat << endl;

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
    actele->LocationVector(*(StructureField()->Discretization()),la,false);

    // get state
    RCP<const Epetra_Vector> disp = (StructureField()->Discretization())->GetState("displacement");
    RCP<const Epetra_Vector> dispmat = (StructureField()->Discretization())->GetState("material displacement");

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
    {
      //cout << "*** particle tracking successfull ***" << endl;
      return;
    }
  }

  // error if element is not found
  if (found == false)
  {
    cout << "coords: " << *XMat1 << "   " << *XMat2 << endl;
    cout << "STRU_ALE::Algorithm::AdvectionMap: Particle tracking not successful." << endl;
  }
  // bye
  return;

}  // STUE_ALE::Algorithm::AdvectionMap()

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     mgit 05/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  SetTimeStep(StructureField()->GetTime(),step);

  return;
}
/*----------------------------------------------------------------------*/
