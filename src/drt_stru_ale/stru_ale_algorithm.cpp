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
#include "../drt_so3/so_hex8.H"
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
    /* START LAGRANGE STEP                                                    */
    /* structural lagrange step with contact                            */
    /********************************************************************/
    // predict and solve structural system
    StructureField()->PrepareTimeStep();
    StructureField()->Solve();
    
    /********************************************************************/
    /* EULERIAN STEP                                                    */  
    /* 1. mesh displacements due to wear from ALE system                */
    /* 2. mapping of results from "old" to "new" mesh                   */
    /********************************************************************/
    
    // 1.----------------------------------------------------------------- 
    // prepare ale output: increase time step
    AleField().PrepareTimeStep();

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
    AleField().Output();

    // 2.-----------------------------------------------------------------     
    // application of mesh displacements to structural field, 
    // update spatial and material displacements
    ApplyMeshDisplacement();

    /********************************************************************/
    /* FINISH LAGRANGE STEP                                             */
    /* Store spatial and material displacements for restart             */
    /* Write Output                                                     */
    /********************************************************************/
    // calculate stresses, strains, energies
    StructureField()->PrepareOutput();

    // update at time step
    StructureField()->Update();

    // write output to screen and files
    StructureField()->Output();

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

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);

  // get vector of contact interfaces
  std::vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  //***********************************************************************
  // we compute the real non-weighted wear vector here
  // 1. store delta_wear to wear
  // 2. compute the non-weighted wear vector
  // FIX: This should be done in a different way
  //***********************************************************************
  cstrategy.StoreNodalQuantities(MORTAR::StrategyBase::wear);
  cstrategy.OutputWear();

  // dimension of the problem
  int dim = strategy.Dim();

  // this currently works only for one interface yet
  if (interface.size()>1)
  {
    cout << "*** Warning: more than 1 contact interface with wear ***" << endl;
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");
  }

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    // get slave row dofs as map
    Teuchos::RCP<Epetra_Map> slavedofs = interface[0]->SlaveRowDofs();

    // additional spatial displacements
    // FIXGIT: transformation of results

    disinterface = Teuchos::rcp(new Epetra_Vector(*slavedofs,true));

    for (int i=0; i<disinterface->MyLength(); ++i)
    {
      if (i%dim == 0)
      {
        // get contact node for the current dofs
        int gid = interface[0]->SlaveRowNodes()->GID(i/dim);
        DRT::Node* node = interface[0]->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

        //calc (w*n)*n --> if <0 then negative wear !!!
        double xp2=(*cstrategy.ContactWear())[i] * cnode->MoData().n()[0];
        double yp2=(*cstrategy.ContactWear())[i+1] * cnode->MoData().n()[1];
        double zp2=0.0;

        if (dim==3)
        {
          zp2=(*cstrategy.ContactWear())[i+2] * cnode->MoData().n()[2];
        }

        double norm = xp2+yp2+zp2;

        if (norm<0.0)
        {
          i=i+dim-1;
        }
        else
        {
          (*disinterface)[i] = -(*cstrategy.ContactWear())[i];
          (*disinterface)[i+1] = -(*cstrategy.ContactWear())[i+1];

          if (dim==3)
          {
            (*disinterface)[i+2] = -(*cstrategy.ContactWear())[i+2];
          }
          i=i+dim-1;
        }
      }
    }
  } // interface loop

  return;
}  // STRU_ALE::Algorithm::InterfaceDisp()

/*----------------------------------------------------------------------*
 | Application of mesh displacement                           mgit 07/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::ApplyMeshDisplacement()
{
  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  RCP<Epetra_Vector> aledisp = AleField().ExtractDispnp();
  RCP<Epetra_Vector> disale = coupalestru_->MasterToSlave(aledisp);

  // vector of current spatial displacements
  RCP<Epetra_Vector> dispn = StructureField()->ExtractDispn();

  // material displacements
  RCP<Epetra_Vector> dismat = Teuchos::rcp(new Epetra_Vector(dispn->Map()),true);
  
  // set state
  (StructureField()->Discretization())->SetState(0,"displacement",dispn);

  // set state
  (StructureField()->Discretization())->SetState(0,"material_displacement",StructureField()->DispMat());
 
  // loop over all row nodes to fill graph
  for (int k=0;k<StructureField()->Discretization()->NumMyRowNodes();++k)
  {
    int gid = StructureField()->Discretization()->NodeRowMap()->GID(k);

    DRT::Node* node = StructureField()->Discretization()->gNode(gid);
    DRT::Element** ElementPtr = node->Elements();
      
    int numelement = node->NumElement();
      
    // create Xmat for 3D problems
    double XMat[3];
    double XMesh[3];

    // local ID
    int locid = 0;

    XMat[0] = node->X()[0];
    XMat[1] = node->X()[1];

    if (ndim==3)
    {
      XMat[2] = node->X()[2];
    }

    // get local id
    if (ndim==2)
    {
      locid=(dispn->Map()).LID(2*gid);
      if (locid==-1) dserror("LID not found on this proc");
    }
    else
    {
      locid=(dispn->Map()).LID(3*gid);
      if (locid==-1) dserror("LID not found on this proc");
    }
    // reference node position + displacement t_n + delta displacement t_n+1
    XMesh[0]=node->X()[0]+(*dispn)[locid]+(*disale)[locid];
    XMesh[1]=node->X()[1]+(*dispn)[locid+1]+(*disale)[locid+1];

    if (ndim==3)
      XMesh[2]=node->X()[2]+(*dispn)[locid+2]+(*disale)[locid+2];

    // create updated  XMat --> via nonlinear interpolation between nodes (like gp projection)
    //if (ndim==2)
    AdvectionMap(&XMat[0],&XMat[1],&XMat[2],&XMesh[0],&XMesh[1],&XMesh[2],ElementPtr,numelement);
    // create delta displacement in material configuration
    (*dismat)[locid] = XMat[0]-node->X()[0];
    (*dismat)[locid+1] = XMat[1]-node->X()[1];

    if (ndim==3)
      (*dismat)[locid+2] = XMat[2]-node->X()[2];
  }

  // update spatial displacements
  dispn->Update(1.0,*disale,1.0);

  // apply material displacements to structural field
  // if advection map is not succesful --> use old xmat
  StructureField()->ApplyDisMat(dismat);

  return;

}  // STRU_ALE::Algorithm::SubtractWear()

/*----------------------------------------------------------------------*
 | material coordinates evaluated from spatial ones          mgit 05/11 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::AdvectionMap(double* XMat1,
                                       double* XMat2,
                                       double* XMat3,
                                       double* XMesh1,
                                       double* XMesh2,
                                       double* XMesh3,
                                       DRT::Element** ElementPtr,
                                       int numelements)
{
  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // found element the spatial coordinate lies in
  bool found = false;

  // loop over adjacent elements
  for(int jele = 0; jele < numelements; jele++)
  {
    // get element
    DRT::Element* actele = ElementPtr[jele];

    // get element location vector, dirichlet flags and ownerships
    // get element location vector
    DRT::Element::LocationArray la(1);
    actele->LocationVector(*(StructureField()->Discretization()),la,false);

    // get state
    RCP<const Epetra_Vector> disp = (StructureField()->Discretization())->GetState("displacement");
    RCP<const Epetra_Vector> dispmat = (StructureField()->Discretization())->GetState("material_displacement");

    if (ndim == 2)
    {
      // structure with ale only for wall element so far
      if (actele->ElementType() != DRT::ELEMENTS::Wall1Type::Instance())
        dserror ("Construction of advection map for 2D problems only for wall element so far.");

      // cast element to wall element
      DRT::ELEMENTS::Wall1* w1ele = static_cast<DRT::ELEMENTS::Wall1*>(actele);

      // checks if the spatial coordinate lies within this element
      // if yes, returns the material displacements
      w1ele->AdvectionMapElement(XMat1,XMat2,XMesh1,XMesh2,disp,dispmat, la,found);
    }
    else
    {
      // structure with ale only for wall element so far
      if (actele->ElementType() != DRT::ELEMENTS::So_hex8Type::Instance())
        dserror ("Construction of advection map for 3D problems only for solid hex8 element so far.");

      // cast element to solid hex8 element
      DRT::ELEMENTS::So_hex8* ele = static_cast<DRT::ELEMENTS::So_hex8*>(actele);

      // checks if the spatial coordinate lies within this element
      // if yes, returns the material displacements
      ele->AdvectionMapElement(XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la,found);
    }

    // leave when element is found
    if (found == true)
    {
      //cout << "coords: " << *XMat1 << "   " << *XMat2 << "   " << *XMat3 << endl;
      //cout << "*** particle tracking successfull ***" << endl;
      return;
    }
  }

  // error if element is not found
  if (found == false)
  {
    cout << "coords: " << *XMat1 << "   " << *XMat2 << "   " << *XMat3 << endl;
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
  AleField().ReadRestart(step);
  SetTimeStep(StructureField()->GetTime(),step);

  return;
}
/*----------------------------------------------------------------------*/
