/*----------------------------------------------------------------------*/
/*!
\file stru_ale_algorithm.cpp

\brief  Basis of all structure approaches with ale
        (Lagrangian step followed by Eulerian step )
<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*
 | headers                                                   mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "stru_ale_algorithm.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"

#include "../drt_mortar/mortar_manager_base.H"

#include "../drt_w1/wall1.H"

#include "../drt_so3/so_hex8.H"

#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_wear_lagrange_strategy.H"
#include "../drt_contact/contact_wear_interface.H"

#include "../drt_structure/stru_aux.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialComm.h"

#include "../linalg/linalg_utils.H"

#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_ale/ale.H"

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
 | time loop                                                farah 10/13 |
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
    /* START LAGRANGE STEP                                              */
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
    Teuchos::RCP<Epetra_Vector> idisale_s, idisale_m, idisale_global;
    InterfaceDisp(idisale_s, idisale_m);

    // system of equation
    AleField().BuildSystemMatrix();

    // merge the both wear vectors for master and slave side to one global vector
    MergeWear(idisale_s,idisale_m,idisale_global);

    // coupling of struct/mortar and ale dofs
    DispCoupling(idisale_global);

    // application of interface displacements as dirichlet conditions
    AleField().ApplyInterfaceDisplacements(idisale_global);

    // solution
    AleField().Solve();

    // 2.-----------------------------------------------------------------     
    // application of mesh displacements to structural field, 
    // update spatial and material displacements
    ApplyMeshDisplacement();

    // update
    AleField().Update();

    // output
    AleField().Output();

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
 | Merge wear from slave and master surface to one           farah 06/13|
 | wear vector                                                          |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::MergeWear(Teuchos::RCP<Epetra_Vector>& disinterface_s ,
                                    Teuchos::RCP<Epetra_Vector>& disinterface_m,
                                    Teuchos::RCP<Epetra_Vector>& disinterface_g)
{
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::CoAbstractStrategy& cstrategy = static_cast<CONTACT::CoAbstractStrategy&>(strategy);
  std::vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();
  Teuchos::RCP<CONTACT::WearInterface> winterface = Teuchos::rcp_dynamic_cast<CONTACT::WearInterface>(interface[0]);
  if (winterface==Teuchos::null) dserror("Casting to WearInterface returned null!");

  disinterface_g = Teuchos::rcp(new Epetra_Vector(*winterface->Discret().DofRowMap()),true);
  Teuchos::RCP<Epetra_Vector>  auxvector = Teuchos::rcp(new Epetra_Vector(*winterface->Discret().DofRowMap()),true);

  LINALG::Export(*disinterface_s,*disinterface_g);
  LINALG::Export(*disinterface_m,*auxvector);

  disinterface_g->Update(1.0,*auxvector,true);

  return;
}

/*----------------------------------------------------------------------*
 | Vector of interface displacements in ALE dofs             farah 05/13|
 | Currently just for 1 interface                                       |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::InterfaceDisp(Teuchos::RCP<Epetra_Vector>& disinterface_s,
                                        Teuchos::RCP<Epetra_Vector>& disinterface_m)
{
  // FIXGIT: From global slave vector
  // FIXGIT: Perhaps master nodes
  
  // get vector of unweighted wear
  Teuchos::RCP<Epetra_Vector> realwear = cmtman_->GetStrategy().ContactWear();
  if(realwear==Teuchos::null) dserror("STRU_ALE::Algorithm::InterfaceDisp: realwear = Teuchos::null");

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::WearLagrangeStrategy& cstrategy = static_cast<CONTACT::WearLagrangeStrategy&>(strategy);


  // get vector of contact interfaces TODO: wear interfaces!!!!
  std::vector<RCP<CONTACT::CoInterface> > interface = cstrategy.ContactInterfaces();

  //***********************************************************************
  // we compute the real non-weighted wear vector here
  // 1. store delta_wear to wear
  // 2. compute the non-weighted wear vector
  // FIX: This should be done in a different way
  //***********************************************************************
  INPAR::CONTACT::WearType wtype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(DRT::Problem::Instance()->WearParams(),"WEARTYPE");

  if(wtype==INPAR::CONTACT::wear_impl)
    cstrategy.OutputWear();
  else
  {
    cstrategy.StoreNodalQuantities(MORTAR::StrategyBase::wear);
    cstrategy.OutputWear();
  }


  // dimension of the problem
  int dim = strategy.Dim();

  // this currently works only for one interface yet
  if (interface.size()>1)
  {
    std::cout << "*** Warning: more than 1 contact interface with wear ***" << std::endl;
    dserror("Error in TSI::Algorithm::ConvertMaps: Only for one interface yet.");
  }

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    Teuchos::RCP<CONTACT::WearInterface> winterface = Teuchos::rcp_dynamic_cast<CONTACT::WearInterface>(interface[m]);
    if (winterface==Teuchos::null) dserror("Casting to WearInterface returned null!");

    // get slave row dofs as map
    Teuchos::RCP<Epetra_Map> slavedofs = winterface->SlaveRowDofs();
    // additional spatial displacements
    // FIXGIT: transformation of results

    disinterface_s = Teuchos::rcp(new Epetra_Vector(*slavedofs,true));

    for (int i=0; i<disinterface_s->MyLength(); ++i)
    {
      if (i%dim == 0)
      {
        // get contact node for the current dofs
        int gid = winterface->SlaveRowNodes()->GID(i/dim);
        DRT::Node* node = winterface->Discret().gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(node);

        if (cstrategy.ContactWear() == Teuchos::null) dserror("unweighted wear vector returned null-pointer");

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
          (*disinterface_s)[i] = -(*cstrategy.ContactWear())[i];
          (*disinterface_s)[i+1] = -(*cstrategy.ContactWear())[i+1];

          if (dim==3)
          {
            (*disinterface_s)[i+2] = -(*cstrategy.ContactWear())[i+2];
          }
          i=i+dim-1;
        }
      }
    }
  } // interface loop


  // *********************************************************************************************
  // and now the same fun for the master surface
  // BUT: in which direction wear should be added? normal to master surface?
  // normal to slave surface? Think about that!
  // *********************************************************************************************
//   loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    Teuchos::RCP<CONTACT::WearInterface> winterface = Teuchos::rcp_dynamic_cast<CONTACT::WearInterface>(interface[m]);
    if (winterface==Teuchos::null) dserror("Casting to WearInterface returned null!");

    // get slave row dofs as map
    Teuchos::RCP<Epetra_Map> masterdofs = winterface->MasterRowDofs();

    disinterface_m = Teuchos::rcp(new Epetra_Vector(*masterdofs,true));

    for (int i=0; i<disinterface_m->MyLength(); ++i)
    {
      if (i%dim == 0)
      {
        (*disinterface_m)[i] = -(*cstrategy.ContactWear2())[i];
        (*disinterface_m)[i+1] = -(*cstrategy.ContactWear2())[i+1];

        if (dim==3)
        {
          (*disinterface_m)[i+2] = -(*cstrategy.ContactWear2())[i+2];
        }
        i=i+dim-1;
      }
    }
  } // interface loop for master side
  return;
}  // STRU_ALE::Algorithm::InterfaceDisp()

/*----------------------------------------------------------------------*
 | Application of mesh displacement                           mgit 07/11|
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::ApplyMeshDisplacement()
{
  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  RCP<const Epetra_Vector> aledispnp = AleField().Dispnp();
  RCP<const Epetra_Vector> aledispn  = AleField().Dispn();

  RCP<Epetra_Vector> disale   = coupalestru_->MasterToSlave(aledispnp);
  RCP<Epetra_Vector> disalen  = coupalestru_->MasterToSlave(aledispn);

  // build increment between aledispn and aledispnp
  disale->Update(-1.0,*disalen,1.0);

  // vector of current spatial displacements
  RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();  // change to ExtractDispn() for overlap

  // material displacements
  RCP<Epetra_Vector> dismat = Teuchos::rcp(new Epetra_Vector(dispnp->Map()),true);
  
  // set state
  (StructureField()->Discretization())->SetState(0,"displacement",dispnp);

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
      locid=(dispnp->Map()).LID(2*gid);
      if (locid==-1) dserror("LID not found on this proc");
    }
    else
    {
      locid=(dispnp->Map()).LID(3*gid);
      if (locid==-1) dserror("LID not found on this proc");
    }
    // reference node position + displacement t_n + delta displacement t_n+1
    XMesh[0]=node->X()[0]+(*dispnp)[locid]+(*disale)[locid];
    XMesh[1]=node->X()[1]+(*dispnp)[locid+1]+(*disale)[locid+1];

    if (ndim==3)
      XMesh[2]=node->X()[2]+(*dispnp)[locid+2]+(*disale)[locid+2];

    // create updated  XMat --> via nonlinear interpolation between nodes (like gp projection)
    AdvectionMap(&XMat[0],&XMat[1],&XMat[2],&XMesh[0],&XMesh[1],&XMesh[2],ElementPtr,numelement);

    // create delta displacement in material configuration
    (*dismat)[locid] = XMat[0]-node->X()[0];
    (*dismat)[locid+1] = XMat[1]-node->X()[1];

    if (ndim==3)
      (*dismat)[locid+2] = XMat[2]-node->X()[2];
  }

  // update spatial displacements
  dispnp->Update(1.0,*disale,1.0);

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
      //std::cout << "coords: " << *XMat1 << "   " << *XMat2 << "   " << *XMat3 << std::endl;
      //std::cout << "*** particle tracking successfull ***" << std::endl;
      return;
    }
  }

  // error if element is not found
  if (found == false)
  {
    //cout << "coords: " << *XMat1 << "   " << *XMat2 << "   " << *XMat3 << endl;
    //cout << "STRU_ALE::Algorithm::AdvectionMap: Particle tracking not successful." << endl;
  }
  // bye
  return;

}  // STUE_ALE::Algorithm::AdvectionMap()

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)    farah 10/13 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  AleField().ReadRestart(step);
  SetTimeStep(StructureField()->GetTime(),step);

  return;
}
/*----------------------------------------------------------------------*/
