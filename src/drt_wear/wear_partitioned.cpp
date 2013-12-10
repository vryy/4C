/*----------------------------------------------------------------------*/
/*!
\file wear_partitioned.cpp

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
 | headers                                                  farah 11/13 |
 *----------------------------------------------------------------------*/
#include "wear_partitioned.H"
#include "wear_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"
#include "../drt_inpar/inpar_ale.H"

#include "../drt_w1/wall1.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_hex20.H"
#include "../drt_so3/so_hex27.H"
#include "../drt_so3/so_tet4.H"
#include "../drt_so3/so_tet10.H"

#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_wear_lagrange_strategy.H"
#include "../drt_contact/contact_wear_interface.H"
#include "../drt_mortar/mortar_manager_base.H"

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
WEAR::Partitioned::Partitioned(const Epetra_Comm& comm)
: Algorithm(comm)
{
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

  // initialize intern variables for wear
  wearnp_i_  = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->Map(AleField().Interface()->cond_ale_wear)),true);
  wearnp_ip_ = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->Map(AleField().Interface()->cond_ale_wear)),true);
  wearincr_  = Teuchos::rcp(new Epetra_Vector(*AleField().Interface()->Map(AleField().Interface()->cond_ale_wear)),true);
  delta_ale_= Teuchos::rcp(new Epetra_Vector(AleField().Dispnp()->Map(),true));

  alepara_   = DRT::Problem::Instance()->AleDynamicParams();
}

/*----------------------------------------------------------------------*
 | destructor (public)                                      farah 11/13 |
 *----------------------------------------------------------------------*/
WEAR::Partitioned::~Partitioned()
{

}

/*----------------------------------------------------------------------*
 | general time loop                                        farah 10/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::TimeLoop()
{
  // get wear paramter list
  const Teuchos::ParameterList& wearpara = DRT::Problem::Instance()->WearParams();
  int timeratio = wearpara.get<double>("WEAR_TIMERATIO");

  int counter = -1;
  bool alestep = false;

  // time loop
  while (NotFinished())
  {
    if ( (int)(Step()/timeratio) >counter )
    {
      counter++;
      alestep=true;
    }

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearCoupAlgo>(wearpara,"WEAR_COUPALGO")
        == INPAR::CONTACT::wear_stagg)
      TimeLoopStagg(alestep);
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearCoupAlgo>(wearpara,"WEAR_COUPALGO")
        == INPAR::CONTACT::wear_iterstagg)
      TimeLoopIterStagg();
    else
      dserror("WEAR::TimeLoop: Algorithm not provided!");

    alestep=false;
  }  // time loop
}

/*----------------------------------------------------------------------*
 | time loop for staggered coupling                         farah 11/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::TimeLoopIterStagg()
{
  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

  //prepare time step for both fields
  PrepareTimeStep();

  bool converged = false;  // converged state?
  bool iterated  = false;  // more than 1 iteration?
  int  iter      = 0;      // iteration counter

  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::WearLagrangeStrategy& cstrategy = static_cast<CONTACT::WearLagrangeStrategy&>(strategy);

  cstrategy.UpdateWearDiscretIterate(false);

  /*************************************************************
   * Nonlinear iterations between Structure and ALE:           *
   * 1. Solve structure + contact to get wear                  *
   * 2. Apply wear increment (i+1 - i) onto ALE (add function) *
   * 3. Employ ALE disp incr (i+1 - i) and spat disp i to get  *
   *    abs mat disp for timestep n+1                          *
   * 4. Upadate spat disp from i to i+1                        *
   * 5. Check for convergence                                  *
   * 6. store ALE disp i = i+1                                 *
   *************************************************************/
  while (converged==false)
  {
    if (iter>0) iterated=true;

    // 1. solution
    StructureField()->Solve();

    // 2. wear as interface displacements in ale dofs
    Teuchos::RCP<Epetra_Vector> idisale_s, idisale_m;
    InterfaceDisp(idisale_s, idisale_m);

    // merge the both wear vectors for master and slave side to one global vector
    MergeWear(idisale_s,idisale_m,wearincr_);

    // coupling of struct/mortar and ale dofs
    DispCoupling(wearincr_);

    // do ale step
    AleStep(wearincr_);

    // 3. application of mesh displacements to structural field,
    // update material displacements
    ApplyMeshDisplacement(iterated);

    // 4. update dispnp
    UpdateDispnp();

    // 5. convergence check fot current iteration
    converged=ConvergenceCheck(iter);

    // store old wear
    cstrategy.UpdateWearDiscretIterate(true);

    ++iter;
  }// end nonlin loop

  // update for structure and ale
  Update();

  // output for structure and ale
  Output();


  return;
}

/*----------------------------------------------------------------------*
 | time loop for oneway coupling                            farah 11/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::TimeLoopStagg(bool alestep)
{
  // stactic cast of mortar strategy to contact strategy
  MORTAR::StrategyBase& strategy = cmtman_->GetStrategy();
  CONTACT::WearLagrangeStrategy& cstrategy = static_cast<CONTACT::WearLagrangeStrategy&>(strategy);

  // counter and print header
  IncrementTimeAndStep();
  PrintHeader();

  //prepare time step for both fields
  PrepareTimeStep();

  /********************************************************************/
  /* START LAGRANGE STEP                                              */
  /* structural lagrange step with contact                            */
  /********************************************************************/

  // solution
  StructureField()->Solve();

  if(alestep)
  {
    if(Comm().MyPID()==0)
      std::cout << "========================= ALE STEP =========================" << std::endl;

    /********************************************************************/
    /* COUPLING                                                         */
    /* Wear from structure solve as dirichlet for ALE                   */
    /********************************************************************/

    // wear as interface displacements in ale dofs
    Teuchos::RCP<Epetra_Vector> idisale_s, idisale_m, idisale_global;
    InterfaceDisp(idisale_s, idisale_m);

    // merge the both wear vectors for master and slave side to one global vector
    MergeWear(idisale_s,idisale_m,idisale_global);

    // coupling of struct/mortar and ale dofs
    DispCoupling(idisale_global);

    /********************************************************************/
    /* EULERIAN STEP                                                    */
    /* 1. mesh displacements due to wear from ALE system                */
    /* 2. mapping of results from "old" to "new" mesh                   */
    /********************************************************************/

    // do all step
    AleStep(idisale_global);

    // application of mesh displacements to structural field,
    // update spatial and material displacements
    ApplyMeshDisplacement();

    /********************************************************************/
    /* FINISH STEP:                                                     */
    /* Update and Write Output                                          */
    /********************************************************************/

    // update dispnp
    UpdateDispnp();

    cstrategy.UpdateWearDiscretIterate(false);
  }
  else
  {
    cstrategy.UpdateWearDiscretAccumulation(true);
  }

  // update for structure and ale
  Update();

  // output for structure and ale
  Output();

  return;
}

/*----------------------------------------------------------------------*
 | prepare time step for ale and structure                  farah 11/13 |
 *----------------------------------------------------------------------*/
bool WEAR::Partitioned::ConvergenceCheck(int iter)
{
  double Wincr = 0.0;
  double ALEincr = 0.0;
  wearincr_->Norm2(&Wincr);
  delta_ale_->Norm2(&ALEincr);

  if (Comm().MyPID()==0)
  {
    std::cout << "-----------------"<< " Step " << iter+1 <<" --------------------" << std::endl;
    std::cout << "Wear incr.= " << Wincr  << "         ALE incr.= " << ALEincr<< std::endl;
    std::cout << "---------------------------------------------" << std::endl;
  }

  if (abs(Wincr)<1e-10 and abs(ALEincr)<1e-10)
    return true;

  if (iter>50)
    dserror("Staggered solution scheme for ale-wear problem unconverged within 50 nonlinear iteration steps!");

  return false;
}

/*----------------------------------------------------------------------*
 | prepare time step for ale and structure                  farah 11/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::PrepareTimeStep()
{
  // predict and solve structural system
  StructureField()->PrepareTimeStep();

  // prepare ale output: increase time step
  AleField().PrepareTimeStep();

  return;
}

/*----------------------------------------------------------------------*
 | update ale and structure                                 farah 11/13 |
 *---------------------------------------------------------------- ------*/
void WEAR::Partitioned::Update()
{
  // update at time step
  StructureField()->Update();

  // update
  AleField().Update();

  return;
}

/*----------------------------------------------------------------------*
 | update spatial displacements                             farah 11/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::UpdateDispnp()
{
  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  RCP<Epetra_Vector> disale   = AleToStructure(AleField().Dispnp());
  RCP<Epetra_Vector> disalen  = AleToStructure(AleField().Dispn());

  // get structure dispnp vector
  RCP<Epetra_Vector> dispnp = StructureField()->WriteAccessDispnp();  // change to ExtractDispn() for overlap

  int aletype = DRT::INPUT::IntegralValue<int>(ParamsAle(),"ALE_TYPE");

  // for incremental lin ale --> in spatial conf.
  if (aletype==INPAR::ALE::incr_lin)
  {
    // update per absolute vector
    dispnp->Update(1.0,*disale,0.0);
  }
  else if (aletype==INPAR::ALE::classic_lin)
  {
    // create increment between n and np
    disale->Update(-1.0,*disalen,1.0);

    // update per increment
    dispnp->Update(1.0,*disale,1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | output ale and structure                                 farah 11/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::Output()
{
  // calculate stresses, strains, energies
  StructureField()->PrepareOutput();

  // write strcture output to screen and files
  StructureField()->Output();

  // output ale
  AleField().Output();

  return;
}

/*----------------------------------------------------------------------*
 | Perform Coupling from struct/mortar to ale dofs          farah 05/13 |
 | This is necessary due to the parallel redistribution                 |
 | of the contact interface                                             |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::DispCoupling(Teuchos::RCP<Epetra_Vector>& disinterface)
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
 | Merge wear from slave and master surface to one          farah 06/13 |
 | wear vector                                                          |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::MergeWear(Teuchos::RCP<Epetra_Vector>& disinterface_s ,
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
 | Vector of interface displacements in ALE dofs            farah 05/13 |
 | Currently just for 1 interface                                       |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::InterfaceDisp(Teuchos::RCP<Epetra_Vector>& disinterface_s,
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
    dserror("Error in WEAR::Partitioned: Only for one interface yet.");
  }

  // loop over all interfaces
  for (int m=0; m<(int)interface.size(); ++m)
  {
    Teuchos::RCP<CONTACT::WearInterface> winterface = Teuchos::rcp_dynamic_cast<CONTACT::WearInterface>(interface[m]);
    if (winterface==Teuchos::null) dserror("Casting to WearInterface returned null!");

    // get slave row dofs as map
    Teuchos::RCP<Epetra_Map> slavedofs = winterface->SlaveRowDofs();
    // additional spatial displacements
    disinterface_s = Teuchos::rcp(new Epetra_Vector(*slavedofs,true));

    // also math. neg. values possible!
    for (int i=0; i<disinterface_s->MyLength(); ++i)
    {
      if (i%dim == 0)
      {
        if (cstrategy.ContactWear() == Teuchos::null) dserror("unweighted wear vector returned null-pointer");

        (*disinterface_s)[i] = -(*cstrategy.ContactWear())[i];
        (*disinterface_s)[i+1] = -(*cstrategy.ContactWear())[i+1];

        if (dim==3)
        {
          (*disinterface_s)[i+2] = -(*cstrategy.ContactWear())[i+2];
        }
        i=i+dim-1;
      }
    }

    // get master row dofs as map
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
  } // interface loop
  return;
}  // STRU_ALE::Algorithm::InterfaceDisp()

/*----------------------------------------------------------------------*
 | Application of mesh displacement                           mgit 07/11|
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::ApplyMeshDisplacement(bool iterated)
{
  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  // mesh displacement from solution of ALE field in structural dofs
  // first perform transformation from ale to structure dofs
  RCP<Epetra_Vector> disale   = AleToStructure(AleField().Dispnp());
  RCP<Epetra_Vector> disalen  = AleToStructure(AleField().Dispn());

  // vector of current spatial displacements
  RCP<const Epetra_Vector> dispnp = StructureField()->Dispnp();  // change to ExtractDispn() for overlap

  int aletype = DRT::INPUT::IntegralValue<int>(ParamsAle(),"ALE_TYPE");

  // for incremental lin ale --> in spatial conf.
  if (aletype==INPAR::ALE::incr_lin)
  {
    disale->Update(-1.0,*dispnp,1.0);
    delta_ale_->Update(1.0,*disale,0.0);
  }
  else if (aletype==INPAR::ALE::classic_lin)
  {
    disale->Update(-1.0,*disalen,1.0);
  }

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
      XMat[2] = node->X()[2];

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
  } // end row node loop

//  double norm2 = 0.0;
//  dismat->Norm2(&norm2);
//  std::cout << "dismat= " << norm2 << std::endl;

  // apply material displacements to structural field
  // if advection map is not succesful --> use old xmat
  StructureField()->ApplyDisMat(dismat);

  return;
}  // STRU_ALE::Algorithm::SubtractWear()

/*----------------------------------------------------------------------*
 | material coordinates evaluated from spatial ones          mgit 05/11 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::AdvectionMap(double* XMat1,
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

  // parameter space coordinates
  double e1 = 0.0;
  double e2 = 0.0;
  double e3 = 0.0;
  double ge1 = 1e12;
  double ge2 = 1e12;
  double ge3 = 1e12;
  int gele = 0;

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
      if (actele->Shape() == DRT::Element::quad4)
      {
        WEAR::UTILS::av<DRT::Element::quad4>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->Shape() == DRT::Element::quad8)
      {
        WEAR::UTILS::av<DRT::Element::quad8>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->Shape() == DRT::Element::quad9)
      {
        WEAR::UTILS::av<DRT::Element::quad8>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->Shape() == DRT::Element::tri3)
      {
        WEAR::UTILS::av<DRT::Element::tri3>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->Shape() == DRT::Element::tri6)
      {
        WEAR::UTILS::av<DRT::Element::tri6>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else
        dserror("shape function not supported!");

      // checks if the spatial coordinate lies within this element
      // if yes, returns the material displacements
//      w1ele->AdvectionMapElement(XMat1,XMat2,XMesh1,XMesh2,disp,dispmat, la,found,e1,e2);

      if (found==false)
      {
        if (abs(ge1)>1.0 and abs(e1)<abs(ge1))
        {
          ge1=e1;
          gele=jele;
        }
        if (abs(ge2)>1.0 and abs(e2)<abs(ge2))
        {
          ge2=e2;
          gele=jele;
        }
      }
    }
    else
    {
      if (actele->ElementType() == DRT::ELEMENTS::So_hex8Type::Instance())
      {
        // cast element to solid hex8 element
        DRT::ELEMENTS::So_hex8* ele = static_cast<DRT::ELEMENTS::So_hex8*>(actele);
        WEAR::UTILS::av<DRT::Element::hex8>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->ElementType() == DRT::ELEMENTS::So_hex20Type::Instance())
      {
        // cast element to solid hex20 element
        DRT::ELEMENTS::So_hex20* ele = static_cast<DRT::ELEMENTS::So_hex20*>(actele);
        WEAR::UTILS::av<DRT::Element::hex20>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->ElementType() == DRT::ELEMENTS::So_hex27Type::Instance())
      {
        // cast element to solid hex27 element
        DRT::ELEMENTS::So_hex27* ele = static_cast<DRT::ELEMENTS::So_hex27*>(actele);
        WEAR::UTILS::av<DRT::Element::hex27>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->ElementType() == DRT::ELEMENTS::So_tet4Type::Instance())
      {
        // cast element to solid tet4 element
        DRT::ELEMENTS::So_tet4* ele = static_cast<DRT::ELEMENTS::So_tet4*>(actele);
        WEAR::UTILS::av<DRT::Element::tet4>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else if (actele->ElementType() == DRT::ELEMENTS::So_tet10Type::Instance())
      {
        // cast element to solid tet10 element
        DRT::ELEMENTS::So_tet10* ele = static_cast<DRT::ELEMENTS::So_tet10*>(actele);
        WEAR::UTILS::av<DRT::Element::tet10>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
      }
      else
        dserror("elementtype not supported!");

      if (found==false)
      {
        if (abs(ge1)>1.0 and abs(e1)<abs(ge1))
        {
          ge1=e1;
          gele=jele;
        }
        if (abs(ge2)>1.0 and abs(e2)<abs(ge2))
        {
          ge2=e2;
          gele=jele;
        }
        if (abs(ge3)>1.0 and abs(e3)<abs(ge3))
        {
          ge3=e3;
          gele=jele;
        }
      }
    }

    // leave when element is found
    if (found == true)
      return;
  } // end loop over adj elements


  // ****************************************
  //  if displ not into elements
  // ****************************************
  DRT::Element* actele = ElementPtr[gele];

  // get element location vector, dirichlet flags and ownerships
  // get element location vector
  DRT::Element::LocationArray la(1);
  actele->LocationVector(*(StructureField()->Discretization()),la,false);

  // get state
  RCP<const Epetra_Vector> disp = (StructureField()->Discretization())->GetState("displacement");
  RCP<const Epetra_Vector> dispmat = (StructureField()->Discretization())->GetState("material_displacement");

  if (ndim == 2)
  {
    if (actele->Shape() == DRT::Element::quad4)
    {
      WEAR::UTILS::av<DRT::Element::quad4>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->Shape() == DRT::Element::quad8)
    {
      WEAR::UTILS::av<DRT::Element::quad8>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->Shape() == DRT::Element::quad9)
    {
      WEAR::UTILS::av<DRT::Element::quad8>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->Shape() == DRT::Element::tri3)
    {
      WEAR::UTILS::av<DRT::Element::tri3>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->Shape() == DRT::Element::tri6)
    {
      WEAR::UTILS::av<DRT::Element::tri6>(actele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else
      dserror("shape function not supported!");
  }
  else
  {
    if (actele->ElementType() == DRT::ELEMENTS::So_hex8Type::Instance())
    {
      // cast element to solid hex8 element
      DRT::ELEMENTS::So_hex8* ele = static_cast<DRT::ELEMENTS::So_hex8*>(actele);
      WEAR::UTILS::av<DRT::Element::hex8>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->ElementType() == DRT::ELEMENTS::So_hex20Type::Instance())
    {
      // cast element to solid hex20 element
      DRT::ELEMENTS::So_hex20* ele = static_cast<DRT::ELEMENTS::So_hex20*>(actele);
      WEAR::UTILS::av<DRT::Element::hex20>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->ElementType() == DRT::ELEMENTS::So_hex27Type::Instance())
    {
      // cast element to solid hex27 element
      DRT::ELEMENTS::So_hex27* ele = static_cast<DRT::ELEMENTS::So_hex27*>(actele);
      WEAR::UTILS::av<DRT::Element::hex27>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->ElementType() == DRT::ELEMENTS::So_tet4Type::Instance())
    {
      // cast element to solid tet4 element
      DRT::ELEMENTS::So_tet4* ele = static_cast<DRT::ELEMENTS::So_tet4*>(actele);
      WEAR::UTILS::av<DRT::Element::tet4>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else if (actele->ElementType() == DRT::ELEMENTS::So_tet10Type::Instance())
    {
      // cast element to solid tet10 element
      DRT::ELEMENTS::So_tet10* ele = static_cast<DRT::ELEMENTS::So_tet10*>(actele);
      WEAR::UTILS::av<DRT::Element::tet10>(ele,XMat1,XMat2,XMat3,XMesh1,XMesh2,XMesh3,disp,dispmat,la[0].lm_,found,e1,e2,e3);
    }
    else
      dserror("elementtype not supported!");
  }

  // bye
  return;

}  // STUE_ALE::Algorithm::AdvectionMap()

/*----------------------------------------------------------------------*
 | Perform ALE step                                         farah 11/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::AleStep(Teuchos::RCP<Epetra_Vector> idisale_global)
{
  int aletype = DRT::INPUT::IntegralValue<int>(ParamsAle(),"ALE_TYPE");

  // for incremental lin ale --> in spatial conf.
  if (aletype==INPAR::ALE::incr_lin)
  {
    // system of equation
    AleField().BuildSystemMatrix();

    RCP<Epetra_Vector> dispnpstru = StructureToAle(StructureField()->Dispnp());
    AleField().WriteAccessDispnp()->Update(1.0,*(dispnpstru),0.0);

    // application of interface displacements as dirichlet conditions
    AleField().ApplyInterfaceDisplacements(idisale_global);

    // solution
    AleField().SolveWear();
  }
  // classical lin in mat. conf --> not correct at all
  else if (aletype==INPAR::ALE::classic_lin)
  {
    // system of equation
    AleField().BuildSystemMatrix();

    // application of interface displacements as dirichlet conditions
    AleField().ApplyInterfaceDisplacements(idisale_global);

    // solution
    AleField().Solve();
  }
  else
    dserror("Choosen ALE type not supported for wear problems");


  return;
}

/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> WEAR::Partitioned::AleToStructure(Teuchos::RCP<Epetra_Vector> vec) const
{
  return coupalestru_->MasterToSlave(vec);
}

/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> WEAR::Partitioned::AleToStructure(Teuchos::RCP<const Epetra_Vector> vec) const
{
  return coupalestru_->MasterToSlave(vec);
}

/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> WEAR::Partitioned::StructureToAle(Teuchos::RCP<Epetra_Vector> vec) const
{
  return coupalestru_->SlaveToMaster(vec);
}

/*----------------------------------------------------------------------*
 | transform from ale to structure map                      farah 11/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> WEAR::Partitioned::StructureToAle(Teuchos::RCP<const Epetra_Vector> vec) const
{
  return coupalestru_->SlaveToMaster(vec);
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)    farah 10/13 |
 *----------------------------------------------------------------------*/
void WEAR::Partitioned::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  AleField().ReadRestart(step);
  SetTimeStep(StructureField()->GetTime(),step);

  return;
}
/*----------------------------------------------------------------------*/
