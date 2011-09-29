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
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                   mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "stru_ale_algorithm.H"
#include "stru_ale_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_w1/wall1.H"
#include "../drt_lib/drt_elementtype.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*
 | constructor (public)                                      mgit 05/11 |
 *----------------------------------------------------------------------*/
STRU_ALE::Algorithm::Algorithm(Epetra_Comm& comm)
 : AlgorithmBase(comm,DRT::Problem::Instance()->StructuralDynamicParams()),
   StructureBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams())
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

    // predict and solve structural system
    StructureField().PrepareTimeStep();
    StructureField().Solve();

    // calculate stresses, strains, energies
    StructureField().PrepareOutput();

    // write output to screen and files
    StructureField().Update();

    // write output to screen and files
    StructureField().Output();

    // subtract wear
    SubtractWear();

  }  // time loop
}  // STRUE_ALE::Algorithm::TimeLoop()

/*----------------------------------------------------------------------*
 | Subtract wear                                              mgit 12/09 |
 *----------------------------------------------------------------------*/
void STRU_ALE::Algorithm::SubtractWear()
{

  if (cmtman_==null)
    dserror("Structure with ale only for contact with wear so far.");

  if (!cmtman_->GetStrategy().Wear())
    return;

  // vector of current spatial displacements
  RCP<Epetra_Vector> dispn = StructureField().ExtractDispn();

  // additional spatial displacements
  RCP<Epetra_Vector> disadditional = rcp(new Epetra_Vector(dispn->Map()),true);

  // new material displacements for slip nodes
  RCP<Epetra_Vector> dismat = rcp(new Epetra_Vector(dispn->Map()),true);

  // set state (spatial displacements)
  (StructureField().Discretization())->SetState(0,"displacement",dispn);

  // set state (material displacements)
  (StructureField().Discretization())->SetState(0,"material displacement",StructureField().DispMat());

  // get vector of unweighted wear
  RCP<Epetra_Vector> realwear = cmtman_->GetStrategy().ContactWear();

  // vector of slipnodes
  RCP<Epetra_Map> sliprownodes = cmtman_->GetStrategy().SlipRowNodes();

  // loop over slipnodes (only they are worn)
  for (int j=0; j<sliprownodes->NumMyElements(); ++j)
  {
    int gid = sliprownodes->GID(j);
    DRT::Node* node = StructureField().Discretization()->gNode(gid);

    // adjacent elemets
    DRT::Element** ElementPtr = node->Elements();
    int numelement = node->NumElement();

    // material displacement
    double XMat[2];
    double XMesh[2];

    int locid = (dispn->Map()).LID(2*gid);
    int locid1 = (realwear->Map()).LID(2*gid);

    if((*realwear)[locid1+1]>0)
    {
      // dserror for negative wear
      //dserror ("Negative wear!");

      printf("REALWEAR1: % e \n",(*realwear)[locid1]);
      printf("REALWEAR2: % e \n",(*realwear)[locid1+1]);
      (*realwear)[locid1+1]=0;
    }

    // new spatial coordinates
    XMesh[0]=node->X()[0]+(*dispn)[locid]-(*realwear)[locid1];
    XMesh[1]=node->X()[1]+(*dispn)[locid+1]-(*realwear)[locid1+1];

    // find material coordinates
    AdvectionMap(&XMat[0],&XMat[1],&XMesh[0],&XMesh[1],ElementPtr,numelement);

    // additional spatial displacements
    (*disadditional)[locid]   = -(*realwear)[locid1];
    (*disadditional)[locid+1] = -(*realwear)[locid1+1];

    // material displacements
    (*dismat)[locid] = XMat[0]-node->X()[0];
    (*dismat)[locid+1] = XMat[1]-node->X()[1];
  }

  // update spatial displacements
  dispn->Update(1.0,*disadditional,1.0);

  // apply material displacements to structural field
  StructureField().ApplyDisMat(dismat);

  return;

//  // vector of current spatial displacements
//  RCP<Epetra_Vector> dispn = StructureField().ExtractDispn();
//
//  // additional spatial displacements
//  RCP<Epetra_Vector> disadditional = rcp(new Epetra_Vector(dispn->Map()),true);
//
//  RCP<Epetra_Vector> dismat = rcp(new Epetra_Vector(dispn->Map()),true);
//
//  // set state
//  (StructureField().Discretization())->SetState(0,"displacement",dispn);
//
//  // set state
//  (StructureField().Discretization())->SetState(0,"material displacement",StructureField().DispMat());
//
//  // loop over all row nodes to fill graph
//    for (int k=0;k<StructureField().Discretization()->NumMyRowNodes();++k)
//    {
//      int gid = StructureField().Discretization()->NodeRowMap()->GID(k);
//      cout << "GID " << gid <<endl;
//
//      DRT::Node* node = StructureField().Discretization()->gNode(gid);
//      if (gid==2)
//      {
//        DRT::Element** ElementPtr = node->Elements();
//
//        int numelement = node->NumElement();
//
//        cout << "NUMELEMENT " << numelement << endl;
//
//        double XMat[2];
//        double XMesh[2];
//
//        int locid = (dispn->Map()).LID(2*gid);
//
//        cout << "STEP " << Step() << endl;
//
//        double ux = 0;
//        double uy =0;
//
//        if (Step() == 1)
//        {
//          ux= -0.2;
//          uy= -0.257196;
//        }
//
//        if (Step() ==2)
//        {
//          ux=0.1;
//          uy = +0.257196;
//        }
//
//
//        XMesh[0]=node->X()[0]+(*dispn)[locid]+ux;
//        XMesh[1]=node->X()[1]+(*dispn)[locid+1]+uy;
//
//
//        StructureField().Discretization()->AdvectionMap(&XMat[0],&XMat[1],&XMesh[0],&XMesh[1],ElementPtr,numelement);
//
//        node->UMat()[0]= XMat[0]-node->X()[0];
//        node->UMat()[1]= XMat[1]-node->X()[1];
//
//        (*disadditional)[locid]   = ux;
//        (*disadditional)[locid+1] = uy;
//
//        (*dismat)[locid] = XMat[0]-node->X()[0];
//        (*dismat)[locid+1] = XMat[1]-node->X()[1];
//
//        cout << "STRUE_ALE " << XMat[0]-node->X()[0] << endl;
//        cout << "STRUE_ALE " << XMat[1]-node->X()[1] << endl;
//      }
//    }
//
//    // update spatial displacements
//    dispn->Update(1.0,*disadditional,1.0);
//
//    // apply material displacements to structural field
//    StructureField().ApplyDisMat(dismat);
//
//    return;

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
    RefCountPtr<const Epetra_Vector> disp = (StructureField().Discretization())->GetState("displacement");
    RefCountPtr<const Epetra_Vector> dispmat = (StructureField().Discretization())->GetState("material displacement");

#ifdef D_WALL1
    // structure with ale only for wall element so far
    if (actele->ElementType() != DRT::ELEMENTS::Wall1Type::Instance())
      dserror ("Construction of advection map only for wall element so far.");

    // cast element to wall element
    DRT::ELEMENTS::Wall1* w1ele = static_cast<DRT::ELEMENTS::Wall1*>(actele);

    // checks if the spatial coordinate lies within this element
    // if yes, returns the material displacements
    w1ele->AdvectionMapElement(XMat1,XMat2,XMesh1,XMesh2,disp,dispmat,la,found);
#endif
    // leave when element is found
    if (found == true)
      return;
  }

  // error if element is not found
  if (found == false)
    dserror("STRU_ALE::Algorithm::AdvectionMap: Particle tracking not successful.");

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
#endif  // CCADISCRET
