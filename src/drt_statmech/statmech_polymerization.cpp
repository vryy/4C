/*!----------------------------------------------------------------------
\file statmech_polymerization.cpp
\brief (de-)polymerization  for StatMech problems


\maintainer Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270

--------------------------------------------------------------------------*/
#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beam3/beam3eb.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_octtree.H"

/*------------------------------------------------------------------------*
 | create barbed and pointed end node sets of filaments                   |
 |                                           (private)     mukherjee 02/15|
 *------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CreateBarbedPointedNodeSets()
{
  //getting a vector consisting of pointers to all filament number conditions set
  std::vector<DRT::Condition*> filamentnumberconditions(0);
  discret_->GetCondition("FilamentNumber",filamentnumberconditions);
  // Create a vector which stores first and last node GIDS of every filaments
  std::vector<int> FirstNode((int)filamentnumberconditions.size());
  std::vector<int> LastNode((int)filamentnumberconditions.size());
  // Temporary vector storing the end node GIDs of a filament
  std::vector<int> EndNodes(2);
  BarbedEnds_ = Teuchos::rcp(new std::vector<double>);
  BarbedEnds_->clear();
  PointedEnds_ = Teuchos::rcp(new std::vector<double>);
  PointedEnds_->clear();
  FirstNode[0]= 0; // Assuming 1st node of the discretization is also 1st node of filament
  //next all the pointers to all the different conditions are looped
  for (int i=0; i<(int)filamentnumberconditions.size(); ++i)
  {
    //get a pointer to nodal cloud covered by the current condition
    const std::vector<int>* nodeids = filamentnumberconditions[i]->Nodes();

    LastNode[i]= FirstNode[i]+(int)nodeids->size()-1;
    // At Last node of filament (LastNodeID +1) doesn't exist
    if (i<(int)filamentnumberconditions.size()-1)
      FirstNode[i+1]=LastNode[i]+1;

    // Assign values of endnodes
    EndNodes[0]=FirstNode[i];
    EndNodes[1]=LastNode[i];

    // Randomly shuffle the vector (Although one could randomly assign
    // barbed ends and poined ends, this information will be different
    // for differnt processor. In that case, one needs to have a seperate
    // treatement)
    // std::random_shuffle(EndNodes.begin(), EndNodes.end());

    PointedEnds_->push_back(EndNodes[0]); // Last node of filament as Barbed end
    BarbedEnds_->push_back(EndNodes[1]);  // First node as pointed end

    /* The last node is marked as Barbed end.
     * Because the new element is added in the direction of the tangent vector,
     * which is pointing outwards (away from the filament). This has to be taken
     * in consideration while modeling  the pointed end */
  }

  return;
}

/*----------------------------------------------------------------------------------*
 | Initialize monomer positions                                                     |
 |(point-like particles and not doubly-bound)              (private) mukherjee 03/15|
 *----------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::MonomerInitiatePos()
{
  int nmonomer = 50;
//  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);

  // create monomer maps
  std::vector<int> gids;
  for (int i=0; i<nmonomer; i++)
    gids.push_back(i);
  // Monomer column map
  MonomerColMap_ = Teuchos::rcp(new Epetra_Map(-1, nmonomer, &gids[0], 0, discret_->Comm()));
  // Monomer row map
  MonomerRowMap_ = Teuchos::rcp(new Epetra_Map(nmonomer, 0, discret_->Comm()));

  // vector holding all polymerizing end gids
  std::vector<int>    PolyEndGIDs;

  // verctor holding GIDs of the row map
  std::vector<int> PolyEndRowGIDs;
   PolyEndRowGIDs.clear();
  for(int i=0; i<(int)BarbedEnds_->size(); i++)
  {
    PolyEndGIDs.push_back((int)PolyEndGIDs.size());
    if(discret_->NodeRowMap()->LID((int)(*BarbedEnds_)[i])!=-1)
      PolyEndRowGIDs.push_back((*BarbedEnds_)[i]);
  }
  // Create row maps
  PolyEndRowMap_ = Teuchos::rcp(new Epetra_Map((int)PolyEndGIDs.size(), (int)PolyEndRowGIDs.size(), &PolyEndRowGIDs[0], 0, discret_->Comm()));

  // create redundant binding spot column map based upon the polymerising ends
  PolyEndColMap_ = Teuchos::rcp(new Epetra_Map(-1, (int)PolyEndGIDs.size(), &PolyEndGIDs[0], 0, discret_->Comm()));

  // initial bond counter is set (no bonds)
  MonomerBindingStatus_ = Teuchos::rcp(new Epetra_Vector(*MonomerColMap_, true));

  std::vector<double> upperbound = *periodlength_;
  // handling both cases: with and without periodic boundary conditions
  if (periodlength_->at(0) == 0.0)
    for(int i=0; i<(int)upperbound.size(); i++)
      upperbound.at(i) = statmechparams_.get<double> ("MaxRandValue", 0.0);

  MonomerPositions_ = Teuchos::rcp(new Epetra_MultiVector(*MonomerColMap_, 3, true));

  Teuchos::RCP<Epetra_MultiVector> MonomerPositionsRow = Teuchos::rcp(new Epetra_MultiVector(*MonomerRowMap_,3,true));

  for (int i=0; i<MonomerPositionsRow->MyLength(); i++)
    for (int j=0; j<MonomerPositionsRow->NumVectors(); j++)
      (*MonomerPositionsRow)[j][i] = upperbound.at(j) * (*uniformgen_)();


  CommunicateMultiVector(MonomerPositionsRow, MonomerPositions_,false,true);

  // initial bonding status is set (no bonds)
  MonomerBond_ = Teuchos::rcp(new Epetra_Vector(*MonomerColMap_, true));
  MonomerBond_->PutScalar(-1.0);

  // initialize the beautiful visuals vector (aka beevee-vector)
  VisualizeMonomerPositions_ = Teuchos::rcp(new Epetra_MultiVector(*MonomerColMap_, 3, false));
  for (int i=0; i<VisualizeMonomerPositions_->MyLength(); i++)
    for (int j=0; j<VisualizeMonomerPositions_->NumVectors(); j++)
      (*VisualizeMonomerPositions_)[j][i] = (*MonomerPositions_)[j][i];



  // calculate  number of total binding spot depending on BSPOTINTERVAL
  NumPolySpots_ = 0;
  NumPolySpots_ = BarbedEnds_->size(); // Check if it returns the correct value

  return;
}//StatMechManager::MonomerInitiatePos

/*---------------------------------------------------------------------------*
| simulation of monomer diffusion                   (private) mukherjee 03/15|
*----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::MonomerDiffusion(double                  mean,
                                                 double           standarddev,
                                                 const double&             dt)
{
  /* Here, the diffusion of monomers is handled.
   * Depending on the number of occupied binding spots of the molecule, its motion
   * is calculated differently.
   */
#ifdef DEBUGCOUT
  Teuchos::RCP<Epetra_MultiVector> crosslinkerdeltatrans = Teuchos::rcp(new Epetra_MultiVector(*MonomerRowMap_, 3, true));
#endif

  Teuchos::RCP<Epetra_MultiVector> MonomerPositionsRow = Teuchos::rcp(new Epetra_MultiVector(*MonomerRowMap_, 3, true));
  CommunicateMultiVector(MonomerPositionsRow, MonomerPositions_, true, false,true);

  // bonding cases
  for (int i=0; i<MonomerPositionsRow->MyLength(); i++)
  {
    int MonomerLID = MonomerColMap_->LID(MonomerRowMap_->GID(i));
    switch ((int)(*MonomerBindingStatus_)[MonomerLID])
    {
      // bonding case 1:  no bonds, diffusion
      case 0:
      {
        for (int j=0; j<MonomerPositionsRow->NumVectors(); j++)
          (*MonomerPositionsRow)[j][i] += standarddev*(*normalgen_)() + mean;

#ifdef DEBUGCOUT
      for (int j=0; j<MonomerPositionsRow->NumVectors(); j++)
        (*MonomerPositionsRow)[j][i] += standarddev*(*normalgen_)() + mean;
#endif
      }
      break;
      // bonding case 2: monomer attached to one filament
      case 1:
      {
      }
      break;
    }
  }

  // check for compliance with periodic boundary conditions if existent
  if (periodlength_->at(0) > 0.0)
    CrosslinkerPeriodicBoundaryShift(MonomerPositionsRow);

  // Update by Broadcast: make this information redundant on all procs
  CommunicateMultiVector(MonomerPositionsRow, MonomerPositions_, false, true);

  return;
}// StatMechManager::MonomerDiffusion

/*----------------------------------------------------------------------*
 | Searches for monomer molecule-filament node pairs and adds actual    |
 | crosslinker elements once certain conditions are met.                |
 | (private)                                           mukherjee (03/15)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SearchAndSetMonomers(const int&                                       istep,
                                                     const double&                                    timen,
                                                     const double&                                    dt,
                                                     const Teuchos::RCP<Epetra_MultiVector>           PolySpotPositions,
                                                     const Teuchos::RCP<Epetra_MultiVector>           PolySpotRotations,
                                                     const Epetra_Vector&                             discol,
                                                     bool                                             printscreen)
{

 return;
}//void StatMechManager::SearchAndSetMonomers

/*-------------------------------------------------------------------------*
 | Assign monomers and nodes to volume partitions                          |
 |                                                (public) mukherjee 03/15 |
 *-------------------------------------------------------------------------*/
void STATMECH::StatMechManager::PartitioningAndSearchMonomer(const Teuchos::RCP<Epetra_MultiVector> PolySpotPositions,
                                                              Teuchos::RCP<Epetra_MultiVector>&     neighbourslid)
{

  return;
}//void StatMechManager::PartitioningAndSearchMonomer

/*----------------------------------------------------------------------*
 | update monomer positions                                             |
 |                                              (public) mukherjee 03/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::MonomerIntermediateUpdate(const Epetra_MultiVector& polyspotpositions,
                                                              const LINALG::SerialDenseMatrix& LID,
                                                              const int& monomernumber)
{
  // case: one-bonded monomer(i.e. two cases: +1 bond (starting from 0 bonds) or -1 bond (molecule is free to diffuse again)
  if (LID.M()==1 && LID.N()==1)
  {
    // set molecule position to node position
    for (int i=0; i < MonomerPositions_->NumVectors(); i++)
      (*MonomerPositions_)[i][monomernumber] = polyspotpositions[i][(int)LID(0,0)];
  }
  else
  {
    // generate vector in random direction of length R_LINK to "reset" monomer position:
    // it may now reenter or leave the bonding proximity
    LINALG::Matrix<3, 1> deltapos;
    for (int i=0; i<(int)deltapos.M(); i++)
      deltapos(i) = (*uniformgen_)();

    deltapos.Scale(statmechparams_.get<double> ("R_LINK", 0.0) / deltapos.Norm2());
    for (int i=0; i<MonomerPositions_->NumVectors(); i++)
      (*MonomerPositions_)[i][monomernumber] += deltapos(i);
  }
  return;
}// StatMechManager::MonomerIntermediateUpdate

/*------------------------------------------------------------------------*
 | add new node at the barbed end of the filament                         |
 |                                           (private)  mukherjee 03/15   |
 *------------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddFilamentNode(int& polyspotgid)
{

  // This method has to be called at the end of each time step
  // If certain criteria are fulfilled
  // After new node is added, the "BarbedEnds_" has to be modified
  // The binding spots are also to be modified.

  // Create a new nodal ID for the new node.
  // This will be essentially the new polymerizing end
  int NewNodalGID= GenerateNewNodeGID(polyspotgid);
  // export row displacement to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);

  /* Create coordinates of the new node */

  DRT::Node* node = discret_->gNode(polyspotgid);

  std::vector<int> DofNode = discret_->Dof(node);

  DRT::Element* Element= node->Elements()[0];

  LINALG::Matrix<3,1> TcurrNode(true);
  LINALG::Matrix<3,1> XrefNewNode(true);
  //if node has also rotational degrees of freedom
  if (discret_->NumDof(node) == 6)
  {
    //Check via dynamic cast, if it's a beam3eb element
    DRT::ELEMENTS::Beam3eb* BeamElement = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(Element);
    //if not, tell the user to use beam3 instead
    if(BeamElement==NULL)
    {
      dserror("The element is not beam element. Can not proceed further");
    }
    else
    {
      double LengthRef= BeamElement->LengthRef();
      LINALG::Matrix<3,1> TrefNode(true);
      TrefNode=BeamElement->Tref()[0];
      for(int j=0; j<3; j++)
      {
        TcurrNode(j)=TrefNode(j)+discol[discret_->DofColMap()->LID(DofNode[j+3])];
        XrefNewNode(j)=node->X()[j]+TcurrNode(j);
      }
    }
  }

  return;
}

/*------------------------------------------------------------------------*
 | add new element at the barbed end of the filament                      |
 |                                           (private)     mukherjee      |
 *------------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddFilamentElement()
{

  // This method has to be called at the end of each time step
  // If certain criteria are fulfilled
  // After new element is added, the "BarbedEnds_" has to be modified

  return;
}

/*------------------------------------------------------------------------*
 | delete element at the pointed end of the filament                      |
 |                                           (private)     mukherjee      |
 *------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DeleteFilamentElement()
{

  // This method has to be called at the end of each time step
  // If certain criteria are fulfilled
  // After element is added, the "PointedEnds_" has to be modified

  return;
}

/*----------------------------------------------------------------------*
 | Generate new nodal  GID for a  new filament element                  |
 |                                                     mukherjee  02/15 |
 *----------------------------------------------------------------------*/
int STATMECH::StatMechManager::GenerateNewNodeGID(int & bspotgid)
{
  // calculate node GID
  int NewNodeGID = (bspotgid + 1)*basisnodes_;

  /* Correction of the nodal GID if there happens to be another node with the same ID.
   * it's a failsafe method. Chances of occurance is very minimal
   * As long as an unused GID cannot be found, the nodal GID keeps getting incremented by 1.*/
  discret_->Comm().Barrier();
  while(1)
  {
    int gidexists = 1;
    // query existance of node on this Proc
    int gidonproc = (int)(discret_->HaveGlobalNode(NewNodeGID));
    // sum over all processors
    discret_->Comm().MaxAll(&gidonproc, &gidexists, 1);
    // calculate new GID if necessary by shifting the initial GID
    if(gidexists>0)
      NewNodeGID++;
    else
      break;
  }
  return NewNodeGID;
}

/*----------------------------------------------------------------------*
 | Generate a new filament element GID                 mukherjee  02/15 |
 *----------------------------------------------------------------------*/
int STATMECH::StatMechManager::GenerateNewFilamentGID(int & bspotgid)
{
  // calculate element GID
  int NewFilamentGID = (bspotgid + 1)*basisnodes_;

  /* Correction of the filament GID if there happens to be another filament with the same ID.
   * This might occur if there is a crosslinker molecule with same GID.
   * As long as an unused GID cannot be found, the filament GID keeps getting incremented by 1.*/
  discret_->Comm().Barrier();
  while(1)
  {
    int gidexists = 1;
    // query existance of node on this Proc
    int gidonproc = (int)(discret_->HaveGlobalElement(NewFilamentGID));
    // sum over all processors
    discret_->Comm().MaxAll(&gidonproc, &gidexists, 1);
    // calculate new GID if necessary by shifting the initial GID
    if(gidexists>0)
      NewFilamentGID++;
    else
      break;
  }
  return NewFilamentGID;
}
