/*!----------------------------------------------------------------------
\file statmech_polymerization.cpp
\brief (de-)polymerization  for StatMech problems

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>
--------------------------------------------------------------------------*/
#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beam3eb/beam3eb.H"
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
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  discret_->GetCondition("FilamentNumber",filamentnumberconditions);
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  // Create a vector which stores first and last node GIDS of every filaments
  std::vector<int> FirstNode((int)filamentnumberconditions.size());
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  std::vector<int> LastNode((int)filamentnumberconditions.size());
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  // Temporary vector storing the end node GIDs of a filament
  std::vector<int> EndNodes(2);
//  std::cout<<__FILE__<<__LINE__ <<std::endl;

  BarbedEnds_ = Teuchos::rcp(new std::vector<double>);
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  BarbedEnds_->clear();
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  PointedEnds_ = Teuchos::rcp(new std::vector<double>);
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  PointedEnds_->clear();
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
  FirstNode[0]= 0; // Assuming 1st node of the discretization is also 1st node of filament
//  std::cout<<__FILE__<<__LINE__ <<std::endl;
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
//  int ncrosslink = statmechparams_.get<int> ("N_crosslink", 0);
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
    {
      PolyEndRowGIDs.push_back((*BarbedEnds_)[i]);
//      for (int i=0; i<(int)PolyEndRowGIDs.size(); i++)
//        std::cout<<"PolyEndRowIDs="<<PolyEndRowGIDs[i]<<std::endl;
    }
//    std::cout<<"polyspotonproc="<<polyspotonproc[i]<<std::endl;
//    std::cout<<"discret_->NodeRowMap()->LID((int)(*BarbedEnds_)[i])="<<discret_->NodeRowMap()->LID((int)(*BarbedEnds_)[i])<<std::endl;
  }
  // Create row maps
  PolyEndRowMap_ = Teuchos::rcp(new Epetra_Map((int)PolyEndGIDs.size(), (int)PolyEndRowGIDs.size(), &PolyEndRowGIDs[0], 0, discret_->Comm()));

  // create redundant binding spot column map based upon the polymerising ends
  PolyEndColMap_ = Teuchos::rcp(new Epetra_Map(-1, (int)PolyEndGIDs.size(), &PolyEndGIDs[0], 0, discret_->Comm()));
//  std::cout<<"discret_->NodeRowMap()="<<std::endl;
//  discret_->NodeRowMap()->Print(std::cout);

//  std::cout<<"PolyEndRowMap_="<<std::endl;
//  PolyEndRowMap_->Print(std::cout);


//  std::cout<<"PolyEndColMap_"<<std::endl;
//  PolyEndColMap_->Print(std::cout);
     //     startindex_ = Teuchos::rcp(new std::vector<double>);


//  bspotcolmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeColMap())));
//  bspotrowmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeRowMap())));
//  bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
//  bspotstatus_->PutScalar(-1.0);

//  // create density-density-correlation-function map with
//  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly ||
//     DRT::INPUT::IntegralValue<int>(statmechparams_, "GMSHOUTPUT"))
//  {
//    std::vector<int> bins;
//    for(int i=0; i<discret_->Comm().NumProc()*numbins; i++)
//      bins.push_back(i);
//    ddcorrcolmap_ = Teuchos::rcp(new Epetra_Map(-1, discret_->Comm().NumProc()*numbins, &bins[0], 0, discret_->Comm()));
//    // create processor-specific density-density-correlation-function map
//    ddcorrrowmap_ = Teuchos::rcp(new Epetra_Map(discret_->Comm().NumProc()*numbins, 0, discret_->Comm()));
//    // create new trafo matrix (for later use in DDCorr Function where we evaluate in layer directions), initialize with identity matrix
//    trafo_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,3,true));
//    for(int i=0; i<trafo_->M(); i++)
//      (*trafo_)(i,i) = 1.0;
//    for(int i=0; i<(int)cog_.M(); i++)
//      cog_(i) = periodlength_->at(i)/2.0;
//
//    // start indices for parallel handling of orientation correlation etc in NumLinkerSpotsAndOrientation()
//    // calculation of start indices for each processor
//    // number of overall independent combinations
//    int numnodes = discret_->NodeColMap()->NumMyElements();
//    int numcombinations = (numnodes*numnodes-numnodes)/2;
//    // combinations on each processor
//    int combinationsperproc = (int)floor((double)numcombinations/(double)discret_->Comm().NumProc());
//    int remainder = numcombinations%combinationsperproc;
//
//    // get starting index tuples for later use
//    startindex_->assign(2*discret_->Comm().NumProc(), 0.0);
//
//    for(int mypid=0; mypid<discret_->Comm().NumProc()-1; mypid++)
//    {
//      std::vector<int> start(2,0);
//      bool continueloop = false;
//      bool quitloop = false;
//      int counter = 0;
//      int appendix = 0;
//      if(mypid==discret_->Comm().NumProc()-1)
//        appendix = remainder;
//
//      // loop over crosslinker pairs
//      for(int i=0; i<discret_->NodeColMap()->NumMyElements(); i++)
//      {
//        for(int j=0; j<discret_->NodeColMap()->NumMyElements(); j++)
//        {
//          if(i==(*startindex_)[2*mypid] && j==(*startindex_)[2*mypid+1])
//            continueloop = true;
//          if(j>i && continueloop)
//          {
//            if(counter<combinationsperproc+appendix)
//              counter++;
//            else
//            {
//              // new start index j
//              if(j==discret_->NodeColMap()->NumMyElements()-1)
//                start[1] = 0;
//              else
//                start[1] = j;
//              quitloop = true;
//              break;
//            }
//          }
//        }
//        if(quitloop)
//        {
//          // new start index i
//          if(start[1]==0)
//            start[0] = i+1;
//          else
//            start[0] = i;
//          // new start tuple
//          (*startindex_)[2*(mypid+1)] = (double)(start[0]);
//          (*startindex_)[2*(mypid+1)+1] = (double)(start[1]);
//          break;
//        }
//      }
//    }
//
//    if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
//    {
//      // initialize unbinding times vector to "-1"
//      crosslinkunbindingtimes_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2));
//      crosslinkunbindingtimes_->PutScalar(-1.0);
//    }
//  }

//    element2crosslink_ = Teuchos::null;
//    unbindingprobability_ = Teuchos::null;
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
    {
      (*MonomerPositionsRow)[j][i] = upperbound.at(j) * (*uniformgen_)();
//      std::cout<<"(*uniformgen_)()="<< (*uniformgen_)()<<"upperbound.at(j)="<<upperbound.at(j)<<std::endl;

    }


  CommunicateMultiVector(MonomerPositionsRow, MonomerPositions_,false,true);

  // initial bonding status is set (no bonds)
  MonomerBond_ = Teuchos::rcp(new Epetra_Vector(*MonomerColMap_, true));
  MonomerBond_->PutScalar(-1.0);

  // crosslinker element IDs of the crosslink molecules
//  crosslink2element_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
//  crosslink2element_->PutScalar(-1.0);

  // initialize the beautiful visuals vector (aka beevee-vector)
  VisualizeMonomerPositions_ = Teuchos::rcp(new Epetra_MultiVector(*MonomerColMap_, 3, false));
  for (int i=0; i<VisualizeMonomerPositions_->MyLength(); i++)
    for (int j=0; j<VisualizeMonomerPositions_->NumVectors(); j++)
      (*VisualizeMonomerPositions_)[j][i] = (*MonomerPositions_)[j][i];



  // calculate  number of total binding spot depending on BSPOTINTERVAL
  NumPolySpots_ = 0;
  NumPolySpots_ = BarbedEnds_->size(); // Check if it returns the correct value

#ifdef DEBUGCOUT
  std::cout<<"Proc "<<discret_->Comm().MyPID()<<": total number of binding spots = "<<numbspots_<<std::endl;
#endif

  return;
}//StatMechManager::MonomerInitiatePos

/*---------------------------------------------------------------------------*
| simulation of monomer diffusion                    (private) mukherjee 03/14|
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
      // bonding case 2: crosslink molecule attached to one filament
      case 1:
      {
      }
      break;
      // bonding case 3: actual crosslinker has been established
      case 2:
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

#ifdef DEBUGCOUT
  Teuchos::RCP<Epetra_MultiVector> crosslinkerdelta = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));
  CommunicateMultiVector(crosslinkerdeltatrans, crosslinkerdelta, false,true);

  for(int i=0; i<crosslinkerpositionstrans->NumVectors(); i++)
    for(int j=0; j<crosslinkerpositionstrans->MyLength(); j++)
      if((*crosslinkerpositionstrans)[i][j]>periodlength_->at(i) ||(*crosslinkerpositionstrans)[i][j]<0.0)
        dserror("Crosslinker %i outside of the periodic box: %d, %d, %d ", j, (*crosslinkerpositionstrans)[0][j], (*crosslinkerpositionstrans)[1][j], (*crosslinkerpositionstrans)[2][j]);

  if(!discret_->Comm().MyPID())
  {
    // check increments of the stochastic process
    std::stringstream crossposfile;
    crossposfile<<"./crosspos.txt";
    FILE* fpcross = NULL;
    fpcross = fopen(crossposfile.str().c_str(), "w");

    std::stringstream crosspos;
    for(int i=0; i<crosslinkerdelta->NumVectors(); i++)
      for(int j=0; j<crosslinkerdelta->MyLength(); j++)
        crosspos<<(*crosslinkerdelta)[i][j]<<" "<<(*crosslinkerdelta)[i][j]<<" "<<(*crosslinkerdelta)[i][j]<<std::endl;
    // move temporary std::stringstream to file and close it
    fprintf(fpcross, crosspos.str().c_str());
    fclose(fpcross);

    // check paths of the stochastic process
    if(crosslinkerpositions_->MyLength()>1)
      dserror("Check only for a single crosslink molecule!");
    std::ostringstream filename;
    filename << outputrootpath_ << "/crosslinkermotion.dat";
    FILE *fp = fopen(filename.str().c_str(), "a");

    std::stringstream linkerposition;
    linkerposition<<(*crosslinkerpositions_)[0][0]<<" "<<(*crosslinkerpositions_)[1][0]<<" "<<(*crosslinkerpositions_)[2][0]<<std::endl;
    fprintf(fp, linkerposition.str().c_str());
    fclose(fp);
  }
#endif
  return;
}// StatMechManager::MonomerDiffusion

/*----------------------------------------------------------------------*
 | Searches for monomer molecule-filament node pairs and adds actual    |
 | crosslinker elements once certain conditions are met.                |
 | (private)                                           mukherjee (03/14)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SearchAndSetMonomers(const int&                                       istep,
                                                     const double&                                    timen,
                                                     const double&                                    dt,
                                                     const Teuchos::RCP<Epetra_MultiVector>           PolySpotPositions,
                                                     const Teuchos::RCP<Epetra_MultiVector>           PolySpotRotations,
                                                     const Epetra_Vector&                             discol,
                                                     bool                                             printscreen)
{

  double t0 = Teuchos::Time::wallTime();
  /*preliminaries*/

  //get current on-rate for crosslinkers
  double kon = 0;
  double starttime = actiontime_->at(1);

  if(timen <= starttime || (timen>starttime && fabs(timen-starttime) < dt/1e4))
    kon = statmechparams_.get<double>("K_ON_monomer",0.0);
  else
    kon = statmechparams_.get<double>("K_ON_end",0.0);

  //probability with which a crosslinker is established between crosslink molecule and neighbour node
  double plink = 1.0 - exp( -dt*kon );

  //Volume partitioning, assignment of nodes and molecules to partitions, search for neighbours
  // map filament (column map, i.e. entire node set on each proc) node positions to volume partitions every SEARCHINTERVAL timesteps

  Teuchos::RCP<Epetra_MultiVector> neighbourslid;

  PartitioningAndSearchMonomer(PolySpotPositions, neighbourslid);

#ifdef MEASURETIME
  double t1 = Teuchos::Time::wallTime();
  double t2 = -1.0e9;
#endif
  /*the following part of the code is executed on processor 0 only and leads to the decision, at which nodes crosslinks shall be set
   *processor 0 goes through all the crosslink molecules and checks whether a crosslink is to be set; this works precisely as follows:
   *(1) the crosslink molecules are looped through in a random order
   *(2) if a node has not yet reached its maximal number of crosslinks, a crosslink may be set
   *(3) a crosslink is set if and only if the node passes a probability check
   *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

  // a vector indicating the crosslink molecule which is going to constitute a crosslinker element
  Teuchos::RCP<Epetra_Vector> addcrosselement = Teuchos::rcp(new Epetra_Vector(*MonomerColMap_, true));

  int numsetelements = 0;
  int numstdlinkers = 0;
  int numactlinkers = 0;

  if(discret_->Comm().MyPID()==0)
  {
    // obtain a random order in which the monomers are addressed
    std::vector<int> order = Permutation(MonomerColMap_->NumMyElements());

    for(int i=0; i<MonomerBindingStatus_->MyLength(); i++)
    {
      int irandom = order[i];
      // only consider crosslink molecule with less than two nodes or if it is actively searching for bonding partners
      if((*MonomerBindingStatus_)[irandom]<1.9)
      {
        // PROBABILITY CHECK & PREPARATION OF VECTORS
        // obtain a random order of neighboursLID indices
        std::vector<int> neighbourorder = Permutation(neighbourslid->NumVectors());

        // loop over neighbouring binding spots
        for(int j=0; j<neighbourslid->NumVectors(); j++)
        {
          // random index
          int index = neighbourorder[j];
          // skip this, if neighbourslid entry is '-2', meaning empty or binding spot is unavailable due to input file specs
          if((*neighbourslid)[index][irandom] > -1.9)
          {
            // check the space between already set linkers and the about to be set linker
            bool setlink = true;
            int bspotinterval = statmechparams_.get<int>("BSPOTINTERVAL", 1);
//            if(bspotinterval>1 && (*neighbourslid)[index][irandom]>-0.9)
//            {
//              int bspotlid = (*neighbourslid)[index][irandom];
//              //getting a vector consisting of pointers to all filament number conditions set
//              std::vector<DRT::Condition*> filaments(0);
//              discret_->GetCondition("FilamentNumber",filaments);
//              int filnum = (*filamentnumber_)[bspotlid];
//              int startlid = bspotlid - int(ceil(double(bspotinterval/2)));
//              int endlid = bspotlid + int(ceil(double(bspotinterval/2)));
//              if(bspotcolmap_->GID(startlid)<(*filaments[filnum]->Nodes())[0])
//                startlid = bspotcolmap_->LID((*filaments[filnum]->Nodes())[0]);
//              if(bspotcolmap_->GID(endlid)>(*filaments[filnum]->Nodes())[(int)filaments[filnum]->Nodes()->size()-1])
//                endlid = bspotcolmap_->LID((*filaments[filnum]->Nodes())[(int)filaments[filnum]->Nodes()->size()-1]);
//
//              for(int k=startlid; k<endlid+1; k++)
//              {
//                if((*bspotstatus_)[k]>-1)
//                {
//                  setlink = false;
//                  break;
//                }
//              }
//            }

            if(setlink)
            {
              // current neighbour LID
              int bspotLID = (int)(*neighbourslid)[index][irandom];

              //skip the rest of this iteration of the i-loop in case the binding spot is occupied
              if((*bspotstatus_)[bspotLID]>-0.1)
                continue;

              // flag indicating loop break after first new bond has been established between i-th crosslink molecule and j-th neighbour node
              bool bondestablished = false;
              // necessary condition to be fulfilled in order to set a crosslinker
              double probability = 0.0;
              // switch between probability for regular inter-filament binding and self-binding crosslinker
              probability = plink;
//              // self binding is activated
//              if(statmechparams_.get<double>("K_ON_SELF", 0.0)>0.0)
//              {
//                if((int)(*numbond_)[irandom]==1)
//                {
//                  int nodelid0 = -1;
//                  int nodelid1 = -1;
//                  nodelid1 = bspotLID;
//                  for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
//                    if((*crosslinkerbond_)[k][irandom]>-0.1)
//                    {
//                      nodelid0 = discret_->NodeColMap()->LID((int)(*crosslinkerbond_)[k][irandom]);
//                      break;
//                    }
//                }
//              }

              if( (*uniformgen_)() < probability )
              {
                int free = 0;
                int occupied = 0;
                bool set = false;

                // check both entries of crosslinkerbond_
                for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
                  if((*crosslinkerbond_)[k][irandom]<-0.9 && !set)
                  {
                    // store free bond position
                    free = k;
                    set = true;
                  }
                  else
                    occupied = k;

                switch((int)(*numbond_)[irandom])
                {
                  // free crosslink molecule creating a bond
                  case 0:
                  {
                    // update of crosslink molecule positions
                    LINALG::SerialDenseMatrix LID(1,1,true);
                    LID(0,0) = bspotLID;
                    (*crosslinkerbond_)[free][irandom] = bspotcolmap_->GID(bspotLID);
                    // update status of binding spot by putting in the crosslinker id
                    ((*bspotstatus_)[bspotLID]) = irandom;
                    // increment the number of bonds of this crosslinker
                    ((*numbond_)[irandom]) = 1.0;
                    MonomerIntermediateUpdate(*PolySpotPositions, LID, irandom);
                    bondestablished = true;
                  }
                  break;
                }// switch((int)(*numbond_)[irandom])
                // for now, break j-loop after a new bond was established, i.e crosslinker elements cannot be established starting from zero bonds
                if(bondestablished)
                  break;
              }//if( (*uniformgen_)() < probability )
            }//
          }// if(neighbourslid...)
        }// for(int j=0; j<(int)neighboursLID.size(); j++)
      }//if((*numbond_)[irandom]<1.9)
    }// for(int i=0; i<numbond_->MyLength(); i++)
#ifdef MEASURETIME
    t2 = Teuchos::Time::wallTime();
#else
    if(printscreen)
      std::cout << "\nsearch + linker admin time: " << Teuchos::Time::wallTime() - t0<< " seconds";
#endif
  }// if(discret_->Comm().MypPID==0)
  /* note: searchforneighbours_ is not communicated
   * to the other Procs because their information is of concern to Proc 0 only.*/
  //synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_,true));
  CommunicateVector(bspotstatusrow, bspotstatus_);

  // transfer vectors
  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
  Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,2,true));
  Teuchos::RCP<Epetra_Vector> numbondtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  Teuchos::RCP<Epetra_Vector> addcrosselementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
  // exports and reimports
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
  CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);
  CommunicateVector(numbondtrans, numbond_);
  CommunicateVector(addcrosselementtrans, addcrosselement);

#ifdef MEASURETIME
  double t3 = Teuchos::Time::wallTime();
#endif

  // ADDING ELEMENTS
  // the node row map is needed in order to adequately assessing ownership of elements to processors
  const Epetra_Map noderowmap(*(discret_->NodeRowMap()));

  int numsetelementsall = -1;
  discret_->Comm().MaxAll(&numsetelements, &numsetelementsall, 1);

  if(numsetelementsall>0)
  {
    // positions from displacement vector
    Teuchos::RCP<Epetra_MultiVector> nodalpositions = Teuchos::null;
    // rotations from displacement vector
    Teuchos::RCP<Epetra_MultiVector> nodalrotations = Teuchos::null;
    // NODAL quaternions from filament elements
    Teuchos::RCP<Epetra_MultiVector> nodalquaternions = Teuchos::null;

    // add elements to problem discretization (processor specific)
    for(int i=0; i<addcrosselement->MyLength(); i++)
    {
      if((*addcrosselement)[i]>0.9)
      {
        /*there is the problem of how to assign to a crosslinker element a GID which is certainly not used by any
         * other element; we know that for a network at the beginning (without crosslinkers) each node has a
         * connectivity of 1 or 2 so that the number of all elements used to discretize the filaments is smaller
         * than basisnodes_. Thus we choose crosslinker GIDs >= basisnodes_ for the additionally added crosslinker
         * elements. To make sure that two crosslinkers never have the same GID we add to basisnodes_ the value
         * basisnodes_*GID1, where GID1 is the GID of the first node of the crosslinker element. Then we add GID2
         * (note: GID2 < basisnodes_), which is the GID of the second node of the crosslinker element.
         * Hence basisnodes_ + GID1*basisnodes_ + GID2 always gives a GID which cannot be used by any other element;
         * the first node of the crosslinker element is always assumed to be the one with the greater GID; the owner
         * of the first node is assumed to be the owner of the crosslinker element*/

        // copied binding spot map (in order to handle the addition of elements based upon an invariant map just in case)
        Epetra_Map polyendcolmap(*PolyEndColMap_);

        // obtain polymerization spot GID
        int polyspotgid =(int)(*MonomerBond_)[i];
        // determine smaller and larger of the GIDs
//        bspotgid.at(1) = std::min((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);
//        bspotgid.at(0) = std::max((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);

//        // different sizes due to different linker elements
//        Teuchos::RCP<std::vector<int> > globalnodeids;
//
//        globalnodeids = Teuchos::rcp(new std::vector<int>(2,0));
//        (*globalnodeids)[0] = bspotgid[0];
//        (*globalnodeids)[1] = bspotgid[1];

        // Create new node




        // calculate element GID
        int newcrosslinkerGID = GenerateNewFilamentGID(polyspotgid);

        /* Create mapping from crosslink molecule to crosslinker element GID
         * Note: No need for the usual procedure of exporting and reimporting to make things redundant
         * because info IS already redundant by design here.*/
        (*crosslink2element_)[i] = newcrosslinkerGID;

        //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
        std::vector<double> rotrefe(6);
        std::vector<double> xrefe(6);
        /* resize in case of interpolated crosslinker element. No need for interpolated Beam3eb element, because rotational
           d.o.f. does not exist */

        // all beam elements
          for(int k=0; k<3; k++)
          {
            xrefe[k]   = (*PolySpotPositions)[k][polyendcolmap.LID(polyspotgid)];
            xrefe[k+3] = (*PolySpotPositions)[k][polyendcolmap.LID(polyspotgid)];
            rotrefe[k] = (*PolySpotRotations)[k][polyendcolmap.LID(polyspotgid)];
            rotrefe[k+3] = (*PolySpotRotations)[k][polyendcolmap.LID(polyspotgid)];
         }


        /*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
         *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
         *crosslinker element will be present as a ghost element only*/
        bool hasrownode = CheckRowNode(noderowmap,polyspotgid);

        if(hasrownode)
        {
          AddFilamentNode(polyspotgid);
//          AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,*discret_,nodalquaternions);
        }



      }
    }

    // synchronization for problem discretization
    discret_->CheckFilledGlobally();
    discret_->FillComplete(true, false, false);

    // synchronization for contact discretization
  }

#ifdef MEASURETIME
  if(!discret_->Comm().MyPID())
  {
    std::cout<<"\n=================Time  Measurement================"<<std::endl;
    std::cout<<"StatMechManager::SearchAndSetCrosslinkers"<<std::endl;
    std::cout<<"Partitioning + Search  :\t"<<std::setprecision(4)<<t1-t0<<"\ts"<<std::endl;
    std::cout<<"Linker administration  :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<std::endl;
    std::cout<<"Vector communication   :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<std::endl;
    std::cout<<"Addition of elements   :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t3<<"\ts"<<std::endl;
    std::cout<<"=================================================="<<std::endl;
    std::cout<<"total time             :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<std::endl;
  }
#endif

  //std::couts
  if(!discret_->Comm().MyPID() && printscreen)
  {
    std::cout<<"\n\n"<<numsetelements<<" crosslinker element(s) added!"<<std::endl;
    if(numstdlinkers+numactlinkers>0)
    {
      std::cout<<" - "<<numstdlinkers<<" standard linkers"<<std::endl;
      std::cout<<" - "<<numactlinkers<<" active   linkers"<<std::endl;
    }
  }
}//void StatMechManager::SearchAndSetCrosslinkers

/*-------------------------------------------------------------------------*
 | Assign monomers and nodes to volume partitions                          |
 |                                                (public) mukherjee 03/15 |
 *-------------------------------------------------------------------------*/
void STATMECH::StatMechManager::PartitioningAndSearchMonomer(const Teuchos::RCP<Epetra_MultiVector> PolySpotPositions,
                                                              Teuchos::RCP<Epetra_MultiVector>&     neighbourslid)
{
  std::vector<double> limit(6,0.0);
  if(periodlength_->at(0)==0.0)
  {
    // initialize
    for(int i=0; i<(int)limit.size(); i++)
    {
      if(i%2==0)
        limit[i] = 1e9;
      else
        limit[i] = -1e9;
    }
    // extreme values
    // binding spots
    for(int i=0; i<(int)PolySpotPositions->NumVectors(); i++)
    {
      for(int j=0; j<(int)PolySpotPositions->MyLength(); j++)
      {
        // max
        if(limit[2*i+1]<(*PolySpotPositions)[i][j])
          limit[2*i+1] = (*PolySpotPositions)[i][j];
        // min
        if(limit[2*i]>(*PolySpotPositions)[i][j])
          limit[2*i] = (*PolySpotPositions)[i][j];
      }
    }
    // linker positions
    for(int i=0; i<(int)MonomerPositions_->NumVectors(); i++)
    {
      for(int j=0; j<(int)MonomerPositions_->MyLength(); j++)
      {
        // max
        if(limit[2*i+1]<(*MonomerPositions_)[i][j])
          limit[2*i+1] = (*MonomerPositions_)[i][j];
        // min
        if(limit[2*i]>(*MonomerPositions_)[i][j])
          limit[2*i] = (*MonomerPositions_)[i][j];
      }
    }
  }
  else
  {
    for(int i=0; i<(int)limit.size(); i++)
      if(i%2==0)
        limit[i] = 0.0;
      else
        limit[i] = (*periodlength_)[(i-1)/2];
  }

  //filament binding spots and monomers are indexed according to their positions within the boundary box volume
  std::vector<std::vector<std::vector<int> > > polyspotinpartition;
  for(int i=0; i<(int)searchres_->size(); i++)
    polyspotinpartition.push_back(std::vector<std::vector<int> >((*searchres_)[i], std::vector<int>()));

  /*nodes*/
  // loop over node positions to map their column map LIDs to partitions
  for(int i=0; i<(int)polyspotinpartition.size(); i++) // note: bspotinpartition.size==3
  {
    for(int j=0;j<PolySpotPositions->MyLength();j++)
    {
      int partition = (int)std::floor(((*PolySpotPositions)[i][j]-limit[2*i])/(limit[2*i+1]-limit[2*i])*(double)(*searchres_)[i]);
      if(partition==(int)(*searchres_)[i])
        partition--;
      polyspotinpartition[i][partition].push_back(PolyEndColMap_->LID(j)); //column lid
    }
  }
  /*Monomers*/
  // Export crosslinkerpositions_ to transfermap_ format (row map format for crosslink molecules)
  Teuchos::RCP<Epetra_MultiVector> monomerpositionsrow = Teuchos::rcp(new Epetra_MultiVector(*MonomerRowMap_, 3, true));
  Teuchos::RCP<Epetra_Vector> MonomerBindingStatusrow = Teuchos::rcp(new Epetra_Vector(*MonomerRowMap_, true));
  Teuchos::RCP<Epetra_MultiVector> monomerspartitionrow = Teuchos::rcp(new Epetra_MultiVector(*MonomerRowMap_, 3, false));

  CommunicateVector(MonomerBindingStatusrow, MonomerBindingStatus_, true, false, false, true);
  CommunicateMultiVector(monomerpositionsrow, MonomerPositions_, true, false, false, true);

  for(int i=0; i<monomerspartitionrow->NumVectors(); i++)
  {
    for(int j=0; j<monomerspartitionrow->MyLength(); j++)
    {
      // mark entries with double-bonded crosslink molecules
      if((*MonomerBindingStatusrow)[j]>1.9)
      {
        (*monomerspartitionrow)[i][j] = -1.0;
        continue;
      }
      else
      {
        int partition = (int)std::floor(((*monomerpositionsrow)[i][j]-limit[2*i])/(limit[2*i+1]-limit[2*i])*(double)(*searchres_)[i]);
        if(partition==(*searchres_)[i])
          partition--;
        (*monomerspartitionrow)[i][j] = partition;
      }
    }
  }

  // detection of nodes within search proximity of the crosslink molecules
  DetectBindingSpots(PolySpotPositions, polyspotinpartition, MonomerBindingStatusrow, monomerpositionsrow, monomerspartitionrow, Teuchos::null, neighbourslid);

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

//   If the new cordinate lies outside the periodic box,
//   the do a periodic boundary shift.




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
