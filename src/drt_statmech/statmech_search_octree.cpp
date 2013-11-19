/*!----------------------------------------------------------------------
\file statmech_search.H
\brief octree search
<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*----------------------------------------------------------------------*/
#include "statmech_search.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                   mueller 09/13|
 *----------------------------------------------------------------------*/
STATMECH::SEARCH::Octree::Octree(Teuchos::RCP<std::vector<double> > periodlength,
                                         Teuchos::RCP<DRT::Discretization>  discret,
                                         Teuchos::RCP<Epetra_Map>           bspotrowmap,
                                         Teuchos::RCP<Epetra_Map>           bspotcolmap,
                                         int                                treedepth,
                                         int                                minbboxesperoct,
                                         double                             bindingradius)
{
	this->bsr_=bindingradius;
	this->bspotrowmap_=bspotrowmap;
	this->bspotcolmap_=bspotcolmap;
	this->discret_=discret;
	this->periodlength_=periodlength;
	this->maxtreedepth_ = treedepth;
	this->minbboxesinoctant_ = minbboxesperoct;
  if((int)periodlength_->size()<3)
    dserror("You only gave %d values for PERIODLENGTH! Check your input file.", (int)periodlength_->size());
  if(periodlength_->at(0)>1e-12)
    periodicBC_ = true;
  else
    periodicBC_ = false;
  //Initialize Maps
  bbox2octant_ = Teuchos::rcp(new Epetra_MultiVector(*(bspotcolmap()),4));
  bbox2octant_->PutScalar(-1.0);

  return;
}

/*----------------------------------------------------------------------*
 |  build octree (public)                                  mueller 09/13|
 *----------------------------------------------------------------------*/
void STATMECH::SEARCH::Octree::BuildOctree(const Teuchos::RCP<Epetra_MultiVector>& bspotpositions)
{
  // initialize class vectors
	allbboxes_ = Teuchos::rcp(new Epetra_MultiVector(*(bspotcolmap_),4, true));
  //Create spherical BoundingBoxes
	//Loop over all Binding Spots
	for (int bspot=0; bspot<bspotpositions->MyLength(); bspot++)
	{
		for(int dof=0; dof<3; dof++)
		  (*allbboxes_)[dof][bspot] = (*bspotpositions)[dof][bspot];
		(*allbboxes_)[allbboxes_->NumVectors()-1][bspot] = bspotcolmap_->GID(bspot);
	}


  // call recursive octtree build
  // clear vector for assigning bounding boxes to octants to be on the safe side before (re)assigning bounding boxes
	LINALG::Matrix<6,1> rootbox;
	for(int i=0; i<(int)rootbox.M(); i++)
	{
		if(i%2==0)
			rootbox(i)= 0.0;
		else
			rootbox(i)=periodlength_->at((i-1)/2);
	}

	//Initialize Tree that enables to find out in which octant a certain point lies (later)
	RootNode_= new Octree::OctreeNode(rootbox);

  // Convert Epetra_MultiVector allbboxes_ to vector(vector<double>)
  std::vector<std::vector<double> > allbboxesstdvec(allbboxes_->MyLength(), std::vector<double>(allbboxes_->NumVectors(),0.0));
  for(int i=0; i < allbboxes_->MyLength(); i++)
    for(int j=0; j<allbboxes_->NumVectors(); j++)
      allbboxesstdvec[i][j] = (*allbboxes_)[j][i];

  //initial tree depth value (will be incremented with each recursive call of locateBox()
  int treedepth = 0;

  // Parameters and Initialization
  std::vector<std::vector<int> > bboxesinoctants;
  bboxesinoctants.clear();
  octreelimits_.clear();
  // Recursively construct octree
  if(!discret_->Comm().MyPID())
    LocateBoundingBox(allbboxesstdvec, rootbox, octreelimits_, bboxesinoctants, treedepth, *RootNode_);

/*
 * CHECK WHETHER NODE TREE WORKS
 *
 *
 *
  LINALG::Matrix<3,1> position;
  position(0)=3.5;
  position(1)=4.5;
  position(2)=1;
  OctreeNode::OctreeNode* pointer = BspotRootNode;
  	bool HIT=false;
  	while(HIT==false)
  	{
  		LINALG::Matrix<6,1> lim;
  		for(int i=0;i<pointer->HaveChildren();i++)
  		{
  			cout<< "IN for schleife HaveChildren()"<<endl;
  			lim=pointer->children[i]->LimitsOfNode();
  			//check whether position lies in this node
  			if(!( (lim(0)>=position(0)) ||(lim(1)<position(0)) ||(position(1)<=lim(2))||(position(1)>lim(3)) ||(position(2)<=lim(4))||(position(2)>lim(5)) ))
  			{
  				pointer= pointer->children[i]; // pointer umbiegen auf dieses KInd. Geht das?
  				cout<< "Number of children (nachumbiegen) "<<pointer->HaveChildren()<< "limits des kinds "<< pointer->LimitsOfNode() <<endl;
  				cout<< " OKTANT ID "<< pointer->GetOctantID()<<endl;
  				if(pointer->HaveChildren()==0) // set flag if arrived and the end of the branch (no more sub-nodes)
  					HIT=true;
  				break; //break iloop and continue in sub-child
  			}
  		}
  	}
  	int ID;
  	ID=pointer->GetOctantID();
  	cout<< "***********************************************************"<<endl;
  	cout<< " OKTANT ID ........................: " << ID << endl;
  	cout<< "***********************************************************"<<endl;
*/

  //determine maximum depth of OctreeMap
  int maxdepthlocal = 0;
  int bboxlengthlocal = 0;
  if(discret_->Comm().MyPID()==0)
  {
    bboxlengthlocal = (int)bboxesinoctants.size();
    for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
      if((int)bboxesinoctants[i].size()>maxdepthlocal)
        maxdepthlocal = (int)bboxesinoctants[i].size();
  }

  int maxdepthglobal = 0;
  int bboxlengthglobal = 0;
  discret_->Comm().MaxAll(&maxdepthlocal, &maxdepthglobal, 1);
  discret_->Comm().MaxAll(&bboxlengthlocal, &bboxlengthglobal, 1);

  /* build temporary, fully overlapping map and row map for octree
   * Note: maxdepthglobal does not occur for a converging Newton iteration. Yet, in some cases, when
   * encountering divergence for the Newton scheme, this might happen.
   * In biopolymer network simulations, this setting is not unlikely and unavoidable. a maximum depth of
   * 0 means, there are no bounding boxes/elements in any octants. Hence, we will not detect any contact and
   * therefore skip the rest of the octree algorithm.*/
  if(maxdepthglobal>0)
  {
    // create octree maps
    std::vector<int> gids;
    for (int i=0 ; i<bboxlengthglobal; i++ )
      gids.push_back(i);
    // crosslinker column and row map
    Epetra_Map octtreerowmap((int)gids.size(), 0, discret_->Comm());
    Epetra_Map octtreemap(-1, (int)gids.size(), &gids[0], 0, discret_->Comm());

    // build Epetra_MultiVectors which hold the BBs of the OctreeMap; for communication
    bboxesinoctants_ = Teuchos::rcp(new Epetra_MultiVector(octtreemap,maxdepthglobal));
    bboxinoctrow_ = Teuchos::rcp(new Epetra_MultiVector(octtreerowmap, maxdepthglobal, true));
    //Epetra_MultiVector bboxinoctrow(octtreerowmap,maxdepthglobal, true);

    // fill bboxinoct for Proc 0
    if(discret_->Comm().MyPID()==0)
    {
      bboxesinoctants_->PutScalar(-9.0);
      for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
        for(int j=0; j<(int)bboxesinoctants[i].size(); j++)
          (*bboxesinoctants_)[j][i] = bboxesinoctants[i][j];
    }

    // Communication
    CommunicateMultiVector(bboxinoctrow_, bboxesinoctants_);
  }
  else
    dserror("(buildbspotoctree()) -> maxdepthglobal=0 ");
}//End of BuildBspotOctree

/*----------------------------------------------------------------------*
 |  Locate bounding box within octree (public)             mueller 09/13|
 *----------------------------------------------------------------------*/
void STATMECH::SEARCH::Octree::LocateBoundingBox(std::vector<std::vector<double> >& allbboxesstdvec,
                                                 LINALG::Matrix<6,1>&               lim,
                                                 std::vector<LINALG::Matrix<6,1> >& OctreeLimits,
                                                 std::vector<std::vector<int> >&    bboxesinoctants,
                                                 int&                               treedepth,
                                                 Octree::OctreeNode&                parentNode)
{
  // Center of octant
  LINALG::Matrix<3,1> center;
  // edge length vector of the suboctants
  LINALG::Matrix<3,1> newedgelength;
  for(int i=0; i<(int)center.M(); i++)
  {
    center(i) = (lim(2*i)+lim(2*i+1))/2.0;
    newedgelength(i) = fabs(lim(2*i+1)-(lim)(2*i))/2.0;
  }
  std::vector<LINALG::Matrix<6,1> > limits;
  limits.clear();
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
      {
        LINALG::Matrix<6,1> sublim;
        sublim(0) = center(0) + (i-1)*newedgelength(0);
        sublim(1) = center(0) +  i   *newedgelength(0);
        sublim(2) = center(1) + (j-1)*newedgelength(1);
        sublim(3) = center(1) +  j   *newedgelength(1);
        sublim(4) = center(2) + (k-1)*newedgelength(2);
        sublim(5) = center(2) +  k   *newedgelength(2);

        limits.push_back(sublim);
      }

  /* Decision to which child box belongs....................
  *
  *           5 ======================== 7
  *           //|                       /||
  *          // |                      //||
  *         //  |                     // ||
  *        //   |                    //  ||
  *       //    |                   //   ||
  *      //     |                  //    ||
  *     //      |                 //     ||
  *    1 ========================= 3     ||
  *    ||       |                ||      ||
  *    ||       |                ||      ||
  *    ||       |      o (center)||      ||
  *    ||      4 ----------------||------ 6
  *    ||      /                 ||     //
  *    ||     /                  ||    //
  *    ||    /                   ||   //
  *    ||   /                    ||  //
  *    ||  /                     || //      y  z
  *    || /                      ||//       | /
  *    ||/                       ||/        |/
  *    0 ========================= 2        ---> x
  *
  */
  Teuchos::RCP<LINALG::Matrix<3,1> > octcenter = Teuchos::null;

  //Goes through all suboctants
  for( int oct=0; oct<8; oct++)
  {
    // Define temporary vector of same size as current allbboxesstdvec
    std::vector<std::vector<double> > bboxsubset;
    bboxsubset.clear();
    //add children to Octree
    OctreeNode* child = new Octree::OctreeNode(limits[oct]); //generates a child = new node with suboctant-limits
    child->SetOctantID(-9);
    parentNode.AddChild(child);

    if(periodlength_->at(0) > 0.0)
    {
    	for( int i=0; i<(int)allbboxesstdvec.size(); i++) //LOOP OVER ALL BSPOTS
      {
				/* Now we check whether bounding box lies in this oktant.
				 * We also shift the bounding box by the periodlength to consider bounding boxes that overlap the periodic boundary
				 * p=0  --> original position of bounding box
				 * p=-1 --> shift positions by -periodlength
				 * p=+1 --> shift positions by +periodlength
				 */
				for(int p=-1;p<2;p++)
				{	if(!((limits[oct](0) >= allbboxesstdvec[i][0] + bsr_ + p*periodlength_->at(0))  || (limits[oct](1) <= allbboxesstdvec[i][0] - bsr_ + p*periodlength_->at(0) ) || (limits[oct](2) >= allbboxesstdvec[i][1] + bsr_ + p*periodlength_->at(1)) || (limits[oct](3) <= allbboxesstdvec[i][1] - bsr_ + p*periodlength_->at(1)) || (limits[oct](4) >= allbboxesstdvec[i][2] + bsr_ + p*periodlength_->at(2)) || (limits[oct](5) <= allbboxesstdvec[i][2] - bsr_ + p*periodlength_->at(2))))
					{
						bboxsubset.push_back(allbboxesstdvec[i]);
						break; //i.e. exit for loop after first hit
					}
				}
      } // end of for-loop which goes through all elements of input
    }
    else // standard procedure without periodic boundary conditions
    {
    	dserror("Octree is only implemented for periodic boundaries");
    }

    // current tree depth
    int currtreedepth = treedepth+1;
    // Check for further recursion by checking number of boxes in octant (first criterion)....................
    int N = (int)bboxsubset.size();

    //If to divide further, let LocateBox call itself with updated inputs
    if (N > minbboxesinoctant_ && currtreedepth < maxtreedepth_-1)
      LocateBoundingBox(bboxsubset, limits[oct], OctreeLimits, bboxesinoctants, currtreedepth,*(child));
    else
    {
      // no further discretization of the volume because either the maximal tree depth or the minimal number of bounding
      // boxes per octant has been reached
      // this vector holds the IDs of the bounding boxes in this octant
      if(N>0)
      {
        std::vector<int> boxids;
        boxids.clear();
        //Push back Limits of suboctants to OctreeLimits
        OctreeLimits.push_back(limits[oct]);

        for (int m = 0; m < (int)bboxsubset.size(); m++)
        {
          // note: the Bounding Box ID is the last column entry of the m-th entry vector bboxsubset
          boxids.push_back((int)bboxsubset[m][(int)bboxsubset[m].size()-1]);
          // assign current octant number to the bounding box
           for(int n=0; n<bbox2octant_->NumVectors(); n++)
            if((*bbox2octant_)[n][bspotcolmap_->LID((int)bboxsubset[m][(int)bboxsubset[m].size()-1])]<-0.9)
            {
              (*bbox2octant_)[n][bspotcolmap_->LID((int)bboxsubset[m][(int)bboxsubset[m].size()-1])] = (double)bboxesinoctants.size();
              break; //leave after finding first empty slot
            }
        }
        //Add Octant ID to OctreeNode
        child->SetOctantID((int)bboxesinoctants.size());
        // add bounding box IDs of this octant to the global vector
        bboxesinoctants.push_back(boxids);
      }
      else
      {
      	// cout<< "let me know if there is no bbox in this octant"<<endl;
      }
     }
  }// end of loop which goes through all suboctants
  return;
} // end of method locateBox

/*----------------------------------------------------------------------*
 |  Locate positions in an octree already built (public)   mueller 09/13|
 *----------------------------------------------------------------------*/
void STATMECH::SEARCH::Octree::LocatePositions(Teuchos::RCP<Epetra_MultiVector>& positions,
                                               Teuchos::RCP<Epetra_Map>          CLmap)
{
  std::vector<std::vector<int> > clinoctants((int)bboxesinoctants_->MyLength());
	LINALG::Matrix<6,1> lim;
	//only on proc 0
  if(discret_->Comm().MyPID()==0)
  {
		for(int j=0;j<positions->MyLength();j++)
		{		//Pointer on rootnode to start at beginnging of tree
			Octree::OctreeNode* pointer = RootNode_;
			//Flag indicating end of branch
			bool HIT=false;

			while(HIT==false)
			{
				for(int i=0;i<pointer->HaveChildren();i++)
				{
					lim=pointer->children[i]->LimitsOfNode();
					//check whether position lies in this node
					if(!( (lim(0)>=(*positions)[0][j]) ||(lim(1)<(*positions)[0][j]) ||(lim(2)>=(*positions)[1][j])||(lim(3)<(*positions)[1][j]) ||(lim(4)>=(*positions)[2][j])||(lim(5)<(*positions)[2][j]) ))
					{
						pointer= pointer->children[i];
						if(pointer->HaveChildren()==0) // set flag if arrived and the end of the branch (no more sub-nodes)
							HIT=true;
						break; //break iloop and continue in sub-child
					}
				}
      }
      //Get octant ID
      if((int)pointer->GetOctantID()>-1)
        clinoctants[(int)pointer->GetOctantID()].push_back(j);
		}//end of for-loop over all crosslinker positions
  }

  //make vector redundant on all procs:
  int maxCrosslinkerInOctantlocal = 0;
  int maxCrosslinkerInOctantglobal = 0;
  if(discret_->Comm().MyPID()==0)
  {
    for (int i=0 ; i<(int)clinoctants.size(); i++ )
      if((int)clinoctants[i].size()>maxCrosslinkerInOctantlocal)
        maxCrosslinkerInOctantlocal = (int)clinoctants[i].size();
  }
  discret_->Comm().MaxAll(&maxCrosslinkerInOctantlocal, &maxCrosslinkerInOctantglobal, 1);
  // build Epetra_MultiVectors which hold the BBs of the OctreeMap; for communication
  crosslinkerinoctants_ = Teuchos::rcp(new Epetra_MultiVector((*bboxesinoctants_).Map(),maxCrosslinkerInOctantglobal));

  // fill bboxinoct for Proc 0
  if(discret_->Comm().MyPID()==0)
  {
    crosslinkerinoctants_->PutScalar(-9.0);
    for (int i=0 ; i<(int)clinoctants.size(); i++ )
      for(int j=0; j<(int)clinoctants[i].size(); j++)
        (*crosslinkerinoctants_)[j][i] = clinoctants[i][j];
  }

  //Communicate this multivector to make it redundant on all Procs
  Teuchos::RCP<Epetra_MultiVector> crosslinkerinoctantsrow = Teuchos::rcp(new Epetra_MultiVector((*bboxinoctrow_).Map(),maxCrosslinkerInOctantglobal,true));
  CommunicateMultiVector(crosslinkerinoctantsrow, crosslinkerinoctants_);

  return;
}//end of LocatePositions()

/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors               mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::SEARCH::Octree::CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                                                      Teuchos::RCP<Epetra_MultiVector> OutVec,
                                                      bool                             doexport,
                                                      bool                             doimport,
                                                      bool                             zerofy,
                                                      bool                             exportinsert)
{
  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec->Map(), InVec->Map());
  Epetra_Import importer(OutVec->Map(), InVec->Map());

  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec->PutScalar(0.0);
    if(exportinsert)
      InVec->Export(*OutVec, exporter, Insert);
    else
      InVec->Export(*OutVec, exporter, Add);
  }
  if(doimport)
    OutVec->Import(*InVec,importer,Insert);
  return;
}

/*----------------------------------------------------------------------*
 |  constructor (public)                                   mueller 09/13|
 *----------------------------------------------------------------------*/
STATMECH::SEARCH::Octree::OctreeNode::OctreeNode(LINALG::Matrix<6,1>& inlimits) //LINALG::Matrix<6,1> limits
{
	this->limits=inlimits;
	this->OctantID=-1;
}

/*----------------------------------------------------------------------*
 |Tell whether this node has children (suboctants)                      |
 |or is at the end of a branch (public)                    mueller 09/13|
 *----------------------------------------------------------------------*/
int STATMECH::SEARCH::Octree::OctreeNode::HaveChildren()
{
	return this->children.size();
}

/*----------------------------------------------------------------------*
 |  return the limits of the octant associated                          |
 |  with this node (public)                                mueller 09/13|
 *----------------------------------------------------------------------*/
LINALG::Matrix<6,1> STATMECH::SEARCH::Octree::OctreeNode::LimitsOfNode()
{
	return this->limits;
}

/*----------------------------------------------------------------------*
 |  Add child   (public)                                   mueller 09/13|
 *----------------------------------------------------------------------*/
void STATMECH::SEARCH::Octree::OctreeNode::AddChild(OctreeNode* new_node)
{
	this->children.push_back(new_node);
}

/*----------------------------------------------------------------------*
 |  Set octant ID of an octree node (public)               mueller 09/13|
 *----------------------------------------------------------------------*/
void STATMECH::SEARCH::Octree::OctreeNode::SetOctantID(int id)
{
	OctantID=id;
}

/*----------------------------------------------------------------------*
 |  Get octant ID of an octree node(public)                mueller 09/13|
 *----------------------------------------------------------------------*/
int STATMECH::SEARCH::Octree::OctreeNode::GetOctantID()
{
	return this->OctantID;
}




