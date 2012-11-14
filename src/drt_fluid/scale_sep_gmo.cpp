#include "../drt_fluid/scale_sep_gmo.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_cut/cut_position.H"


/*----------------------------------------------------------------------*
 | Destructor (public)                                   rasthofer 07/11|
 *----------------------------------------------------------------------*/
FLD::LESScaleSeparation::~LESScaleSeparation()
{
  return;
}


/*----------------------------------------------------------------------*
 | Constructor (public)                                  rasthofer 07/11|
 *----------------------------------------------------------------------*/
FLD::LESScaleSeparation::LESScaleSeparation(
  INPAR::FLUID::ScaleSeparation scale_sep,
  RCP<DRT::Discretization> discret):
  scale_sep_(scale_sep),
  discret_(discret),
  sepmat_build_(false)
{

 switch(scale_sep_)
 {
 case INPAR::FLUID::box_filter:
 {
   break;
 }
 case INPAR::FLUID::algebraic_multigrid_operator:
 {
   break;
 }
 case INPAR::FLUID::geometric_multigrid_operator:
 {
   ConstructSepMatGeoMultigrid();
   break;
 }
 default:
 {
   dserror("Unkonwn filter type for les!");
   break;
 }
 }

  return;
}



/*----------------------------------------------------------------------*
 | Build separation matrix for geometric multigrid       rasthofer 07/11|
 *----------------------------------------------------------------------*/
void FLD::LESScaleSeparation::ConstructSepMatGeoMultigrid()
{
  if (sepmat_build_)
   dserror("Separation-matrix has already been build!");

  /*
   * This function builds the separation matrix based on a geometric multigrid operator.
   * At the moment, it is tested for 3d-problems, only. 2d might be possible, too. Moreover,
   * as the coarse scale grid is build out of the complete grid, which is given in the input
   * file, the number of elements in each direction has to be even. PBCs are included and
   * cause no problems.
   */

  //----------------------------------------------------------------------
  // create sets of coordinates
  //----------------------------------------------------------------------

  // the criterion allows differences in coordinates by 1e-9
  set<double,LineSortCriterion> x1coords;
  set<double,LineSortCriterion> x2coords;
  set<double,LineSortCriterion> x3coords;

  // loop nodes and build sets of lines in x1-,x2- and x3-direction
  // accessible on this proc
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);

      x1coords.insert(node->X()[0]);
      x2coords.insert(node->X()[1]);
      x3coords.insert(node->X()[2]);
  }

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates to all procs
  //--------------------------------------------------------------------

  {
#ifdef PARALLEL
    int myrank  =discret_->Comm().MyPID();
#endif
    int numprocs=discret_->Comm().NumProc();

    vector<char> sblock;
    vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point communication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (set<double,LineSortCriterion>::iterator x1line=x1coords.begin();
           x1line!=x1coords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      data.StartPacking();
      for (set<double,LineSortCriterion>::iterator x1line=x1coords.begin();
           x1line!=x1coords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid  = (myrank+1)%numprocs;

      int length = sblock.size();

      exporter.ISend(frompid,topid,&(sblock[0]),sblock.size(),tag,request);

      rblock.clear();

      // receive from predecessor
      frompid = (myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblock,length);

      if(tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);

      {
        // for safety
        exporter.Comm().Barrier();
      }
#else
      // dummy communication
      rblock.clear();
      rblock=sblock;
#endif

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        vector<double> coordsvec;

        coordsvec.clear();

        vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x1coords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (set<double,LineSortCriterion>::iterator x2line=x2coords.begin();
           x2line!=x2coords.end();
           ++x2line)
      {
        DRT::ParObject::AddtoPack(data,*x2line);
      }
      data.StartPacking();
      for (set<double,LineSortCriterion>::iterator x2line=x2coords.begin();
           x2line!=x2coords.end();
           ++x2line)
      {
        DRT::ParObject::AddtoPack(data,*x2line);
      }
      swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank+1)%numprocs;

      int length = sblock.size();

      exporter.ISend(frompid,topid,&(sblock[0]),sblock.size(),tag,request);

      rblock.clear();

      // receive from predecessor
      frompid=(myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblock,length);

      if(tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);

      {
        // for safety
        exporter.Comm().Barrier();
      }
#else
      // dummy communication
      rblock.clear();
      rblock=sblock;
#endif

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        vector<double> coordsvec;

        coordsvec.clear();

        vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2coords.insert(onecoord);
        }
      }
    }

    // third, communicate coordinates in x3-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (set<double,LineSortCriterion>::iterator x3line=x3coords.begin();
           x3line!=x3coords.end();
           ++x3line)
      {
        DRT::ParObject::AddtoPack(data,*x3line);
      }
      data.StartPacking();
      for (set<double,LineSortCriterion>::iterator x3line=x3coords.begin();
           x3line!=x3coords.end();
           ++x3line)
      {
        DRT::ParObject::AddtoPack(data,*x3line);
      }
      swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank+1)%numprocs;

      int length = sblock.size();

      exporter.ISend(frompid,topid,&(sblock[0]),sblock.size(),tag,request);

      rblock.clear();

      // receive from predecessor
      frompid=(myrank+numprocs-1)%numprocs;
      exporter.ReceiveAny(frompid,tag,rblock,length);

      if(tag!=(myrank+numprocs-1)%numprocs)
      {
        dserror("received wrong message (ReceiveAny)");
      }

      exporter.Wait(request);

      {
        // for safety
        exporter.Comm().Barrier();
      }
#else
      // dummy communication
      rblock.clear();
      rblock=sblock;
#endif

      //--------------------------------------------------
      // Unpack received block into set of all planes.
      {
        vector<double> coordsvec;

        coordsvec.clear();

        vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x3coords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------

  RCP<vector<double> > x1coordinates;
  x1coordinates = Teuchos::rcp(new vector<double> );
  RCP<vector<double> > x2coordinates;
  x2coordinates = Teuchos::rcp(new vector<double> );
  RCP<vector<double> > x3coordinates;
  x3coordinates = Teuchos::rcp(new vector<double> );

  for(set<double,LineSortCriterion>::iterator coord1=x1coords.begin();
      coord1!=x1coords.end();
      ++coord1)
  {
    x1coordinates->push_back(*coord1);
    //std::cout << *coord1 << std::endl;
  }

  for(set<double,LineSortCriterion>::iterator coord2=x2coords.begin();
      coord2!=x2coords.end();
      ++coord2)
  {
    x2coordinates->push_back(*coord2);
    //std::cout << *coord2 << std::endl;
  }

  for(set<double,LineSortCriterion>::iterator coord3=x3coords.begin();
      coord3!=x3coords.end();
      ++coord3)
  {
    x3coordinates->push_back(*coord3);
    //std::cout << *coord3 << std::endl;
  }

  //----------------------------------------------------------------------
  // extract coarse node coordinates and push them in vectors
  //----------------------------------------------------------------------

  // start with first coordinate and take ever other
  RCP<vector<double> > x1coarsecoordinates;
  x1coarsecoordinates = Teuchos::rcp(new vector<double> );
  RCP<vector<double> > x2coarsecoordinates;
  x2coarsecoordinates = Teuchos::rcp(new vector<double> );
  RCP<vector<double> > x3coarsecoordinates;
  x3coarsecoordinates = Teuchos::rcp(new vector<double> );

  if ((x1coordinates->size()%2==0) or (x3coordinates->size()%2==0) or (x3coordinates->size()%2==0))
    dserror("Even number of elements expected!");

  for (size_t coords=0; coords<x1coordinates->size();coords++)
  {
    if(coords%2==0)
     x1coarsecoordinates->push_back(x1coordinates->at(coords));
  }

  for (size_t coords=0; coords<x2coordinates->size();coords++)
  {
    if(coords%2==0)
     x2coarsecoordinates->push_back(x2coordinates->at(coords));
  }

  for (size_t coords=0; coords<x3coordinates->size();coords++)
  {
    if(coords%2==0)
     x3coarsecoordinates->push_back(x3coordinates->at(coords));
  }

  //----------------------------------------------------------------------
  // build separation matrix
  //----------------------------------------------------------------------

  // separation matrix
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  RCP<Epetra_CrsMatrix> crsPRmat;
  crsPRmat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*dofrowmap,8,true));

  // loop all nodes on this proc
  for (int n=0;n<discret_->NumMyRowNodes();++n)
  {
    // get the node
    DRT::Node* node = discret_->lRowNode(n);

    // check pbc: the filtering is only done on the master side
    // check whether we have a pbc condition for this node
     vector<DRT::Condition*> mypbc;
     // lines are not considered at the moment as a 3d-domain are expected
     node->GetCondition("SurfacePeriodic",mypbc);

     bool isslave = false;
     // yes, we have surface pbcs for this nodes
     if (mypbc.size()>0)
     {
       // loop them and check, whether this is a pbc pure master node
       unsigned nummaster = 0;
       for (unsigned numcond=0;numcond<mypbc.size();++numcond)
       {
         const string* mymasterslavetoggle
         = mypbc[numcond]->Get<string>("Is slave periodic boundary condition");

         // yes, there is a master condition
         if(*mymasterslavetoggle=="Master")
         {
           ++nummaster;
         }
       }

       // do we have only master conditions, here?
       if(nummaster!=mypbc.size())
       {
         //no, this node is a slave
         isslave = true;
       }
     }

    if (isslave == false)
    {

    bool foundx = false;
    bool foundy = false;
    bool foundz = false;

    double x1coarse_min = -10e9;
    double x1coarse_max = 10e9;
    double x2coarse_min = -10e9;
    double x2coarse_max = 10e9;
    double x3coarse_min = -10e9;
    double x3coarse_max = 10e9;

    for (size_t coord=0; coord<x1coarsecoordinates->size(); coord++)
    {
      if ((node->X()[0] <= (x1coarsecoordinates->at(coord)+1e-9)) and (node->X()[0] >= (x1coarsecoordinates->at(coord)-1e-9)))
      {
        foundx = true;
        x1coarse_min = -10e9;
        x1coarse_max = -10e9;
      }
      else if (node->X()[0] > x1coarsecoordinates->at(coord))
      {
        if (x1coarsecoordinates->at(coord) > x1coarse_min)
        {
          x1coarse_min = x1coarsecoordinates->at(coord);
        }
      }
      else
      {
        if (x1coarsecoordinates->at(coord) < x1coarse_max)
        {
          x1coarse_max = x1coarsecoordinates->at(coord);
        }
      }
    }

    for (size_t coord=0; coord<x2coarsecoordinates->size(); coord++)
    {
      if ((node->X()[1] <= (x2coarsecoordinates->at(coord)+1e-9)) and (node->X()[1] >= (x2coarsecoordinates->at(coord)-1e-9)))
      {
        foundy = true;
        x2coarse_min = -10e9;
        x2coarse_max = -10e9;
      }
      else if (node->X()[1] > x2coarsecoordinates->at(coord))
      {
        if (x2coarsecoordinates->at(coord) > x2coarse_min)
        {
          x2coarse_min = x2coarsecoordinates->at(coord);
        }
      }
      else
      {
        if (x2coarsecoordinates->at(coord) < x2coarse_max)
        {
          x2coarse_max = x2coarsecoordinates->at(coord);
        }
      }
    }

    for (size_t coord=0; coord<x3coarsecoordinates->size(); coord++)
    {
      if ((node->X()[2] <= (x3coarsecoordinates->at(coord)+1e-9)) and (node->X()[2] >= (x3coarsecoordinates->at(coord)-1e-9)))
      {
        foundz = true;
        x3coarse_min = -10e9;
        x3coarse_max = -10e9;
      }
      else if (node->X()[2] > x3coarsecoordinates->at(coord))
      {
        if (x3coarsecoordinates->at(coord) > x3coarse_min)
        {
          x3coarse_min = x3coarsecoordinates->at(coord);
        }
      }
      else
      {
        if (x3coarsecoordinates->at(coord) < x3coarse_max)
        {
          x3coarse_max = x3coarsecoordinates->at(coord);
        }
      }
    }

    // vector to insert neighboring nodes
    std::vector<DRT::Node*> all_nodes;
    // vector to insert weights
    std::vector<double> weights;
    std::vector<int> colids;

    // calculate weights and insert them in vector
    if ((foundx == true) and (foundy == true) and (foundz == true))
    {
      // node is coarse node

      // weight is 1.0
      std::vector<double> funct;
      funct.push_back(1.0);
      // store node
      all_nodes.push_back(node);
      // store weight and id
      StoreInVector(all_nodes,funct,weights,colids);
    }
    else if (((foundx == true) and (foundy == true) and (foundz == false))
           or((foundx == true) and (foundy == false) and (foundz == true))
           or((foundx == false) and (foundy == true) and (foundz == true)))
    {
      // node on line between to coarse nodes

      // get coordinates of current node
      LINALG::Matrix<3,1> xyz;
      xyz(0,0) = node->X()[0];
      xyz(1,0) = node->X()[1];
      xyz(2,0) = node->X()[2];

      // interpolate between ...
      LINALG::Matrix<3,2> ele_xyz(true);
      if ((foundx == true) and (foundy == true) and (foundz == false))
      {
        // ... z-coordinates
        ele_xyz(0,0) = node->X()[0];
        ele_xyz(0,1) = node->X()[0];
        ele_xyz(1,0) = node->X()[1];
        ele_xyz(1,1) = node->X()[1];
        ele_xyz(2,0) = x3coarse_min;
        ele_xyz(2,1) = x3coarse_max;
      }
      else if ((foundx == true) and (foundy == false) and (foundz == true))
      {
        // ... y-coordinates
        ele_xyz(0,0) = node->X()[0];
        ele_xyz(0,1) = node->X()[0];
        ele_xyz(1,0) = x2coarse_min;
        ele_xyz(1,1) = x2coarse_max;
        ele_xyz(2,0) = node->X()[2];
        ele_xyz(2,1) = node->X()[2];
      }
      else if ((foundx == false) and (foundy == true) and (foundz == true))
      {
        // ... x-coordinates
        ele_xyz(0,0) = x1coarse_min;
        ele_xyz(0,1) = x1coarse_max;
        ele_xyz(1,0) = node->X()[1];
        ele_xyz(1,1) = node->X()[1];
        ele_xyz(2,0) = node->X()[2];
        ele_xyz(2,1) = node->X()[2];
      }
      else
       dserror("Something went wrong!");

      // transfer physical node coordinate to element coordinate
      LINALG::Matrix<1,1> eta;
      GEO::CurrentToLineElementCoordinates(DRT::Element::line2,ele_xyz,xyz,eta);

      // evaluate shape functions at current node
      LINALG::Matrix<2,1> shpfunct;
      DRT::UTILS::shape_function<DRT::Element::line2>(eta,shpfunct);

       // change format
      std::vector<double> funct;
      for (int nen=0; nen<2; nen++)
        funct.push_back(shpfunct(nen,0));

      // get nodes which contribute to current node and store this nodes
      all_nodes.push_back(FindNode(ele_xyz(0,0),ele_xyz(1,0),ele_xyz(2,0)));
      all_nodes.push_back(FindNode(ele_xyz(0,1),ele_xyz(1,1),ele_xyz(2,1)));

      // store weights and (col)ids of contributing nodes
      StoreInVector(all_nodes,funct,weights,colids);
    }
    else if (((foundx == true) and (foundy == false) and (foundz == false))
          or ((foundx == false) and (foundy == true) and (foundz == false))
          or ((foundx == false) and (foundy == false) and (foundz == true)))
    {
      // node on surface between to coarse nodes

      // get coordinates of current node
      LINALG::Matrix<3,1> xyz;
      xyz(0,0) = node->X()[0];
      xyz(1,0) = node->X()[1];
      xyz(2,0) = node->X()[2];

      // interpolate between ...
      LINALG::Matrix<3,4> ele_xyz(true);
      if ((foundx == false) and (foundy == false) and (foundz == true))
      {
        // ... xy-coordinates
        ele_xyz(0,0) = x1coarse_min;
        ele_xyz(1,0) = x2coarse_min;
        ele_xyz(2,0) = node->X()[2];
        ele_xyz(0,1) = x1coarse_max;
        ele_xyz(1,1) = x2coarse_min;
        ele_xyz(2,1) = node->X()[2];
        ele_xyz(0,2) = x1coarse_max;
        ele_xyz(1,2) = x2coarse_max;
        ele_xyz(2,2) = node->X()[2];
        ele_xyz(0,3) = x1coarse_min;
        ele_xyz(1,3) = x2coarse_max;
        ele_xyz(2,3) = node->X()[2];
      }
      else if ((foundx == false) and (foundy == true) and (foundz == false))
      {
        // ... zx-coordinates
        ele_xyz(0,0) = x1coarse_min;
        ele_xyz(1,0) = node->X()[1];
        ele_xyz(2,0) = x3coarse_min;
        ele_xyz(0,1) = x1coarse_min;
        ele_xyz(1,1) = node->X()[1];
        ele_xyz(2,1) = x3coarse_max;
        ele_xyz(0,2) = x1coarse_max;
        ele_xyz(1,2) = node->X()[1];
        ele_xyz(2,2) = x3coarse_max;
        ele_xyz(0,3) = x1coarse_max;
        ele_xyz(1,3) = node->X()[1];
        ele_xyz(2,3) = x3coarse_min;
      }
      else if ((foundx == true) and (foundy == false) and (foundz == false))
      {
        // ... yz-coordinates
        ele_xyz(0,0) = node->X()[0];
        ele_xyz(1,0) = x2coarse_min;
        ele_xyz(2,0) = x3coarse_min;
        ele_xyz(0,1) = node->X()[0];
        ele_xyz(1,1) = x2coarse_max;
        ele_xyz(2,1) = x3coarse_min;
        ele_xyz(0,2) = node->X()[0];
        ele_xyz(1,2) = x2coarse_max;
        ele_xyz(2,2) = x3coarse_max;
        ele_xyz(0,3) = node->X()[0];
        ele_xyz(1,3) = x2coarse_min;
        ele_xyz(2,3) = x3coarse_max;
      }
      else
       dserror("Something went wrong!");

      // transfer physical node coordinate to element coordinate
      LINALG::Matrix<2,1> eta;
      GEO::CurrentToSurfaceElementCoordinates(DRT::Element::quad4,ele_xyz,xyz,eta);

      // evaluate shape functions at current node
      LINALG::Matrix<4,1> shpfunct;
      DRT::UTILS::shape_function<DRT::Element::quad4>(eta,shpfunct);

      // change format
     std::vector<double> funct;
     for (int nen=0; nen<4; nen++)
       funct.push_back(shpfunct(nen,0));

     // get nodes which contribute to current node and store this nodes
     all_nodes.push_back(FindNode(ele_xyz(0,0),ele_xyz(1,0),ele_xyz(2,0)));
     all_nodes.push_back(FindNode(ele_xyz(0,1),ele_xyz(1,1),ele_xyz(2,1)));
     all_nodes.push_back(FindNode(ele_xyz(0,2),ele_xyz(1,2),ele_xyz(2,2)));
     all_nodes.push_back(FindNode(ele_xyz(0,3),ele_xyz(1,3),ele_xyz(2,3)));

     // store weights and (col)ids of contributing nodes
     StoreInVector(all_nodes,funct,weights,colids);
    }
    else if ((foundx == false) and (foundy == false) and (foundz == false))
    {
      // node in volume between to coarse nodes

      // get coordinates of current node
      LINALG::Matrix<3,1> xyz;
      xyz(0,0) = node->X()[0];
      xyz(1,0) = node->X()[1];
      xyz(2,0) = node->X()[2];

      // interpolate between ...
      LINALG::Matrix<3,8> ele_xyz(true);
      // ... xyz-coordinates
      ele_xyz(0,0) = x1coarse_min;
      ele_xyz(1,0) = x2coarse_min;
      ele_xyz(2,0) = x3coarse_max;

      ele_xyz(0,1) = x1coarse_max;
      ele_xyz(1,1) = x2coarse_min;
      ele_xyz(2,1) = x3coarse_max;;

      ele_xyz(0,2) = x1coarse_max;
      ele_xyz(1,2) = x2coarse_min;
      ele_xyz(2,2) = x3coarse_min;

      ele_xyz(0,3) = x1coarse_min;
      ele_xyz(1,3) = x2coarse_min;
      ele_xyz(2,3) = x3coarse_min;

      ele_xyz(0,4) = x1coarse_min;
      ele_xyz(1,4) = x2coarse_max;
      ele_xyz(2,4) = x3coarse_max;

      ele_xyz(0,5) = x1coarse_max;
      ele_xyz(1,5) = x2coarse_max;
      ele_xyz(2,5) = x3coarse_max;;

      ele_xyz(0,6) = x1coarse_max;
      ele_xyz(1,6) = x2coarse_max;
      ele_xyz(2,6) = x3coarse_min;

      ele_xyz(0,7) = x1coarse_min;
      ele_xyz(1,7) = x2coarse_max;
      ele_xyz(2,7) = x3coarse_min;

      // transfer physical node coordinate to element coordinate
      LINALG::Matrix<3,1> eta;
      //GEO::currentToVolumeElementCoordinates(DRT::Element::hex8,ele_xyz,xyz,eta);
      std::cout << "Warning: GEO::CUT::Position seems to work wrong on the cluster!" << std::endl;
      GEO::CUT::Position<DRT::Element::hex8>pos(ele_xyz,xyz);

      // evaluate shape functions at node
      LINALG::Matrix<8,1> shpfunct;
      DRT::UTILS::shape_function<DRT::Element::hex8>(eta,shpfunct);

      // change format
     std::vector<double> funct;
     for (int nen=0; nen<8; nen++)
       funct.push_back(shpfunct(nen,0));

     // get nodes which contribute to current node and store this nodes
     all_nodes.push_back(FindNode(ele_xyz(0,0),ele_xyz(1,0),ele_xyz(2,0)));
     all_nodes.push_back(FindNode(ele_xyz(0,1),ele_xyz(1,1),ele_xyz(2,1)));
     all_nodes.push_back(FindNode(ele_xyz(0,2),ele_xyz(1,2),ele_xyz(2,2)));
     all_nodes.push_back(FindNode(ele_xyz(0,3),ele_xyz(1,3),ele_xyz(2,3)));
     all_nodes.push_back(FindNode(ele_xyz(0,4),ele_xyz(1,4),ele_xyz(2,4)));
     all_nodes.push_back(FindNode(ele_xyz(0,5),ele_xyz(1,5),ele_xyz(2,5)));
     all_nodes.push_back(FindNode(ele_xyz(0,6),ele_xyz(1,6),ele_xyz(2,6)));
     all_nodes.push_back(FindNode(ele_xyz(0,7),ele_xyz(1,7),ele_xyz(2,7)));

     // store weights and (col)ids of contributing nodes
     StoreInVector(all_nodes,funct,weights,colids);
    }

    // store vector in matrix
    vector<int> dofs = discret_->Dof(node);
    // store value in vector
    for(int d=0;d<discret_->NumDof(node)-1;++d)
    {
      int id =dofs[d];

      std::vector<int> dofcolids;
      for (size_t numnodes=0; numnodes<all_nodes.size(); numnodes++)
      {
        vector<int> coldofs = discret_->Dof(all_nodes[numnodes]);
        int dofcolid = coldofs[d];
        dofcolids.push_back(dofcolid);
      }

      crsPRmat->InsertGlobalValues(id,weights.size(),&weights[0],&dofcolids[0]);
    }

  }

  }// end loop all nodes expect the ones with slave conditions

  crsPRmat->FillComplete();

  // store matrix
  Sep_ = Teuchos::rcp(new LINALG::SparseMatrix(crsPRmat));
  //complete scale-separation matrix and check maps
  Sep_->Complete(Sep_->DomainMap(),Sep_->RangeMap());

  sepmat_build_ = true;

  return;
}


void FLD::LESScaleSeparation::ApplyScaleSeparation(
  Teuchos::RCP<Epetra_Vector> vel,
  Teuchos::RCP<Epetra_Vector> fsvel)
{
  if (!sepmat_build_)
   dserror("Separation-matrix has not been build!");

  switch(scale_sep_)
  {
  case INPAR::FLUID::box_filter:
  {
    dserror("Not yet implemented!");
    break;
  }
  case INPAR::FLUID::algebraic_multigrid_operator:
  {
    dserror("Not yet implemented!");
    break;
  }
  case INPAR::FLUID::geometric_multigrid_operator:
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    RCP<Epetra_Vector> tmp;
    tmp = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));
    // get the coarse velocity
    Sep_->Multiply(false,*vel,*tmp);
    // and calculate the fine-scale velocity
    fsvel->Update(1.0,*tmp,0.0);
    fsvel->Update(1.0,*vel,-1.0);
    break;
  }
  default:
  {
    dserror("Unknown filter type for les!");
    break;
  }
  }

  return;
}


//---------------------------
// some useful functions
//---------------------------

DRT::Node* FLD::LESScaleSeparation::FindNode(
  double x,
  double y,
  double z)
{
  DRT::Node* actnode = 0;

  // loop all nodes on this proc
  for (int n=0;n<discret_->NumMyColNodes();++n)
  {
    // get the node
    DRT::Node* node = discret_->lColNode(n);

    if (((node->X()[0] <= (x+1E-9)) and (node->X()[0] >= (x-1E-9)))
    and ((node->X()[1] <= (y+1E-9)) and (node->X()[1] >= (y-1E-9)))
    and ((node->X()[2] <= (z+1E-9)) and (node->X()[2] >= (z-1E-9))))
    {
      actnode = node;
      break;
    }
  }

  if (actnode == 0)
   dserror("Could not find node!");

  return actnode;
}


void FLD::LESScaleSeparation::StoreInVector(
  const std::vector<DRT::Node*>& nodes,
  const std::vector<double>& funct,
  std::vector<double>& weights,
  std::vector<int>& colids)
{
  for (size_t n=0; n<nodes.size(); n++)
  {
     //same weight for all dofs of a node
      int id = nodes[n]->Id();
      weights.push_back(funct[n]);
      colids.push_back(id);
  }

  return;
}

