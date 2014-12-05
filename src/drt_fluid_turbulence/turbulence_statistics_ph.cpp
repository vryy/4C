/*!----------------------------------------------------------------------
\file turbulence_statistics_Ph.cpp

\brief calculate pressures, mean velocity values and fluctuations for
turbulent flow over a backward-facing step

<pre>
o Create sets for various evaluation lines in domain
  (Construction based on a round robin communication pattern):
  - 21 lines in x2-direction
  - lines along upper and lower wall

o loop nodes closest to / on lines

o values on lines are averaged in time over all steps between two
  outputs

Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#include "turbulence_statistics_ph.H"

//#define COMBINE_SAMPLES

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for lines

  o Allocate distributed vector for squares
  </pre>

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsPh::TurbulenceStatisticsPh(
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::ParameterList&          params):
  discret_      (actdis),
  params_       (params),
  x1statlocations_(true)
//  x2statlocations_(true)
//  inflowchannel_(DRT::INPUT::IntegralValue<int>(params_.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")),
//  inflowmax_    (params_.sublist("TURBULENT INFLOW").get<double>("INFLOW_CHA_SIDE",0.0))
{
  if (discret_->Comm().MyPID()==0)
  {
    std::cout << "This is the turbulence statistics manager of periodic hill problem" << std::endl;
    std::cout << "based on the geometry of ERCOFTAC" << std::endl;
    std::cout << "If additional output in front of the step is required, it has to be set manually (look for numx1supplocations_)." << std::endl;
  }

  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");

  // type of fluid flow solver: incompressible, Boussinesq approximation, varying density
//  const INPAR::FLUID::PhysicalType physicaltype = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type");


  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  squaredvelnp_ = LINALG::CreateVector(*dofrowmap,true);
  squaredscanp_ = LINALG::CreateVector(*dofrowmap,true);
  invscanp_ = LINALG::CreateVector(*dofrowmap,true);
  squaredinvscanp_ = LINALG::CreateVector(*dofrowmap,true);

  toggleu_      = LINALG::CreateVector(*dofrowmap,true);
  togglev_      = LINALG::CreateVector(*dofrowmap,true);
  togglew_      = LINALG::CreateVector(*dofrowmap,true);
  togglep_      = LINALG::CreateVector(*dofrowmap,true);

  // bounds for extension of flow domain in x2-direction
  x2min_ = +10e+19;
  x2max_ = -10e+19;
  // bounds for extension of flow domain in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates
  //----------------------------------------------------------------------
  x1coordinates_ = Teuchos::rcp(new std::vector<double> );
  x2coordinates_ = Teuchos::rcp(new std::vector<double> );

  // the criterion allows differences in coordinates by 1e-9
  std::set<double,LineSortCriterion> x1avcoords;
  std::set<double,LineSortCriterion> x2avcoords;
  std::set<double,LineSortCriterion> x2loccc;
  std::set<double,LineSortCriterion> x2statlocat;


  // loop nodes and build sets of lines in x1- and x2-direction
  // accessible on this proc
  // For x1-direction: consider horizontal line at x2=0
  // and assume no change in discretization behind the step
  // For x2-direction: consider vertical line at x1=0
  // and assume no change in discretization behind the step
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {

    DRT::Node* node = discret_->lRowNode(i);
//    std::cout << *node << std::endl;

    if (node->X()[1]<85.008+2e-9 && node->X()[1]>85.008-2e-9)
      x1avcoords.insert(node->X()[0]);

    if (node->X()[0]<2e-9 && node->X()[0]>-2e-9)
      x2avcoords.insert(node->X()[1]);

    // find mins and maxs

    if (x1min_>node->X()[0]) x1min_=node->X()[0];
    if (x1max_<node->X()[0]) x1max_=node->X()[0];


    if (x2min_>node->X()[1]) x2min_=node->X()[1];
    if (x2max_<node->X()[1]) x2max_=node->X()[1];


    if (x3min_>node->X()[2]) x3min_=node->X()[2];
    if (x3max_<node->X()[2]) x3max_=node->X()[2];
  }
//  std::cout << "x1avcoords: " << x1avcoords.size() << std::endl;
//  std::cout << "x2avcoords: " << x2avcoords.size() << std::endl;
  // communicate x1mins and x1maxs
  double min1;
  discret_->Comm().MinAll(&x1min_,&min1,1);
  x1min_=min1;

  double max1;
  discret_->Comm().MaxAll(&x1max_,&max1,1);
  x1max_=max1;

  // communicate x2mins and x2maxs
  double min2;
  discret_->Comm().MinAll(&x2min_,&min2,1);
  x2min_=min2;

  double max2;
  discret_->Comm().MaxAll(&x2max_,&max2,1);
  x2max_=max2;

  // communicate x3mins and x3maxs
  double min3;
  discret_->Comm().MinAll(&x3min_,&min3,1);
  x3min_=min3;

  double max3;
  discret_->Comm().MaxAll(&x3max_,&max3,1);
  x3max_=max3;
  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates to all procs
  //--------------------------------------------------------------------
  {
#ifdef PARALLEL
    int myrank  =discret_->Comm().MyPID();
#endif
    int numprocs=discret_->Comm().NumProc();

    std::vector<char> sblock;
    std::vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point communication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

//      for (std::set<double,LineSortCriterion>::iterator x1line=x1avcoords.begin();
//           x1line!=x1avcoords.end();
//           ++x1line)
//      {
//        DRT::ParObject::AddtoPack(data,*x1line);
//      }
      data.StartPacking();
      for (std::set<double,LineSortCriterion>::iterator x1line=x1avcoords.begin();
           x1line!=x1avcoords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      std::swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblock.size();

      exporter.ISend(frompid,topid,
                     &(sblock[0]),sblock.size(),
                     tag,request);

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
        std::vector<double> coordsvec;

        coordsvec.clear();

        std::vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x1avcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (std::set<double,LineSortCriterion>::iterator x2line=x2avcoords.begin();
           x2line!=x2avcoords.end();
           ++x2line)
      {
        DRT::ParObject::AddtoPack(data,*x2line);
      }
      data.StartPacking();
      for (std::set<double,LineSortCriterion>::iterator x2line=x2avcoords.begin();
           x2line!=x2avcoords.end();
           ++x2line)
      {
        DRT::ParObject::AddtoPack(data,*x2line);
      }
      std::swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblock.size();

      exporter.ISend(frompid,topid,
                     &(sblock[0]),sblock.size(),
                     tag,request);

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
        std::vector<double> coordsvec;

        coordsvec.clear();

        std::vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2avcoords.insert(onecoord);
        }
      }
    }
  }
  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------
  {
    for(std::set<double,LineSortCriterion>::iterator coord1=x1avcoords.begin();
        coord1!=x1avcoords.end();
        ++coord1)
    {
      x1coordinates_->push_back(*coord1);
//      std::cout << x1coordinates_ << std::endl;
    }

    for(std::set<double,LineSortCriterion>::iterator coord2=x2avcoords.begin();
        coord2!=x2avcoords.end();
        ++coord2)
    {
      x2coordinates_->push_back(*coord2);
//      std::cout << x2coordinates_ << std::endl;
    }
  }
  //----------------------------------------------------------------------
  // number of coordinates in x1- and x2-direction
  //----------------------------------------------------------------------
  numx1coor_ = x1coordinates_->size();
  numx2coor_ = x2coordinates_->size();

//  std::cout << numx1coor_ << std::endl;
//  std::cout << numx2coor_ << std::endl;

//  x1statlocations_(1,1)=4.5;
//  std::cout << "x1statlocations_(1,1): " << x1statlocations_(1,1) << std::endl;


  //----------------------------------------------------------------------
  // define locations in x1-direction for statistical evaluation
  //----------------------------------------------------------------------
  //define vector containing x1-coords for sampling

  x1setstatlocations_ = Teuchos::rcp(new std::vector<double> );

  //! x1-coordinates for statistical sampling
//  double x1setstatloc []={10,30,120,180};

  x1setstatlocations_->push_back(1.4);
  x1setstatlocations_->push_back(56);
  x1setstatlocations_->push_back(168);
  x1setstatlocations_->push_back(224);
//
//  for (unsigned x1stat=0; x1stat < 4; x1stat++)
//      x1setstatlocations_->push_back(x1setstatloc[x1stat]);

  numx1statlocations_ = x1setstatlocations_->size();
  std::cout << "numx1statlocations: " << numx1statlocations_ << std::endl;
//  x1dist_ = x1max_ - x1min_;
//  numx1elem_ = numx1coor_-1;

  x1elemlengthhalf = (x1max_ - x1min_)/(numx1coor_-1);
  for (int x1line =0 ; x1line < numx1statlocations_; ++x1line)
  {
    int pos = 0;
    mindist = x1elemlengthhalf;
    for (int x1linecoords = 0; x1linecoords < numx1coor_; ++x1linecoords)
    {
      dist = abs(x1setstatlocations_->at(x1line) - x1coordinates_->at(x1linecoords));
      if (dist < mindist)
      {
        mindist = dist;
        pos = x1linecoords;
      }
    }
    x1statlocations_(x1line,0) = x1coordinates_->at(pos);
    std::cout << "x1statlocations_(x1line,1): " << x1statlocations_(x1line,0) << std::endl;
  }

  //----------------------------------------------------------------------
  // define locations in x2-direction for statistical evaluation
  //----------------------------------------------------------------------
  //define vector containing x2-coords for sampling

  for (int x1stat = 0 ; x1stat < numx1statlocations_; ++x1stat)
  {
    x2statlocat.clear();

    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      if (node->X()[2] < +2e-9 && node->X()[2] > -2e-9
          && node->X()[0] < x1statlocations_(x1stat,0)+2e-9 && node->X()[0] > x1statlocations_(x1stat,0)-2e-9)
      {
        std::cout << "node X1: " << node->X()[1] << std::endl;
        x2statlocat.insert(node->X()[1]);
      }
    }

    for (std::set<double,LineSortCriterion>::iterator x2=x2statlocat.begin();
         x2!=x2statlocat.end();
         ++x2)
    {
      std::cout << " proc  " << discret_-> Comm().MyPID()<< "  " << *x2 << std::endl;
    }
    std::cout << "sizeofthisproc  "<< x2statlocat.size() << std::endl;
//////////////////////////////////////////////
//  communication of x2-ccordinates
//////////////////////////////////////////////
  {
#ifdef PARALLEL
    int myrank  =discret_->Comm().MyPID();
#endif
    int numprocs=discret_->Comm().NumProc();

    std::vector<char> sblock;
    std::vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point communication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (std::set<double,LineSortCriterion>::iterator x2=x2statlocat.begin();
           x2!=x2statlocat.end();
           ++x2)
      {
        DRT::ParObject::AddtoPack(data,*x2);
      }

      data.StartPacking();

      for (std::set<double,LineSortCriterion>::iterator x2=x2statlocat.begin();
           x2!=x2statlocat.end();
           ++x2)
      {
        DRT::ParObject::AddtoPack(data,*x2);
      }

      std::swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int         tag    =myrank;

      int         frompid=myrank;
      int         topid  =(myrank+1)%numprocs;

      int         length=sblock.size();

      exporter.ISend(frompid,topid,
                     &(sblock[0]),sblock.size(),
                     tag,request);

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
        std::cout <<  " proc "<< discret_-> Comm().MyPID() << std::endl;
        std::vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          std::cout << "onecoord: "  << onecoord << " proc "<< discret_-> Comm().MyPID() << std::endl;
          x2statlocat.insert(onecoord);
        }
//        int x2it = -1;
//        for (std::set<double,LineSortCriterion>::iterator x2locc=x2loccc.begin();
//             x2locc!=x2loccc.end();
//             ++x2locc)
//        {
//          x2it += 1;
//          //std::cout << "PID, x1stat, x2it, x2locc: " << discret_-> Comm().MyPID() << ", " << x1stat << ", " << x2it << ", " << *x2locc << std::endl;
//          x2statlocations(x1stat,x2it) = *x2locc;
//        }
      }

    }

    }
//  int x2it = -1;
//  for (std::set<double,LineSortCriterion>::iterator x2write=x2loccc.begin();
//       x2write!=x2loccc.end(); ++x2write)
//    {
//    x2it += 1;
//    x2statlocations(x1stat,x2it) = *x2write;
//    }
  int x2it = -1;
  for (std::set<double,LineSortCriterion>::iterator x2locc=x2statlocat.begin();
       x2locc!=x2statlocat.end();
       ++x2locc)
  {
    x2it += 1;
    //std::cout << "PID, x1stat, x2it, x2locc: " << discret_-> Comm().MyPID() << ", " << x1stat << ", " << x2it << ", " << *x2locc << std::endl;
    x2statlocations(x1stat,x2it) = *x2locc;
  }
std::cout << "x2statlocations: " << x2statlocations << std::endl;
  }


//----------------------------------------------------------------------
// define locations in x2-direction for statistical evaluation
//----------------------------------------------------------------------
///////////////////////////////////////////
//    end of communication in x2-direction
///////////////////////////////////////////


//  std::cout << x2statlocations_ << std::endl;
//  numx2statlocations_=17;

  //----------------------------------------------------------------------
  // define locations in x2-direction for statistical evaluation
  // (lower and upper wall)
  //----------------------------------------------------------------------
//  if (geotype_ == TurbulenceStatisticsPh::geometry_LES_flow_with_heating)
//  {
//    // num2statlocations_ also defines number of supplocations
//    numx2statlocations_ = 2;
//
//    x2statlocations_(0) = x2min_;
//    x2statlocations_(1) = x2max_;
//    //----------------------------------------------------------------------
//    // define supplementary locations in x2-direction for statistical
//    // evaluation of velocity derivative at wall
//    // (first nodes off lower and upper wall, respectively)
//    //----------------------------------------------------------------------
//    x2supplocations_(0) = x2coordinates_->at(1);
//    x2supplocations_(1) = x2coordinates_->at(x2coordinates_->size()-2);
//  }
//  else if (geotype_ == TurbulenceStatisticsPh::geometry_DNS_incomp_flow)
//  {
//    // num2statlocations_ also defines number of supplocations
//    numx2statlocations_ = 1;
//
//    x2statlocations_(0) = x2min_;
//    x2statlocations_(1) = x2max_; //not needed here, upper wall is slip wall
//    //----------------------------------------------------------------------
//    // define supplementary locations in x2-direction for statistical
//    // evaluation of velocity derivative at wall
//    // (first nodes off lower and upper wall, respectively)
//    //----------------------------------------------------------------------
//    x2supplocations_(0) = x2coordinates_->at(1);
//    x2supplocations_(1) = x2coordinates_->at(x2coordinates_->size()-2); //not needed here, upper wall is slip wall
//  }
//  else
//    dserror("Unknown geometry of backward facing step!");

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  // x1-direction
  x1sump_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sump_->Reshape(numx2statlocations_,numx1coor_);

  x1sumu_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumu_->Reshape(numx2statlocations_,numx1coor_);


  // x2-direction
  // first-order moments
  x2sumu_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumu_->Reshape(numx1statlocations_,numx2coor_);

  x2sumv_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumv_->Reshape(numx1statlocations_,numx2coor_);

  x2sumw_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumw_->Reshape(numx1statlocations_,numx2coor_);

  x2sump_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sump_->Reshape(numx1statlocations_,numx2coor_);

  // second-order moments
  x2sumsqu_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqu_->Reshape(numx1statlocations_,numx2coor_);

  x2sumsqv_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqv_->Reshape(numx1statlocations_,numx2coor_);

  x2sumsqw_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqw_->Reshape(numx1statlocations_,numx2coor_);

  x2sumsqp_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqp_->Reshape(numx1statlocations_,numx2coor_);

  x2sumuv_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumuv_->Reshape(numx1statlocations_,numx2coor_);

  x2sumuw_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumuw_->Reshape(numx1statlocations_,numx2coor_);

  x2sumvw_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumvw_->Reshape(numx1statlocations_,numx2coor_);

  // set number of samples to zero
  numsamp_ = 0;

  //----------------------------------------------------------------------
  // define homogeneous direction to compute averages of Smagorinsky constant

  Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));
  // check if we want to compute averages of Smagorinsky constant
  if (modelparams->get<std::string>("PHYSICAL_MODEL","no_model") == "Dynamic_Smagorinsky")
  {
    // store them in parameterlist for access on the element
    modelparams->set<Teuchos::RCP<std::vector<double> > >("dir1coords_",x1coordinates_);
    modelparams->set<Teuchos::RCP<std::vector<double> > >("dir2coords_",x2coordinates_);
  }

  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file

  Teuchos::RCP<std::ofstream> log;

//  std::cout << "Comm().MyPID: " << discret_->Comm().MyPID() << std::endl;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<std::string>("statistics outfile");
//    std::cout << s << std::endl;

      s.append(".flow_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
      (*log) << "# Statistics for turbulent incompressible flow over periodic hill (first- and second-order moments)\n\n";
//    }

    log->flush();
  }

  return;
}// TurbulenceStatisticsPh::TurbulenceStatisticsPh

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsPh::~TurbulenceStatisticsPh()
{
  return;
}// TurbulenceStatisticsPh::~TurbulenceStatisticsPh()


//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsPh::DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp)
{
  if(discret_->Comm().MyPID()==0)
    std::cout << "------------Time Sampling Routine begins---------" <<std::endl;
  // compute squared values of velocity
  squaredvelnp_->Multiply(1.0,*velnp,*velnp,0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

//  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
//  for (std::vector<double>::iterator x1line=x1coordinates_->begin();
//       x1line!=x1coordinates_->end();
//       ++x1line)
//  {
//    x1nodnum++;
//
//    for (int x2nodnum=0;x2nodnum<numx2statlocations_;++x2nodnum)
//    {
//      // current x2-coordinate of respective wall
//      double x2cwall = x2statlocations_(x2nodnum);
//
//      // current x2-coordinate of supplementary location to respective wall
//      double x2csupp = x2supplocations_(x2nodnum);
//
//      // toggle vectors are one in the position of a dof of this node,
//      // else 0
//      toggleu_->PutScalar(0.0);
//      togglep_->PutScalar(0.0);
//
//      // count the number of nodes in x3-direction contributing to this nodal value
//      int countnodes=0;
//
//      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
//      {
//        DRT::Node* node = discret_->lRowNode(nn);
//
//        // this is the wall node
//        if ((node->X()[0]<(*x1line+2e-9) and node->X()[0]>(*x1line-2e-9)) and
//            (node->X()[1]<(x2cwall+2e-5) and node->X()[1]>(x2cwall-2e-5)))
//        {
//          std::vector<int> dof = discret_->Dof(node);
//          double           one = 1.0;
//
//          togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));
//
//          countnodes++;
//        }
//        // this is the supplementary node
//        else if ((node->X()[0]<(*x1line+2e-9) and node->X()[0]>(*x1line-2e-9)) and
//                 (node->X()[1]<(x2csupp+2e-5) and node->X()[1]>(x2csupp-2e-5)))
//        {
//          std::vector<int> dof = discret_->Dof(node);
//          double           one = 1.0;
//
//          toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
//        }
//      }
//
//      int countnodesonallprocs=0;
//
//      discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);
//
//      // reduce by 1 due to periodic boundary condition
//      countnodesonallprocs-=1;
//
//      if (countnodesonallprocs)
//      {
//        //----------------------------------------------------------------------
//        // get values for velocity derivative and pressure
//        //----------------------------------------------------------------------
//        double u;
//        velnp->Dot(*toggleu_,&u);
//        double p;
//        velnp->Dot(*togglep_,&p);
//
//        //----------------------------------------------------------------------
//        // calculate spatial means
//        //----------------------------------------------------------------------
//        double usm=u/countnodesonallprocs;
//        double psm=p/countnodesonallprocs;
//
//        //----------------------------------------------------------------------
//        // add spatial mean values to statistical sample
//        //----------------------------------------------------------------------
//        (*x1sumu_)(x2nodnum,x1nodnum)+=usm;
//        (*x1sump_)(x2nodnum,x1nodnum)+=psm;
//      }
//    }
//  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------

  for (int x1nodnum=0;x1nodnum<numx1statlocations_;++x1nodnum)
  {
    // current x1-coordinate
//    double x1c = 1.0e20;
//    x1c = x1statlocations_(x1nodnum,1);

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction
    //----------------------------------------------------------------------
    for (int x2line = 0; x2line < numx2coor_; ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes=0;

      for (int nn=0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the node
        if (node->X()[0] < x1statlocations_(x1nodnum,0)+2e-5 && node->X()[0] > x1statlocations_(x1nodnum,0)-2e-5 &&
            node->X()[1] < x2statlocations(x1nodnum,x2line)+2e-9 && node->X()[1] > x2statlocations(x1nodnum,x2line)-2e-9 &&
            node->X()[2] < x3max_+2e-5)  //the last node in x3 direction is not sampled because of the periodic boundary condition
        {
          std::vector<int> dof = discret_->Dof(node);
          double           one = 1.0;

          toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
          togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
          togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
          togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs=0;

      discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs-=1;

//      std::cout << "countnodesonallprocs: " << countnodesonallprocs << std::endl;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity and pressure on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp->Dot(*toggleu_,&u);
        velnp->Dot(*togglev_,&v);
        velnp->Dot(*togglew_,&w);
        velnp->Dot(*togglep_,&p);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->Dot(*toggleu_,&uu);
        squaredvelnp_->Dot(*togglev_,&vv);
        squaredvelnp_->Dot(*togglew_,&ww);
        squaredvelnp_->Dot(*togglep_,&pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr=1;rr<velnp->MyLength();++rr)
        {
          locuv += ((*velnp)[rr-1]*(*toggleu_)[rr-1]) * ((*velnp)[rr]*(*togglev_)[rr]);
        }
        discret_->Comm().SumAll(&locuv,&uv,1);
        for (int rr=2;rr<velnp->MyLength();++rr)
        {
          locuw += ((*velnp)[rr-2]*(*toggleu_)[rr-2]) * ((*velnp)[rr]*(*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locuw,&uw,1);
        for (int rr=2;rr<velnp->MyLength();++rr)
        {
          locvw += ((*velnp)[rr-1]*(*togglev_)[rr-1]) * ((*velnp)[rr]*(*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locvw,&vw,1);

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum,x2nodnum)+=u/countnodesonallprocs;
        (*x2sumv_)(x1nodnum,x2nodnum)+=v/countnodesonallprocs;
        (*x2sumw_)(x1nodnum,x2nodnum)+=w/countnodesonallprocs;
        (*x2sump_)(x1nodnum,x2nodnum)+=p/countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum,x2nodnum)+=uu/countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum,x2nodnum)+=vv/countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum,x2nodnum)+=ww/countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum,x2nodnum)+=pp/countnodesonallprocs;

        (*x2sumuv_)(x1nodnum,x2nodnum)+=uv/countnodesonallprocs;
        (*x2sumuw_)(x1nodnum,x2nodnum)+=uw/countnodesonallprocs;
        (*x2sumvw_)(x1nodnum,x2nodnum)+=vw/countnodesonallprocs;
      }
    }
  }

  return;
}// TurbulenceStatisticsPh::DoTimeSample


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsPh::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<std::string>("statistics outfile");
    s.append(".flow_statistics");
    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow for a periodic hill (first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
//    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n\n\n";
//    (*log) << std::scientific;
//    (*log) << "\n\n\n";
//    (*log) << "# lower wall behind step\n";
//    (*log) << "#     x1";
//    (*log) << "           duxdy         pmean\n";
//    // distance from wall to first node off wall
//    double dist =  x2statlocations_(0);
//    for (unsigned i=0; i<x1coordinates_->size(); ++i)
//    {
////      std::cout <<"x1coordinates: "<< (*x1coordinates_)[i] << std:endl;
//      if (x1coordinates_->at(i) > 54.-2e-9)
//      {
//        double lwx1u     = (*x1sumu_)(17,i)/numsamp_;
//        double lwx1duxdy = lwx1u/dist;
//        double lwx1p     = (*x1sump_)(17,i)/numsamp_;
//        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << x1coordinates_->at(i);
//        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1duxdy;
//        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1p;
//        (*log) << "\n";
//      }
//    }

//    if (geotype_ == TurbulenceStatisticsPh::geometry_LES_flow_with_heating)
//    {
//      (*log) << "\n\n\n";
//      (*log) << "# upper wall\n";
//      (*log) << "#     x1";
//      (*log) << "           duxdy         pmean\n";
//
//      // distance from wall to first node off wall
//      dist = x2statlocations_(1) - x2supplocations_(1);
//
//      for (unsigned i=0; i<x1coordinates_->size(); ++i)
//      {
//        double uwx1u     = (*x1sumu_)(1,i)/numsamp_;
//        double uwx1duxdy = uwx1u/dist;
//        double uwx1p     = (*x1sump_)(1,i)/numsamp_;
//
//        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
//        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1duxdy;
//        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1p;
//        (*log) << "\n";
//      }
//    }
    for (int i=0; i<numx1statlocations_; ++i)
    {
      // current x1-coordinate
      double x1 = 1.0e20;
      x1 = x1statlocations_(i);

      (*log) << "\n\n\n";
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(4) << x1 << "\n";
      (*log) << "#     x2";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "         urms          vrms          wrms          prms";
      (*log) << "          u'v'          u'w'          v'w'\n";
      for (unsigned j=0; j<x2coordinates_->size(); ++j)
      {
        double x2u  = (*x2sumu_)(i,j)/numsamp_;
        double x2v  = (*x2sumv_)(i,j)/numsamp_;
        double x2w  = (*x2sumw_)(i,j)/numsamp_;
        double x2p  = (*x2sump_)(i,j)/numsamp_;
        double x2urms  = 0.0;
        double x2vrms  = 0.0;
        double x2wrms  = 0.0;
        double x2prms  = 0.0;
        if (((*x2sumsqu_)(i,j)/numsamp_-x2u*x2u) > 0.0)
          x2urms  = std::sqrt((*x2sumsqu_)(i,j)/numsamp_-x2u*x2u);
        if (((*x2sumsqv_)(i,j)/numsamp_-x2v*x2v) > 0.0)
          x2vrms  = std::sqrt((*x2sumsqv_)(i,j)/numsamp_-x2v*x2v);
        if (((*x2sumsqw_)(i,j)/numsamp_-x2w*x2w) > 0.0)
          x2wrms  = std::sqrt((*x2sumsqw_)(i,j)/numsamp_-x2w*x2w);
        if (((*x2sumsqp_)(i,j)/numsamp_-x2p*x2p) > 0.0)
          x2prms  = std::sqrt((*x2sumsqp_)(i,j)/numsamp_-x2p*x2p);
        double x2uv   = (*x2sumuv_)(i,j)/numsamp_-x2u*x2v;
        double x2uw   = (*x2sumuw_)(i,j)/numsamp_-x2u*x2w;
        double x2vw   = (*x2sumvw_)(i,j)/numsamp_-x2v*x2w;
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << x2statlocations(i,j);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2u;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2v;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2p;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2urms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vrms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2wrms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2prms;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uv;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uw;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vw;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsPh::DumpStatistics





