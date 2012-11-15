/*!----------------------------------------------------------------------
\file turbulence_statistics_ldc.cpp

\brief calculate mean values and fluctuations for turbulent flow in a
lid-driven cavity.

<pre>
o Create sets for centerlines in x1- and x2-direction
  (Construction based on a round robin communication pattern)

o loop nodes closest to centerlines

o values on centerlines are averaged in time over all steps between two
  outputs

Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#include "turbulence_statistics_ldc.H"

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for centerlines in x1- and x2-direction

  o Allocate distributed vector for squares
  </pre>

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsLdc::TurbulenceStatisticsLdc(
  RCP<DRT::Discretization> actdis,
  ParameterList&                   params)
  :
  discret_(actdis),
  params_ (params)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");

  INPAR::FLUID::PhysicalType physicaltype = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type");

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  toggleu_ = LINALG::CreateVector(*dofrowmap,true);
  togglev_ = LINALG::CreateVector(*dofrowmap,true);
  togglew_ = LINALG::CreateVector(*dofrowmap,true);
  togglep_ = LINALG::CreateVector(*dofrowmap,true);

  // bounds for extension of cavity in x3-direction
  x1min_ = +10e+19;
  x1max_ = -10e+19;
  // bounds for extension of cavity in x3-direction
  x2min_ = +10e+19;
  x2max_ = -10e+19;
  // bounds for extension of cavity in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates for centerlines in x1-, x2- and x3-direction
  //----------------------------------------------------------------------
  x1coordinates_ = Teuchos::rcp(new vector<double> );
  x2coordinates_ = Teuchos::rcp(new vector<double> );
  x3coordinates_ = Teuchos::rcp(new vector<double> );

  // the criterion allows differences in coordinates by 1e-9
  set<double,LineSortCriterion> x1avcoords;
  set<double,LineSortCriterion> x2avcoords;
  set<double,LineSortCriterion> x3avcoords;

  // loop nodes, build sets of centerlines accessible on this proc and
  // calculate extension of cavity in x3-direction
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);
    x1avcoords.insert(node->X()[0]);
    x2avcoords.insert(node->X()[1]);
    x3avcoords.insert(node->X()[2]);

    if (x1min_>node->X()[2])
    {
      x1min_=node->X()[2];
    }
    if (x1max_<node->X()[2])
    {
      x1max_=node->X()[2];
    }
    if (x2min_>node->X()[2])
    {
      x2min_=node->X()[2];
    }
    if (x2max_<node->X()[2])
    {
      x2max_=node->X()[2];
    }
    if (x3min_>node->X()[2])
    {
      x3min_=node->X()[2];
    }
    if (x3max_<node->X()[2])
    {
      x3max_=node->X()[2];
    }
  }

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
  // round robin loop to communicate coordinates in x1-, x2- and x3-
  // direction to all procs
  //--------------------------------------------------------------------
  {
#ifdef PARALLEL
    int myrank  =discret_->Comm().MyPID();
#endif
    int numprocs=discret_->Comm().NumProc();

    vector<char> sblock;
    vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (std::set<double,LineSortCriterion>::iterator x1line=x1avcoords.begin();
           x1line!=x1avcoords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      data.StartPacking();
      for (std::set<double,LineSortCriterion>::iterator x1line=x1avcoords.begin();
           x1line!=x1avcoords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      swap( sblock, data() );

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
        vector<double> coordsvec;

        coordsvec.clear();

        vector<char>::size_type index = 0;
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
      swap( sblock, data() );
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
        vector<double> coordsvec;

        coordsvec.clear();

        vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2avcoords.insert(onecoord);
        }
      }
    }

    // third, communicate coordinates in x3-direction
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (std::set<double,LineSortCriterion>::iterator x3line=x3avcoords.begin();
           x3line!=x3avcoords.end();
           ++x3line)
      {
        DRT::ParObject::AddtoPack(data,*x3line);
      }
      data.StartPacking();
      for (std::set<double,LineSortCriterion>::iterator x3line=x3avcoords.begin();
           x3line!=x3avcoords.end();
           ++x3line)
      {
        DRT::ParObject::AddtoPack(data,*x3line);
      }
      swap( sblock, data() );

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
        vector<double> coordsvec;

        coordsvec.clear();

        vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x3avcoords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in x1-, x2- and x3-direction in a vector
  //----------------------------------------------------------------------
  {
    x1coordinates_ = Teuchos::rcp(new vector<double> );
    x2coordinates_ = Teuchos::rcp(new vector<double> );
    x3coordinates_ = Teuchos::rcp(new vector<double> );

    for(std::set<double,LineSortCriterion>::iterator coord1=x1avcoords.begin();
        coord1!=x1avcoords.end();
        ++coord1)
    {
      x1coordinates_->push_back(*coord1);
    }

    for(std::set<double,LineSortCriterion>::iterator coord2=x2avcoords.begin();
        coord2!=x2avcoords.end();
        ++coord2)
    {
      x2coordinates_->push_back(*coord2);
    }

    for(std::set<double,LineSortCriterion>::iterator coord3=x3avcoords.begin();
        coord3!=x3avcoords.end();
        ++coord3)
    {
      x3coordinates_->push_back(*coord3);
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  int size1 = x1coordinates_->size();
  int size2 = x2coordinates_->size();
  int size3 = x3coordinates_->size();

  // first-order moments
  x1sumu_ =  Teuchos::rcp(new vector<double> );
  x1sumu_->resize(size1,0.0);
  x2sumu_ =  Teuchos::rcp(new vector<double> );
  x2sumu_->resize(size2,0.0);
  x3sumu_ =  Teuchos::rcp(new vector<double> );
  x3sumu_->resize(size3,0.0);

  x1sumv_ =  Teuchos::rcp(new vector<double> );
  x1sumv_->resize(size1,0.0);
  x2sumv_ =  Teuchos::rcp(new vector<double> );
  x2sumv_->resize(size2,0.0);
  x3sumv_ =  Teuchos::rcp(new vector<double> );
  x3sumv_->resize(size3,0.0);

  x1sumw_ =  Teuchos::rcp(new vector<double> );
  x1sumw_->resize(size1,0.0);
  x2sumw_ =  Teuchos::rcp(new vector<double> );
  x2sumw_->resize(size2,0.0);
  x3sumw_ =  Teuchos::rcp(new vector<double> );
  x3sumw_->resize(size3,0.0);

  x1sump_ =  Teuchos::rcp(new vector<double> );
  x1sump_->resize(size1,0.0);
  x2sump_ =  Teuchos::rcp(new vector<double> );
  x2sump_->resize(size2,0.0);
  x3sump_ =  Teuchos::rcp(new vector<double> );
  x3sump_->resize(size3,0.0);

  // second-order moments
  x1sumsqu_ =  Teuchos::rcp(new vector<double> );
  x1sumsqu_->resize(size1,0.0);
  x2sumsqu_ =  Teuchos::rcp(new vector<double> );
  x2sumsqu_->resize(size2,0.0);
  x3sumsqu_ =  Teuchos::rcp(new vector<double> );
  x3sumsqu_->resize(size3,0.0);

  x1sumsqv_ =  Teuchos::rcp(new vector<double> );
  x1sumsqv_->resize(size1,0.0);
  x2sumsqv_ =  Teuchos::rcp(new vector<double> );
  x2sumsqv_->resize(size2,0.0);
  x3sumsqv_ =  Teuchos::rcp(new vector<double> );
  x3sumsqv_->resize(size3,0.0);

  x1sumsqw_ =  Teuchos::rcp(new vector<double> );
  x1sumsqw_->resize(size1,0.0);
  x2sumsqw_ =  Teuchos::rcp(new vector<double> );
  x2sumsqw_->resize(size2,0.0);
  x3sumsqw_ =  Teuchos::rcp(new vector<double> );
  x3sumsqw_->resize(size3,0.0);

  x1sumsqp_ =  Teuchos::rcp(new vector<double> );
  x1sumsqp_->resize(size1,0.0);
  x2sumsqp_ =  Teuchos::rcp(new vector<double> );
  x2sumsqp_->resize(size2,0.0);
  x3sumsqp_ =  Teuchos::rcp(new vector<double> );
  x3sumsqp_->resize(size3,0.0);

  x1sumuv_ =  Teuchos::rcp(new vector<double> );
  x1sumuv_->resize(size1,0.0);
  x2sumuv_ =  Teuchos::rcp(new vector<double> );
  x2sumuv_->resize(size2,0.0);
  x3sumuv_ =  Teuchos::rcp(new vector<double> );
  x3sumuv_->resize(size3,0.0);

  x1sumuw_ =  Teuchos::rcp(new vector<double> );
  x1sumuw_->resize(size1,0.0);
  x2sumuw_ =  Teuchos::rcp(new vector<double> );
  x2sumuw_->resize(size2,0.0);
  x3sumuw_ =  Teuchos::rcp(new vector<double> );
  x3sumuw_->resize(size3,0.0);

  x1sumvw_ =  Teuchos::rcp(new vector<double> );
  x1sumvw_->resize(size1,0.0);
  x2sumvw_ =  Teuchos::rcp(new vector<double> );
  x2sumvw_->resize(size2,0.0);
  x3sumvw_ =  Teuchos::rcp(new vector<double> );
  x3sumvw_->resize(size3,0.0);

  // the following vectors are only necessary for low-Mach-number flow
  // first-order moments
  x1sumrho_ =  Teuchos::rcp(new vector<double> );
  x1sumrho_->resize(size1,0.0);
  x2sumrho_ =  Teuchos::rcp(new vector<double> );
  x2sumrho_->resize(size2,0.0);
  x3sumrho_ =  Teuchos::rcp(new vector<double> );
  x3sumrho_->resize(size3,0.0);

  x1sumT_ =  Teuchos::rcp(new vector<double> );
  x1sumT_->resize(size1,0.0);
  x2sumT_ =  Teuchos::rcp(new vector<double> );
  x2sumT_->resize(size2,0.0);
  x3sumT_ =  Teuchos::rcp(new vector<double> );
  x3sumT_->resize(size3,0.0);

  // second-order moments
  x1sumsqrho_ =  Teuchos::rcp(new vector<double> );
  x1sumsqrho_->resize(size1,0.0);
  x2sumsqrho_ =  Teuchos::rcp(new vector<double> );
  x2sumsqrho_->resize(size2,0.0);
  x3sumsqrho_ =  Teuchos::rcp(new vector<double> );
  x3sumsqrho_->resize(size3,0.0);

  x1sumsqT_ =  Teuchos::rcp(new vector<double> );
  x1sumsqT_->resize(size1,0.0);
  x2sumsqT_ =  Teuchos::rcp(new vector<double> );
  x2sumsqT_->resize(size2,0.0);
  x3sumsqT_ =  Teuchos::rcp(new vector<double> );
  x3sumsqT_->resize(size3,0.0);

  x1sumuT_ =  Teuchos::rcp(new vector<double> );
  x1sumuT_->resize(size1,0.0);
  x2sumuT_ =  Teuchos::rcp(new vector<double> );
  x2sumuT_->resize(size2,0.0);
  x3sumuT_ =  Teuchos::rcp(new vector<double> );
  x3sumuT_->resize(size3,0.0);

  x1sumvT_ =  Teuchos::rcp(new vector<double> );
  x1sumvT_->resize(size1,0.0);
  x2sumvT_ =  Teuchos::rcp(new vector<double> );
  x2sumvT_->resize(size2,0.0);
  x3sumvT_ =  Teuchos::rcp(new vector<double> );
  x3sumvT_->resize(size3,0.0);

  x1sumwT_ =  Teuchos::rcp(new vector<double> );
  x1sumwT_->resize(size1,0.0);
  x2sumwT_ =  Teuchos::rcp(new vector<double> );
  x2sumwT_->resize(size2,0.0);
  x3sumwT_ =  Teuchos::rcp(new vector<double> );
  x3sumwT_->resize(size3,0.0);

  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file

  Teuchos::RCP<std::ofstream> log;

  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");

    if (physicaltype == INPAR::FLUID::loma)
    {
      s.append(".loma_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
      (*log) << "# Statistics for turbulent variable-density flow in a lid-driven cavity at low Mach number (first- and second-order moments)\n\n";
    }
    else
    {
      s.append(".flow_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
      (*log) << "# Statistics for turbulent incompressible flow in a lid-driven cavity (first- and second-order moments)\n\n";
    }

    log->flush();
  }

  // clear statistics
  ClearStatistics();

  return;
}// TurbulenceStatisticsLdc::TurbulenceStatisticsLdc

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsLdc::~TurbulenceStatisticsLdc()
{
  return;
}// TurbulenceStatisticsLdc::~TurbulenceStatisticsLdc()

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsLdc::DoTimeSample(
Teuchos::RCP<Epetra_Vector> velnp
  )
{

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x1-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x1line=x1coordinates_->begin();
       x1line!=x1coordinates_->end();
       ++x1line)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on x1-centerline
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[0]<(*x1line+2e-9) && node->X()[0]>(*x1line-2e-9)) &&
          (node->X()[1]<((x2max_+x2min_)/2.0)+2e-9 &&
           node->X()[1]>((x2max_+x2min_)/2.0)-2e-9) &&
          (node->X()[2]<((x3max_+x3min_)/2.0)+2e-9 &&
           node->X()[2]>((x3max_+x3min_)/2.0)-2e-9))
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        countnodes++;
      }
    }

    int countnodesonallprocs=0;

    discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

    if (countnodesonallprocs)
    {
      //----------------------------------------------------------------------
      // get values for velocity and pressure on this centerline
      //----------------------------------------------------------------------
      double u;
      double v;
      double w;
      double p;
      velnp->Dot(*toggleu_,&u);
      velnp->Dot(*togglev_,&v);
      velnp->Dot(*togglew_,&w);
      velnp->Dot(*togglep_,&p);

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this centerline
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x1sumu_)[x1nodnum]+=usm;
      (*x1sumv_)[x1nodnum]+=vsm;
      (*x1sumw_)[x1nodnum]+=wsm;
      (*x1sump_)[x1nodnum]+=psm;

      (*x1sumsqu_)[x1nodnum]+=usm*usm;
      (*x1sumsqv_)[x1nodnum]+=vsm*vsm;
      (*x1sumsqw_)[x1nodnum]+=wsm*wsm;
      (*x1sumsqp_)[x1nodnum]+=psm*psm;

      (*x1sumuv_)[x1nodnum]+=usm*vsm;
      (*x1sumuw_)[x1nodnum]+=usm*wsm;
      (*x1sumvw_)[x1nodnum]+=vsm*wsm;
    }
    x1nodnum++;

  }

  int x2nodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x2-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x2line=x2coordinates_->begin();
       x2line!=x2coordinates_->end();
       ++x2line)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on x2-centerline
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
     if ((node->X()[1]<(*x2line+2e-9) && node->X()[1]>(*x2line-2e-9)) &&
         (node->X()[0]<((x1max_+x1min_)/2.0)+2e-9 &&
          node->X()[0]>((x1max_+x1min_)/2.0)-2e-9) &&
         (node->X()[2]<((x3max_+x3min_)/2.0)+2e-9 &&
          node->X()[2]>((x3max_+x3min_)/2.0)-2e-9))
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        countnodes++;
      }
    }

    int countnodesonallprocs=0;

    discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

    if (countnodesonallprocs)
    {
      //----------------------------------------------------------------------
      // get values for velocity and pressure on this centerline
      //----------------------------------------------------------------------
      double u;
      double v;
      double w;
      double p;
      velnp->Dot(*toggleu_,&u);
      velnp->Dot(*togglev_,&v);
      velnp->Dot(*togglew_,&w);
      velnp->Dot(*togglep_,&p);

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this centerline
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x2sumu_)[x2nodnum]+=usm;
      (*x2sumv_)[x2nodnum]+=vsm;
      (*x2sumw_)[x2nodnum]+=wsm;
      (*x2sump_)[x2nodnum]+=psm;

      (*x2sumsqu_)[x2nodnum]+=usm*usm;
      (*x2sumsqv_)[x2nodnum]+=vsm*vsm;
      (*x2sumsqw_)[x2nodnum]+=wsm*wsm;
      (*x2sumsqp_)[x2nodnum]+=psm*psm;

      (*x2sumuv_)[x2nodnum]+=usm*vsm;
      (*x2sumuw_)[x2nodnum]+=usm*wsm;
      (*x2sumvw_)[x2nodnum]+=vsm*wsm;
    }
    x2nodnum++;

  }

  int x3nodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x3-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x3line=x3coordinates_->begin();
       x3line!=x3coordinates_->end();
       ++x3line)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on x3-centerline
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x3-direction
      if ((node->X()[2]<(*x3line+2e-9) && node->X()[2]>(*x3line-2e-9)) &&
          (node->X()[0]<((x1max_+x1min_)/2.0)+2e-9 &&
           node->X()[0]>((x1max_+x1min_)/2.0)-2e-9) &&
          (node->X()[1]<((x2max_+x2min_)/2.0)+2e-9 &&
           node->X()[1]>((x2max_+x2min_)/2.0)-2e-9))
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        countnodes++;
      }
    }

    int countnodesonallprocs=0;

    discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

    if (countnodesonallprocs)
    {
      //----------------------------------------------------------------------
      // get values for velocity and pressure on this centerline
      //----------------------------------------------------------------------
      double u;
      double v;
      double w;
      double p;
      velnp->Dot(*toggleu_,&u);
      velnp->Dot(*togglev_,&v);
      velnp->Dot(*togglew_,&w);
      velnp->Dot(*togglep_,&p);

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this centerline
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x3sumu_)[x3nodnum]+=usm;
      (*x3sumv_)[x3nodnum]+=vsm;
      (*x3sumw_)[x3nodnum]+=wsm;
      (*x3sump_)[x3nodnum]+=psm;

      (*x3sumsqu_)[x3nodnum]+=usm*usm;
      (*x3sumsqv_)[x3nodnum]+=vsm*vsm;
      (*x3sumsqw_)[x3nodnum]+=wsm*wsm;
      (*x3sumsqp_)[x3nodnum]+=psm*psm;

      (*x3sumuv_)[x3nodnum]+=usm*vsm;
      (*x3sumuw_)[x3nodnum]+=usm*wsm;
      (*x3sumvw_)[x3nodnum]+=vsm*wsm;
    }
    x3nodnum++;

  }

  return;
}// TurbulenceStatisticsLdc::DoTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsLdc::DoLomaTimeSample(
Teuchos::RCP<Epetra_Vector> velnp,
Teuchos::RCP<Epetra_Vector> scanp,
Epetra_Vector &                     force,
const double                        eosfac)
{
  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x1-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x1line=x1coordinates_->begin();
       x1line!=x1coordinates_->end();
       ++x1line)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on x1-centerline
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[0]<(*x1line+2e-9) && node->X()[0]>(*x1line-2e-9)) &&
          (node->X()[1]<(0.5+2e-9) && node->X()[1]>(0.5-2e-9)) &&
          (node->X()[2]<((x3max_-x3min_)/2.0)+2e-9 &&
           node->X()[2]>((x3max_-x3min_)/2.0)-2e-9))
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        countnodes++;
      }
    }

    int countnodesonallprocs=0;

    discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

    if (countnodesonallprocs)
    {
      //----------------------------------------------------------------------
      // get values for velocity, pressure, density and temperature
      // on this centerline
      //----------------------------------------------------------------------
      double u;
      double v;
      double w;
      double p;
      velnp->Dot(*toggleu_,&u);
      velnp->Dot(*togglev_,&v);
      velnp->Dot(*togglew_,&w);
      velnp->Dot(*togglep_,&p);

      double T;
      scanp->Dot(*togglep_,&T);

      //----------------------------------------------------------------------
      // calculate spatial means for vel., press., dens. and temp.
      // on this centerline
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;
      double Tsm=T/countnodesonallprocs;
      // compute density: rho = eosfac/T
      double rsm = eosfac/Tsm;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x1sumu_)[x1nodnum]  +=usm;
      (*x1sumv_)[x1nodnum]  +=vsm;
      (*x1sumw_)[x1nodnum]  +=wsm;
      (*x1sump_)[x1nodnum]  +=psm;
      (*x1sumrho_)[x1nodnum]+=rsm;
      (*x1sumT_)[x1nodnum]  +=Tsm;

      (*x1sumsqu_)[x1nodnum]  +=usm*usm;
      (*x1sumsqv_)[x1nodnum]  +=vsm*vsm;
      (*x1sumsqw_)[x1nodnum]  +=wsm*wsm;
      (*x1sumsqp_)[x1nodnum]  +=psm*psm;
      (*x1sumsqrho_)[x1nodnum]+=rsm*rsm;
      (*x1sumsqT_)[x1nodnum]  +=Tsm*Tsm;

      (*x1sumuv_)[x1nodnum]+=usm*vsm;
      (*x1sumuw_)[x1nodnum]+=usm*wsm;
      (*x1sumvw_)[x1nodnum]+=vsm*wsm;
      (*x1sumuT_)[x1nodnum]+=usm*Tsm;
      (*x1sumvT_)[x1nodnum]+=vsm*Tsm;
      (*x1sumwT_)[x1nodnum]+=wsm*Tsm;
    }
    x1nodnum++;

  }

  int x2nodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x2-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x2line=x2coordinates_->begin();
       x2line!=x2coordinates_->end();
       ++x2line)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on x2-centerline
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
      if ((node->X()[1]<(*x2line+2e-9) && node->X()[1]>(*x2line-2e-9)) &&
          (node->X()[0]<(0.5+2e-9) && node->X()[0]>(0.5-2e-9)) &&
          (node->X()[2]<((x3max_-x3min_)/2.0)+2e-9 &&
           node->X()[2]>((x3max_-x3min_)/2.0)-2e-9))
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        countnodes++;
      }
    }

    int countnodesonallprocs=0;

    discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

    if (countnodesonallprocs)
    {
      //----------------------------------------------------------------------
      // get values for velocity, pressure, density and temperature
      // on this centerline
      //----------------------------------------------------------------------
      double u;
      double v;
      double w;
      double p;
      velnp->Dot(*toggleu_,&u);
      velnp->Dot(*togglev_,&v);
      velnp->Dot(*togglew_,&w);
      velnp->Dot(*togglep_,&p);

      double T;
      scanp->Dot(*togglep_,&T);

      //----------------------------------------------------------------------
      // calculate spatial means for vel., press., dens. and temp.
      // on this centerline
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;
      double Tsm=T/countnodesonallprocs;
      // compute density: rho = eosfac/T
      double rsm = eosfac/Tsm;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x2sumu_)[x2nodnum]  +=usm;
      (*x2sumv_)[x2nodnum]  +=vsm;
      (*x2sumw_)[x2nodnum]  +=wsm;
      (*x2sump_)[x2nodnum]  +=psm;
      (*x2sumrho_)[x2nodnum]+=rsm;
      (*x2sumT_)[x2nodnum]  +=Tsm;

      (*x2sumsqu_)[x2nodnum]  +=usm*usm;
      (*x2sumsqv_)[x2nodnum]  +=vsm*vsm;
      (*x2sumsqw_)[x2nodnum]  +=wsm*wsm;
      (*x2sumsqp_)[x2nodnum]  +=psm*psm;
      (*x2sumsqrho_)[x2nodnum]+=rsm*rsm;
      (*x2sumsqT_)[x2nodnum]  +=Tsm*Tsm;

      (*x2sumuv_)[x2nodnum]+=usm*vsm;
      (*x2sumuw_)[x2nodnum]+=usm*wsm;
      (*x2sumvw_)[x2nodnum]+=vsm*wsm;
      (*x2sumuT_)[x2nodnum]+=usm*Tsm;
      (*x2sumvT_)[x2nodnum]+=vsm*Tsm;
      (*x2sumwT_)[x2nodnum]+=wsm*Tsm;
    }
    x2nodnum++;

  }

  int x3nodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x3-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x3line=x3coordinates_->begin();
       x3line!=x3coordinates_->end();
       ++x3line)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on x3-centerline
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x3-direction
      if ((node->X()[2]<(*x3line+2e-9) && node->X()[2]>(*x3line-2e-9)) &&
          (node->X()[0]<(0.5+2e-9) && node->X()[0]>(0.5-2e-9)) &&
          (node->X()[1]<(0.5+2e-9) && node->X()[1]>(0.5-2e-9)))
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        countnodes++;
      }
    }

    int countnodesonallprocs=0;

    discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

    if (countnodesonallprocs)
    {
      //----------------------------------------------------------------------
      // get values for velocity, pressure, density and temperature
      // on this centerline
      //----------------------------------------------------------------------
      double u;
      double v;
      double w;
      double p;
      velnp->Dot(*toggleu_,&u);
      velnp->Dot(*togglev_,&v);
      velnp->Dot(*togglew_,&w);
      velnp->Dot(*togglep_,&p);

      double T;
      scanp->Dot(*togglep_,&T);

      //----------------------------------------------------------------------
      // calculate spatial means for vel., press., dens. and temp.
      // on this centerline
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;
      double Tsm=T/countnodesonallprocs;
      // compute density: rho = eosfac/T
      double rsm = eosfac/Tsm;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x3sumu_)[x3nodnum]  +=usm;
      (*x3sumv_)[x3nodnum]  +=vsm;
      (*x3sumw_)[x3nodnum]  +=wsm;
      (*x3sump_)[x3nodnum]  +=psm;
      (*x3sumrho_)[x3nodnum]+=rsm;
      (*x3sumT_)[x3nodnum]  +=Tsm;

      (*x3sumsqu_)[x3nodnum]  +=usm*usm;
      (*x3sumsqv_)[x3nodnum]  +=vsm*vsm;
      (*x3sumsqw_)[x3nodnum]  +=wsm*wsm;
      (*x3sumsqp_)[x3nodnum]  +=psm*psm;
      (*x3sumsqrho_)[x3nodnum]+=rsm*rsm;
      (*x3sumsqT_)[x3nodnum]  +=Tsm*Tsm;

      (*x3sumuv_)[x3nodnum]+=usm*vsm;
      (*x3sumuw_)[x3nodnum]+=usm*wsm;
      (*x3sumvw_)[x3nodnum]+=vsm*wsm;
      (*x3sumuT_)[x3nodnum]+=usm*Tsm;
      (*x3sumvT_)[x3nodnum]+=vsm*Tsm;
      (*x3sumwT_)[x3nodnum]+=wsm*Tsm;
    }
    x3nodnum++;

  }

  return;
}// TurbulenceStatisticsLdc::DoLomaTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsLdc::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow in a lid-driven cavity (first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "          u'v'          u'w'          v'w'          prms\n";

    (*log) << std::scientific;
    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      double x1u    = (*x1sumu_)[i]/numsamp_;
      double x1v    = (*x1sumv_)[i]/numsamp_;
      double x1w    = (*x1sumw_)[i]/numsamp_;
      double x1p    = (*x1sump_)[i]/numsamp_;

      double x1urms = sqrt((*x1sumsqu_)[i]/numsamp_-x1u*x1u);
      double x1vrms = sqrt((*x1sumsqv_)[i]/numsamp_-x1v*x1v);
      double x1wrms = sqrt((*x1sumsqw_)[i]/numsamp_-x1w*x1w);
      double x1prms = sqrt((*x1sumsqp_)[i]/numsamp_-x1p*x1p);

      double x1uv   = (*x1sumuv_)[i]/numsamp_-x1u*x1v;
      double x1uw   = (*x1sumuw_)[i]/numsamp_-x1u*x1w;
      double x1vw   = (*x1sumvw_)[i]/numsamp_-x1v*x1w;

      (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1prms;
      (*log) << "\n";
    }

    (*log) << "\n\n\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "          u'v'          u'w'          v'w'          prms\n";

    (*log) << std::scientific;
    for(unsigned i=0; i<x2coordinates_->size(); ++i)
    {
      double x2u    = (*x2sumu_)[i]/numsamp_;
      double x2v    = (*x2sumv_)[i]/numsamp_;
      double x2w    = (*x2sumw_)[i]/numsamp_;
      double x2p    = (*x2sump_)[i]/numsamp_;

      double x2urms = sqrt((*x2sumsqu_)[i]/numsamp_-x2u*x2u);
      double x2vrms = sqrt((*x2sumsqv_)[i]/numsamp_-x2v*x2v);
      double x2wrms = sqrt((*x2sumsqw_)[i]/numsamp_-x2w*x2w);
      double x2prms = sqrt((*x2sumsqp_)[i]/numsamp_-x2p*x2p);

      double x2uv   = (*x2sumuv_)[i]/numsamp_-x2u*x2v;
      double x2uw   = (*x2sumuw_)[i]/numsamp_-x2u*x2w;
      double x2vw   = (*x2sumvw_)[i]/numsamp_-x2v*x2w;

      (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x2coordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2prms;
      (*log) << "\n";
    }

    (*log) << "\n\n\n";
    (*log) << "#     x3";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "          u'v'          u'w'          v'w'          prms\n";

    (*log) << std::scientific;
    for(unsigned i=0; i<x3coordinates_->size(); ++i)
    {
      double x3u    = (*x3sumu_)[i]/numsamp_;
      double x3v    = (*x3sumv_)[i]/numsamp_;
      double x3w    = (*x3sumw_)[i]/numsamp_;
      double x3p    = (*x3sump_)[i]/numsamp_;

      double x3urms = sqrt((*x3sumsqu_)[i]/numsamp_-x3u*x3u);
      double x3vrms = sqrt((*x3sumsqv_)[i]/numsamp_-x3v*x3v);
      double x3wrms = sqrt((*x3sumsqw_)[i]/numsamp_-x3w*x3w);
      double x3prms = sqrt((*x3sumsqp_)[i]/numsamp_-x3p*x3p);

      double x3uv   = (*x3sumuv_)[i]/numsamp_-x3u*x3v;
      double x3uw   = (*x3sumuw_)[i]/numsamp_-x3u*x3w;
      double x3vw   = (*x3sumvw_)[i]/numsamp_-x3v*x3w;

      (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x3coordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3prms;
      (*log) << "\n";
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsLdc::DumpStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsLdc::DumpLomaStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".loma_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
    (*log) << "# Statistics for turbulent variable-density flow in a lid-driven cavity at low Mach number (first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean       rhomean         Tmean";
    (*log) << "         urms          vrms          wrms          prms        rhorms          Trms";
    (*log) << "          u'v'          u'w'          v'w'          u'T'          v'T'          w'T'\n";

    (*log) << std::scientific;
    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      double x1u    = (*x1sumu_)[i]/numsamp_;
      double x1v    = (*x1sumv_)[i]/numsamp_;
      double x1w    = (*x1sumw_)[i]/numsamp_;
      double x1p    = (*x1sump_)[i]/numsamp_;
      double x1rho  = (*x1sumrho_)[i]/numsamp_;
      double x1T    = (*x1sumT_)[i]/numsamp_;

      double x1urms   = sqrt((*x1sumsqu_)[i]/numsamp_-x1u*x1u);
      double x1vrms   = sqrt((*x1sumsqv_)[i]/numsamp_-x1v*x1v);
      double x1wrms   = sqrt((*x1sumsqw_)[i]/numsamp_-x1w*x1w);
      double x1prms   = sqrt((*x1sumsqp_)[i]/numsamp_-x1p*x1p);
      double x1rhorms = sqrt((*x1sumsqrho_)[i]/numsamp_-x1rho*x1rho);
      double x1Trms   = sqrt((*x1sumsqT_)[i]/numsamp_-x1T*x1T);

      double x1uv   = (*x1sumuv_)[i]/numsamp_-x1u*x1v;
      double x1uw   = (*x1sumuw_)[i]/numsamp_-x1u*x1w;
      double x1vw   = (*x1sumvw_)[i]/numsamp_-x1v*x1w;
      double x1uT   = (*x1sumuT_)[i]/numsamp_-x1u*x1T;
      double x1vT   = (*x1sumvT_)[i]/numsamp_-x1v*x1T;
      double x1wT   = (*x1sumwT_)[i]/numsamp_-x1w*x1T;

      (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1rho;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1T;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1prms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1rhorms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1Trms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1uT;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1vT;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1wT;
      (*log) << "\n";
    }

    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean       rhomean         Tmean";
    (*log) << "         urms          vrms          wrms          prms        rhorms          Trms";
    (*log) << "          u'v'          u'w'          v'w'          u'T'          v'T'          w'T'\n";

    (*log) << std::scientific;
    for(unsigned i=0; i<x2coordinates_->size(); ++i)
    {
      double x2u    = (*x2sumu_)[i]/numsamp_;
      double x2v    = (*x2sumv_)[i]/numsamp_;
      double x2w    = (*x2sumw_)[i]/numsamp_;
      double x2p    = (*x2sump_)[i]/numsamp_;
      double x2rho  = (*x2sumrho_)[i]/numsamp_;
      double x2T    = (*x2sumT_)[i]/numsamp_;

      double x2urms   = sqrt((*x2sumsqu_)[i]/numsamp_-x2u*x2u);
      double x2vrms   = sqrt((*x2sumsqv_)[i]/numsamp_-x2v*x2v);
      double x2wrms   = sqrt((*x2sumsqw_)[i]/numsamp_-x2w*x2w);
      double x2prms   = sqrt((*x2sumsqp_)[i]/numsamp_-x2p*x2p);
      double x2rhorms = sqrt((*x2sumsqrho_)[i]/numsamp_-x2rho*x2rho);
      double x2Trms   = sqrt((*x2sumsqT_)[i]/numsamp_-x2T*x2T);

      double x2uv   = (*x2sumuv_)[i]/numsamp_-x2u*x2v;
      double x2uw   = (*x2sumuw_)[i]/numsamp_-x2u*x2w;
      double x2vw   = (*x2sumvw_)[i]/numsamp_-x2v*x2w;
      double x2uT   = (*x2sumuT_)[i]/numsamp_-x2u*x2T;
      double x2vT   = (*x2sumvT_)[i]/numsamp_-x2v*x2T;
      double x2wT   = (*x2sumwT_)[i]/numsamp_-x2w*x2T;

      (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x2coordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2rho;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2T;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2prms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2rhorms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2Trms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2uT;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2vT;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2wT;
      (*log) << "\n";
    }

    (*log) << "#     x3";
    (*log) << "           umean         vmean         wmean         pmean       rhomean         Tmean";
    (*log) << "         urms          vrms          wrms          prms        rhorms          Trms";
    (*log) << "          u'v'          u'w'          v'w'          u'T'          v'T'          w'T'\n";

    (*log) << std::scientific;
    for(unsigned i=0; i<x3coordinates_->size(); ++i)
    {
      double x3u    = (*x3sumu_)[i]/numsamp_;
      double x3v    = (*x3sumv_)[i]/numsamp_;
      double x3w    = (*x3sumw_)[i]/numsamp_;
      double x3p    = (*x3sump_)[i]/numsamp_;
      double x3rho  = (*x3sumrho_)[i]/numsamp_;
      double x3T    = (*x3sumT_)[i]/numsamp_;

      double x3urms   = sqrt((*x3sumsqu_)[i]/numsamp_-x3u*x3u);
      double x3vrms   = sqrt((*x3sumsqv_)[i]/numsamp_-x3v*x3v);
      double x3wrms   = sqrt((*x3sumsqw_)[i]/numsamp_-x3w*x3w);
      double x3prms   = sqrt((*x3sumsqp_)[i]/numsamp_-x3p*x3p);
      double x3rhorms = sqrt((*x3sumsqrho_)[i]/numsamp_-x3rho*x3rho);
      double x3Trms   = sqrt((*x3sumsqT_)[i]/numsamp_-x3T*x3T);

      double x3uv   = (*x3sumuv_)[i]/numsamp_-x3u*x3v;
      double x3uw   = (*x3sumuw_)[i]/numsamp_-x3u*x3w;
      double x3vw   = (*x3sumvw_)[i]/numsamp_-x3v*x3w;
      double x3uT   = (*x3sumuT_)[i]/numsamp_-x3u*x3T;
      double x3vT   = (*x3sumvT_)[i]/numsamp_-x3v*x3T;
      double x3wT   = (*x3sumwT_)[i]/numsamp_-x3w*x3T;

      (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x3coordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3rho;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3T;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3prms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3rhorms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3Trms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3uT;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3vT;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x3wT;
      (*log) << "\n";
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsLdc::DumpLomaStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsLdc::ClearStatistics()
{
  numsamp_ = 0;

  for(unsigned i=0; i<x1coordinates_->size(); ++i)
  {
    (*x1sumu_)[i]=0;
    (*x1sumv_)[i]=0;
    (*x1sumw_)[i]=0;
    (*x1sump_)[i]=0;
    (*x1sumrho_)[i]=0;
    (*x1sumT_)[i]=0;

    (*x1sumsqu_)[i]=0;
    (*x1sumsqv_)[i]=0;
    (*x1sumsqw_)[i]=0;
    (*x1sumsqp_)[i]=0;
    (*x1sumsqrho_)[i]=0;
    (*x1sumsqT_)[i]=0;

    (*x1sumuv_)[i] =0;
    (*x1sumuw_)[i] =0;
    (*x1sumvw_)[i] =0;
    (*x1sumuT_)[i] =0;
    (*x1sumvT_)[i] =0;
    (*x1sumwT_)[i] =0;
  }

  for(unsigned i=0; i<x2coordinates_->size(); ++i)
  {
    (*x2sumu_)[i]=0;
    (*x2sumv_)[i]=0;
    (*x2sumw_)[i]=0;
    (*x2sump_)[i]=0;
    (*x2sumrho_)[i]=0;
    (*x2sumT_)[i]=0;

    (*x2sumsqu_)[i]=0;
    (*x2sumsqv_)[i]=0;
    (*x2sumsqw_)[i]=0;
    (*x2sumsqp_)[i]=0;
    (*x2sumsqrho_)[i]=0;
    (*x2sumsqT_)[i]=0;

    (*x2sumuv_)[i] =0;
    (*x2sumuw_)[i] =0;
    (*x2sumvw_)[i] =0;
    (*x2sumuT_)[i] =0;
    (*x2sumvT_)[i] =0;
    (*x2sumwT_)[i] =0;
  }

  for(unsigned i=0; i<x3coordinates_->size(); ++i)
  {
    (*x3sumu_)[i]=0;
    (*x3sumv_)[i]=0;
    (*x3sumw_)[i]=0;
    (*x3sump_)[i]=0;
    (*x3sumrho_)[i]=0;
    (*x3sumT_)[i]=0;

    (*x3sumsqu_)[i]=0;
    (*x3sumsqv_)[i]=0;
    (*x3sumsqw_)[i]=0;
    (*x3sumsqp_)[i]=0;
    (*x3sumsqrho_)[i]=0;
    (*x3sumsqT_)[i]=0;

    (*x3sumuv_)[i] =0;
    (*x3sumuw_)[i] =0;
    (*x3sumvw_)[i] =0;
    (*x3sumuT_)[i] =0;
    (*x3sumvT_)[i] =0;
    (*x3sumwT_)[i] =0;
  }

  return;
}// TurbulenceStatisticsLdc::ClearStatistics


