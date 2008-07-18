/*!----------------------------------------------------------------------
\file turbulence_statistics_sqc.cpp

\brief calculate pressures, mean velocity values and fluctuations for
turbulent flow past a square cylinder.

<pre>
o Create sets for various evaluation lines in domain
  (Construction based on a round robin communication pattern):
  - centerline in x1-direction
  - centerline (with respect to cylinder center) in x2-direction
  - lines in wake at x1=7.5 and x1=11.5 in x2-direction
  - lines around cylinder

o loop nodes closest to / on lines

o values on lines are averaged in time over all steps between two
  outputs

Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "turbulence_statistics_sqc.H"

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

  <pre>
  o Create sets for lines in x1- and x2-direction

  o Allocate distributed vector for squares
  </pre>

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsSqc::TurbulenceStatisticsSqc(
  RefCountPtr<DRT::Discretization> actdis,
  ParameterList&                   params)
  :
  discret_(actdis),
  params_ (params)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
  {
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");
  }

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  squaredvelnp_ = LINALG::CreateVector(*dofrowmap,true);

  toggleu_      = LINALG::CreateVector(*dofrowmap,true);
  togglev_      = LINALG::CreateVector(*dofrowmap,true);
  togglew_      = LINALG::CreateVector(*dofrowmap,true);
  togglep_      = LINALG::CreateVector(*dofrowmap,true);

  // bounds for extension of flow domain in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates for required evaluation lines
  //----------------------------------------------------------------------
  x1ccoordinates_ = rcp(new vector<double> );
  x2ccoordinates_ = rcp(new vector<double> );
  x2wcoordinates_ = rcp(new vector<double> );
  clrcoordinates_ = rcp(new vector<double> );
  ctbcoordinates_ = rcp(new vector<double> );

  // the criterion allows differences in coordinates by 1e-9
  set<double,LineSortCriterion> x1cavcoords;
  set<double,LineSortCriterion> x2cavcoords;
  set<double,LineSortCriterion> x2wavcoords;
  set<double,LineSortCriterion> clravcoords;
  set<double,LineSortCriterion> ctbavcoords;

  // loop nodes and build sets of lines accessible on this proc
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);

    if (node->X()[1]<(7.0+2e-9) && node->X()[1]>(7.0-2e-9))
      x1cavcoords.insert(node->X()[0]);
    if (node->X()[0]<(5.0+2e-9) && node->X()[0]>(5.0-2e-9))
      x2cavcoords.insert(node->X()[1]);
    if (node->X()[0]<(7.5+2e-9) && node->X()[0]>(7.5-2e-9))
      x2wavcoords.insert(node->X()[1]);
    if ((node->X()[0]<(4.5+2e-9) && node->X()[0]>(4.5-2e-9)) &&
        (node->X()[1]<(7.5+2e-9) && node->X()[1]>(6.5-2e-9)))
      clravcoords.insert(node->X()[1]);
    if ((node->X()[0]<(5.5+2e-9) && node->X()[0]>(4.5-2e-9)) &&
        (node->X()[1]<(7.5+2e-9) && node->X()[1]>(7.5-2e-9)))
      ctbavcoords.insert(node->X()[0]);

    if (x3min_>node->X()[2])
    {
      x3min_=node->X()[2];
    }
    if (x3max_<node->X()[2])
    {
      x3max_=node->X()[2];
    }
  }

  // communicate x3mins and x3maxs
  double min;
  discret_->Comm().MinAll(&x3min_,&min,1);
  x3min_=min;

  double max;
  discret_->Comm().MaxAll(&x3max_,&max,1);
  x3max_=max;

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates on respective lines to
  // all procs
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

    // first, communicate coordinates in x1-centerline-direction
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x1cline=x1cavcoords.begin();
           x1cline!=x1cavcoords.end();
           ++x1cline)
      {
        DRT::ParObject::AddtoPack(sblock,*x1cline);
      }
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

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x1cavcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-centerline-direction
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x2cline=x2cavcoords.begin();
           x2cline!=x2cavcoords.end();
           ++x2cline)
      {
        DRT::ParObject::AddtoPack(sblock,*x2cline);
      }
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

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2cavcoords.insert(onecoord);
        }
      }
    }

    // third, communicate coordinates in x2-wakeline-direction
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x2wline=x2wavcoords.begin();
           x2wline!=x2wavcoords.end();
           ++x2wline)
      {
        DRT::ParObject::AddtoPack(sblock,*x2wline);
      }
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

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2wavcoords.insert(onecoord);
        }
      }
    }

    // fourth, communicate coordinates on left/right cylinder boundary
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator clrline=clravcoords.begin();
           clrline!=clravcoords.end();
           ++clrline)
      {
        DRT::ParObject::AddtoPack(sblock,*clrline);
      }
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

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          clravcoords.insert(onecoord);
        }
      }
    }

    // fifth, communicate coordinates on top/bottom cylinder boundary
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator ctbline=ctbavcoords.begin();
           ctbline!=ctbavcoords.end();
           ++ctbline)
      {
        DRT::ParObject::AddtoPack(sblock,*ctbline);
      }
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

        int index = 0;
        while (index < (int)rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          ctbavcoords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------
  {
    x1ccoordinates_ = rcp(new vector<double> );
    x2ccoordinates_ = rcp(new vector<double> );
    x2wcoordinates_ = rcp(new vector<double> );
    clrcoordinates_ = rcp(new vector<double> );
    ctbcoordinates_ = rcp(new vector<double> );

    for(set<double,LineSortCriterion>::iterator coord1c=x1cavcoords.begin();
        coord1c!=x1cavcoords.end();
        ++coord1c)
    {
      x1ccoordinates_->push_back(*coord1c);
    }

    for(set<double,LineSortCriterion>::iterator coord2c=x2cavcoords.begin();
        coord2c!=x2cavcoords.end();
        ++coord2c)
    {
      x2ccoordinates_->push_back(*coord2c);
    }

    for(set<double,LineSortCriterion>::iterator coord2w=x2wavcoords.begin();
        coord2w!=x2wavcoords.end();
        ++coord2w)
    {
      x2wcoordinates_->push_back(*coord2w);
    }

    for(set<double,LineSortCriterion>::iterator coordlr=clravcoords.begin();
        coordlr!=clravcoords.end();
        ++coordlr)
    {
      clrcoordinates_->push_back(*coordlr);
    }

    for(set<double,LineSortCriterion>::iterator coordtb=ctbavcoords.begin();
        coordtb!=ctbavcoords.end();
        ++coordtb)
    {
      ctbcoordinates_->push_back(*coordtb);
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  int size1c = x1ccoordinates_->size();
  int size2c = x2ccoordinates_->size();
  int size2w = x2wcoordinates_->size();
  int sizelr = clrcoordinates_->size();
  int sizetb = ctbcoordinates_->size();

  // first-order moments
  x1csumu_ =  rcp(new vector<double> );
  x1csumu_->resize(size1c,0.0);
  x2csumu_ =  rcp(new vector<double> );
  x2csumu_->resize(size2c,0.0);
  x2w1sumu_ =  rcp(new vector<double> );
  x2w1sumu_->resize(size2w,0.0);
  x2w2sumu_ =  rcp(new vector<double> );
  x2w2sumu_->resize(size2w,0.0);
  cyllsumu_ =  rcp(new vector<double> );
  cyllsumu_->resize(sizelr,0.0);
  cyltsumu_ =  rcp(new vector<double> );
  cyltsumu_->resize(sizetb,0.0);
  cylrsumu_ =  rcp(new vector<double> );
  cylrsumu_->resize(sizelr,0.0);
  cylbsumu_ =  rcp(new vector<double> );
  cylbsumu_->resize(sizetb,0.0);

  x1csumv_ =  rcp(new vector<double> );
  x1csumv_->resize(size1c,0.0);
  x2csumv_ =  rcp(new vector<double> );
  x2csumv_->resize(size2c,0.0);
  x2w1sumv_ =  rcp(new vector<double> );
  x2w1sumv_->resize(size2w,0.0);
  x2w2sumv_ =  rcp(new vector<double> );
  x2w2sumv_->resize(size2w,0.0);
  cyllsumv_ =  rcp(new vector<double> );
  cyllsumv_->resize(sizelr,0.0);
  cyltsumv_ =  rcp(new vector<double> );
  cyltsumv_->resize(sizetb,0.0);
  cylrsumv_ =  rcp(new vector<double> );
  cylrsumv_->resize(sizelr,0.0);
  cylbsumv_ =  rcp(new vector<double> );
  cylbsumv_->resize(sizetb,0.0);

  x1csumw_ =  rcp(new vector<double> );
  x1csumw_->resize(size1c,0.0);
  x2csumw_ =  rcp(new vector<double> );
  x2csumw_->resize(size2c,0.0);
  x2w1sumw_ =  rcp(new vector<double> );
  x2w1sumw_->resize(size2w,0.0);
  x2w2sumw_ =  rcp(new vector<double> );
  x2w2sumw_->resize(size2w,0.0);
  cyllsumw_ =  rcp(new vector<double> );
  cyllsumw_->resize(sizelr,0.0);
  cyltsumw_ =  rcp(new vector<double> );
  cyltsumw_->resize(sizetb,0.0);
  cylrsumw_ =  rcp(new vector<double> );
  cylrsumw_->resize(sizelr,0.0);
  cylbsumw_ =  rcp(new vector<double> );
  cylbsumw_->resize(sizetb,0.0);

  x1csump_ =  rcp(new vector<double> );
  x1csump_->resize(size1c,0.0);
  x2csump_ =  rcp(new vector<double> );
  x2csump_->resize(size2c,0.0);
  x2w1sump_ =  rcp(new vector<double> );
  x2w1sump_->resize(size2w,0.0);
  x2w2sump_ =  rcp(new vector<double> );
  x2w2sump_->resize(size2w,0.0);
  cyllsump_ =  rcp(new vector<double> );
  cyllsump_->resize(sizelr,0.0);
  cyltsump_ =  rcp(new vector<double> );
  cyltsump_->resize(sizetb,0.0);
  cylrsump_ =  rcp(new vector<double> );
  cylrsump_->resize(sizelr,0.0);
  cylbsump_ =  rcp(new vector<double> );
  cylbsump_->resize(sizetb,0.0);

  // second-order moments
  x1csumsqu_ =  rcp(new vector<double> );
  x1csumsqu_->resize(size1c,0.0);
  x2csumsqu_ =  rcp(new vector<double> );
  x2csumsqu_->resize(size2c,0.0);
  x2w1sumsqu_ =  rcp(new vector<double> );
  x2w1sumsqu_->resize(size2w,0.0);
  x2w2sumsqu_ =  rcp(new vector<double> );
  x2w2sumsqu_->resize(size2w,0.0);
  cyllsumsqu_ =  rcp(new vector<double> );
  cyllsumsqu_->resize(sizelr,0.0);
  cyltsumsqu_ =  rcp(new vector<double> );
  cyltsumsqu_->resize(sizetb,0.0);
  cylrsumsqu_ =  rcp(new vector<double> );
  cylrsumsqu_->resize(sizelr,0.0);
  cylbsumsqu_ =  rcp(new vector<double> );
  cylbsumsqu_->resize(sizetb,0.0);

  x1csumsqv_ =  rcp(new vector<double> );
  x1csumsqv_->resize(size1c,0.0);
  x2csumsqv_ =  rcp(new vector<double> );
  x2csumsqv_->resize(size2c,0.0);
  x2w1sumsqv_ =  rcp(new vector<double> );
  x2w1sumsqv_->resize(size2w,0.0);
  x2w2sumsqv_ =  rcp(new vector<double> );
  x2w2sumsqv_->resize(size2w,0.0);
  cyllsumsqv_ =  rcp(new vector<double> );
  cyllsumsqv_->resize(sizelr,0.0);
  cyltsumsqv_ =  rcp(new vector<double> );
  cyltsumsqv_->resize(sizetb,0.0);
  cylrsumsqv_ =  rcp(new vector<double> );
  cylrsumsqv_->resize(sizelr,0.0);
  cylbsumsqv_ =  rcp(new vector<double> );
  cylbsumsqv_->resize(sizetb,0.0);

  x1csumsqw_ =  rcp(new vector<double> );
  x1csumsqw_->resize(size1c,0.0);
  x2csumsqw_ =  rcp(new vector<double> );
  x2csumsqw_->resize(size2c,0.0);
  x2w1sumsqw_ =  rcp(new vector<double> );
  x2w1sumsqw_->resize(size2w,0.0);
  x2w2sumsqw_ =  rcp(new vector<double> );
  x2w2sumsqw_->resize(size2w,0.0);
  cyllsumsqw_ =  rcp(new vector<double> );
  cyllsumsqw_->resize(sizelr,0.0);
  cyltsumsqw_ =  rcp(new vector<double> );
  cyltsumsqw_->resize(sizetb,0.0);
  cylrsumsqw_ =  rcp(new vector<double> );
  cylrsumsqw_->resize(sizelr,0.0);
  cylbsumsqw_ =  rcp(new vector<double> );
  cylbsumsqw_->resize(sizetb,0.0);

  x1csumuv_ =  rcp(new vector<double> );
  x1csumuv_->resize(size1c,0.0);
  x2csumuv_ =  rcp(new vector<double> );
  x2csumuv_->resize(size2c,0.0);
  x2w1sumuv_ =  rcp(new vector<double> );
  x2w1sumuv_->resize(size2w,0.0);
  x2w2sumuv_ =  rcp(new vector<double> );
  x2w2sumuv_->resize(size2w,0.0);
  cyllsumuv_ =  rcp(new vector<double> );
  cyllsumuv_->resize(sizelr,0.0);
  cyltsumuv_ =  rcp(new vector<double> );
  cyltsumuv_->resize(sizetb,0.0);
  cylrsumuv_ =  rcp(new vector<double> );
  cylrsumuv_->resize(sizelr,0.0);
  cylbsumuv_ =  rcp(new vector<double> );
  cylbsumuv_->resize(sizetb,0.0);

  x1csumuw_ =  rcp(new vector<double> );
  x1csumuw_->resize(size1c,0.0);
  x2csumuw_ =  rcp(new vector<double> );
  x2csumuw_->resize(size2c,0.0);
  x2w1sumuw_ =  rcp(new vector<double> );
  x2w1sumuw_->resize(size2w,0.0);
  x2w2sumuw_ =  rcp(new vector<double> );
  x2w2sumuw_->resize(size2w,0.0);
  cyllsumuw_ =  rcp(new vector<double> );
  cyllsumuw_->resize(sizelr,0.0);
  cyltsumuw_ =  rcp(new vector<double> );
  cyltsumuw_->resize(sizetb,0.0);
  cylrsumuw_ =  rcp(new vector<double> );
  cylrsumuw_->resize(sizelr,0.0);
  cylbsumuw_ =  rcp(new vector<double> );
  cylbsumuw_->resize(sizetb,0.0);

  x1csumvw_ =  rcp(new vector<double> );
  x1csumvw_->resize(size1c,0.0);
  x2csumvw_ =  rcp(new vector<double> );
  x2csumvw_->resize(size2c,0.0);
  x2w1sumvw_ =  rcp(new vector<double> );
  x2w1sumvw_->resize(size2w,0.0);
  x2w2sumvw_ =  rcp(new vector<double> );
  x2w2sumvw_->resize(size2w,0.0);
  cyllsumvw_ =  rcp(new vector<double> );
  cyllsumvw_->resize(sizelr,0.0);
  cyltsumvw_ =  rcp(new vector<double> );
  cyltsumvw_->resize(sizetb,0.0);
  cylrsumvw_ =  rcp(new vector<double> );
  cylrsumvw_->resize(sizelr,0.0);
  cylbsumvw_ =  rcp(new vector<double> );
  cylbsumvw_->resize(sizetb,0.0);

  x1csumsqp_ =  rcp(new vector<double> );
  x1csumsqp_->resize(size1c,0.0);
  x2csumsqp_ =  rcp(new vector<double> );
  x2csumsqp_->resize(size2c,0.0);
  x2w1sumsqp_ =  rcp(new vector<double> );
  x2w1sumsqp_->resize(size2w,0.0);
  x2w2sumsqp_ =  rcp(new vector<double> );
  x2w2sumsqp_->resize(size2w,0.0);
  cyllsumsqp_ =  rcp(new vector<double> );
  cyllsumsqp_->resize(sizelr,0.0);
  cyltsumsqp_ =  rcp(new vector<double> );
  cyltsumsqp_->resize(sizetb,0.0);
  cylrsumsqp_ =  rcp(new vector<double> );
  cylrsumsqp_->resize(sizelr,0.0);
  cylbsumsqp_ =  rcp(new vector<double> );
  cylbsumsqp_->resize(sizetb,0.0);

  // clear statistics
  this->ClearStatistics();

  return;
}// TurbulenceStatisticsSqc::TurbulenceStatisticsSqc

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsSqc::~TurbulenceStatisticsSqc()
{
  return;
}// TurbulenceStatisticsSqc::~TurbulenceStatisticsSqc()


//----------------------------------------------------------------------
// sampling of lift/drag values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsSqc::DoLiftDragTimeSample(double dragforce,
                                                   double liftforce)
{

  double cdrag =2.0*dragforce/(x3max_-x3min_);
  double clift =2.0*liftforce/(x3max_-x3min_);

  drag_   +=cdrag;
  lift_   +=clift;
  dragsq_ +=cdrag*cdrag;
  liftsq_ +=clift*clift;

  return;
}// TurbulenceStatisticsSqc::DoLiftDragTimeSample


//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsSqc::DoTimeSample(
Teuchos::RefCountPtr<Epetra_Vector> velnp
  )
{

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1cnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x1-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x1cline=x1ccoordinates_->begin();
       x1cline!=x1ccoordinates_->end();
       ++x1cline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[0]<(*x1cline+2e-9) && node->X()[0]>(*x1cline-2e-9)) &&
          (node->X()[1]<(7.0+2e-9) && node->X()[1]>(7.0-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x1csumu_)[x1cnodnum]+=usm;
      (*x1csumv_)[x1cnodnum]+=vsm;
      (*x1csumw_)[x1cnodnum]+=wsm;
      (*x1csump_)[x1cnodnum]+=psm;

      (*x1csumsqu_)[x1cnodnum]+=usm*usm;
      (*x1csumsqv_)[x1cnodnum]+=vsm*vsm;
      (*x1csumsqw_)[x1cnodnum]+=wsm*wsm;
      (*x1csumsqp_)[x1cnodnum]+=psm*psm;

      (*x1csumuv_)[x1cnodnum]+=usm*vsm;
      (*x1csumuw_)[x1cnodnum]+=usm*wsm;
      (*x1csumvw_)[x1cnodnum]+=vsm*wsm;
    }
    x1cnodnum++;

  }

  int x2cnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x2-direction (with respect to cylinder
  // center) and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x2cline=x2ccoordinates_->begin();
       x2cline!=x2ccoordinates_->end();
       ++x2cline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
      if ((node->X()[1]<(*x2cline+2e-9) && node->X()[1]>(*x2cline-2e-9)) &&
          (node->X()[0]<(5.0+2e-9) && node->X()[0]>(5.0-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x2csumu_)[x2cnodnum]+=usm;
      (*x2csumv_)[x2cnodnum]+=vsm;
      (*x2csumw_)[x2cnodnum]+=wsm;
      (*x2csump_)[x2cnodnum]+=psm;

      (*x2csumsqu_)[x2cnodnum]+=usm*usm;
      (*x2csumsqv_)[x2cnodnum]+=vsm*vsm;
      (*x2csumsqw_)[x2cnodnum]+=wsm*wsm;
      (*x2csumsqp_)[x2cnodnum]+=psm*psm;

      (*x2csumuv_)[x2cnodnum]+=usm*vsm;
      (*x2csumuw_)[x2cnodnum]+=usm*wsm;
      (*x2csumvw_)[x2cnodnum]+=vsm*wsm;
    }
    x2cnodnum++;

  }

  int x2wnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on first wakeline in x2-direction (x1=7.5) and calculate
  // pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x2wline=x2wcoordinates_->begin();
       x2wline!=x2wcoordinates_->end();
       ++x2wline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
      if ((node->X()[1]<(*x2wline+2e-9) && node->X()[1]>(*x2wline-2e-9)) &&
          (node->X()[0]<(7.5+2e-9) && node->X()[0]>(7.5-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x2w1sumu_)[x2wnodnum]+=usm;
      (*x2w1sumv_)[x2wnodnum]+=vsm;
      (*x2w1sumw_)[x2wnodnum]+=wsm;
      (*x2w1sump_)[x2wnodnum]+=psm;

      (*x2w1sumsqu_)[x2wnodnum]+=usm*usm;
      (*x2w1sumsqv_)[x2wnodnum]+=vsm*vsm;
      (*x2w1sumsqw_)[x2wnodnum]+=wsm*wsm;
      (*x2w1sumsqp_)[x2wnodnum]+=psm*psm;

      (*x2w1sumuv_)[x2wnodnum]+=usm*vsm;
      (*x2w1sumuw_)[x2wnodnum]+=usm*wsm;
      (*x2w1sumvw_)[x2wnodnum]+=vsm*wsm;
    }
    x2wnodnum++;

  }

  x2wnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on second wakeline in x2-direction (x1=11.5) and calculate
  // pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator x2wline=x2wcoordinates_->begin();
       x2wline!=x2wcoordinates_->end();
       ++x2wline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
      if ((node->X()[1]<(*x2wline+2e-9) && node->X()[1]>(*x2wline-2e-9)) &&
          (node->X()[0]<(11.5+2e-9) && node->X()[0]>(11.5-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*x2w2sumu_)[x2wnodnum]+=usm;
      (*x2w2sumv_)[x2wnodnum]+=vsm;
      (*x2w2sumw_)[x2wnodnum]+=wsm;
      (*x2w2sump_)[x2wnodnum]+=psm;

      (*x2w2sumsqu_)[x2wnodnum]+=usm*usm;
      (*x2w2sumsqv_)[x2wnodnum]+=vsm*vsm;
      (*x2w2sumsqw_)[x2wnodnum]+=wsm*wsm;
      (*x2w2sumsqp_)[x2wnodnum]+=psm*psm;

      (*x2w2sumuv_)[x2wnodnum]+=usm*vsm;
      (*x2w2sumuw_)[x2wnodnum]+=usm*wsm;
      (*x2w2sumvw_)[x2wnodnum]+=vsm*wsm;
    }
    x2wnodnum++;

  }

  int clrnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on left cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator clrline=clrcoordinates_->begin();
       clrline!=clrcoordinates_->end();
       ++clrline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[1]<(*clrline+2e-9) && node->X()[1]>(*clrline-2e-9)) &&
          (node->X()[0]<(4.5+2e-9) && node->X()[0]>(4.5-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*cyllsumu_)[clrnodnum]+=usm;
      (*cyllsumv_)[clrnodnum]+=vsm;
      (*cyllsumw_)[clrnodnum]+=wsm;
      (*cyllsump_)[clrnodnum]+=psm;

      (*cyllsumsqu_)[clrnodnum]+=usm*usm;
      (*cyllsumsqv_)[clrnodnum]+=vsm*vsm;
      (*cyllsumsqw_)[clrnodnum]+=wsm*wsm;
      (*cyllsumsqp_)[clrnodnum]+=psm*psm;

      (*cyllsumuv_)[clrnodnum]+=usm*vsm;
      (*cyllsumuw_)[clrnodnum]+=usm*wsm;
      (*cyllsumvw_)[clrnodnum]+=vsm*wsm;
    }
    clrnodnum++;

  }

  int ctbnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on top cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator ctbline=ctbcoordinates_->begin();
       ctbline!=ctbcoordinates_->end();
       ++ctbline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[0]<(*ctbline+2e-9) && node->X()[0]>(*ctbline-2e-9)) &&
          (node->X()[1]<(7.5+2e-9) && node->X()[1]>(7.5-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*cyltsumu_)[ctbnodnum]+=usm;
      (*cyltsumv_)[ctbnodnum]+=vsm;
      (*cyltsumw_)[ctbnodnum]+=wsm;
      (*cyltsump_)[ctbnodnum]+=psm;

      (*cyltsumsqu_)[ctbnodnum]+=usm*usm;
      (*cyltsumsqv_)[ctbnodnum]+=vsm*vsm;
      (*cyltsumsqw_)[ctbnodnum]+=wsm*wsm;
      (*cyltsumsqp_)[ctbnodnum]+=psm*psm;

      (*cyltsumuv_)[ctbnodnum]+=usm*vsm;
      (*cyltsumuw_)[ctbnodnum]+=usm*wsm;
      (*cyltsumvw_)[ctbnodnum]+=vsm*wsm;
    }
    ctbnodnum++;

  }

  clrnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on right cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator clrline=clrcoordinates_->begin();
       clrline!=clrcoordinates_->end();
       ++clrline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[1]<(*clrline+2e-9) && node->X()[1]>(*clrline-2e-9)) &&
          (node->X()[0]<(5.5+2e-9) && node->X()[0]>(5.5-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this centerline node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*cylrsumu_)[clrnodnum]+=usm;
      (*cylrsumv_)[clrnodnum]+=vsm;
      (*cylrsumw_)[clrnodnum]+=wsm;
      (*cylrsump_)[clrnodnum]+=psm;

      (*cylrsumsqu_)[clrnodnum]+=usm*usm;
      (*cylrsumsqv_)[clrnodnum]+=vsm*vsm;
      (*cylrsumsqw_)[clrnodnum]+=wsm*wsm;
      (*cylrsumsqp_)[clrnodnum]+=psm*psm;

      (*cylrsumuv_)[clrnodnum]+=usm*vsm;
      (*cylrsumuw_)[clrnodnum]+=usm*wsm;
      (*cylrsumvw_)[clrnodnum]+=vsm*wsm;
    }
    clrnodnum++;

  }

  ctbnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on bottom cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (vector<double>::iterator ctbline=ctbcoordinates_->begin();
       ctbline!=ctbcoordinates_->end();
       ++ctbline)
  {

    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes=0;

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if ((node->X()[0]<(*ctbline+2e-9) && node->X()[0]>(*ctbline-2e-9)) &&
          (node->X()[1]<(6.5+2e-9) && node->X()[1]>(6.5-2e-9)) &&
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

      //----------------------------------------------------------------------
      // calculate spatial means for velocity and pressure on this line
      // (if more than one node contributing to this line node)
      //----------------------------------------------------------------------
      double usm=u/countnodesonallprocs;
      double vsm=v/countnodesonallprocs;
      double wsm=w/countnodesonallprocs;
      double psm=p/countnodesonallprocs;

      //----------------------------------------------------------------------
      // add spatial mean values to statistical sample
      //----------------------------------------------------------------------
      (*cylbsumu_)[ctbnodnum]+=usm;
      (*cylbsumv_)[ctbnodnum]+=vsm;
      (*cylbsumw_)[ctbnodnum]+=wsm;
      (*cylbsump_)[ctbnodnum]+=psm;

      (*cylbsumsqu_)[ctbnodnum]+=usm*usm;
      (*cylbsumsqv_)[ctbnodnum]+=vsm*vsm;
      (*cylbsumsqw_)[ctbnodnum]+=wsm*wsm;
      (*cylbsumsqp_)[ctbnodnum]+=psm*psm;

      (*cylbsumuv_)[ctbnodnum]+=usm*vsm;
      (*cylbsumuw_)[ctbnodnum]+=usm*wsm;
      (*cylbsumvw_)[ctbnodnum]+=vsm*wsm;
    }
    ctbnodnum++;

  }

  return;
}// TurbulenceStatisticsSqc::DoTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsSqc::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistic");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Flow statistics for turbulent flow past a square-section cylinder (first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";
    (*log) << "\n";

    double liftmean = lift_/numsamp_;
    double dragmean = drag_/numsamp_;
    double liftrms  = sqrt(liftsq_/numsamp_-liftmean*liftmean);
    double dragrms  = sqrt(dragsq_/numsamp_-dragmean*dragmean);

    (*log) << "# lift and drag values\n";
    (*log) << "# mean lift" << setw(11) << setprecision(4) << liftmean;
    (*log) << "\n";
    (*log) << "# mean drag" << setw(11) << setprecision(4) << dragmean;
    (*log) << "\n";
    (*log) << "# rms lift " << setw(11) << setprecision(4) << liftrms;
    (*log) << "\n";
    (*log) << "# rms drag " << setw(11) << setprecision(4) << dragrms;
    (*log) << "\n";
    (*log) << "\n";

    (*log) << "# centerline in x1-direction\n";
    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<x1ccoordinates_->size(); ++i)
    {
      double x1cu    = (*x1csumu_)[i]/numsamp_;
      double x1cv    = (*x1csumv_)[i]/numsamp_;
      double x1cw    = (*x1csumw_)[i]/numsamp_;
      double x1cp    = (*x1csump_)[i]/numsamp_;

      double x1curms = sqrt((*x1csumsqu_)[i]/numsamp_-x1cu*x1cu);
      double x1cvrms = sqrt((*x1csumsqv_)[i]/numsamp_-x1cv*x1cv);
      double x1cwrms = sqrt((*x1csumsqw_)[i]/numsamp_-x1cw*x1cw);
      double x1cprms = sqrt((*x1csumsqp_)[i]/numsamp_-x1cp*x1cp);

      double x1cuv   = (*x1csumuv_)[i]/numsamp_-x1cu*x1cv;
      double x1cuw   = (*x1csumuw_)[i]/numsamp_-x1cu*x1cw;
      double x1cvw   = (*x1csumvw_)[i]/numsamp_-x1cv*x1cw;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x1ccoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << x1cu;
      (*log) << "   " << setw(11) << setprecision(4) << x1cv;
      (*log) << "   " << setw(11) << setprecision(4) << x1cw;
      (*log) << "   " << setw(11) << setprecision(4) << x1cp;
      (*log) << "   " << setw(11) << setprecision(4) << x1curms;
      (*log) << "   " << setw(11) << setprecision(4) << x1cvrms;
      (*log) << "   " << setw(11) << setprecision(4) << x1cwrms;
      (*log) << "   " << setw(11) << setprecision(4) << x1cuv;
      (*log) << "   " << setw(11) << setprecision(4) << x1cuw;
      (*log) << "   " << setw(11) << setprecision(4) << x1cvw;
      (*log) << "   " << setw(11) << setprecision(4) << x1cprms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# centerline in x2-direction (with respect to cylinder center)\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<x2ccoordinates_->size(); ++i)
    {
      double x2cu    = (*x2csumu_)[i]/numsamp_;
      double x2cv    = (*x2csumv_)[i]/numsamp_;
      double x2cw    = (*x2csumw_)[i]/numsamp_;
      double x2cp    = (*x2csump_)[i]/numsamp_;

      double x2curms = sqrt((*x2csumsqu_)[i]/numsamp_-x2cu*x2cu);
      double x2cvrms = sqrt((*x2csumsqv_)[i]/numsamp_-x2cv*x2cv);
      double x2cwrms = sqrt((*x2csumsqw_)[i]/numsamp_-x2cw*x2cw);
      double x2cprms = sqrt((*x2csumsqp_)[i]/numsamp_-x2cp*x2cp);

      double x2cuv   = (*x2csumuv_)[i]/numsamp_-x2cu*x2cv;
      double x2cuw   = (*x2csumuw_)[i]/numsamp_-x2cu*x2cw;
      double x2cvw   = (*x2csumvw_)[i]/numsamp_-x2cv*x2cw;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x2ccoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << x2cu;
      (*log) << "   " << setw(11) << setprecision(4) << x2cv;
      (*log) << "   " << setw(11) << setprecision(4) << x2cw;
      (*log) << "   " << setw(11) << setprecision(4) << x2cp;
      (*log) << "   " << setw(11) << setprecision(4) << x2curms;
      (*log) << "   " << setw(11) << setprecision(4) << x2cvrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2cwrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2cuv;
      (*log) << "   " << setw(11) << setprecision(4) << x2cuw;
      (*log) << "   " << setw(11) << setprecision(4) << x2cvw;
      (*log) << "   " << setw(11) << setprecision(4) << x2cprms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# first wakeline in x2-direction (at x1=7.5)\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<x2wcoordinates_->size(); ++i)
    {
      double x2w1u    = (*x2w1sumu_)[i]/numsamp_;
      double x2w1v    = (*x2w1sumv_)[i]/numsamp_;
      double x2w1w    = (*x2w1sumw_)[i]/numsamp_;
      double x2w1p    = (*x2w1sump_)[i]/numsamp_;

      double x2w1urms = sqrt((*x2w1sumsqu_)[i]/numsamp_-x2w1u*x2w1u);
      double x2w1vrms = sqrt((*x2w1sumsqv_)[i]/numsamp_-x2w1v*x2w1v);
      double x2w1wrms = sqrt((*x2w1sumsqw_)[i]/numsamp_-x2w1w*x2w1w);
      double x2w1prms = sqrt((*x2w1sumsqp_)[i]/numsamp_-x2w1p*x2w1p);

      double x2w1uv   = (*x2w1sumuv_)[i]/numsamp_-x2w1u*x2w1v;
      double x2w1uw   = (*x2w1sumuw_)[i]/numsamp_-x2w1u*x2w1w;
      double x2w1vw   = (*x2w1sumvw_)[i]/numsamp_-x2w1v*x2w1w;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x2wcoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << x2w1u;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1v;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1w;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1p;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1urms;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1vrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1wrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1uv;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1uw;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1vw;
      (*log) << "   " << setw(11) << setprecision(4) << x2w1prms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# second wakeline in x2-direction (at x1=11.5)\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<x2wcoordinates_->size(); ++i)
    {
      double x2w2u    = (*x2w2sumu_)[i]/numsamp_;
      double x2w2v    = (*x2w2sumv_)[i]/numsamp_;
      double x2w2w    = (*x2w2sumw_)[i]/numsamp_;
      double x2w2p    = (*x2w2sump_)[i]/numsamp_;

      double x2w2urms = sqrt((*x2w2sumsqu_)[i]/numsamp_-x2w2u*x2w2u);
      double x2w2vrms = sqrt((*x2w2sumsqv_)[i]/numsamp_-x2w2v*x2w2v);
      double x2w2wrms = sqrt((*x2w2sumsqw_)[i]/numsamp_-x2w2w*x2w2w);
      double x2w2prms = sqrt((*x2w2sumsqp_)[i]/numsamp_-x2w2p*x2w2p);

      double x2w2uv   = (*x2w2sumuv_)[i]/numsamp_-x2w2u*x2w2v;
      double x2w2uw   = (*x2w2sumuw_)[i]/numsamp_-x2w2u*x2w2w;
      double x2w2vw   = (*x2w2sumvw_)[i]/numsamp_-x2w2v*x2w2w;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x2wcoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << x2w2u;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2v;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2w;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2p;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2urms;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2vrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2wrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2uv;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2uw;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2vw;
      (*log) << "   " << setw(11) << setprecision(4) << x2w2prms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# left cylinder boundary\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<clrcoordinates_->size(); ++i)
    {
      double cyllu    = (*cyllsumu_)[i]/numsamp_;
      double cyllv    = (*cyllsumv_)[i]/numsamp_;
      double cyllw    = (*cyllsumw_)[i]/numsamp_;
      double cyllp    = (*cyllsump_)[i]/numsamp_;

      double cyllurms = sqrt((*cyllsumsqu_)[i]/numsamp_-cyllu*cyllu);
      double cyllvrms = sqrt((*cyllsumsqv_)[i]/numsamp_-cyllv*cyllv);
      double cyllwrms = sqrt((*cyllsumsqw_)[i]/numsamp_-cyllw*cyllw);
      double cyllprms = sqrt((*cyllsumsqp_)[i]/numsamp_-cyllp*cyllp);

      double cylluv   = (*cyllsumuv_)[i]/numsamp_-cyllu*cyllv;
      double cylluw   = (*cyllsumuw_)[i]/numsamp_-cyllu*cyllw;
      double cyllvw   = (*cyllsumvw_)[i]/numsamp_-cyllv*cyllw;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*clrcoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << cyllu;
      (*log) << "   " << setw(11) << setprecision(4) << cyllv;
      (*log) << "   " << setw(11) << setprecision(4) << cyllw;
      (*log) << "   " << setw(11) << setprecision(4) << cyllp;
      (*log) << "   " << setw(11) << setprecision(4) << cyllurms;
      (*log) << "   " << setw(11) << setprecision(4) << cyllvrms;
      (*log) << "   " << setw(11) << setprecision(4) << cyllwrms;
      (*log) << "   " << setw(11) << setprecision(4) << cylluv;
      (*log) << "   " << setw(11) << setprecision(4) << cylluw;
      (*log) << "   " << setw(11) << setprecision(4) << cyllvw;
      (*log) << "   " << setw(11) << setprecision(4) << cyllprms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# top cylinder boundary\n";
    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<ctbcoordinates_->size(); ++i)
    {
      double cyltu    = (*cyltsumu_)[i]/numsamp_;
      double cyltv    = (*cyltsumv_)[i]/numsamp_;
      double cyltw    = (*cyltsumw_)[i]/numsamp_;
      double cyltp    = (*cyltsump_)[i]/numsamp_;

      double cylturms = sqrt((*cyltsumsqu_)[i]/numsamp_-cyltu*cyltu);
      double cyltvrms = sqrt((*cyltsumsqv_)[i]/numsamp_-cyltv*cyltv);
      double cyltwrms = sqrt((*cyltsumsqw_)[i]/numsamp_-cyltw*cyltw);
      double cyltprms = sqrt((*cyltsumsqp_)[i]/numsamp_-cyltp*cyltp);

      double cyltuv   = (*cyltsumuv_)[i]/numsamp_-cyltu*cyltv;
      double cyltuw   = (*cyltsumuw_)[i]/numsamp_-cyltu*cyltw;
      double cyltvw   = (*cyltsumvw_)[i]/numsamp_-cyltv*cyltw;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*ctbcoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << cyltu;
      (*log) << "   " << setw(11) << setprecision(4) << cyltv;
      (*log) << "   " << setw(11) << setprecision(4) << cyltw;
      (*log) << "   " << setw(11) << setprecision(4) << cyltp;
      (*log) << "   " << setw(11) << setprecision(4) << cylturms;
      (*log) << "   " << setw(11) << setprecision(4) << cyltvrms;
      (*log) << "   " << setw(11) << setprecision(4) << cyltwrms;
      (*log) << "   " << setw(11) << setprecision(4) << cyltuv;
      (*log) << "   " << setw(11) << setprecision(4) << cyltuw;
      (*log) << "   " << setw(11) << setprecision(4) << cyltvw;
      (*log) << "   " << setw(11) << setprecision(4) << cyltprms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# right cylinder boundary\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<clrcoordinates_->size(); ++i)
    {
      double cylru    = (*cylrsumu_)[i]/numsamp_;
      double cylrv    = (*cylrsumv_)[i]/numsamp_;
      double cylrw    = (*cylrsumw_)[i]/numsamp_;
      double cylrp    = (*cylrsump_)[i]/numsamp_;

      double cylrurms = sqrt((*cylrsumsqu_)[i]/numsamp_-cylru*cylru);
      double cylrvrms = sqrt((*cylrsumsqv_)[i]/numsamp_-cylrv*cylrv);
      double cylrwrms = sqrt((*cylrsumsqw_)[i]/numsamp_-cylrw*cylrw);
      double cylrprms = sqrt((*cylrsumsqp_)[i]/numsamp_-cylrp*cylrp);

      double cylruv   = (*cylrsumuv_)[i]/numsamp_-cylru*cylrv;
      double cylruw   = (*cylrsumuw_)[i]/numsamp_-cylru*cylrw;
      double cylrvw   = (*cylrsumvw_)[i]/numsamp_-cylrv*cylrw;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*clrcoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << cylru;
      (*log) << "   " << setw(11) << setprecision(4) << cylrv;
      (*log) << "   " << setw(11) << setprecision(4) << cylrw;
      (*log) << "   " << setw(11) << setprecision(4) << cylrp;
      (*log) << "   " << setw(11) << setprecision(4) << cylrurms;
      (*log) << "   " << setw(11) << setprecision(4) << cylrvrms;
      (*log) << "   " << setw(11) << setprecision(4) << cylrwrms;
      (*log) << "   " << setw(11) << setprecision(4) << cylruv;
      (*log) << "   " << setw(11) << setprecision(4) << cylruw;
      (*log) << "   " << setw(11) << setprecision(4) << cylrvw;
      (*log) << "   " << setw(11) << setprecision(4) << cylrprms;
      (*log) << "   \n";
    }

    (*log) << "\n";
    (*log) << "# bottom cylinder boundary\n";
    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<ctbcoordinates_->size(); ++i)
    {
      double cylbu    = (*cylbsumu_)[i]/numsamp_;
      double cylbv    = (*cylbsumv_)[i]/numsamp_;
      double cylbw    = (*cylbsumw_)[i]/numsamp_;
      double cylbp    = (*cylbsump_)[i]/numsamp_;

      double cylburms = sqrt((*cylbsumsqu_)[i]/numsamp_-cylbu*cylbu);
      double cylbvrms = sqrt((*cylbsumsqv_)[i]/numsamp_-cylbv*cylbv);
      double cylbwrms = sqrt((*cylbsumsqw_)[i]/numsamp_-cylbw*cylbw);
      double cylbprms = sqrt((*cylbsumsqp_)[i]/numsamp_-cylbp*cylbp);

      double cylbuv   = (*cylbsumuv_)[i]/numsamp_-cylbu*cylbv;
      double cylbuw   = (*cylbsumuw_)[i]/numsamp_-cylbu*cylbw;
      double cylbvw   = (*cylbsumvw_)[i]/numsamp_-cylbv*cylbw;

      (*log) <<  " "  << setw(11) << setprecision(4) << (*ctbcoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << cylbu;
      (*log) << "   " << setw(11) << setprecision(4) << cylbv;
      (*log) << "   " << setw(11) << setprecision(4) << cylbw;
      (*log) << "   " << setw(11) << setprecision(4) << cylbp;
      (*log) << "   " << setw(11) << setprecision(4) << cylburms;
      (*log) << "   " << setw(11) << setprecision(4) << cylbvrms;
      (*log) << "   " << setw(11) << setprecision(4) << cylbwrms;
      (*log) << "   " << setw(11) << setprecision(4) << cylbuv;
      (*log) << "   " << setw(11) << setprecision(4) << cylbuw;
      (*log) << "   " << setw(11) << setprecision(4) << cylbvw;
      (*log) << "   " << setw(11) << setprecision(4) << cylbprms;
      (*log) << "   \n";
    }

    log->flush();
  }

  return;

}// TurbulenceStatisticsSqc::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsSqc::ClearStatistics()
{
  numsamp_ =0;

  lift_   =0;
  drag_   =0;
  liftsq_ =0;
  dragsq_ =0;

  for(unsigned i=0; i<x1ccoordinates_->size(); ++i)
  {
    (*x1csumu_)[i]=0;
    (*x1csumv_)[i]=0;
    (*x1csumw_)[i]=0;
    (*x1csump_)[i]=0;

    (*x1csumsqu_)[i]=0;
    (*x1csumsqv_)[i]=0;
    (*x1csumsqw_)[i]=0;
    (*x1csumuv_)[i] =0;
    (*x1csumuw_)[i] =0;
    (*x1csumvw_)[i] =0;
    (*x1csumsqp_)[i]=0;
  }

  for(unsigned i=0; i<x2ccoordinates_->size(); ++i)
  {
    (*x2csumu_)[i]=0;
    (*x2csumv_)[i]=0;
    (*x2csumw_)[i]=0;
    (*x2csump_)[i]=0;

    (*x2csumsqu_)[i]=0;
    (*x2csumsqv_)[i]=0;
    (*x2csumsqw_)[i]=0;
    (*x2csumuv_)[i] =0;
    (*x2csumuw_)[i] =0;
    (*x2csumvw_)[i] =0;
    (*x2csumsqp_)[i]=0;
  }

  for(unsigned i=0; i<x2wcoordinates_->size(); ++i)
  {
    (*x2w1sumu_)[i]=0;
    (*x2w1sumv_)[i]=0;
    (*x2w1sumw_)[i]=0;
    (*x2w1sump_)[i]=0;

    (*x2w1sumsqu_)[i]=0;
    (*x2w1sumsqv_)[i]=0;
    (*x2w1sumsqw_)[i]=0;
    (*x2w1sumuv_)[i] =0;
    (*x2w1sumuw_)[i] =0;
    (*x2w1sumvw_)[i] =0;
    (*x2w1sumsqp_)[i]=0;

    (*x2w2sumu_)[i]=0;
    (*x2w2sumv_)[i]=0;
    (*x2w2sumw_)[i]=0;
    (*x2w2sump_)[i]=0;

    (*x2w2sumsqu_)[i]=0;
    (*x2w2sumsqv_)[i]=0;
    (*x2w2sumsqw_)[i]=0;
    (*x2w2sumuv_)[i] =0;
    (*x2w2sumuw_)[i] =0;
    (*x2w2sumvw_)[i] =0;
    (*x2w2sumsqp_)[i]=0;
  }

  for(unsigned i=0; i<clrcoordinates_->size(); ++i)
  {
    (*cyllsumu_)[i]=0;
    (*cyllsumv_)[i]=0;
    (*cyllsumw_)[i]=0;
    (*cyllsump_)[i]=0;

    (*cyllsumsqu_)[i]=0;
    (*cyllsumsqv_)[i]=0;
    (*cyllsumsqw_)[i]=0;
    (*cyllsumuv_)[i] =0;
    (*cyllsumuw_)[i] =0;
    (*cyllsumvw_)[i] =0;
    (*cyllsumsqp_)[i]=0;

    (*cylrsumu_)[i]=0;
    (*cylrsumv_)[i]=0;
    (*cylrsumw_)[i]=0;
    (*cylrsump_)[i]=0;

    (*cylrsumsqu_)[i]=0;
    (*cylrsumsqv_)[i]=0;
    (*cylrsumsqw_)[i]=0;
    (*cylrsumuv_)[i] =0;
    (*cylrsumuw_)[i] =0;
    (*cylrsumvw_)[i] =0;
    (*cylrsumsqp_)[i]=0;
  }

  for(unsigned i=0; i<ctbcoordinates_->size(); ++i)
  {
    (*cyltsumu_)[i]=0;
    (*cyltsumv_)[i]=0;
    (*cyltsumw_)[i]=0;
    (*cyltsump_)[i]=0;

    (*cyltsumsqu_)[i]=0;
    (*cyltsumsqv_)[i]=0;
    (*cyltsumsqw_)[i]=0;
    (*cyltsumuv_)[i] =0;
    (*cyltsumuw_)[i] =0;
    (*cyltsumvw_)[i] =0;
    (*cyltsumsqp_)[i]=0;

    (*cylbsumu_)[i]=0;
    (*cylbsumv_)[i]=0;
    (*cylbsumw_)[i]=0;
    (*cylbsump_)[i]=0;

    (*cylbsumsqu_)[i]=0;
    (*cylbsumsqv_)[i]=0;
    (*cylbsumsqw_)[i]=0;
    (*cylbsumuv_)[i] =0;
    (*cylbsumuw_)[i] =0;
    (*cylbsumvw_)[i] =0;
    (*cylbsumsqp_)[i]=0;
  }

  return;
}// TurbulenceStatisticsSqc::ClearStatistics



#endif /* CCADISCRET       */
