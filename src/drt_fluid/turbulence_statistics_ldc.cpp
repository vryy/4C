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
#ifdef CCADISCRET

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
TurbulenceStatisticsLdc::TurbulenceStatisticsLdc(
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

  // bounds for extension of cavity in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates for centerlines in x1- and x2-direction
  //----------------------------------------------------------------------
  x1coordinates_ = rcp(new vector<double> );
  x2coordinates_ = rcp(new vector<double> );

  // the criterion allows differences in coordinates by 1e-9
  set<double,LineSortCriterion> x1avcoords;
  set<double,LineSortCriterion> x2avcoords;

  // loop nodes, build sets of centerlines accessible on this proc and
  // calculate extension of cavity in x3-direction
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);
    x1avcoords.insert(node->X()[0]);
    x2avcoords.insert(node->X()[1]);

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
  // round robin loop to communicate coordinates in x1- and x2-direction
  // to all procs
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
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x1line=x1avcoords.begin();
           x1line!=x1avcoords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(sblock,*x1line);
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
          x1avcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-direction
    for (int np=0;np<numprocs;++np)
    {
      // export set to sendbuffer
      sblock.clear();

      for (set<double,LineSortCriterion>::iterator x2line=x2avcoords.begin();
           x2line!=x2avcoords.end();
           ++x2line)
      {
        DRT::ParObject::AddtoPack(sblock,*x2line);
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
          x2avcoords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in x1- and x2-direction in a vector
  //----------------------------------------------------------------------
  {
    x1coordinates_ = rcp(new vector<double> );
    x2coordinates_ = rcp(new vector<double> );

    for(set<double,LineSortCriterion>::iterator coord1=x1avcoords.begin();
        coord1!=x1avcoords.end();
        ++coord1)
    {
      x1coordinates_->push_back(*coord1);
    }

    for(set<double,LineSortCriterion>::iterator coord2=x2avcoords.begin();
        coord2!=x2avcoords.end();
        ++coord2)
    {
      x2coordinates_->push_back(*coord2);
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  int size1 = x1coordinates_->size();
  int size2 = x2coordinates_->size();

  // first-order moments
  x1sumu_ =  rcp(new vector<double> );
  x1sumu_->resize(size1,0.0);
  x2sumu_ =  rcp(new vector<double> );
  x2sumu_->resize(size2,0.0);

  x1sumv_ =  rcp(new vector<double> );
  x1sumv_->resize(size1,0.0);
  x2sumv_ =  rcp(new vector<double> );
  x2sumv_->resize(size2,0.0);

  x1sumw_ =  rcp(new vector<double> );
  x1sumw_->resize(size1,0.0);
  x2sumw_ =  rcp(new vector<double> );
  x2sumw_->resize(size2,0.0);

  x1sump_ =  rcp(new vector<double> );
  x1sump_->resize(size1,0.0);
  x2sump_ =  rcp(new vector<double> );
  x2sump_->resize(size2,0.0);

  // second-order moments
  x1sumsqu_ =  rcp(new vector<double> );
  x1sumsqu_->resize(size1,0.0);
  x2sumsqu_ =  rcp(new vector<double> );
  x2sumsqu_->resize(size2,0.0);

  x1sumsqv_ =  rcp(new vector<double> );
  x1sumsqv_->resize(size1,0.0);
  x2sumsqv_ =  rcp(new vector<double> );
  x2sumsqv_->resize(size2,0.0);

  x1sumsqw_ =  rcp(new vector<double> );
  x1sumsqw_->resize(size1,0.0);
  x2sumsqw_ =  rcp(new vector<double> );
  x2sumsqw_->resize(size2,0.0);

  x1sumuv_ =  rcp(new vector<double> );
  x1sumuv_->resize(size1,0.0);
  x2sumuv_ =  rcp(new vector<double> );
  x2sumuv_->resize(size2,0.0);

  x1sumuw_ =  rcp(new vector<double> );
  x1sumuw_->resize(size1,0.0);
  x2sumuw_ =  rcp(new vector<double> );
  x2sumuw_->resize(size2,0.0);

  x1sumvw_ =  rcp(new vector<double> );
  x1sumvw_->resize(size1,0.0);
  x2sumvw_ =  rcp(new vector<double> );
  x2sumvw_->resize(size2,0.0);

  x1sumsqp_ =  rcp(new vector<double> );
  x1sumsqp_->resize(size1,0.0);
  x2sumsqp_ =  rcp(new vector<double> );
  x2sumsqp_->resize(size2,0.0);

  // clear statistics
  this->ClearStatistics();

  return;
}// TurbulenceStatisticsLdc::TurbulenceStatisticsLdc

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
TurbulenceStatisticsLdc::~TurbulenceStatisticsLdc()
{
  return;
}// TurbulenceStatisticsLdc::~TurbulenceStatisticsLdc()

void TurbulenceStatisticsLdc::DoTimeSample(
Teuchos::RefCountPtr<Epetra_Vector> velnp
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

  return;
}// TurbulenceStatisticsLdc::DoTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void TurbulenceStatisticsLdc::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistic");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Flow statistics for turbulent flow in a lid-driven cavity (first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "           urms          vrms          wrms";
    (*log) << "           u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      double x1u    = (*x1sumu_)[i]/numsamp_;
      double x1v    = (*x1sumv_)[i]/numsamp_;
      double x1w    = (*x1sumw_)[i]/numsamp_;
      double x1p    = (*x1sump_)[i]/numsamp_;

      double x1urms = 10*sqrt((*x1sumsqu_)[i]/numsamp_-x1u*x1u);
      double x1vrms = 10*sqrt((*x1sumsqv_)[i]/numsamp_-x1v*x1v);
      double x1wrms = 10*sqrt((*x1sumsqw_)[i]/numsamp_-x1w*x1w);
      double x1prms = 10*sqrt((*x1sumsqp_)[i]/numsamp_-x1p*x1p);

      double x1uv   = 500*((*x1sumuv_)[i]/numsamp_-x1u*x1v);
      double x1uw   = 500*((*x1sumuw_)[i]/numsamp_-x1u*x1w);
      double x1vw   = 500*((*x1sumvw_)[i]/numsamp_-x1v*x1w);

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x1coordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << x1u;
      (*log) << "   " << setw(11) << setprecision(4) << x1v;
      (*log) << "   " << setw(11) << setprecision(4) << x1w;
      (*log) << "   " << setw(11) << setprecision(4) << x1p;
      (*log) << "   " << setw(11) << setprecision(4) << x1urms;
      (*log) << "   " << setw(11) << setprecision(4) << x1vrms;
      (*log) << "   " << setw(11) << setprecision(4) << x1wrms;
      (*log) << "   " << setw(11) << setprecision(4) << x1uv;
      (*log) << "   " << setw(11) << setprecision(4) << x1uw;
      (*log) << "   " << setw(11) << setprecision(4) << x1vw;
      (*log) << "   " << setw(11) << setprecision(4) << x1prms;
      (*log) << "   \n";
    }

    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "           urms          vrms          wrms";
    (*log) << "           u'v'          u'w'          v'w'          prms   \n";

    (*log) << scientific;
    for(unsigned i=0; i<x2coordinates_->size(); ++i)
    {
      double x2u    = (*x2sumu_)[i]/numsamp_;
      double x2v    = (*x2sumv_)[i]/numsamp_;
      double x2w    = (*x2sumw_)[i]/numsamp_;
      double x2p    = (*x2sump_)[i]/numsamp_;

      double x2urms = 10*sqrt((*x2sumsqu_)[i]/numsamp_-x2u*x2u);
      double x2vrms = 10*sqrt((*x2sumsqv_)[i]/numsamp_-x2v*x2v);
      double x2wrms = 10*sqrt((*x2sumsqw_)[i]/numsamp_-x2w*x2w);
      double x2prms = 10*sqrt((*x2sumsqp_)[i]/numsamp_-x2p*x2p);

      double x2uv   = 500*((*x2sumuv_)[i]/numsamp_-x2u*x2v);
      double x2uw   = 500*((*x2sumuw_)[i]/numsamp_-x2u*x2w);
      double x2vw   = 500*((*x2sumvw_)[i]/numsamp_-x2v*x2w);

      (*log) <<  " "  << setw(11) << setprecision(4) << (*x2coordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << x2u;
      (*log) << "   " << setw(11) << setprecision(4) << x2v;
      (*log) << "   " << setw(11) << setprecision(4) << x2w;
      (*log) << "   " << setw(11) << setprecision(4) << x2p;
      (*log) << "   " << setw(11) << setprecision(4) << x2urms;
      (*log) << "   " << setw(11) << setprecision(4) << x2vrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2wrms;
      (*log) << "   " << setw(11) << setprecision(4) << x2uv;
      (*log) << "   " << setw(11) << setprecision(4) << x2uw;
      (*log) << "   " << setw(11) << setprecision(4) << x2vw;
      (*log) << "   " << setw(11) << setprecision(4) << x2prms;
      (*log) << "   \n";
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsLdc::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void TurbulenceStatisticsLdc::ClearStatistics()
{
  numsamp_ =0;

  for(unsigned i=0; i<x1coordinates_->size(); ++i)
  {
    (*x1sumu_)[i]=0;
    (*x1sumv_)[i]=0;
    (*x1sumw_)[i]=0;
    (*x1sump_)[i]=0;

    (*x1sumsqu_)[i]=0;
    (*x1sumsqv_)[i]=0;
    (*x1sumsqw_)[i]=0;
    (*x1sumuv_)[i] =0;
    (*x1sumuw_)[i] =0;
    (*x1sumvw_)[i] =0;
    (*x1sumsqp_)[i]=0;
  }

  for(unsigned i=0; i<x2coordinates_->size(); ++i)
  {
    (*x2sumu_)[i]=0;
    (*x2sumv_)[i]=0;
    (*x2sumw_)[i]=0;
    (*x2sump_)[i]=0;

    (*x2sumsqu_)[i]=0;
    (*x2sumsqv_)[i]=0;
    (*x2sumsqw_)[i]=0;
    (*x2sumuv_)[i] =0;
    (*x2sumuw_)[i] =0;
    (*x2sumvw_)[i] =0;
    (*x2sumsqp_)[i]=0;
  }

  return;
}// TurbulenceStatisticsLdc::ClearStatistics



#endif /* CCADISCRET       */
