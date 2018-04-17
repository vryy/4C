/*!----------------------------------------------------------------------
\file turbulence_statistics_ph.cpp

\brief calculate pressures, mean velocity values and fluctuations for
turbulent flow over a periodic constrictions
literature see Froehlich Mellen Rodi Temmerman Leschziner (2004) and Rapp Breuer Manhart Peller (2010)
created in a Semesterarbeit by Stephan Jaeger 2015

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/

#include "turbulence_statistics_ph.H"

#define SAMP_ALL

//----------------------------------------------------------------------
// constructor
//----------------------------------------------------------------------
FLD::TurbulenceStatisticsPh::TurbulenceStatisticsPh(
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::ParameterList&           params,
  const std::string&                statistics_outfilename):
  discret_      (actdis),
  params_       (params),
  statistics_outfilename_(statistics_outfilename)
{
  if (discret_->Comm().MyPID()==0)
  {
    std::cout << "This is the turbulence statistics manager of periodic hill problem" << std::endl;
    std::cout << "based on the geometry of ERCOFTAC" << std::endl;
  }

  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");


  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  squaredvelnp_ = LINALG::CreateVector(*dofrowmap,true);

  toggleu_      = LINALG::CreateVector(*dofrowmap,true);
  togglev_      = LINALG::CreateVector(*dofrowmap,true);
  togglew_      = LINALG::CreateVector(*dofrowmap,true);
  togglep_      = LINALG::CreateVector(*dofrowmap,true);

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


  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);

    if (node->X()[1]<85.008+2e-9 && node->X()[1]>85.008-2e-9)
      x1avcoords.insert(node->X()[0]);

    if (node->X()[0]<2e-9 && node->X()[0]>-2e-9)
      x2avcoords.insert(node->X()[1]);
  }

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
    }

    for(std::set<double,LineSortCriterion>::iterator coord2=x2avcoords.begin();
        coord2!=x2avcoords.end();
        ++coord2)
    {
      x2coordinates_->push_back(*coord2);
    }
  }
  //----------------------------------------------------------------------
  // number of coordinates in x1- and x2-direction
  //----------------------------------------------------------------------
  numx1coor_ = x1coordinates_->size();
  numx2coor_ = x2coordinates_->size();

  //----------------------------------------------------------------------
  // define locations in x1-direction for statistical evaluation
  //----------------------------------------------------------------------

#ifndef SAMP_ALL
  //define vector containing x1-coords for sampling

  x1setstatlocations_ = Teuchos::rcp(new std::vector<double> );

  //! x1-coordinates for statistical sampling  !!! in sampling routine the offwall distance has to be modified for different statlocations
  x1setstatlocations_->push_back(0);
  x1setstatlocations_->push_back(14);
  x1setstatlocations_->push_back(28);
  x1setstatlocations_->push_back(56);
  x1setstatlocations_->push_back(84);
  x1setstatlocations_->push_back(112);
  x1setstatlocations_->push_back(140);
  x1setstatlocations_->push_back(168);
  x1setstatlocations_->push_back(196);
  x1setstatlocations_->push_back(224);
//
#endif
#ifdef SAMP_ALL
  x1setstatlocations_ = x1coordinates_;
#endif
  numx1statlocations_ = x1setstatlocations_->size();

  if(discret_->Comm().MyPID()==0)
  {
    std::cout << "Sample at " << numx1statlocations_ << " numx1statlocations_ \n" << std::endl;
  }

// Resize matrix for x1statlocations and x2statlocations coordinates
  x1statlocations.reshape(numx1statlocations_,1);
  x2statlocations.reshape(numx1statlocations_,numx2coor_);

#ifndef SAMP_ALL
  x1elemlengthhalf_ = (x1max_ - x1min_)/(numx1coor_-1);
  for (int x1line =0 ; x1line < numx1statlocations_; ++x1line)
  {
    int pos = 0;
    mindist_ = x1elemlengthhalf_;
    for (int x1linecoords = 0; x1linecoords < numx1coor_; ++x1linecoords)
    {
      dist_ = abs(x1setstatlocations_->at(x1line) - x1coordinates_->at(x1linecoords));
      if (dist_ < mindist_)
      {
        mindist_ = dist_;
        pos = x1linecoords;
      }
    }
    x1statlocations(x1line,0) = x1coordinates_->at(pos);
  }
#endif
#ifdef SAMP_ALL
  for (int x1line =0 ; x1line < numx1statlocations_; ++x1line)
  {
    x1statlocations(x1line,0) = x1setstatlocations_->at(x1line);
  }
#endif

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
          && node->X()[0] < x1statlocations(x1stat,0)+2e-9 && node->X()[0] > x1statlocations(x1stat,0)-2e-9)
      {
        x2statlocat.insert(node->X()[1]);
      }
    }

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
        std::vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          x2statlocat.insert(onecoord);
        }
      }
    }
    }

  int x2it = -1;
  for (std::set<double,LineSortCriterion>::iterator x2locc=x2statlocat.begin();
       x2locc!=x2statlocat.end();
       ++x2locc)
  {
    x2it += 1;
    x2statlocations(x1stat,x2it) = *x2locc;
  }
  }

///////////////////////////////////////////
//    end of communication in x2-direction
///////////////////////////////////////////

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------

  // x1-direction
  x1sump_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sump_->Reshape(numx1statlocations_,2);

  x1sumu_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumu_->Reshape(numx1statlocations_,1);

  x1sumf_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumf_->Reshape(numx1statlocations_,3);


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

  if (discret_->Comm().MyPID()==0)
  {
    std::string s(statistics_outfilename_);

      s.append(".flow_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
      (*log) << "# Statistics for turbulent incompressible flow over periodic hill (first- and second-order moments)\n\n";

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
void FLD::TurbulenceStatisticsPh::DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp,Teuchos::RCP<Epetra_Vector> stresses)
{
  if(discret_->Comm().MyPID()==0)
    std::cout << "------------Time Sampling Routine begins---------" <<std::endl;


  // compute squared values of velocity
  squaredvelnp_->Multiply(1.0,*velnp,*velnp,0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

#ifdef SAMP_ALL

    for (int k=0; k <numx1statlocations_; ++k)
    {
      // current x2-coordinate of respective wall
      //constant offset from wall at the middle bottom; has to be adapted for different statlocations
      double x2cwall = x2statlocations(k,1);
      double x1cwall = x1statlocations(k,0);
      double x2cwall_puw = x2statlocations(k,numx2coor_-1);
      double x2cwall_plw = x2statlocations(k,0);

      // toggle vectors are one in the position of a dof of this node,
      toggleu_->PutScalar(0.0);

      //du1/dx
      {
        // count the number of nodes in x3-direction contributing to this nodal value
        int countnodes=0;

        for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
        {
          DRT::Node* node = discret_->lRowNode(nn);
          // this is the wall node
          if (node->X()[0]<(x1cwall+2e-9) and node->X()[0]>(x1cwall-2e-9) and
              node->X()[1]<(x2cwall+2e-9) and node->X()[1]>(x2cwall-2e-9))
          {
            std::vector<int> dof = discret_->Dof(node);
            double           one = 1.0;

            toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
            countnodes++;
          }
        }
        int countnodesonallprocs=0;

        discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

        // reduce by 1 due to periodic boundary condition
        countnodesonallprocs-=1;
        if (countnodesonallprocs)
        {
          //----------------------------------------------------------------------
          // get values for velocity derivative and pressure
          //----------------------------------------------------------------------
          double u;
          velnp->Dot(*toggleu_,&u);

          //----------------------------------------------------------------------
          // calculate spatial means
          //----------------------------------------------------------------------
          double usm=u/countnodesonallprocs;
          //----------------------------------------------------------------------
          // add spatial mean values to statistical sample
          //----------------------------------------------------------------------
          (*x1sumu_)(k,0)+=usm;
        }
      }

      togglep_->PutScalar(0.0);
      //p, lower wall
      {
        // count the number of nodes in x3-direction contributing to this nodal value
        int countnodes=0;

        for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
        {
          DRT::Node* node = discret_->lRowNode(nn);
          // this is the wall node
          if (node->X()[0]<(x1cwall+2e-9) and node->X()[0]>(x1cwall-2e-9) and
              node->X()[1]<(x2cwall_plw+2e-9) and node->X()[1]>(x2cwall_plw-2e-9))
          {
            std::vector<int> dof = discret_->Dof(node);
            double           one = 1.0;

            togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));
            countnodes++;
          }
        }
        int countnodesonallprocs=0;

        discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

        // reduce by 1 due to periodic boundary condition
        countnodesonallprocs-=1;
        if (countnodesonallprocs)
        {
          //----------------------------------------------------------------------
          // get values for velocity derivative and pressure
          //----------------------------------------------------------------------
          double p;
          velnp->Dot(*togglep_,&p);

          //----------------------------------------------------------------------
          // calculate spatial means
          //----------------------------------------------------------------------
          double lwpsm=p/countnodesonallprocs;
          //----------------------------------------------------------------------
          // add spatial mean values to statistical sample
          //----------------------------------------------------------------------
          (*x1sump_)(k,0)+=lwpsm;
        }
      }

      togglep_->PutScalar(0.0);
      // p, upper wall
      {
        // count the number of nodes in x3-direction contributing to this nodal value
        int countnodes=0;

        for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
        {
          DRT::Node* node = discret_->lRowNode(nn);
          // this is the wall node
          if (node->X()[0]<(x1cwall+2e-9) and node->X()[0]>(x1cwall-2e-9) and
              node->X()[1]<(x2cwall_puw+2e-9) and node->X()[1]>(x2cwall_puw-2e-9))
          {
            std::vector<int> dof = discret_->Dof(node);
            double           one = 1.0;

            togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));
            countnodes++;
          }
        }
        int countnodesonallprocs=0;

        discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

        // reduce by 1 due to periodic boundary condition
        countnodesonallprocs-=1;
        if (countnodesonallprocs)
        {
          //----------------------------------------------------------------------
          // get values for velocity derivative and pressure
          //----------------------------------------------------------------------
          double p;
          velnp->Dot(*togglep_,&p);

          //----------------------------------------------------------------------
          // calculate spatial means
          //----------------------------------------------------------------------
          double uwpsm=p/countnodesonallprocs;
          //----------------------------------------------------------------------
          // add spatial mean values to statistical sample
          //----------------------------------------------------------------------
          (*x1sump_)(k,1)+=uwpsm;
        }
      }
    }

    //wall stresses
    for (int k=0; k <numx1statlocations_; ++k)
    {
      // current x2-coordinate of respective wall
      //constant offset from wall at the middle bottom; has to be adapted for different statlocations
      double x2cwall = x2statlocations(k,0);
      double x1cwall = x1statlocations(k,0);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes=0;

      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);
        // this is the wall node
        if (node->X()[0]<(x1cwall+2e-9) and node->X()[0]>(x1cwall-2e-9) and
            node->X()[1]<(x2cwall+2e-9) and node->X()[1]>(x2cwall-2e-9))
        {
          std::vector<int> dof = discret_->Dof(node);
          double           one = 1.0;

          toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
          togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
          togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
          countnodes++;
        }
      }
      int countnodesonallprocs=0;

      discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs-=1;
      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative and pressure
        //----------------------------------------------------------------------
        double fx;
        stresses->Dot(*toggleu_,&fx);
        double fy;
        stresses->Dot(*togglev_,&fy);
        double fz;
        stresses->Dot(*togglew_,&fz);

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumf_)(k,0)+=fx/(double)countnodesonallprocs;
        (*x1sumf_)(k,1)+=fy/(double)countnodesonallprocs;
        (*x1sumf_)(k,2)+=fz/(double)countnodesonallprocs;
      }
    }

    //----------------------------------------------------------------------
    // end of loop velocity gradient at lower wall
    //----------------------------------------------------------------------
#endif

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------

  for (int x1nodnum=0; x1nodnum<numx1statlocations_; ++x1nodnum)
  {

    //----------------------------------------------------------------------
    // loop nodes in x2-direction
    //----------------------------------------------------------------------
    for (int x2line = 0; x2line < numx2coor_; ++x2line)
    {

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
        if (node->X()[0] < x1statlocations(x1nodnum,0)+2e-5 && node->X()[0] > x1statlocations(x1nodnum,0)-2e-5 &&
            node->X()[1] < x2statlocations(x1nodnum,x2line)+2e-9 && node->X()[1] > x2statlocations(x1nodnum,x2line)-2e-9) //the last node in x3 direction has to be sampled
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
        (*x2sumu_)(x1nodnum,x2line)+=u/countnodesonallprocs;
        (*x2sumv_)(x1nodnum,x2line)+=v/countnodesonallprocs;
        (*x2sumw_)(x1nodnum,x2line)+=w/countnodesonallprocs;
        (*x2sump_)(x1nodnum,x2line)+=p/countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum,x2line)+=uu/countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum,x2line)+=vv/countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum,x2line)+=ww/countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum,x2line)+=pp/countnodesonallprocs;

        (*x2sumuv_)(x1nodnum,x2line)+=uv/countnodesonallprocs;
        (*x2sumuw_)(x1nodnum,x2line)+=uw/countnodesonallprocs;
        (*x2sumvw_)(x1nodnum,x2line)+=vw/countnodesonallprocs;
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
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow for a periodic hill (first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n\n\n";
    (*log) << std::scientific;
    (*log) << "\n\n\n";

#ifdef SAMP_ALL
    (*log) << "# lower wall behind step\n";
    (*log) << "#     x1";
    (*log) << "           duxdy       pmean (lw)    pmean (uw)     tau_w          dist\n";
    // distance from wall to first node off wall
    for (int i=0; i<numx1statlocations_; ++i)
    {
        double lwx1u     = (*x1sumu_)(i,0)/numsamp_;
        double dist =  x2statlocations(i,1) - x2statlocations(i,0);
        double lwx1duxdy = lwx1u/dist;
        double lwx1p     = (*x1sump_)(i,0)/numsamp_;
        double uwx1p     = (*x1sump_)(i,1)/numsamp_;
        double fx     = (*x1sumf_)(i,0)/numsamp_;
        double fy     = (*x1sumf_)(i,1)/numsamp_;
        double fz     = (*x1sumf_)(i,2)/numsamp_;
        double f = std::sqrt(fx*fx+fy*fy+fz*fz);
        if(fx<0.0)
          f*=-1.0; //recover sign
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << x1statlocations(i,0);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1duxdy;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1p;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1p;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << f;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << dist;
        (*log) << "\n";
    }
#endif
    for (int i=0; i<numx1statlocations_; ++i)
    {
      // current x1-coordinate
      double x1 = 1.0e20;
      x1 = x1statlocations(i,0);

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





