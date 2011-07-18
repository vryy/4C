/*!-----------------------------------------------------------------------------------------------*
\file turbulence_statistics_oracles.cpp

\brief statistical data processing for ORACLES problem

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "turbulence_statistics_oracles.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_globalproblem.H"

#define NODETOL 1e-9


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::TurbulenceStatisticsORACLES::TurbulenceStatisticsORACLES(
  RefCountPtr<DRT::Discretization> discret,
  ParameterList&                   params,
  const string&                    geotype)
  :
  discret_    (discret),
  params_     (params),
  h_          (0.0299),
  x1min_      (-9.0*h_), // (-5.0*h_),
  x1max_      (16.0*h_), //(-6.0*h_),
  x2min_      (-0.0653),//(0.010/2.),//
  x2max_      (+0.0653),//(0.0708/2.),//
  x2inflowmin_(-0.0708/2.),
  x2inflowmax_(+0.0708/2.),
  x2inflowchannelmin_(0.010/2.),
  x2inflowchannelmax_(0.0708/2.),
  x3min_      (-0.1505/2.),//(-0.07525/2.),
  x3max_      (+0.1505/2.),//(+0.07525/2.),
  midupchan_  (0.0855),
  midlowchan_ (0.0451),
  midchamber_ (0.0653),
  x2first_    (0.00004),
  x2max_first_(x2max_-x2first_),
  x2min_first_(x2min_+x2first_),
  x2inflowchannelmax_first_(x2inflowchannelmax_-x2first_),
  x2inflowchannelmin_first_(x2inflowchannelmin_+x2first_),
  x3position_ (0.0)
{
  // assumptions
  // - nodes positioned on the mid lines of both inflow channels
  // - nodes positioned on the mid line of the chamber

  //----------------------------
  // initial plausibility checks
  //----------------------------
  {
    if (geotype != "geometry_ORACLES")
      dserror("statistics manager did not receive the ORACLES geometry");

    // 3D problem expected
    const int numdim = params_.get<int>("number of velocity degrees of freedom");
    if (numdim!=3)
      dserror("evaluation of turbulence statistics for ORACLES available only for 3D");

    // type of fluid flow solver
    const INPAR::FLUID::PhysicalType physicaltype = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type");
    if (physicaltype != INPAR::FLUID::incompressible)
      dserror("incompressible flow expected");

    // find min and max coordinates for every dimension on this proc
    double myx1min = +1.0E12;
    double myx1max = -1.0E12;
    double myx2min = +1.0E12;
    double myx2max = -1.0E12;
    double myx3min = +1.0E12;
    double myx3max = -1.0E12;

    // loop all nodes on this proc and check if they are within the given bounds
    for (int inode=0; inode<discret_->NumMyRowNodes(); ++inode)
    {
      DRT::Node* node = discret_->lRowNode(inode);

      static LINALG::Matrix<3,1> x;
      x(0) = node->X()[0];
      x(1) = node->X()[1];
      x(2) = node->X()[2];

      //-----------------------------------------------------------
      // initial check if all nodes are within the ORACLES geometry
      //-----------------------------------------------------------
      if ( !(x(0)>x1min_-NODETOL and x(0)<x1max_+NODETOL and
             x(1)>x2min_-NODETOL and x(1)<x2max_+NODETOL and
             x(2)>x3min_-NODETOL and x(2)<x3max_+NODETOL) )
        dserror("node %d not within bounds of ORACLES geometry",node->Id() );

      if (myx1min>x(0)) myx1min=x(0);
      if (myx1max<x(0)) myx1max=x(0);

      if (myx2min>x(1)) myx2min=x(1);
      if (myx2max<x(1)) myx2max=x(1);

      if (myx3min>x(2)) myx3min=x(2);
      if (myx3max<x(2)) myx3max=x(2);
    }

    //----------------------------------------------------------------------
    // initial check if parallel discretization matches the ORACLES geometry
    //----------------------------------------------------------------------
    // communicate min and max coordinates of every proc to find the global min and max coordinates
    double globx1min=+1.0E12;
    discret_->Comm().MinAll(&myx1min,&globx1min,1);
    if( abs(globx1min-x1min_)>NODETOL ) dserror("global x1 min coordinate does not correspond to ORACLES geometry");

    double globx1max=-1.0E12;
    discret_->Comm().MaxAll(&myx1max,&globx1max,1);
    if( abs(globx1max-x1max_)>NODETOL ) dserror("global x1 max coordinate does not correspond to ORACLES geometry");

    double globx2min=+1.0E12;
    discret_->Comm().MinAll(&myx2min,&globx2min,1);
    if( abs(globx2min-x2min_)>NODETOL ) dserror("global x2 min coordinate does not correspond to ORACLES geometry");

    double globx2max=-1.0E12;
    discret_->Comm().MaxAll(&myx2max,&globx2max,1);
    if( abs(globx2max-x2max_)>NODETOL ) dserror("global x2 max coordinate does not correspond to ORACLES geometry");

    double globx3min=+1.0E12;
    discret_->Comm().MinAll(&myx3min,&globx3min,1);
    if( abs(globx3min-x3min_)>NODETOL ) dserror("global x3 min coordinate does not correspond to ORACLES geometry");

    double globx3max=-1.0E12;
    discret_->Comm().MaxAll(&myx3max,&globx3max,1);
    if( abs(globx3max-x3max_)>NODETOL ) dserror("global x3 max coordinate does not correspond to ORACLES geometry");
  }

  {
    // get fluid viscosity from material definition --- for computation of ltau
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id<0) dserror("could not find Newtonian fluid material");
    else
    {
      const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
      // we need the kinematic viscosity here
      dens_ = actmat->density_;
      visc_ = actmat->viscosity_/actmat->density_;
    }
  }

  // initialize counters
  numsamp_ = 0;
  countrecord_ = 0;

  //--------------------------------------------
  // define positions for statistical evaluation
  //--------------------------------------------
  {
    x1positions_.Clear();

    // get x2-coordinates from x2_inflow_
    x1positions_(0)  = x1min_;           // -5.*h_
    x1positions_(1)  = -(0.0704+0.0397); // begin ramp
    // -3.*h_ located on splitter ramp -> not interesting
    // get x2-coordinates from x2_0h_
    x1positions_(2)  = -2.*h_;
    x1positions_(3)  = -1.*h_;
    x1positions_(4)  = 0.*h_;
    // get x2-coordinates from x2_1h_
    x1positions_(5)  = 1.*h_;
    // get x2-coordinates from x2_2h_
    x1positions_(6)  = 2.*h_;
    // get x2-coordinates from x2_3h_
    x1positions_(7)  = 3.*h_;
    x1positions_(8)  = 4.*h_;
    x1positions_(9)  = 5.*h_;
    x1positions_(10) = 6.*h_;
    x1positions_(11) = 7.*h_;
    x1positions_(12) = 8.*h_;
    x1positions_(13) = 9.*h_;
    x1positions_(14) = 10.*h_;
    x1positions_(15) = 11.*h_;
    x1positions_(16) = 12.*h_;
    x1positions_(17) = 13.*h_;
    x1positions_(18) = 14.*h_;
    x1positions_(19) = 15.*h_;
    x1positions_(20) = 16.*h_;
  }

  //------------------------
  // initialize some vectors
  //------------------------
  {
    // allocate some toggle vectors used to compute sums via scalar products
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    toggleu_ = LINALG::CreateVector(*dofrowmap,true);
    togglev_ = LINALG::CreateVector(*dofrowmap,true);
    togglew_ = LINALG::CreateVector(*dofrowmap,true);
    togglep_ = LINALG::CreateVector(*dofrowmap,true);

    vel_   = LINALG::CreateVector(*dofrowmap,true);
    force_ = LINALG::CreateVector(*dofrowmap,true);
  }

  //----------------------------------------------------------------------
  // create vectors of coordinates at positions for statistical evaluation
  //----------------------------------------------------------------------

  // the sort criterion allows differences in coordinates by 1e-9
  set<double,LineSortCriterion> x1_midupperchannel;
  set<double,LineSortCriterion> x1_midlowerchannel;
  set<double,LineSortCriterion> x1_midchamber;
  set<double,LineSortCriterion> x1_walltopchamber;
  set<double,LineSortCriterion> x1_wallbottomchamber;
  set<double,LineSortCriterion> x1_wallinflowchannel;

  set<double,LineSortCriterion> x2_inflow;
  set<double,LineSortCriterion> x2_0h;
  set<double,LineSortCriterion> x2_1h;
  set<double,LineSortCriterion> x2_2h;
  set<double,LineSortCriterion> x2_3h;

  for (int inode=0; inode<discret_->NumMyRowNodes(); ++inode)
  {
    DRT::Node* node = discret_->lRowNode(inode);

    static LINALG::Matrix<3,1> x;
    x(0) = node->X()[0];
    x(1) = node->X()[1];
    x(2) = node->X()[2];

    if (x(1) > midupchan_-NODETOL  and x(1) < midupchan_ +NODETOL and x(0) < NODETOL)
      x1_midupperchannel.insert(x(0));
    if (x(1) > midlowchan_-NODETOL and x(1) < midlowchan_+NODETOL and x(0) < NODETOL)
      x1_midlowerchannel.insert(x(0));
    if (x(1) > midchamber_-NODETOL and x(1) < midchamber_+NODETOL and x(0) > -NODETOL)
      x1_midchamber.insert(x(0));
    if (x(1) > x2max_-NODETOL and x(1) < x2max_+NODETOL and x(0) > -NODETOL)
      x1_walltopchamber.insert(x(0));
    if (x(1) > x2min_-NODETOL and x(1) < x2min_+NODETOL and x(0) > -NODETOL)
      x1_wallbottomchamber.insert(x(0));
    if (x(1) > x2inflowchannelmax_-NODETOL and x(1) < x2inflowchannelmax_+NODETOL and x(0) < -6.0*h_+NODETOL)
      x1_wallinflowchannel.insert(x(0));


    if (x(0) > x1min_-NODETOL and x(0) < x1min_+NODETOL) // x1 = -5h
      x2_inflow.insert(x(1));
    if (x(0) > -NODETOL and x(0)<NODETOL and x(1) > x2inflowmin_-NODETOL and x(1) < x2inflowmax_+NODETOL ) // x1 = 0h
      x2_0h.insert(x(1));
    if (x(0) > h_-NODETOL and x(0) < h_+NODETOL) // x1 = h
      x2_1h.insert(x(1));
    if (x(0) > 2.*h_-NODETOL and x(0) < 2.*h_+NODETOL) // x1 = 2h
      x2_2h.insert(x(1));
    if (x(0) > 3.*h_-NODETOL and x(0) < 3.*h_+NODETOL) // x1 = 3h
      x2_3h.insert(x(1));
  }

  // send averaging locations around until they are known on all procs
  ExportLocation(x1_midupperchannel);
  ExportLocation(x1_midlowerchannel);
  ExportLocation(x1_midchamber);
  ExportLocation(x1_walltopchamber);
  ExportLocation(x1_wallbottomchamber);
  ExportLocation(x1_wallinflowchannel);

  ExportLocation(x2_inflow);
  ExportLocation(x2_0h);
  ExportLocation(x2_1h);
  ExportLocation(x2_2h);
  ExportLocation(x2_3h);

  //----------------------------
  // fill vectors of coordinates
  //----------------------------
  {
    // horizontal line vectors
    x1_midupperchannel_ = rcp(new vector<double> ); // x1-coordinates of nodes at y = 85.5mm = 2.86h (mid line upper channel)
    x1_midlowerchannel_ = rcp(new vector<double> ); // x1-coordinates of nodes at y = 45.1mm = 1.51h (mid line lower channel)
    x1_midchamber_      = rcp(new vector<double> ); // x1-coordinates of nodes at y = 65.3mm = 2.18h (mid line chamber)
    x1_wallchamber_     = rcp(new vector<double> ); // x1-coordinates of nodes at the walls
    x1_wallinflowchannel_ = rcp(new vector<double> ); // x1-coordinates of nodes at the walls

    // vertical line vectors
    x2_inflow_ = rcp(new vector<double> ); // x2-coordinates of nodes at x =-5h (both inflow channels)
    x2_0h_     = rcp(new vector<double> ); // x2-coordinates of nodes at x = 0h (step/expansion)
    x2_1h_     = rcp(new vector<double> ); // x2-coordinates of nodes at x = 1h
    x2_2h_     = rcp(new vector<double> ); // x2-coordinates of nodes at x = 2h
    x2_3h_     = rcp(new vector<double> ); // x2-coordinates of nodes at x > 3h

    for(set<double,LineSortCriterion>::const_iterator coorditer=x1_midupperchannel.begin(); coorditer!=x1_midupperchannel.end(); ++coorditer)
      x1_midupperchannel_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x1_midlowerchannel.begin(); coorditer!=x1_midlowerchannel.end(); ++coorditer)
      x1_midlowerchannel_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x1_midchamber.begin(); coorditer!=x1_midchamber.end(); ++coorditer)
      x1_midchamber_->push_back(*coorditer);

    //if (x1_walltopchamber.size() != x1_wallbottomchamber.size()) dserror("grid mismatch between top and bottom wall");
    for(set<double,LineSortCriterion>::const_iterator coorditer=x1_walltopchamber.begin(); coorditer!=x1_walltopchamber.end(); ++coorditer)
      x1_wallchamber_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x1_wallinflowchannel.begin(); coorditer!=x1_wallinflowchannel.end(); ++coorditer)
      x1_wallinflowchannel_->push_back(*coorditer);

    Teuchos::RCP<vector<double> > x1_wallchamberbottom;
    x1_wallchamberbottom = rcp(new vector<double> );
    for(set<double,LineSortCriterion>::const_iterator coorditer=x1_wallbottomchamber.begin(); coorditer!=x1_wallbottomchamber.end(); ++coorditer)
      x1_wallchamberbottom->push_back(*coorditer);
    if ((*x1_wallchamber_).size() != (*x1_wallchamberbottom).size()) dserror("grid mismatch between top and bottom wall");
    for (size_t ipos=0; ipos<(*x1_wallchamber_).size(); ++ipos)
    {
      if ( ((*x1_wallchamber_)[ipos]-(*x1_wallchamberbottom)[ipos])>1.0E-12 ) dserror("x1-coordinates at walls do not match");
    }


    for(set<double,LineSortCriterion>::const_iterator coorditer=x2_inflow.begin(); coorditer!=x2_inflow.end(); ++coorditer)
      x2_inflow_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x2_0h.begin(); coorditer!=x2_0h.end(); ++coorditer)
      x2_0h_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x2_1h.begin(); coorditer!=x2_1h.end(); ++coorditer)
      x2_1h_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x2_2h.begin(); coorditer!=x2_2h.end(); ++coorditer)
      x2_2h_->push_back(*coorditer);
    for(set<double,LineSortCriterion>::const_iterator coorditer=x2_3h.begin(); coorditer!=x2_3h.end(); ++coorditer)
      x2_3h_->push_back(*coorditer);
  }

  //---------------------------------------------------
  // check number of coordinates along lines to average
  //---------------------------------------------------
  {
    const size_t numx1_midupperchannel = x1_midupperchannel_->size();
    const size_t numx1_midlowerchannel = x1_midlowerchannel_->size();
    if (numx1_midupperchannel != numx1_midlowerchannel) dserror("grid mismatch in upper and lower inflow channel");
    numx1_midinflowchannel_ = numx1_midupperchannel;
    numx1_midchamber_ = x1_midchamber_->size();
    numx1_wallchamber_ = x1_wallchamber_->size();
    numx1_wallinflowchannel_ = x1_wallinflowchannel_->size();

    numx2_inflow_ = x2_inflow_->size();
    numx2_0h_ = x2_0h_->size();
    // number of nodes at inflows -1 = number of nodes at expansion, since two nodes coincide at edge of splitter plate
//TODO put check back in
    //    if (numx2_inflow_-1 != numx2_0h_) dserror("structured grid expected between inflow and expansion");

    const size_t numx2_1h = x2_1h_->size();
    const size_t numx2_2h = x2_2h_->size();
    const size_t numx2_3h = x2_3h_->size();
    if ( !( numx2_1h == numx2_2h and numx2_1h == numx2_3h) ) dserror("structured grid expected in chamber");
    numx2_chamber_ = numx2_3h;
  }

  //----------------------------------------------
  // allocate arrays holding time mean Cs profiles
  //----------------------------------------------
  wallforceu_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallforcev_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallforcew_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));

  wallvelu_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallvelv_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallvelw_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));

  wallinflowchannelforceu_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelforcev_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelforcew_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));

  wallinflowchannelvelu_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelvelv_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelvelw_ = rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));

  //------------------------------------------------------------------
  // allocate arrays holding time mean profiles of first order moments
  //------------------------------------------------------------------
  // time average of vertical profiles in inflow zone (-5h, -4h)
  vertinflowu_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowv_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinfloww_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowp_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));

  vertmixingu_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingv_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingw_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingp_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));

  vert1hu_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hv_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hw_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hp_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));

  vert2hu_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hv_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hw_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hp_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));

  vertchamberu_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchamberv_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchamberw_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchamberp_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));

  //-------------------------------------------------------------------
  // allocate arrays holding time mean profiles of second order moments
  //-------------------------------------------------------------------
  vertinflowuu_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowvv_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowww_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowpp_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));

  vertmixinguu_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingvv_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingww_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingpp_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));

  vert1huu_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hvv_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hww_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hpp_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));

  vert2huu_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hvv_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hww_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hpp_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));

  vertchamberuu_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchambervv_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchamberww_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchamberpp_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));

  //-------------------------------------------------------------------
  // allocate arrays holding time mean profiles of second order moments
  //-------------------------------------------------------------------
  vertinflowuv_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowuw_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));
  vertinflowvw_ = rcp(new LINALG::SerialDenseMatrix(2,x2_inflow_->size(),true));

  vertmixinguv_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixinguw_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));
  vertmixingvw_ = rcp(new LINALG::SerialDenseMatrix(3,x2_0h_->size(),true));

  vert1huv_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1huw_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));
  vert1hvw_ = rcp(new LINALG::SerialDenseMatrix(1,x2_1h_->size(),true));

  vert2huv_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2huw_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));
  vert2hvw_ = rcp(new LINALG::SerialDenseMatrix(1,x2_2h_->size(),true));

  vertchamberuv_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchamberuw_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));
  vertchambervw_ = rcp(new LINALG::SerialDenseMatrix(14,x2_3h_->size(),true));

  return;
}


/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::TurbulenceStatisticsORACLES::~TurbulenceStatisticsORACLES()
{
  return;
}


/*------------------------------------------------------------------------------------------------*
 | take a sample from the flow field (first and second order momentum profiles, Cs)   henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::DoTimeSample(
  Teuchos::RCP<Epetra_Vector> velnp,
  Teuchos::RCP<Epetra_Vector> residual
  )
{
  if (discret_->Comm().MyPID()==0)
    std::cout << "---  sampling for ORACLES flow statistics ... " << std::flush;

  //------------------------
  // increase sample counter
  //------------------------
  numsamp_++;

  // meanvelnp is a refcount copy of velnp
  vel_->Update(1.0,*velnp,0.0);
  force_->Update(1.0,*residual,0.0);

  //----------------------------------------------------------------------
  // loop planes and calculate integral means in each plane

  this->ExtractProfiles();

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane

  //this->EvaluatePointwiseMeanValuesInPlanes();

  //--------------------------------------------------
  // compute forces on top and bottom walls of chamber
  //--------------------------------------------------
  {
    // top and bottom planes in chamber
    const size_t numplanes = 2;

    vector<double> forceplanes(2);
    forceplanes[0] = x2min_;
    forceplanes[1] = x2max_;

    vector<double> velplanes(2);
    velplanes[0] = x2min_first_;
    velplanes[1] = x2max_first_;

    const int dim = 1; // x2-coordinate

    //// top and bottom planes in inflow zone
    //vector<double> xzplanesinflow(2);
    //xzplanesinflow[0] = x2inflowmin_;
    //xzplanesinflow[1] = x2inflowmax_;
    //int dim = 1; // x2-coordinate

    vector<double> locforceu(numx1_wallchamber_);
    vector<double> locforcev(numx1_wallchamber_);
    vector<double> locforcew(numx1_wallchamber_);

    vector<double> globforceu(numx1_wallchamber_);
    vector<double> globforcev(numx1_wallchamber_);
    vector<double> globforcew(numx1_wallchamber_);

    vector<double> locvelu(numx1_wallchamber_);
    vector<double> locvelv(numx1_wallchamber_);
    vector<double> locvelw(numx1_wallchamber_);

    vector<double> globvelu(numx1_wallchamber_);
    vector<double> globvelv(numx1_wallchamber_);
    vector<double> globvelw(numx1_wallchamber_);

    for (size_t iplane=0; iplane<numplanes; ++iplane)
    {
      // toggle vectors are one in the position of a dof in this plane, else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);

      for (unsigned i=0;i<numx1_wallchamber_;++i)
      {
        locforceu[i] = 0.0;
        locforcev[i] = 0.0;
        locforcew[i] = 0.0;
        globforceu[i]= 0.0;
        globforcev[i]= 0.0;
        globforcew[i]= 0.0;

        locvelu[i] = 0.0;
        locvelv[i] = 0.0;
        locvelw[i] = 0.0;
        globvelu[i]= 0.0;
        globvelv[i]= 0.0;
        globvelw[i]= 0.0;
      }

      for (size_t ipos=0; ipos<numx1_wallchamber_; ++ipos)
      {
        for (int inode=0; inode<discret_->NumMyRowNodes(); ++inode)
        {
          DRT::Node* node = discret_->lRowNode(inode);

          if ( node->X()[2] < x3position_+NODETOL and node->X()[2] > x3position_-NODETOL )
          {
            if (node->X()[0]>(*x1_wallchamber_)[ipos]-NODETOL and node->X()[0]<(*x1_wallchamber_)[ipos]+NODETOL)
            {
              // this node belongs to the plane under consideration
              if ( node->X()[dim] > forceplanes[iplane]-NODETOL and node->X()[dim]<forceplanes[iplane]+NODETOL )
              {
                vector<int> dof = discret_->Dof(node);

                // extract local values from the global vector
                vector<double> myforce(dof.size());
                DRT::UTILS::ExtractMyValues(*force_, myforce, dof);

                locforceu[ipos]+=myforce[0];
                locforcev[ipos]+=myforce[1];
                locforcew[ipos]+=myforce[2];

                double      one = 1.0;
                toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
                //togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
                //togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
              }
              else if ( node->X()[dim] > velplanes[iplane]-NODETOL and node->X()[dim]<velplanes[iplane]+NODETOL )
              {
                vector<int> dof = discret_->Dof(node);

                // extract local values from the global vector
                vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locvelu[ipos]+=myvel[0];
                locvelv[ipos]+=myvel[1];
                locvelw[ipos]+=myvel[2];
              }
            }
          }
        }
      }

      // TODO comment out safety check
//      double inc=0.0;
//      {
//        double local_inc=0.0;
//        for(int rr=0;rr<(*toggleu_).MyLength();++rr)
//        {
//          local_inc+=(*toggleu_)[rr]*(*toggleu_)[rr];
//        }
//        discret_->Comm().SumAll(&local_inc,&inc,1);
//
//        if(abs(inc)<1.0E-9)
//          dserror("there are no forced nodes on the boundary\n");
//      }

      discret_->Comm().SumAll(&(locforceu[0]),&(globforceu[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locforcev[0]),&(globforcev[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locforcew[0]),&(globforcew[0]),numx1_wallchamber_);

      discret_->Comm().SumAll(&(locvelu[0]),&(globvelu[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locvelv[0]),&(globvelv[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locvelw[0]),&(globvelw[0]),numx1_wallchamber_);

      for (size_t ipos=0; ipos<numx1_wallchamber_; ++ipos)
      {
        (*wallforceu_)(iplane,ipos) += globforceu[ipos];
        (*wallforcev_)(iplane,ipos) += globforcev[ipos];
        (*wallforcew_)(iplane,ipos) += globforcew[ipos];

        (*wallvelu_)(iplane,ipos) += globvelu[ipos];
        (*wallvelv_)(iplane,ipos) += globvelv[ipos];
        (*wallvelw_)(iplane,ipos) += globvelw[ipos];
      }
    }
  }

  //----------------------------------------------------------
  // compute forces on top and bottom walls of inflow channels
  //----------------------------------------------------------
  {
    // top and bottom planes in chamber
    const size_t numplanes = 2;

    vector<double> forceplanes(2);
    forceplanes[0] = x2inflowchannelmin_;
    forceplanes[1] = x2inflowchannelmax_;

    vector<double> velplanes(2);
    velplanes[0] = x2inflowchannelmin_first_;
    velplanes[1] = x2inflowchannelmax_first_;

    const int dim = 1; // x2-coordinate

    //// top and bottom planes in inflow zone
    //vector<double> xzplanesinflow(2);
    //xzplanesinflow[0] = x2inflowmin_;
    //xzplanesinflow[1] = x2inflowmax_;
    //int dim = 1; // x2-coordinate

    vector<double> locforceu(numx1_wallinflowchannel_);
    vector<double> locforcev(numx1_wallinflowchannel_);
    vector<double> locforcew(numx1_wallinflowchannel_);

    vector<double> globforceu(numx1_wallinflowchannel_);
    vector<double> globforcev(numx1_wallinflowchannel_);
    vector<double> globforcew(numx1_wallinflowchannel_);

    vector<double> locvelu(numx1_wallinflowchannel_);
    vector<double> locvelv(numx1_wallinflowchannel_);
    vector<double> locvelw(numx1_wallinflowchannel_);

    vector<double> globvelu(numx1_wallinflowchannel_);
    vector<double> globvelv(numx1_wallinflowchannel_);
    vector<double> globvelw(numx1_wallinflowchannel_);

    for (size_t iplane=0; iplane<numplanes; ++iplane)
    {
      // toggle vectors are one in the position of a dof in this plane, else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);

      for (unsigned i=0;i<numx1_wallinflowchannel_;++i)
      {
        locforceu[i] = 0.0;
        locforcev[i] = 0.0;
        locforcew[i] = 0.0;
        globforceu[i]= 0.0;
        globforcev[i]= 0.0;
        globforcew[i]= 0.0;

        locvelu[i] = 0.0;
        locvelv[i] = 0.0;
        locvelw[i] = 0.0;
        globvelu[i]= 0.0;
        globvelv[i]= 0.0;
        globvelw[i]= 0.0;
      }

      for (size_t ipos=0; ipos<numx1_wallinflowchannel_; ++ipos)
      {
        for (int inode=0; inode<discret_->NumMyRowNodes(); ++inode)
        {
          DRT::Node* node = discret_->lRowNode(inode);

          if ( node->X()[2] < x3position_+NODETOL and node->X()[2] > x3position_-NODETOL )
          {
            if (node->X()[0]>(*x1_wallinflowchannel_)[ipos]-NODETOL and node->X()[0]<(*x1_wallinflowchannel_)[ipos]+NODETOL)
            {
              // this node belongs to the plane under consideration
              if ( node->X()[dim] > forceplanes[iplane]-NODETOL and node->X()[dim]<forceplanes[iplane]+NODETOL )
              {
                vector<int> dof = discret_->Dof(node);

                // extract local values from the global vector
                vector<double> myforce(dof.size());
                DRT::UTILS::ExtractMyValues(*force_, myforce, dof);

                locforceu[ipos]+=myforce[0];
                locforcev[ipos]+=myforce[1];
                locforcew[ipos]+=myforce[2];

                double      one = 1.0;
                toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
                //togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
                //togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
              }
              else if ( node->X()[dim] > velplanes[iplane]-NODETOL and node->X()[dim]<velplanes[iplane]+NODETOL )
              {
                vector<int> dof = discret_->Dof(node);

                // extract local values from the global vector
                vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locvelu[ipos]+=myvel[0];
                locvelv[ipos]+=myvel[1];
                locvelw[ipos]+=myvel[2];
              }
            }
          }
        }
      }

      // TODO comment out safety check
//      double inc=0.0;
//      {
//        double local_inc=0.0;
//        for(int rr=0;rr<(*toggleu_).MyLength();++rr)
//        {
//          local_inc+=(*toggleu_)[rr]*(*toggleu_)[rr];
//        }
//        discret_->Comm().SumAll(&local_inc,&inc,1);
//
//        if(abs(inc)<1.0E-9)
//          dserror("there are no forced nodes on the boundary\n");
//      }

      discret_->Comm().SumAll(&(locforceu[0]),&(globforceu[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locforcev[0]),&(globforcev[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locforcew[0]),&(globforcew[0]),numx1_wallinflowchannel_);

      discret_->Comm().SumAll(&(locvelu[0]),&(globvelu[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locvelv[0]),&(globvelv[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locvelw[0]),&(globvelw[0]),numx1_wallinflowchannel_);

      for (size_t ipos=0; ipos<numx1_wallinflowchannel_; ++ipos)
      {
        (*wallinflowchannelforceu_)(iplane,ipos) += globforceu[ipos];
        (*wallinflowchannelforcev_)(iplane,ipos) += globforcev[ipos];
        (*wallinflowchannelforcew_)(iplane,ipos) += globforcew[ipos];

        (*wallinflowchannelvelu_)(iplane,ipos) += globvelu[ipos];
        (*wallinflowchannelvelv_)(iplane,ipos) += globvelv[ipos];
        (*wallinflowchannelvelw_)(iplane,ipos) += globvelw[ipos];
      }
    }
  }

  if (discret_->Comm().MyPID()==0)
    std::cout << "done" << std::endl;
}


/*------------------------------------------------------------------------------------------------*
 | extract profiles at specified locations                                            henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::ExtractProfiles()
{
  //----------------------------------------
  // select vertical profiles in inflow zone
  //----------------------------------------
  {
    vector<double> x1locations(2);
    x1locations[0] = x1positions_(0); // x1=-5h
    x1locations[1] = x1positions_(1); // x1=-4h

    ExtractSetOfProfiles(x1locations, *x2_inflow_, *vertinflowu_,  *vertinflowv_,  *vertinfloww_, *vertinflowp_,
                                                   *vertinflowuu_, *vertinflowvv_, *vertinflowww_, *vertinflowpp_,
                                                   *vertinflowuv_, *vertinflowuw_, *vertinflowvw_);
  }

  //----------------------------------------
  // select vertical profiles in mixing zone
  //----------------------------------------
  {
    vector<double> x1locations(3);
    x1locations[0] = x1positions_(2); // x1=-2h
    x1locations[1] = x1positions_(3); // x1=-1h
    x1locations[2] = x1positions_(4); // x1= 0h

    ExtractSetOfProfiles(x1locations, *x2_0h_, *vertmixingu_,  *vertmixingv_,  *vertmixingw_, *vertmixingp_,
                                               *vertmixinguu_, *vertmixingvv_, *vertmixingww_, *vertmixingpp_,
                                               *vertmixinguv_, *vertmixinguw_, *vertmixingvw_);
  }

  //-----------------------------------------------------
  // select vertical profiles at 1h in recirculation zone
  //-----------------------------------------------------
  {
    vector<double> x1locations(1);
    x1locations[0] = x1positions_(5); // x1= 1h

    ExtractSetOfProfiles(x1locations, *x2_1h_, *vert1hu_,  *vert1hv_,  *vert1hw_, *vert1hp_,
                                               *vert1huu_, *vert1hvv_, *vert1hww_, *vert1hpp_,
                                               *vert1huv_, *vert1huw_, *vert1hvw_);
  }

  //-----------------------------------------------------
  // select vertical profiles at 2h in recirculation zone
  //-----------------------------------------------------
  {
    vector<double> x1locations(1);
    x1locations[0] = x1positions_(6); // x1= 2h

    ExtractSetOfProfiles(x1locations, *x2_2h_, *vert2hu_,  *vert2hv_,  *vert2hw_, *vert2hp_,
                                               *vert2huu_, *vert2hvv_, *vert2hww_, *vert2hpp_,
                                               *vert2huv_, *vert2huw_, *vert2hvw_);
  }

  //------------------------------------
  // select vertical profiles in chamber
  //------------------------------------
  {
    vector<double> x1locations(14);
    x1locations[0]  = x1positions_(7);  // x1= 3h;
    x1locations[1]  = x1positions_(8);  // x1= 4h;
    x1locations[2]  = x1positions_(9);  // x1= 5h;
    x1locations[3]  = x1positions_(10); // x1= 6h;
    x1locations[4]  = x1positions_(11); // x1= 7h;
    x1locations[5]  = x1positions_(12); // x1= 8h;
    x1locations[6]  = x1positions_(13); // x1= 9h;
    x1locations[7]  = x1positions_(14); // x1= 10h;
    x1locations[8]  = x1positions_(15); // x1= 11h;
    x1locations[9]  = x1positions_(16); // x1= 12h;
    x1locations[10] = x1positions_(17); // x1= 13h;
    x1locations[11] = x1positions_(18); // x1= 14h;
    x1locations[12] = x1positions_(19); // x1= 15h;
    x1locations[13] = x1positions_(20); // x1= 16h;

    ExtractSetOfProfiles(x1locations, *x2_3h_, *vertchamberu_,  *vertchamberv_,  *vertchamberw_, *vertchamberp_,
                                               *vertchamberuu_, *vertchambervv_, *vertchamberww_, *vertchamberpp_,
                                               *vertchamberuv_, *vertchamberuw_, *vertchambervw_);
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | extract a set of profiles at specified locations within a homogeneous grid zone    henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::ExtractSetOfProfiles(
    const vector<double>&      x1locations,
    const vector<double>&      x2locations,
    LINALG::SerialDenseMatrix& profilesu,
    LINALG::SerialDenseMatrix& profilesv,
    LINALG::SerialDenseMatrix& profilesw,
    LINALG::SerialDenseMatrix& profilesp,
    LINALG::SerialDenseMatrix& profilesuu,
    LINALG::SerialDenseMatrix& profilesvv,
    LINALG::SerialDenseMatrix& profilesww,
    LINALG::SerialDenseMatrix& profilespp,
    LINALG::SerialDenseMatrix& profilesuv,
    LINALG::SerialDenseMatrix& profilesuw,
    LINALG::SerialDenseMatrix& profilesvw
)
{
  // store contributions of each processor (local) to profile
  vector<vector<double> > locu(x1locations.size(),vector<double>(x2locations.size()));
  vector<vector<double> > locv(x1locations.size(),vector<double>(x2locations.size()));
  vector<vector<double> > locw(x1locations.size(),vector<double>(x2locations.size()));
  vector<vector<double> > locp(x1locations.size(),vector<double>(x2locations.size()));

  vector<vector<double> > globu(x1locations.size(),vector<double>(x2locations.size()));
  vector<vector<double> > globv(x1locations.size(),vector<double>(x2locations.size()));
  vector<vector<double> > globw(x1locations.size(),vector<double>(x2locations.size()));
  vector<vector<double> > globp(x1locations.size(),vector<double>(x2locations.size()));

  //vector<double> locu(x1locations.size()*x2locations.size());
  //std::fill(locu.begin(),locu.end(),0);


  for(size_t ix1pos=0; ix1pos<x1locations.size(); ++ix1pos)
  {
    for (size_t irow=0; irow<x1locations.size(); ++irow)
    {
      for (size_t icol=0; icol<x2locations.size(); ++icol)
      {
        locu[irow][icol] = 0.0;
        locv[irow][icol] = 0.0;
        locw[irow][icol] = 0.0;
        locp[irow][icol] = 0.0;

        globu[irow][icol] = 0.0;
        globv[irow][icol] = 0.0;
        globw[irow][icol] = 0.0;
        globp[irow][icol] = 0.0;
      }
    }

    for(size_t ix2pos=0; ix2pos<x2locations.size(); ++ix2pos)
    {
      for (int inode=0; inode<discret_->NumMyRowNodes(); ++inode)
      {
        DRT::Node* node = discret_->lRowNode(inode);

        static LINALG::Matrix<3,1> x;
        x(0) = node->X()[0];
        x(1) = node->X()[1];
        x(2) = node->X()[2];

        if ( x(2)>(x3position_-NODETOL) and x(2)<(x3position_+NODETOL) ) // x3==0;
        {
          if ( x(0)>(x1locations[ix1pos]-NODETOL) and x(0)<(x1locations[ix1pos]+NODETOL) )
          {
            if ( x(1)>(x2locations[ix2pos]-NODETOL) and x(1)<(x2locations[ix2pos]+NODETOL) )
            {

              vector<int> dof = discret_->Dof(node);

              // extract local values from the global vector
              vector<double> myvel(dof.size());
              DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

              locu[ix1pos][ix2pos]=myvel[0];
              locv[ix1pos][ix2pos]=myvel[1];
              locw[ix1pos][ix2pos]=myvel[2];
              locp[ix1pos][ix2pos]=myvel[3];
            }
          }
        }
      }
    }
    discret_->Comm().SumAll( &(locu[ix1pos][0]), &(globu[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locv[ix1pos][0]), &(globv[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locw[ix1pos][0]), &(globw[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locp[ix1pos][0]), &(globp[ix1pos][0]),x2locations.size() );

    for(size_t ix2pos=0; ix2pos<x2locations.size(); ++ix2pos)
    {
      profilesu(ix1pos,ix2pos) += globu[ix1pos][ix2pos];
      profilesv(ix1pos,ix2pos) += globv[ix1pos][ix2pos];
      profilesw(ix1pos,ix2pos) += globw[ix1pos][ix2pos];
      profilesp(ix1pos,ix2pos) += globp[ix1pos][ix2pos];

      profilesuu(ix1pos,ix2pos) += globu[ix1pos][ix2pos]*globu[ix1pos][ix2pos];
      profilesvv(ix1pos,ix2pos) += globv[ix1pos][ix2pos]*globv[ix1pos][ix2pos];
      profilesww(ix1pos,ix2pos) += globw[ix1pos][ix2pos]*globw[ix1pos][ix2pos];
      profilespp(ix1pos,ix2pos) += globp[ix1pos][ix2pos]*globp[ix1pos][ix2pos];

      profilesuv(ix1pos,ix2pos) += globu[ix1pos][ix2pos]*globv[ix1pos][ix2pos];
      profilesuw(ix1pos,ix2pos) += globu[ix1pos][ix2pos]*globw[ix1pos][ix2pos];
      profilesvw(ix1pos,ix2pos) += globv[ix1pos][ix2pos]*globw[ix1pos][ix2pos];
    }

    //std::copy(globu[ix1pos].begin(), globu[ix1pos].end(), profilesu(ix1pos,0));
    //std::copy(globv[ix1pos].begin(), globv[ix1pos].end(), profilesv(ix1pos,0));
    //std::copy(globw[ix1pos].begin(), globw[ix1pos].end(), profilesw(ix1pos,0));
  }

}


/*------------------------------------------------------------------------------------------------*
 | compute time average from all samples belonging to the last record                 henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::TimeAverageStatistics()
{
  if (numsamp_ <= 0)
    dserror("No samples to do time average");

  const double avfac = 1./numsamp_;

  (*wallforceu_).Scale(avfac);
  (*wallforcev_).Scale(avfac);
  (*wallforcew_).Scale(avfac);

  (*wallvelu_).Scale(avfac);
  (*wallvelv_).Scale(avfac);
  (*wallvelw_).Scale(avfac);

  (*wallinflowchannelforceu_).Scale(avfac);
  (*wallinflowchannelforcev_).Scale(avfac);
  (*wallinflowchannelforcew_).Scale(avfac);

  (*wallinflowchannelvelu_).Scale(avfac);
  (*wallinflowchannelvelv_).Scale(avfac);
  (*wallinflowchannelvelw_).Scale(avfac);

  //! first order momentum
  (*vertinflowu_).Scale(avfac);
  (*vertinflowv_).Scale(avfac);
  (*vertinfloww_).Scale(avfac);
  (*vertinflowp_).Scale(avfac);

  (*vertmixingu_).Scale(avfac);
  (*vertmixingv_).Scale(avfac);
  (*vertmixingw_).Scale(avfac);
  (*vertmixingp_).Scale(avfac);

  (*vert1hu_).Scale(avfac);
  (*vert1hv_).Scale(avfac);
  (*vert1hw_).Scale(avfac);
  (*vert1hp_).Scale(avfac);

  (*vert2hu_).Scale(avfac);
  (*vert2hv_).Scale(avfac);
  (*vert2hw_).Scale(avfac);
  (*vert2hp_).Scale(avfac);

  (*vertchamberu_).Scale(avfac);
  (*vertchamberv_).Scale(avfac);
  (*vertchamberw_).Scale(avfac);
  (*vertchamberp_).Scale(avfac);

  //! second order momentum
  (*vertinflowuu_).Scale(avfac);
  (*vertinflowvv_).Scale(avfac);
  (*vertinflowww_).Scale(avfac);
  (*vertinflowpp_).Scale(avfac);

  (*vertmixinguu_).Scale(avfac);
  (*vertmixingvv_).Scale(avfac);
  (*vertmixingww_).Scale(avfac);
  (*vertmixingpp_).Scale(avfac);

  (*vert1huu_).Scale(avfac);
  (*vert1hvv_).Scale(avfac);
  (*vert1hww_).Scale(avfac);
  (*vert1hpp_).Scale(avfac);

  (*vert2huu_).Scale(avfac);
  (*vert2hvv_).Scale(avfac);
  (*vert2hww_).Scale(avfac);
  (*vert2hpp_).Scale(avfac);

  (*vertchamberuu_).Scale(avfac);
  (*vertchambervv_).Scale(avfac);
  (*vertchamberww_).Scale(avfac);
  (*vertchamberpp_).Scale(avfac);

  //! second order momentum
  (*vertinflowuv_).Scale(avfac);
  (*vertinflowuw_).Scale(avfac);
  (*vertinflowvw_).Scale(avfac);

  (*vertmixinguv_).Scale(avfac);
  (*vertmixinguw_).Scale(avfac);
  (*vertmixingvw_).Scale(avfac);

  (*vert1huv_).Scale(avfac);
  (*vert1huw_).Scale(avfac);
  (*vert1hvw_).Scale(avfac);

  (*vert2huv_).Scale(avfac);
  (*vert2huw_).Scale(avfac);
  (*vert2hvw_).Scale(avfac);

  (*vertchamberuv_).Scale(avfac);
  (*vertchamberuw_).Scale(avfac);
  (*vertchambervw_).Scale(avfac);
}


/*------------------------------------------------------------------------------------------------*
 | write statistics of last record into a file                                        henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::OutputStatistics(int step)
{
  // only first proc writes output
  if (discret_->Comm().MyPID()==0)
  {
    //------------------------------------
    // compute u_tau and l_tau (and tau_W)
    //------------------------------------
    // evaluate area
//    LINALG::SerialDenseMatrix area(1,numx1_wallchamber_,true);
//
//    for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
//    {
//      size_t leftpos = ix1pos-1;
//      size_t rightpos = ix1pos+1;
//      double meanfac = 0.5;
//      if (ix1pos==0) // farther left node
//      {
//        leftpos = 0;
//        meanfac = 1.;
//      }
//      if (ix1pos==numx1_wallchamber_) // farther right node
//      {
//        rightpos = numx1_wallchamber_;
//        meanfac = 1.;
//      }
//      // use average length in x1-direction (mean of left and right position)
//      area(0,ix1pos) = (x3max_-x3min_) * meanfac*((*x1_wallchamber_)[rightpos]-(*x1_wallchamber_)[leftpos]);
//    }
    const double area = 1.0;//(x3max_-x3min_) * (x1max_-0.0);

    // nonzero forces (tractions) only expected in the streamwise (x1) direction

    // compute ltau (used to compute y+)
    LINALG::SerialDenseMatrix ltau(2,numx1_wallchamber_,true);
    // reference ltau (at x=10h)
    double ltau_ref = 1.0;

    //TODO clean
    for(size_t iplane=0; iplane<2; ++iplane)
    {
      for (size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        //cout << abs((*wallforceu_)(iplane,ix1pos)) << endl;
        //cout << abs((*wallforcev_)(iplane,ix1pos)) << endl;
        //cout << abs((*wallforcew_)(iplane,ix1pos)) << endl;

        //if ( (abs((*wallforceu_)(iplane,ix1pos)) > abs((*wallforcev_)(iplane,ix1pos))) and
        //     (abs((*wallforceu_)(iplane,ix1pos)) > abs((*wallforcew_)(iplane,ix1pos))) )
        //{
//          if( abs((*wallforceu_)(iplane,ix1pos)) < 1.0E-12 )
//            dserror("zero force during computation of wall shear stress\n");

          //ltau(iplane,ix1pos) = visc_ / (sqrt( (*wallforceu_)(iplane,ix1pos) / (area(0,ix1pos)*dens_) ));
          ltau(iplane,ix1pos) = visc_ / (sqrt( abs((*wallforceu_)(iplane,ix1pos)) / (area*dens_) ));
//cout << "wallforce " << (*wallforceu_)(iplane,ix1pos) << endl;
//cout << "ltau " << ltau(iplane,ix1pos) << endl;

          if ( (x1positions_(ix1pos) > 10*h_-NODETOL) and (x1positions_(ix1pos) < 10*h_+NODETOL) )
          {
            ltau_ref = ltau(iplane,ix1pos);
          }
        //}
        //else
        //{
        //  dserror("Cannot determine flow direction by traction (seems to be not unique)");
        //}
      }
    }

    //---------------------------------------
    // write horizontal ltau profiles chamber
    //---------------------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      // define file name suffix
      s.append(".flow_statistics_horiz");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at bottom wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallforceu_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallforcev_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallforcew_)(0,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at top wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallforceu_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallforcev_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallforcew_)(1,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from bottom wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallvelu_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallvelv_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallvelw_)(0,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from top wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallvelu_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallvelv_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallvelw_)(1,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();
    }
    //-----------------------------------------------
    // write horizontal ltau profiles inflow channels
    //-----------------------------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      // define file name suffix
      s.append(".flow_statistics_inflowchannel_horiz");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at bottom wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelforceu_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelforcev_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelforcew_)(0,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at top wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelforceu_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelforcev_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelforcew_)(1,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from bottom wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelvelu_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelvelv_)(0,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelvelw_)(0,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from top wall \n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelvelu_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelvelv_)(1,ix1pos)/area/dens_;
        (*log) << "   " << setw(11) << setprecision(4) << (*wallinflowchannelvelw_)(1,ix1pos)/area/dens_;
        (*log) << &endl;
      }
      log->flush();
    }
    //-------------------------------
    // write vertical inflow profiles
    //-------------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      s.append(".flow_statistics_vert_inflow");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "#     y            y+";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "        mean u^2      mean v^2      mean w^2";
      (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2 \n";
      (*log) << scientific;

      for(int ix1pos=0; ix1pos<(*vertinflowu_).M(); ++ix1pos)
      {
        (*log) << "\n";
        (*log) << "# profile at x-position " << setprecision(2) << x1positions_(0+ix1pos)/h_ << " h \n";

        for(int ix2pos=(*vertinflowu_).N()-1; ix2pos>=0; --ix2pos)
        {
          // y and y+
          (*log) <<  " "  << setw(11) << setprecision(4) << (*x2_inflow_)[ix2pos];
          (*log) << "   " << setw(11) << setprecision(4) << (*x2_inflow_)[ix2pos]/ltau_ref;

          // time mean values
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowu_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowv_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinfloww_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowp_)(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowuu_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowvv_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowww_)(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowuv_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowuw_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowvw_)(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertinflowpp_)(ix1pos,ix2pos);

          (*log) << "\n";
        }
        log->flush();
      }
    }
    //-------------------------------
    // write vertical mixing profiles
    //-------------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      s.append(".flow_statistics_vert_mixing");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "#     y            y+";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "        mean u^2      mean v^2      mean w^2";
      (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2 \n";
      (*log) << scientific;

      for(int ix1pos=0; ix1pos<(*vertmixingu_).M(); ++ix1pos)
      {
        (*log) << "\n";
        (*log) << "# profile at x-position " << setprecision(2) << x1positions_(2+ix1pos)/h_ << " h \n";

        for(int ix2pos=(*vertmixingu_).N()-1; ix2pos>=0; --ix2pos)
        {
          // y and y+
          (*log) <<  " "  << setw(11) << setprecision(4) << (*x2_0h_)[ix2pos];
          (*log) << "   " << setw(11) << setprecision(4) << (*x2_0h_)[ix2pos]/ltau_ref;

          // time mean values
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingu_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingv_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingw_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingp_)(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixinguu_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingvv_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingww_)(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixinguv_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixinguw_)(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingvw_)(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertmixingpp_)(ix1pos,ix2pos);

          (*log) << "\n";
        }
        log->flush();
      }
    }
    //-----------------------------
    // write vertical profile at 1h
    //-----------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      s.append(".flow_statistics_vert_1h");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "#     y            y+";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "        mean u^2      mean v^2      mean w^2";
      (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2 \n";
      (*log) << scientific;

      for(int ix1pos=0; ix1pos<(*vert1hu_).M(); ++ix1pos)
      {
        (*log) << "\n";
        (*log) << "# profile at x-position " << setprecision(2) << x1positions_(5+ix1pos)/h_ << " h \n";

        for(int ix2pos=(*vert1hu_).N()-1; ix2pos>=0; --ix2pos)
        {
          // y and y+
          (*log) <<  " "  << setw(11) << setprecision(4) << (*x2_1h_)[ix2pos];
          (*log) << "   " << setw(11) << setprecision(4) << (*x2_1h_)[ix2pos]/ltau_ref;

          // time mean values
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hu_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hw_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hp_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vert1huu_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hvv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hww_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vert1huv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1huw_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hvw_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vert1hpp_ )(ix1pos,ix2pos);

          (*log) << "\n";
        }
        log->flush();
      }
    }
    //-----------------------------
    // write vertical profile at 2h
    //-----------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      s.append(".flow_statistics_vert_2h");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "#     y            y+";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "        mean u^2      mean v^2      mean w^2";
      (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2 \n";
      (*log) << scientific;

      for(int ix1pos=0; ix1pos<(*vert2hu_).M(); ++ix1pos)
      {
        (*log) << "\n";
        (*log) << "# profile at x-position " << setprecision(2) << x1positions_(6+ix1pos)/h_ << " h \n";

        for(int ix2pos=(*vert2hu_).N()-1; ix2pos>=0; --ix2pos)
        {
          // y and y+
          (*log) <<  " "  << setw(11) << setprecision(4) << (*x2_2h_)[ix2pos];
          (*log) << "   " << setw(11) << setprecision(4) << (*x2_2h_)[ix2pos]/ltau_ref;

          // time mean values
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hu_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hw_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hp_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vert2huu_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hvv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hww_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vert2huv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2huw_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hvw_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vert2hpp_ )(ix1pos,ix2pos);

          (*log) << "\n";
        }
        log->flush();
      }
    }
    //---------------------------------
    // write vertical chamber profiles
    //---------------------------------
    {
      std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
      s.append(".flow_statistics_vert_chamber");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
      (*log) << "\n\n\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "#     y            y+";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "        mean u^2      mean v^2      mean w^2";
      (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2 \n";
      (*log) << scientific;

      for(int ix1pos=0; ix1pos<(*vertchamberu_).M(); ++ix1pos)
      {
        (*log) << "\n";
        (*log) << "# profile at x-position " << setprecision(2) << x1positions_(7+ix1pos)/h_ << " h \n";

        for(int ix2pos=(*vertchamberu_).N()-1; ix2pos>=0; --ix2pos)
        {
          // y and y+
          (*log) <<  " "  << setw(11) << setprecision(4) << (*x2_3h_)[ix2pos];
          (*log) << "   " << setw(11) << setprecision(4) << (*x2_3h_)[ix2pos]/ltau_ref;

          // time mean values
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberu_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberw_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberp_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberuu_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchambervv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberww_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberuv_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberuw_ )(ix1pos,ix2pos);
          (*log) << "   " << setw(11) << setprecision(4) << (*vertchambervw_ )(ix1pos,ix2pos);

          (*log) << "   " << setw(11) << setprecision(4) << (*vertchamberpp_ )(ix1pos,ix2pos);

          (*log) << "\n";
        }
        log->flush();
      }
    }

  }

  // increase counter of written records
  countrecord_++;

  return;
}


/*------------------------------------------------------------------------------------------------*
 | reset data structures for a new record                                             henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::ClearStatistics()
{
  //------------------------
  // reset sample counter
  //------------------------
  numsamp_ = 0;

  //-------------
  // reset arrays
  //-------------
  (*wallforceu_).Zero();
  (*wallforcev_).Zero();
  (*wallforcew_).Zero();

  (*wallvelu_).Zero();
  (*wallvelv_).Zero();
  (*wallvelw_).Zero();

  (*wallinflowchannelforceu_).Zero();
  (*wallinflowchannelforcev_).Zero();
  (*wallinflowchannelforcew_).Zero();

  (*wallinflowchannelvelu_).Zero();
  (*wallinflowchannelvelv_).Zero();
  (*wallinflowchannelvelw_).Zero();

  //! first order momentum
  (*vertinflowu_).Zero();
  (*vertinflowv_).Zero();
  (*vertinfloww_).Zero();
  (*vertinflowp_).Zero();

  (*vertmixingu_).Zero();
  (*vertmixingv_).Zero();
  (*vertmixingw_).Zero();
  (*vertmixingp_).Zero();

  (*vert1hu_).Zero();
  (*vert1hv_).Zero();
  (*vert1hw_).Zero();
  (*vert1hp_).Zero();

  (*vert2hu_).Zero();
  (*vert2hv_).Zero();
  (*vert2hw_).Zero();
  (*vert2hp_).Zero();

  (*vertchamberu_).Zero();
  (*vertchamberv_).Zero();
  (*vertchamberw_).Zero();
  (*vertchamberp_).Zero();

  //! second order momentum
  (*vertinflowuu_).Zero();
  (*vertinflowvv_).Zero();
  (*vertinflowww_).Zero();
  (*vertinflowpp_).Zero();

  (*vertmixinguu_).Zero();
  (*vertmixingvv_).Zero();
  (*vertmixingww_).Zero();
  (*vertmixingpp_).Zero();

  (*vert1huu_).Zero();
  (*vert1hvv_).Zero();
  (*vert1hww_).Zero();
  (*vert1hpp_).Zero();

  (*vert2huu_).Zero();
  (*vert2hvv_).Zero();
  (*vert2hww_).Zero();
  (*vert2hpp_).Zero();

  (*vertchamberuu_).Zero();
  (*vertchambervv_).Zero();
  (*vertchamberww_).Zero();
  (*vertchamberpp_).Zero();

  //! second order momentum
  (*vertinflowuv_).Zero();
  (*vertinflowuw_).Zero();
  (*vertinflowvw_).Zero();

  (*vertmixinguv_).Zero();
  (*vertmixinguw_).Zero();
  (*vertmixingvw_).Zero();

  (*vert1huv_).Zero();
  (*vert1huw_).Zero();
  (*vert1hvw_).Zero();

  (*vert2huv_).Zero();
  (*vert2huw_).Zero();
  (*vert2hvw_).Zero();

  (*vertchamberuv_).Zero();
  (*vertchamberuw_).Zero();
  (*vertchambervw_).Zero();

  return;
}

/*------------------------------------------------------------------------------------------------*
 | write statistic data to log file                                                   henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::DumpStatistics(const int step)
{
#if 0
  Teuchos::RCP<std::ofstream> log;

  // only proc 0 writes
  if (discret_->Comm().MyPID()==0)
  {
    // get name of the output *.control file
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# statistics for incompressible flow for turbulent premixed combustion benchmark problem ORACLES (first- and second-order moments) \n";
    (*log) << "\n";
    (*log) << "# statistics record for steps " << step-numsamp_+1 << "--" << step <<"\n";
    (*log) << "\n\n";

    (*log) << scientific;
    (*log) << "# lower wall behind step \n";
    (*log) << "#     x1           duxdy         pmean \n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i=0; i<x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u     = (*x1sumu_)(0,i)/numsamp_;
        double lwx1duxdy = lwx1u/dist;
        double lwx1p     = (*x1sump_)(0,i)/numsamp_;

        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1coordinates_)[i];
        (*log) << "   " << setw(11) << setprecision(4) << lwx1duxdy;
        (*log) << "   " << setw(11) << setprecision(4) << lwx1p;
        (*log) << "\n";
      }
    }

//    if (geotype_ == TurbulenceStatisticsORACLES::geometry_LES_flow_with_heating)
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
//        (*log) <<  " "  << setw(11) << setprecision(4) << (*x1coordinates_)[i];
//        (*log) << "   " << setw(11) << setprecision(4) << uwx1duxdy;
//        (*log) << "   " << setw(11) << setprecision(4) << uwx1p;
//        (*log) << "\n";
//      }
//    }

    for (int i=0; i<numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first (only DNS geometry)
      double x1 = 1.0e20;
      if (i < numx1supplocations_)
      {
        x1 = x1supplocations_(i);
      }
      else
      {
        x1 = x1statlocations_(i - numx1supplocations_);
      }

      (*log) << "\n\n\n";
      (*log) << "# line in x2-direction at x1 = " << setw(11) << setprecision(4) << x1 << "\n";
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

        double x2urms  = sqrt((*x2sumsqu_)(i,j)/numsamp_-x2u*x2u);
        double x2vrms  = sqrt((*x2sumsqv_)(i,j)/numsamp_-x2v*x2v);
        double x2wrms  = sqrt((*x2sumsqw_)(i,j)/numsamp_-x2w*x2w);
        double x2prms  = sqrt((*x2sumsqp_)(i,j)/numsamp_-x2p*x2p);

        double x2uv   = (*x2sumuv_)(i,j)/numsamp_-x2u*x2v;
        double x2uw   = (*x2sumuw_)(i,j)/numsamp_-x2u*x2w;
        double x2vw   = (*x2sumvw_)(i,j)/numsamp_-x2v*x2w;

        (*log) <<  " "  << setw(11) << setprecision(4) << (*x2coordinates_)[j];
        (*log) << "   " << setw(11) << setprecision(4) << x2u;
        (*log) << "   " << setw(11) << setprecision(4) << x2v;
        (*log) << "   " << setw(11) << setprecision(4) << x2w;
        (*log) << "   " << setw(11) << setprecision(4) << x2p;
        (*log) << "   " << setw(11) << setprecision(4) << x2urms;
        (*log) << "   " << setw(11) << setprecision(4) << x2vrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2wrms;
        (*log) << "   " << setw(11) << setprecision(4) << x2prms;
        (*log) << "   " << setw(11) << setprecision(4) << x2uv;
        (*log) << "   " << setw(11) << setprecision(4) << x2uw;
        (*log) << "   " << setw(11) << setprecision(4) << x2vw;
        (*log) << "\n";
      }
    }

    log->flush();
  }
#endif
  return;

}


/*------------------------------------------------------------------------------------------------*
 | export location from this proc to all procs (round robin)                          henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::ExportLocation(std::set<double,LineSortCriterion>& locations)
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

  for (int np=0;np<numprocs;++np)
  {
    DRT::PackBuffer data;

    for (set<double,LineSortCriterion>::iterator line=locations.begin();
        line!=locations.end();
        ++line)
    {
      DRT::ParObject::AddtoPack(data,*line);
    }
    data.StartPacking();
    for (set<double,LineSortCriterion>::iterator line=locations.begin();
        line!=locations.end();
        ++line)
    {
      DRT::ParObject::AddtoPack(data,*line);
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

    // unpack received block into set of all planes.
    {
      vector<double> coordsvec;

      coordsvec.clear();

      vector<char>::size_type index = 0;
      while (index < rblock.size())
      {
        double onecoord;
        DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
        locations.insert(onecoord);
      }
    }
  }

}

#endif
