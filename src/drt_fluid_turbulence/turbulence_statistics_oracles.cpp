/*!-----------------------------------------------------------------------------------------------*
\file turbulence_statistics_oracles.cpp

\brief statistical data processing for ORACLES problem

\level 2

<pre>
\maintainer Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *------------------------------------------------------------------------------------------------*/



#include "turbulence_statistics_oracles.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_utils.H"

#define NODETOL 1e-9


/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::TurbulenceStatisticsORACLES::TurbulenceStatisticsORACLES(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::ParameterList&          params,
  const std::string&               statistics_outfilename,
  const std::string&               geotype,
  const bool                       withscatra)
  :
  discret_    (discret),
  params_     (params),
  statistics_outfilename_(statistics_outfilename),
  withscatra_ (withscatra),
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
  utau_       (0.48),
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

  if (discret_->Comm().MyPID()==0)
  {
    std::cout << "\n" << "parameters used for ORACLES output: "       << std::endl;
    std::cout         << "  utau: "               << utau_            << std::endl;
    std::cout         << "  first y-coordinate: " << x2first_ << "\n" << std::endl;
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
    x1positions_(2)  = 0.0; // not used, since the mesh changesin y-direction here
    x1positions_(3)  = 0.0; // not used, since the mesh changesin y-direction here
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

    vel_   = LINALG::CreateVector(*dofrowmap,true);
    force_ = LINALG::CreateVector(*dofrowmap,true);
  }

  //----------------------------------------------------------------------
  // create vectors of coordinates at positions for statistical evaluation
  //----------------------------------------------------------------------

  // the sort criterion allows differences in coordinates by 1e-9
  std::set<double,LineSortCriterion> x1_midupperchannel;
  std::set<double,LineSortCriterion> x1_midlowerchannel;
  std::set<double,LineSortCriterion> x1_midchamber;
  std::set<double,LineSortCriterion> x1_walltopchamber;
  std::set<double,LineSortCriterion> x1_wallbottomchamber;
  std::set<double,LineSortCriterion> x1_wallinflowchannel;

  std::set<double,LineSortCriterion> x2_inflow;
  std::set<double,LineSortCriterion> x2_0h;
  std::set<double,LineSortCriterion> x2_1h;
  std::set<double,LineSortCriterion> x2_2h;
  std::set<double,LineSortCriterion> x2_3h;

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
    x1_midupperchannel_ = Teuchos::rcp(new std::vector<double> ); // x1-coordinates of nodes at y = 85.5mm = 2.86h (mid line upper channel)
    x1_midlowerchannel_ = Teuchos::rcp(new std::vector<double> ); // x1-coordinates of nodes at y = 45.1mm = 1.51h (mid line lower channel)
    x1_midchamber_      = Teuchos::rcp(new std::vector<double> ); // x1-coordinates of nodes at y = 65.3mm = 2.18h (mid line chamber)
    x1_wallchamber_     = Teuchos::rcp(new std::vector<double> ); // x1-coordinates of nodes at the walls
    x1_wallinflowchannel_ = Teuchos::rcp(new std::vector<double> ); // x1-coordinates of nodes at the walls

    // vertical line vectors
    x2_inflow_ = Teuchos::rcp(new std::vector<double> ); // x2-coordinates of nodes at x =-5h (both inflow channels)
    x2_0h_     = Teuchos::rcp(new std::vector<double> ); // x2-coordinates of nodes at x = 0h (step/expansion)
    x2_1h_     = Teuchos::rcp(new std::vector<double> ); // x2-coordinates of nodes at x = 1h
    x2_2h_     = Teuchos::rcp(new std::vector<double> ); // x2-coordinates of nodes at x = 2h
    x2_3h_     = Teuchos::rcp(new std::vector<double> ); // x2-coordinates of nodes at x > 3h

    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x1_midupperchannel.begin(); coorditer!=x1_midupperchannel.end(); ++coorditer)
      x1_midupperchannel_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x1_midlowerchannel.begin(); coorditer!=x1_midlowerchannel.end(); ++coorditer)
      x1_midlowerchannel_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x1_midchamber.begin(); coorditer!=x1_midchamber.end(); ++coorditer)
      x1_midchamber_->push_back(*coorditer);

    //if (x1_walltopchamber.size() != x1_wallbottomchamber.size()) dserror("grid mismatch between top and bottom wall");
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x1_walltopchamber.begin(); coorditer!=x1_walltopchamber.end(); ++coorditer)
      x1_wallchamber_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x1_wallinflowchannel.begin(); coorditer!=x1_wallinflowchannel.end(); ++coorditer)
      x1_wallinflowchannel_->push_back(*coorditer);

    Teuchos::RCP<std::vector<double> > x1_wallchamberbottom;
    x1_wallchamberbottom = Teuchos::rcp(new std::vector<double> );
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x1_wallbottomchamber.begin(); coorditer!=x1_wallbottomchamber.end(); ++coorditer)
      x1_wallchamberbottom->push_back(*coorditer);
    if ((*x1_wallchamber_).size() != (*x1_wallchamberbottom).size()) dserror("grid mismatch between top and bottom wall");
    for (size_t ipos=0; ipos<(*x1_wallchamber_).size(); ++ipos)
    {
      if ( ((*x1_wallchamber_)[ipos]-(*x1_wallchamberbottom)[ipos])>1.0E-12 ) dserror("x1-coordinates at walls do not match");
    }


    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x2_inflow.begin(); coorditer!=x2_inflow.end(); ++coorditer)
      x2_inflow_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x2_0h.begin(); coorditer!=x2_0h.end(); ++coorditer)
      x2_0h_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x2_1h.begin(); coorditer!=x2_1h.end(); ++coorditer)
      x2_1h_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x2_2h.begin(); coorditer!=x2_2h.end(); ++coorditer)
      x2_2h_->push_back(*coorditer);
    for(std::set<double,LineSortCriterion>::const_iterator coorditer=x2_3h.begin(); coorditer!=x2_3h.end(); ++coorditer)
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

    const size_t numx2_1h = x2_1h_->size();
    const size_t numx2_2h = x2_2h_->size();
    const size_t numx2_3h = x2_3h_->size();
    if ( !( numx2_1h == numx2_2h and numx2_1h == numx2_3h) ) dserror("structured grid expected in chamber");
  }

  //----------------------------------------------
  // allocate arrays holding time mean Cs profiles
  //----------------------------------------------
  wallforceu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallforcev_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallforcew_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallp_      = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));

  wallvelu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallvelv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));
  wallvelw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallchamber_->size(),true));

  wallinflowchannelforceu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelforcev_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelforcew_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelp_      = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));

  wallinflowchannelvelu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelvelv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));
  wallinflowchannelvelw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(2,x1_wallinflowchannel_->size(),true));

  //------------------------------------------------------------------
  // allocate arrays holding time mean profiles of first order moments
  //------------------------------------------------------------------
  // time average of vertical profiles in inflow zone (-5h, -4h)
  vertinflowu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinfloww_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowg_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));

  vertmixingu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingg_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));

  vert1hu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hg_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));

  vert2hu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hg_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));

  vertchamberu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberg_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));

  //-------------------------------------------------------------------
  // allocate arrays holding time mean profiles of second order moments
  //-------------------------------------------------------------------
  vertinflowuu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowvv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowww_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowpp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));

  vertmixinguu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingvv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingww_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingpp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));

  vert1huu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hvv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hww_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hpp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));

  vert2huu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hvv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hww_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hpp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));

  vertchamberuu_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchambervv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberww_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberpp_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));

  //-------------------------------------------------------------------
  // allocate arrays holding time mean profiles of second order moments
  //-------------------------------------------------------------------
  vertinflowuv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowuw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));
  vertinflowvw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_inflow_->size(),2,true));

  vertmixinguv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixinguw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));
  vertmixingvw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_0h_->size(),3,true));

  vert1huv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1huw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));
  vert1hvw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_1h_->size(),1,true));

  vert2huv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2huw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));
  vert2hvw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_2h_->size(),1,true));

  vertchamberuv_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchamberuw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));
  vertchambervw_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(x2_3h_->size(),14,true));

  if (discret_->Comm().MyPID()==0)
  {
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.-5h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.ramp.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.00h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.01h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.02h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.03h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.04h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.05h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.06h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.07h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.08h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.09h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.10h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.11h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.12h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.13h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.14h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.15h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.16h.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES vertical flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.horiz.chamber.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES horizontal chamber flow statistics file\n";
      title.flush();
    }
    {
      std::string outfile(statistics_outfilename_);
      outfile.append(".oracles.horiz.inflow.flow_statistics");
      std::ofstream title(outfile.c_str(),std::ios::out);
      title << "# ORACLES horizontal inflow flow statistics file\n";
      title.flush();
    }
  }

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
  Teuchos::RCP<Epetra_Vector>     velnp,
  Teuchos::RCP<Epetra_Vector>     residual,
  Teuchos::RCP<Epetra_Vector>     phi,
  Teuchos::RCP<const DRT::DofSet> stddofset
  )
{
  if (discret_->Comm().MyPID()==0)
    std::cout << "---  sampling for ORACLES flow statistics ... " << std::flush;

  //------------------------
  // increase sample counter
  //------------------------
  numsamp_++;

  // call from combustion fluid (XFEM) with G-function field
  if (withscatra_)
  {
    if (discret_->Comm().MyPID()==0)
      std::cout << "sampling with scatra... " << std::flush;
    if (velnp != Teuchos::null)
      vel_=velnp;
    else
      dserror("turbulence statistics ORACLES did not receive velnp vector");

    if (residual != Teuchos::null)
      force_ = residual;
    else
      dserror("turbulence statistics ORACLES did not receive residual vector");

    if (phi != Teuchos::null)
      phi_ = phi;
    else
      dserror("turbulence statistics ORACLES did not receive phi vector");

    if (stddofset != Teuchos::null)
      stddofset_ = stddofset;
  }
  // call from standard fluid
  else
  {
    if (discret_->Comm().MyPID()==0)
      std::cout << "sampling without scatra... " << std::flush;
    vel_  ->Update(1.0,*velnp,0.0);
    force_->Update(1.0,*residual,0.0);
    phi_       = Teuchos::null;
    stddofset_ = Teuchos::null;
  }

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

    std::vector<double> forceplanes(2);
    forceplanes[0] = x2min_;
    forceplanes[1] = x2max_;

    std::vector<double> velplanes(2);
    velplanes[0] = x2min_first_;
    velplanes[1] = x2max_first_;

    const int dim = 1; // x2-coordinate

    //// top and bottom planes in inflow zone
    //vector<double> xzplanesinflow(2);
    //xzplanesinflow[0] = x2inflowmin_;
    //xzplanesinflow[1] = x2inflowmax_;
    //int dim = 1; // x2-coordinate

    std::vector<double> locforceu(numx1_wallchamber_);
    std::vector<double> locforcev(numx1_wallchamber_);
    std::vector<double> locforcew(numx1_wallchamber_);
    std::vector<double> locp     (numx1_wallchamber_);

    std::vector<double> globforceu(numx1_wallchamber_);
    std::vector<double> globforcev(numx1_wallchamber_);
    std::vector<double> globforcew(numx1_wallchamber_);
    std::vector<double> globp     (numx1_wallchamber_);

    std::vector<double> locvelu(numx1_wallchamber_);
    std::vector<double> locvelv(numx1_wallchamber_);
    std::vector<double> locvelw(numx1_wallchamber_);

    std::vector<double> globvelu(numx1_wallchamber_);
    std::vector<double> globvelv(numx1_wallchamber_);
    std::vector<double> globvelw(numx1_wallchamber_);

    for (size_t iplane=0; iplane<numplanes; ++iplane)
    {
      for (unsigned i=0;i<numx1_wallchamber_;++i)
      {
        locforceu[i] = 0.0;
        locforcev[i] = 0.0;
        locforcew[i] = 0.0;
        locp[i] = 0.0;
        globforceu[i]= 0.0;
        globforcev[i]= 0.0;
        globforcew[i]= 0.0;
        globp[i] = 0.0;

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
                std::vector<int> dof;
                if (withscatra_)
                  dof = stddofset_->Dof(node);
                else
                  dof = discret_->Dof(node);

                // extract local values from the global vector
                std::vector<double> myforce(dof.size());
                DRT::UTILS::ExtractMyValues(*force_, myforce, dof);
                // extract pressure
                std::vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locforceu[ipos]+=myforce[0];
                locforcev[ipos]+=myforce[1];
                locforcew[ipos]+=myforce[2];
                locp     [ipos]+=myvel[3];
              }
              else if ( node->X()[dim] > velplanes[iplane]-NODETOL and node->X()[dim]<velplanes[iplane]+NODETOL )
              {
                std::vector<int> dof = discret_->Dof(node);

                // extract local values from the global vector
                std::vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locvelu[ipos]+=myvel[0];
                locvelv[ipos]+=myvel[1];
                locvelw[ipos]+=myvel[2];
              }
            }
          }
        }
      }

      discret_->Comm().SumAll(&(locforceu[0]),&(globforceu[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locforcev[0]),&(globforcev[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locforcew[0]),&(globforcew[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locp     [0]),&(globp     [0]),numx1_wallchamber_);

      discret_->Comm().SumAll(&(locvelu[0]),&(globvelu[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locvelv[0]),&(globvelv[0]),numx1_wallchamber_);
      discret_->Comm().SumAll(&(locvelw[0]),&(globvelw[0]),numx1_wallchamber_);

      for (size_t ipos=0; ipos<numx1_wallchamber_; ++ipos)
      {
        (*wallforceu_)(iplane,ipos) += globforceu[ipos];
        (*wallforcev_)(iplane,ipos) += globforcev[ipos];
        (*wallforcew_)(iplane,ipos) += globforcew[ipos];
        (*wallp_)     (iplane,ipos) += globp     [ipos];

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

    std::vector<double> forceplanes(2);
    forceplanes[0] = x2inflowchannelmin_;
    forceplanes[1] = x2inflowchannelmax_;

    std::vector<double> velplanes(2);
    velplanes[0] = x2inflowchannelmin_first_;
    velplanes[1] = x2inflowchannelmax_first_;

    const int dim = 1; // x2-coordinate

    //// top and bottom planes in inflow zone
    //vector<double> xzplanesinflow(2);
    //xzplanesinflow[0] = x2inflowmin_;
    //xzplanesinflow[1] = x2inflowmax_;
    //int dim = 1; // x2-coordinate

    std::vector<double> locforceu(numx1_wallinflowchannel_);
    std::vector<double> locforcev(numx1_wallinflowchannel_);
    std::vector<double> locforcew(numx1_wallinflowchannel_);
    std::vector<double> locp     (numx1_wallinflowchannel_);

    std::vector<double> globforceu(numx1_wallinflowchannel_);
    std::vector<double> globforcev(numx1_wallinflowchannel_);
    std::vector<double> globforcew(numx1_wallinflowchannel_);
    std::vector<double> globp     (numx1_wallinflowchannel_);

    std::vector<double> locvelu(numx1_wallinflowchannel_);
    std::vector<double> locvelv(numx1_wallinflowchannel_);
    std::vector<double> locvelw(numx1_wallinflowchannel_);

    std::vector<double> globvelu(numx1_wallinflowchannel_);
    std::vector<double> globvelv(numx1_wallinflowchannel_);
    std::vector<double> globvelw(numx1_wallinflowchannel_);

    for (size_t iplane=0; iplane<numplanes; ++iplane)
    {
      for (unsigned i=0;i<numx1_wallinflowchannel_;++i)
      {
        locforceu[i] = 0.0;
        locforcev[i] = 0.0;
        locforcew[i] = 0.0;
        locp[i] = 0.0;
        globforceu[i]= 0.0;
        globforcev[i]= 0.0;
        globforcew[i]= 0.0;
        globp[i] = 0.0;

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
                std::vector<int> dof;
                if (withscatra_)
                  dof = stddofset_->Dof(node);
                else
                  dof = discret_->Dof(node);

                // extract local values from the global vector
                std::vector<double> myforce(dof.size());
                DRT::UTILS::ExtractMyValues(*force_, myforce, dof);
                // extract pressure
                std::vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locforceu[ipos]+=myforce[0];
                locforcev[ipos]+=myforce[1];
                locforcew[ipos]+=myforce[2];
                locp     [ipos]+=myvel[3];
              }
              else if ( node->X()[dim] > velplanes[iplane]-NODETOL and node->X()[dim]<velplanes[iplane]+NODETOL )
              {
                std::vector<int> dof = discret_->Dof(node);

                // extract local values from the global vector
                std::vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locvelu[ipos]+=myvel[0];
                locvelv[ipos]+=myvel[1];
                locvelw[ipos]+=myvel[2];
              }
            }
          }
        }
      }

      discret_->Comm().SumAll(&(locforceu[0]),&(globforceu[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locforcev[0]),&(globforcev[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locforcew[0]),&(globforcew[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locp     [0]),&(globp     [0]),numx1_wallinflowchannel_);

      discret_->Comm().SumAll(&(locvelu[0]),&(globvelu[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locvelv[0]),&(globvelv[0]),numx1_wallinflowchannel_);
      discret_->Comm().SumAll(&(locvelw[0]),&(globvelw[0]),numx1_wallinflowchannel_);

      for (size_t ipos=0; ipos<numx1_wallinflowchannel_; ++ipos)
      {
        (*wallinflowchannelforceu_)(iplane,ipos) += globforceu[ipos];
        (*wallinflowchannelforcev_)(iplane,ipos) += globforcev[ipos];
        (*wallinflowchannelforcew_)(iplane,ipos) += globforcew[ipos];
        (*wallinflowchannelp_)     (iplane,ipos) += globp     [ipos];

        (*wallinflowchannelvelu_)(iplane,ipos) += globvelu[ipos];
        (*wallinflowchannelvelv_)(iplane,ipos) += globvelv[ipos];
        (*wallinflowchannelvelw_)(iplane,ipos) += globvelw[ipos];
      }
    }
  }

  // clean up pointers
  if (withscatra_)
  {
    vel_       = Teuchos::null;
    force_     = Teuchos::null;
    phi_       = Teuchos::null;
    stddofset_ = Teuchos::null;
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
    std::vector<double> x1locations(2);
    x1locations[0] = x1positions_(0); // x1=-5h
    x1locations[1] = x1positions_(1); // x1=-4h

    ExtractSetOfProfiles(x1locations, *x2_inflow_, *vertinflowu_,  *vertinflowv_,  *vertinfloww_, *vertinflowp_,
                                                   *vertinflowuu_, *vertinflowvv_, *vertinflowww_, *vertinflowpp_,
                                                   *vertinflowuv_, *vertinflowuw_, *vertinflowvw_, *vertinflowg_);
  }

  //----------------------------------------
  // select vertical profiles in mixing zone
  //----------------------------------------
  {
    std::vector<double> x1locations(3);
    x1locations[0] = x1positions_(2); // x1=-2h
    x1locations[1] = x1positions_(3); // x1=-1h
    x1locations[2] = x1positions_(4); // x1= 0h

    ExtractSetOfProfiles(x1locations, *x2_0h_, *vertmixingu_,  *vertmixingv_,  *vertmixingw_, *vertmixingp_,
                                               *vertmixinguu_, *vertmixingvv_, *vertmixingww_, *vertmixingpp_,
                                               *vertmixinguv_, *vertmixinguw_, *vertmixingvw_, *vertmixingg_);
  }

  //-----------------------------------------------------
  // select vertical profiles at 1h in recirculation zone
  //-----------------------------------------------------
  {
    std::vector<double> x1locations(1);
    x1locations[0] = x1positions_(5); // x1= 1h

    ExtractSetOfProfiles(x1locations, *x2_1h_, *vert1hu_,  *vert1hv_,  *vert1hw_, *vert1hp_,
                                               *vert1huu_, *vert1hvv_, *vert1hww_, *vert1hpp_,
                                               *vert1huv_, *vert1huw_, *vert1hvw_, *vert1hg_);
  }

  //-----------------------------------------------------
  // select vertical profiles at 2h in recirculation zone
  //-----------------------------------------------------
  {
    std::vector<double> x1locations(1);
    x1locations[0] = x1positions_(6); // x1= 2h

    ExtractSetOfProfiles(x1locations, *x2_2h_, *vert2hu_,  *vert2hv_,  *vert2hw_, *vert2hp_,
                                               *vert2huu_, *vert2hvv_, *vert2hww_, *vert2hpp_,
                                               *vert2huv_, *vert2huw_, *vert2hvw_, *vert2hg_);
  }

  //------------------------------------
  // select vertical profiles in chamber
  //------------------------------------
  {
    std::vector<double> x1locations(14);
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
                                               *vertchamberuv_, *vertchamberuw_, *vertchambervw_, *vertchamberg_);
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 | extract a set of profiles at specified locations within a homogeneous grid zone    henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::ExtractSetOfProfiles(
    const std::vector<double>&  x1locations,
    const std::vector<double>&  x2locations,
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
    LINALG::SerialDenseMatrix& profilesvw,
    LINALG::SerialDenseMatrix& profilesg
)
{
  // store contributions of each processor (local) to profile
  std::vector<std::vector<double> > locu(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > locv(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > locw(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > locp(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > locg(x1locations.size(),std::vector<double>(x2locations.size()));

  std::vector<std::vector<double> > globu(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > globv(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > globw(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > globp(x1locations.size(),std::vector<double>(x2locations.size()));
  std::vector<std::vector<double> > globg(x1locations.size(),std::vector<double>(x2locations.size()));

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
        locg[irow][icol] = 0.0;

        globu[irow][icol] = 0.0;
        globv[irow][icol] = 0.0;
        globw[irow][icol] = 0.0;
        globp[irow][icol] = 0.0;
        globg[irow][icol] = 0.0;
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

              if (withscatra_)
              {
                std::vector<int> dof = stddofset_->Dof(node);
                // extract local values from the global vector
                std::vector<double> myvel(dof.size());
                DRT::UTILS::ExtractMyValues(*vel_, myvel, dof);

                locu[ix1pos][ix2pos]=myvel[0];
                locv[ix1pos][ix2pos]=myvel[1];
                locw[ix1pos][ix2pos]=myvel[2];
                locp[ix1pos][ix2pos]=myvel[3];

                const int nodegid = node->Id();
                const int lid = phi_->Map().LID(nodegid);

                locg[ix1pos][ix2pos] = (*phi_)[lid];
              }
              else
              {
                std::vector<int> dof = discret_->Dof(node);
                // extract local values from the global vector
                std::vector<double> myvel(dof.size());
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
    }
    discret_->Comm().SumAll( &(locu[ix1pos][0]), &(globu[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locv[ix1pos][0]), &(globv[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locw[ix1pos][0]), &(globw[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locp[ix1pos][0]), &(globp[ix1pos][0]),x2locations.size() );
    discret_->Comm().SumAll( &(locg[ix1pos][0]), &(globg[ix1pos][0]),x2locations.size() );

    for(size_t ix2pos=0; ix2pos<x2locations.size(); ++ix2pos)
    {
      profilesu(ix2pos,ix1pos) += globu[ix1pos][ix2pos];
      profilesv(ix2pos,ix1pos) += globv[ix1pos][ix2pos];
      profilesw(ix2pos,ix1pos) += globw[ix1pos][ix2pos];
      profilesp(ix2pos,ix1pos) += globp[ix1pos][ix2pos];

      profilesuu(ix2pos,ix1pos) += globu[ix1pos][ix2pos]*globu[ix1pos][ix2pos];
      profilesvv(ix2pos,ix1pos) += globv[ix1pos][ix2pos]*globv[ix1pos][ix2pos];
      profilesww(ix2pos,ix1pos) += globw[ix1pos][ix2pos]*globw[ix1pos][ix2pos];
      profilespp(ix2pos,ix1pos) += globp[ix1pos][ix2pos]*globp[ix1pos][ix2pos];

      profilesuv(ix2pos,ix1pos) += globu[ix1pos][ix2pos]*globv[ix1pos][ix2pos];
      profilesuw(ix2pos,ix1pos) += globu[ix1pos][ix2pos]*globw[ix1pos][ix2pos];
      profilesvw(ix2pos,ix1pos) += globv[ix1pos][ix2pos]*globw[ix1pos][ix2pos];

      profilesg(ix2pos,ix1pos) += globg[ix1pos][ix2pos];
    }

    //std::copy(globu[ix1pos].begin(), globu[ix1pos].end(), profilesu(0,ix1pos));
    //std::copy(globv[ix1pos].begin(), globv[ix1pos].end(), profilesv(0,ix1pos));
    //std::copy(globw[ix1pos].begin(), globw[ix1pos].end(), profilesw(0,ix1pos));
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
  (*wallp_).Scale(avfac);

  (*wallvelu_).Scale(avfac);
  (*wallvelv_).Scale(avfac);
  (*wallvelw_).Scale(avfac);

  (*wallinflowchannelforceu_).Scale(avfac);
  (*wallinflowchannelforcev_).Scale(avfac);
  (*wallinflowchannelforcew_).Scale(avfac);
  (*wallinflowchannelp_).Scale(avfac);

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
    const double area = 1.0;//(x3max_-x3min_) * (x1max_-0.0);

    // nonzero forces (tractions) only expected in the streamwise (x1) direction

    // compute ltau (used to compute y+)
    LINALG::SerialDenseMatrix ltau(2,numx1_wallchamber_,true);

    for(size_t iplane=0; iplane<2; ++iplane)
    {
      for (size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        //std::cout << abs((*wallforceu_)(iplane,ix1pos)) << std::endl;
        //std::cout << abs((*wallforcev_)(iplane,ix1pos)) << std::endl;
        //std::cout << abs((*wallforcew_)(iplane,ix1pos)) << std::endl;

        //ltau(iplane,ix1pos) = visc_ / (sqrt( abs((*wallforceu_)(iplane,ix1pos)) / (area(0,ix1pos)*dens_) ));
        ltau(iplane,ix1pos) = visc_ / (sqrt( abs((*wallforceu_)(iplane,ix1pos)) / (area*dens_) ));
      }
    }

    //---------------------------------------
    // write horizontal ltau profiles chamber
    //---------------------------------------
    {
      std::string s(statistics_outfilename_);
      // define file name suffix
      s.append(".oracles.horiz.chamber.flow_statistics");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::app));
      (*log) << "\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at bottom wall \n";
      (*log) << "#   x               force u        force v        force w        p\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallforceu_)(0,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallforcev_)(0,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallforcew_)(0,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallp_)     (0,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at top wall \n";
      (*log) << "#   x               force u        force v        force w        p\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallforceu_)(1,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallforcev_)(1,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallforcew_)(1,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallp_)     (1,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from bottom wall \n";
      (*log) << "#   x               u              v              w\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallvelu_)(0,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallvelv_)(0,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallvelw_)(0,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from top wall \n";
      (*log) << "#   x               u              v              w\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallchamber_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallchamber_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallvelu_)(1,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallvelv_)(1,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallvelw_)(1,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();
    }
    //-----------------------------------------------
    // write horizontal ltau profiles inflow channels
    //-----------------------------------------------
    {
      std::string s(statistics_outfilename_);
      // define file name suffix
      s.append(".oracles.horiz.inflow.flow_statistics");

      // output to log-file
      Teuchos::RCP<std::ofstream> log;
      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::app));
      (*log) << "\n";
      (*log) << "# Statistics record " << countrecord_;
      (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at bottom wall \n";
      (*log) << "#   x               force u        force v        force w        p\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelforceu_)(0,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelforcev_)(0,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelforcew_)(0,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelp_)     (0,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# (u_tau)^2 = tau_W/rho at top wall \n";
      (*log) << "#   x               force u        force v        force w        p\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelforceu_)(1,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelforcev_)(1,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelforcew_)(1,ix1pos)/area/dens_;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelp_)     (1,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from bottom wall \n";
      (*log) << "#   x               u              v              w\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelvelu_)(0,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelvelv_)(0,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelvelw_)(0,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();

      (*log) << "\n";
      (*log) << "# u first node away from top wall \n";
      (*log) << "#   x               u              v              w\n";
      for(size_t ix1pos=0; ix1pos<numx1_wallinflowchannel_; ++ix1pos)
      {
        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*x1_wallinflowchannel_)[ix1pos];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelvelu_)(1,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelvelv_)(1,ix1pos);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*wallinflowchannelvelw_)(1,ix1pos);
        (*log) << &std::endl;
      }
      log->flush();
    }
    //-------------------------------
    // write vertical inflow profiles
    //-------------------------------
    {
      WriteStatisticsFile(step, "oracles.-5h", 0, *x2_inflow_,
          (*vertinflowu_) [0], (*vertinflowv_) [0], (*vertinfloww_) [0], (*vertinflowp_) [0],
          (*vertinflowuu_)[0], (*vertinflowvv_)[0], (*vertinflowww_)[0], (*vertinflowpp_)[0],
          (*vertinflowuv_)[0], (*vertinflowuw_)[0], (*vertinflowvw_)[0], (*vertinflowg_) [0]);

      WriteStatisticsFile(step, "oracles.ramp", 1, *x2_inflow_,
          (*vertinflowu_) [1], (*vertinflowv_) [1], (*vertinfloww_) [1], (*vertinflowp_) [1],
          (*vertinflowuu_)[1], (*vertinflowvv_)[1], (*vertinflowww_)[1], (*vertinflowpp_)[1],
          (*vertinflowuv_)[1], (*vertinflowuw_)[1], (*vertinflowvw_)[1], (*vertinflowg_) [1]);
    }
    //-------------------------------
    // write vertical mixing profiles
    //-------------------------------
    {
      WriteStatisticsFile(step, "oracles.00h", 4, *x2_0h_,
          (*vertmixingu_) [2], (*vertmixingv_) [2], (*vertmixingw_) [2], (*vertmixingp_) [2],
          (*vertmixinguu_)[2], (*vertmixingvv_)[2], (*vertmixingww_)[2], (*vertmixingpp_)[2],
          (*vertmixinguv_)[2], (*vertmixinguw_)[2], (*vertmixingvw_)[2], (*vertmixingg_) [2]);
    }
    //-----------------------------
    // write vertical profile at 1h
    //-----------------------------
    {
      WriteStatisticsFile(step, "oracles.01h", 5, *x2_1h_,
          (*vert1hu_) [0], (*vert1hv_) [0], (*vert1hw_) [0], (*vert1hp_) [0],
          (*vert1huu_)[0], (*vert1hvv_)[0], (*vert1hww_)[0], (*vert1hpp_)[0],
          (*vert1huv_)[0], (*vert1huw_)[0], (*vert1hvw_)[0], (*vert1hg_) [0]);
    }
    //-----------------------------
    // write vertical profile at 2h
    //-----------------------------
    {
      WriteStatisticsFile(step, "oracles.02h", 6, *x2_2h_,
          (*vert2hu_) [0], (*vert2hv_) [0], (*vert2hw_) [0], (*vert2hp_) [0],
          (*vert2huu_)[0], (*vert2hvv_)[0], (*vert2hww_)[0], (*vert2hpp_)[0],
          (*vert2huv_)[0], (*vert2huw_)[0], (*vert2hvw_)[0], (*vert2hg_) [0]);
    }
    //---------------------------------
    // write vertical chamber profiles
    //---------------------------------
    {
      WriteStatisticsFile(step, "oracles.03h", 7, *x2_3h_,
          (*vertchamberu_) [0], (*vertchamberv_) [0], (*vertchamberw_) [0], (*vertchamberp_) [0],
          (*vertchamberuu_)[0], (*vertchambervv_)[0], (*vertchamberww_)[0], (*vertchamberpp_)[0],
          (*vertchamberuv_)[0], (*vertchamberuw_)[0], (*vertchambervw_)[0], (*vertchamberg_) [0]);

      WriteStatisticsFile(step, "oracles.04h", 8, *x2_3h_,
          (*vertchamberu_) [1], (*vertchamberv_) [1], (*vertchamberw_) [1], (*vertchamberp_) [1],
          (*vertchamberuu_)[1], (*vertchambervv_)[1], (*vertchamberww_)[1], (*vertchamberpp_)[1],
          (*vertchamberuv_)[1], (*vertchamberuw_)[1], (*vertchambervw_)[1], (*vertchamberg_) [1]);

      WriteStatisticsFile(step, "oracles.05h", 9, *x2_3h_,
          (*vertchamberu_) [2], (*vertchamberv_) [2], (*vertchamberw_) [2], (*vertchamberp_) [2],
          (*vertchamberuu_)[2], (*vertchambervv_)[2], (*vertchamberww_)[2], (*vertchamberpp_)[2],
          (*vertchamberuv_)[2], (*vertchamberuw_)[2], (*vertchambervw_)[2], (*vertchamberg_) [2]);

      WriteStatisticsFile(step, "oracles.06h", 10, *x2_3h_,
          (*vertchamberu_) [3], (*vertchamberv_) [3], (*vertchamberw_) [3], (*vertchamberp_) [3],
          (*vertchamberuu_)[3], (*vertchambervv_)[3], (*vertchamberww_)[3], (*vertchamberpp_)[3],
          (*vertchamberuv_)[3], (*vertchamberuw_)[3], (*vertchambervw_)[3], (*vertchamberg_) [3]);

      WriteStatisticsFile(step, "oracles.07h", 11, *x2_3h_,
          (*vertchamberu_) [4], (*vertchamberv_) [4], (*vertchamberw_) [4], (*vertchamberp_) [4],
          (*vertchamberuu_)[4], (*vertchambervv_)[4], (*vertchamberww_)[4], (*vertchamberpp_)[4],
          (*vertchamberuv_)[4], (*vertchamberuw_)[4], (*vertchambervw_)[4], (*vertchamberg_) [4]);

      WriteStatisticsFile(step, "oracles.08h", 12, *x2_3h_,
          (*vertchamberu_) [5], (*vertchamberv_) [5], (*vertchamberw_) [5], (*vertchamberp_) [5],
          (*vertchamberuu_)[5], (*vertchambervv_)[5], (*vertchamberww_)[5], (*vertchamberpp_)[5],
          (*vertchamberuv_)[5], (*vertchamberuw_)[5], (*vertchambervw_)[5], (*vertchamberg_) [5]);

      WriteStatisticsFile(step, "oracles.09h", 13, *x2_3h_,
          (*vertchamberu_) [6], (*vertchamberv_) [6], (*vertchamberw_) [6], (*vertchamberp_) [6],
          (*vertchamberuu_)[6], (*vertchambervv_)[6], (*vertchamberww_)[6], (*vertchamberpp_)[6],
          (*vertchamberuv_)[6], (*vertchamberuw_)[6], (*vertchambervw_)[6], (*vertchamberg_) [6]);

      WriteStatisticsFile(step, "oracles.10h", 14, *x2_3h_,
          (*vertchamberu_) [7], (*vertchamberv_) [7], (*vertchamberw_) [7], (*vertchamberp_) [7],
          (*vertchamberuu_)[7], (*vertchambervv_)[7], (*vertchamberww_)[7], (*vertchamberpp_)[7],
          (*vertchamberuv_)[7], (*vertchamberuw_)[7], (*vertchambervw_)[7], (*vertchamberg_) [7]);

      WriteStatisticsFile(step, "oracles.11h", 15, *x2_3h_,
          (*vertchamberu_) [8], (*vertchamberv_) [8], (*vertchamberw_) [8], (*vertchamberp_) [8],
          (*vertchamberuu_)[8], (*vertchambervv_)[8], (*vertchamberww_)[8], (*vertchamberpp_)[8],
          (*vertchamberuv_)[8], (*vertchamberuw_)[8], (*vertchambervw_)[8], (*vertchamberg_) [8]);

      WriteStatisticsFile(step, "oracles.12h", 16, *x2_3h_,
          (*vertchamberu_) [9], (*vertchamberv_) [9], (*vertchamberw_) [9], (*vertchamberp_) [9],
          (*vertchamberuu_)[9], (*vertchambervv_)[9], (*vertchamberww_)[9], (*vertchamberpp_)[9],
          (*vertchamberuv_)[9], (*vertchamberuw_)[9], (*vertchambervw_)[9], (*vertchamberg_) [9]);

      WriteStatisticsFile(step, "oracles.13h", 17, *x2_3h_,
          (*vertchamberu_) [10], (*vertchamberv_) [10], (*vertchamberw_) [10], (*vertchamberp_) [10],
          (*vertchamberuu_)[10], (*vertchambervv_)[10], (*vertchamberww_)[10], (*vertchamberpp_)[10],
          (*vertchamberuv_)[10], (*vertchamberuw_)[10], (*vertchambervw_)[10], (*vertchamberg_)[10]);

      WriteStatisticsFile(step, "oracles.14h", 18, *x2_3h_,
          (*vertchamberu_) [11], (*vertchamberv_) [11], (*vertchamberw_) [11], (*vertchamberp_) [11],
          (*vertchamberuu_)[11], (*vertchambervv_)[11], (*vertchamberww_)[11], (*vertchamberpp_)[11],
          (*vertchamberuv_)[11], (*vertchamberuw_)[11], (*vertchambervw_)[11], (*vertchamberg_) [11]);

      WriteStatisticsFile(step, "oracles.15h", 19, *x2_3h_,
          (*vertchamberu_) [12], (*vertchamberv_) [12], (*vertchamberw_) [12], (*vertchamberp_) [12],
          (*vertchamberuu_)[12], (*vertchambervv_)[12], (*vertchamberww_)[12], (*vertchamberpp_)[12],
          (*vertchamberuv_)[12], (*vertchamberuw_)[12], (*vertchambervw_)[12], (*vertchamberg_) [12]);

      WriteStatisticsFile(step, "oracles.16h", 20, *x2_3h_,
          (*vertchamberu_) [13], (*vertchamberv_) [13], (*vertchamberw_) [13], (*vertchamberp_) [13],
          (*vertchamberuu_)[13], (*vertchambervv_)[13], (*vertchamberww_)[13], (*vertchamberpp_)[13],
          (*vertchamberuv_)[13], (*vertchamberuw_)[13], (*vertchambervw_)[13], (*vertchamberg_) [13]);
    }
  }

  // increase counter of written records
  countrecord_++;

  return;
}

/*------------------------------------------------------------------------------------------------*
 | reset data structures for a new record                                             henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::WriteStatisticsFile(
      const int                  step,
      const std::string&              suffix,
      const int                  x1pos,
      const std::vector<double>&  x2locations,
      double* profileu,
      double* profilev,
      double* profilew,
      double* profilep,
      double* profileuu,
      double* profilevv,
      double* profileww,
      double* profilepp,
      double* profileuv,
      double* profileuw,
      double* profilevw,
      double* profileg
)
{
  std::string s(statistics_outfilename_);

  std::ostringstream filename;
  filename << "." << suffix << ".flow_statistics";
  s.append(filename.str());

  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::app));

  (*log) << "\n";
  (*log) << "# Statistics record " << countrecord_;
  (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

  (*log) << "# profile at x-position " << std::setw(5) << std::setprecision(2) << x1positions_(x1pos)/h_ << " h \n";
  (*log) << "\n";

  (*log) << "#     y             y+";
  (*log) << "          umean         vmean         wmean         pmean";
  (*log) << "        mean u^2      mean v^2      mean w^2";
  (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2      mean G\n";
  (*log) << std::scientific;

  for(int ix2pos=x2locations.size()-1; ix2pos>=0; --ix2pos)
  {
    // y and y+
    (*log) <<  " "  << std::setw(11) << std::setprecision(4) << x2locations[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << x2locations[ix2pos]*utau_/visc_;

    // time mean values
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profileu[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profilev[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profilew[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profilep[ix2pos];

    (*log) << "   " << std::setw(11) << std::setprecision(4) << profileuu[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profilevv[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profileww[ix2pos];

    (*log) << "   " << std::setw(11) << std::setprecision(4) << profileuv[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profileuw[ix2pos];
    (*log) << "   " << std::setw(11) << std::setprecision(4) << profilevw[ix2pos];

    (*log) << "   " << std::setw(11) << std::setprecision(4) << profilepp[ix2pos];

    (*log) << "   " << std::setw(11) << std::setprecision(4) << profileg[ix2pos];

    (*log) << "\n";
  }
  log->flush();

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
  (*wallp_).Zero();

  (*wallvelu_).Zero();
  (*wallvelv_).Zero();
  (*wallvelw_).Zero();

  (*wallinflowchannelforceu_).Zero();
  (*wallinflowchannelforcev_).Zero();
  (*wallinflowchannelforcew_).Zero();
  (*wallinflowchannelp_).Zero();

  (*wallinflowchannelvelu_).Zero();
  (*wallinflowchannelvelv_).Zero();
  (*wallinflowchannelvelw_).Zero();

  //! first order momentum
  (*vertinflowu_).Zero();
  (*vertinflowv_).Zero();
  (*vertinfloww_).Zero();
  (*vertinflowp_).Zero();
  (*vertinflowg_).Zero();

  (*vertmixingu_).Zero();
  (*vertmixingv_).Zero();
  (*vertmixingw_).Zero();
  (*vertmixingp_).Zero();
  (*vertmixingg_).Zero();

  (*vert1hu_).Zero();
  (*vert1hv_).Zero();
  (*vert1hw_).Zero();
  (*vert1hp_).Zero();
  (*vert1hg_).Zero();

  (*vert2hu_).Zero();
  (*vert2hv_).Zero();
  (*vert2hw_).Zero();
  (*vert2hp_).Zero();
  (*vert2hg_).Zero();

  (*vertchamberu_).Zero();
  (*vertchamberv_).Zero();
  (*vertchamberw_).Zero();
  (*vertchamberp_).Zero();
  (*vertchamberg_).Zero();

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
 | export location from this proc to all procs (round robin)                          henke 06/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsORACLES::ExportLocation(std::set<double,LineSortCriterion>& locations)
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

  for (int np=0;np<numprocs;++np)
  {
    DRT::PackBuffer data;

    for (std::set<double,LineSortCriterion>::iterator line=locations.begin();
        line!=locations.end();
        ++line)
    {
      DRT::ParObject::AddtoPack(data,*line);
    }
    data.StartPacking();
    for (std::set<double,LineSortCriterion>::iterator line=locations.begin();
        line!=locations.end();
        ++line)
    {
      DRT::ParObject::AddtoPack(data,*line);
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

    // unpack received block into set of all planes.
    {
      std::vector<double> coordsvec;

      coordsvec.clear();

      std::vector<char>::size_type index = 0;
      while (index < rblock.size())
      {
        double onecoord;
        DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
        locations.insert(onecoord);
      }
    }
  }

}

