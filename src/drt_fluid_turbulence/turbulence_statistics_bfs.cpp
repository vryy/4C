/*----------------------------------------------------------------------*/
/*! \file

\brief calculate pressures, mean velocity values and fluctuations for
turbulent flow over a backward-facing step

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/
#include <fstream>

#include "turbulence_statistics_bfs.H"

//#define COMBINE_SAMPLES

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

    o Create sets for lines

  o Allocate distributed vector for squares

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfs::TurbulenceStatisticsBfs(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::ParameterList& params, const std::string& statistics_outfilename,
    const std::string& geotype)
    : discret_(actdis),
      params_(params),
      geotype_(TurbulenceStatisticsBfs::none),
      inflowchannel_(
          DRT::INPUT::IntegralValue<int>(params_.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW")),
      inflowmax_(params_.sublist("TURBULENT INFLOW").get<double>("INFLOW_CHA_SIDE", 0.0)),
      statistics_outfilename_(statistics_outfilename)
{
  if (discret_->Comm().MyPID() == 0)
  {
    std::cout << "This is the turbulence statistics manager of backward-facing step problems:"
              << std::endl;
    std::cout << "It can deal with the following two geometries:" << std::endl;
    std::cout << "- geometry of DNS by Le,Moin,Kim (expansion ratio 1.2) and" << std::endl;
    std::cout << "- geometry of experiment by Kasagi,Matsunaga (expansion ratio 1.5) and"
              << std::endl;
    std::cout << "- geometry of experiment by Vogel, Eaton (expansion ratio 1.25) at Re 28,000."
              << std::endl;
    std::cout << "If additional output in front of the step is required, it has to be set manually "
                 "(look for numx1supplocations_)."
              << std::endl;
  }

  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3) dserror("Evaluation of turbulence statistics only for 3d flow problems!");

  // type of fluid flow solver: incompressible, Boussinesq approximation, varying density, loma
  const INPAR::FLUID::PhysicalType physicaltype =
      DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type");

  // geometry of bfs
  convertStringToGeoType(geotype);

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  squaredvelnp_ = LINALG::CreateVector(*dofrowmap, true);
  squaredscanp_ = LINALG::CreateVector(*dofrowmap, true);
  invscanp_ = LINALG::CreateVector(*dofrowmap, true);
  squaredinvscanp_ = LINALG::CreateVector(*dofrowmap, true);

  toggleu_ = LINALG::CreateVector(*dofrowmap, true);
  togglev_ = LINALG::CreateVector(*dofrowmap, true);
  togglew_ = LINALG::CreateVector(*dofrowmap, true);
  togglep_ = LINALG::CreateVector(*dofrowmap, true);

  // bounds for extension of flow domain in x2-direction
  x2min_ = +10e+19;
  x2max_ = -10e+19;
  // bounds for extension of flow domain in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates
  //----------------------------------------------------------------------
  x1coordinates_ = Teuchos::rcp(new std::vector<double>);
  x2coordinates_ = Teuchos::rcp(new std::vector<double>);

  // the criterion allows differences in coordinates by 1e-9
  std::set<double, LineSortCriterion> x1avcoords;
  std::set<double, LineSortCriterion> x2avcoords;

  // loop nodes and build sets of lines in x1- and x2-direction
  // accessible on this proc
  // For x1-direction: consider horizontal line at x2=0
  // and assume no change in discretization behind the step
  // For x2-direction: consider vertical line at x1=0
  // and assume no change in discretization behind the step
  for (int i = 0; i < discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);

    if (inflowchannel_ == true)
    {
      // store also x-coordinate of outflow of inflow channel
      if ((node->X()[1] < 2e-9 && node->X()[1] > -2e-9) && (node->X()[0] > (inflowmax_ - 2e-9)))
        x1avcoords.insert(node->X()[0]);
    }
    else
    {
      if (node->X()[1] < 2e-9 && node->X()[1] > -2e-9) x1avcoords.insert(node->X()[0]);
    }
    if (node->X()[0] < 2e-9 && node->X()[0] > -2e-9) x2avcoords.insert(node->X()[1]);

    // find mins and maxs
    // we do not look for x1min and x1max as they depend
    // on the inflow generation technique
    if (x2min_ > node->X()[1]) x2min_ = node->X()[1];
    if (x2max_ < node->X()[1]) x2max_ = node->X()[1];

    if (x3min_ > node->X()[2]) x3min_ = node->X()[2];
    if (x3max_ < node->X()[2]) x3max_ = node->X()[2];
  }

  // communicate x2mins and x2maxs
  double min2;
  discret_->Comm().MinAll(&x2min_, &min2, 1);
  x2min_ = min2;

  double max2;
  discret_->Comm().MaxAll(&x2max_, &max2, 1);
  x2max_ = max2;

  // communicate x3mins and x3maxs
  double min3;
  discret_->Comm().MinAll(&x3min_, &min3, 1);
  x3min_ = min3;

  double max3;
  discret_->Comm().MaxAll(&x3max_, &max3, 1);
  x3max_ = max3;

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates to all procs
  //--------------------------------------------------------------------
  {
#ifdef PARALLEL
    int myrank = discret_->Comm().MyPID();
#endif
    int numprocs = discret_->Comm().NumProc();

    std::vector<char> sblock;
    std::vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point communication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-direction
    for (int np = 0; np < numprocs; ++np)
    {
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator x1line = x1avcoords.begin();
           x1line != x1avcoords.end(); ++x1line)
      {
        DRT::ParObject::AddtoPack(data, *x1line);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator x1line = x1avcoords.begin();
           x1line != x1avcoords.end(); ++x1line)
      {
        DRT::ParObject::AddtoPack(data, *x1line);
      }
      std::swap(sblock, data());

#ifdef PARALLEL
      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank + 1) % numprocs;

      int length = sblock.size();

      exporter.ISend(frompid, topid, &(sblock[0]), sblock.size(), tag, request);

      rblock.clear();

      // receive from predecessor
      frompid = (myrank + numprocs - 1) % numprocs;
      exporter.ReceiveAny(frompid, tag, rblock, length);

      if (tag != (myrank + numprocs - 1) % numprocs)
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
      rblock = sblock;
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
          DRT::ParObject::ExtractfromPack(index, rblock, onecoord);
          x1avcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-direction
    for (int np = 0; np < numprocs; ++np)
    {
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator x2line = x2avcoords.begin();
           x2line != x2avcoords.end(); ++x2line)
      {
        DRT::ParObject::AddtoPack(data, *x2line);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator x2line = x2avcoords.begin();
           x2line != x2avcoords.end(); ++x2line)
      {
        DRT::ParObject::AddtoPack(data, *x2line);
      }
      std::swap(sblock, data());

#ifdef PARALLEL
      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid = (myrank + 1) % numprocs;

      int length = sblock.size();

      exporter.ISend(frompid, topid, &(sblock[0]), sblock.size(), tag, request);

      rblock.clear();

      // receive from predecessor
      frompid = (myrank + numprocs - 1) % numprocs;
      exporter.ReceiveAny(frompid, tag, rblock, length);

      if (tag != (myrank + numprocs - 1) % numprocs)
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
      rblock = sblock;
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
          DRT::ParObject::ExtractfromPack(index, rblock, onecoord);
          x2avcoords.insert(onecoord);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------
  {
    x1coordinates_ = Teuchos::rcp(new std::vector<double>);
    x2coordinates_ = Teuchos::rcp(new std::vector<double>);

    for (std::set<double, LineSortCriterion>::iterator coord1 = x1avcoords.begin();
         coord1 != x1avcoords.end(); ++coord1)
    {
      x1coordinates_->push_back(*coord1);
      // std::cout << *coord1 << std::endl;
    }

    for (std::set<double, LineSortCriterion>::iterator coord2 = x2avcoords.begin();
         coord2 != x2avcoords.end(); ++coord2)
    {
      x2coordinates_->push_back(*coord2);
      // std::cout << *coord2 << std::endl;
    }
  }

  //----------------------------------------------------------------------
  // number of coordinates in x1- and x2-direction
  //----------------------------------------------------------------------
  numx1coor_ = x1coordinates_->size();
  numx2coor_ = x2coordinates_->size();

  //----------------------------------------------------------------------
  // number of locations in x1-direction for statistical evaluation
  //----------------------------------------------------------------------
  numx1statlocations_ = 21;
  numx1supplocations_ = 0;

  //----------------------------------------------------------------------
  // define locations in x1-direction for statistical evaluation
  //----------------------------------------------------------------------

  // compute step height
  double h = 0.0;
  if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating)
  {
    h = (x2max_ - x2min_) / 3.0;
  }
  else if (geotype_ == TurbulenceStatisticsBfs::geometry_DNS_incomp_flow)
  {
    h = (x2max_ - x2min_) / 6.0;
  }
  else if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
  {
    h = (x2max_ - x2min_) / 5.0;
    numx1statlocations_ = 10;
  }

  if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
  {
    // locactions given by Vogel&Eaton
    double givenpos[10] = {2.2, 3.0, 3.73, 4.47, 5.2, 5.93, 6.67, 7.4, 8.13, 8.87};
    if (numx1statlocations_ != 10) dserror("wrong number of numx1statlocations_");
    // find locations
    for (int rr = 0; rr < numx1statlocations_; rr++)
    {
      double actpos = givenpos[rr];
      double dist = 30.0 * h;
      double mindist = 30.0 * h;
      int pos = 0;

      for (int ll = 0; ll < numx1coor_; ll++)
      {
        dist = abs(actpos * h - x1coordinates_->at(ll));
        if (dist < mindist)
        {
          mindist = dist;
          pos = ll;
        }
      }

      x1statlocations_(rr) = x1coordinates_->at(pos);
    }
  }
  else
  {
    // find locations x/h=0 ... x/h=20
    for (int rr = 0; rr < numx1statlocations_; rr++)
    {
      int actpos = rr;
      double dist = 30 * h;
      double mindist = 30.0 * h;
      int pos = 0;

      for (int ll = 0; ll < numx1coor_; ll++)
      {
        dist = abs(actpos * h - x1coordinates_->at(ll));
        if (dist < mindist)
        {
          mindist = dist;
          pos = ll;
        }
      }

      x1statlocations_(rr) = x1coordinates_->at(pos);
    }
  }

  // find supplementary locations x/h=-2 and -1
  // remark1: number of supplementary location depends on length of inlet section
  // remark2: as separate channel may be located in front of the step, we have to
  //         use this rather complicated way
  numx1supplocations_ = 2;
  numx1statlocations_ += numx1supplocations_;
  for (int rr = 0; rr < numx1supplocations_; rr++)
  {
    int actpos = rr - numx1supplocations_;
    double dist = 30 * h;
    double mindist = 30.0 * h;
    int pos = 0;

    for (int ll = 0; ll < numx1coor_; ll++)
    {
      dist = abs(actpos * h - x1coordinates_->at(ll));
      if (dist < mindist)
      {
        mindist = dist;
        pos = ll;
      }
    }

    x1supplocations_(rr) = x1coordinates_->at(pos);
    // std::cout << x1supplocations_(rr) << std::endl;
  }

  //----------------------------------------------------------------------
  // define locations in x2-direction for statistical evaluation
  // (lower and upper wall)
  //----------------------------------------------------------------------
  if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating ||
      geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
  {
    // num2statlocations_ also defines number of supplocations
    numx2statlocations_ = 2;

    x2statlocations_(0) = x2min_;
    x2statlocations_(1) = x2max_;
    //----------------------------------------------------------------------
    // define supplementary locations in x2-direction for statistical
    // evaluation of velocity derivative at wall
    // (first nodes off lower and upper wall, respectively)
    //----------------------------------------------------------------------
    x2supplocations_(0) = x2coordinates_->at(1);
    x2supplocations_(1) = x2coordinates_->at(x2coordinates_->size() - 2);
  }
  else if (geotype_ == TurbulenceStatisticsBfs::geometry_DNS_incomp_flow)
  {
    // num2statlocations_ also defines number of supplocations
    numx2statlocations_ = 1;

    x2statlocations_(0) = x2min_;
    x2statlocations_(1) = x2max_;  // not needed here, upper wall is slip wall
    //----------------------------------------------------------------------
    // define supplementary locations in x2-direction for statistical
    // evaluation of velocity derivative at wall
    // (first nodes off lower and upper wall, respectively)
    //----------------------------------------------------------------------
    x2supplocations_(0) = x2coordinates_->at(1);
    x2supplocations_(1) =
        x2coordinates_->at(x2coordinates_->size() - 2);  // not needed here, upper wall is slip wall
  }
  else
    dserror("Unknown geometry of backward facing step!");

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
  // x1-direction
  x1sump_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sump_->Reshape(numx2statlocations_, numx1coor_);

  x1sumu_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumu_->Reshape(numx2statlocations_, numx1coor_);

  // the following vectors are only necessary for low-Mach-number flow
  x1sumrho_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumrho_->Reshape(numx2statlocations_, numx1coor_);

  x1sumT_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumT_->Reshape(numx2statlocations_, numx1coor_);

  x1sumtauw_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x1sumtauw_->Reshape(numx2statlocations_, numx1coor_);

  // x2-direction
  // first-order moments
  x2sumu_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumu_->Reshape(numx1statlocations_, numx2coor_);

  x2sumv_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumv_->Reshape(numx1statlocations_, numx2coor_);

  x2sumw_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumw_->Reshape(numx1statlocations_, numx2coor_);

  x2sump_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sump_->Reshape(numx1statlocations_, numx2coor_);

  // second-order moments
  x2sumsqu_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqu_->Reshape(numx1statlocations_, numx2coor_);

  x2sumsqv_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqv_->Reshape(numx1statlocations_, numx2coor_);

  x2sumsqw_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqw_->Reshape(numx1statlocations_, numx2coor_);

  x2sumsqp_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqp_->Reshape(numx1statlocations_, numx2coor_);

  x2sumuv_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumuv_->Reshape(numx1statlocations_, numx2coor_);

  x2sumuw_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumuw_->Reshape(numx1statlocations_, numx2coor_);

  x2sumvw_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumvw_->Reshape(numx1statlocations_, numx2coor_);

  // the following vectors are only necessary for low-Mach-number flow
  // first-order moments
  x2sumrho_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumrho_->Reshape(numx1statlocations_, numx2coor_);

  x2sumT_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumT_->Reshape(numx1statlocations_, numx2coor_);

  // second-order moments
  x2sumsqrho_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqrho_->Reshape(numx1statlocations_, numx2coor_);

  x2sumsqT_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumsqT_->Reshape(numx1statlocations_, numx2coor_);

  x2sumrhou_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumrhou_->Reshape(numx1statlocations_, numx2coor_);

  x2sumuT_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumuT_->Reshape(numx1statlocations_, numx2coor_);

  x2sumrhov_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumrhov_->Reshape(numx1statlocations_, numx2coor_);

  x2sumvT_ = Teuchos::rcp(new Epetra_SerialDenseMatrix);
  x2sumvT_->Reshape(numx1statlocations_, numx2coor_);

  // set number of samples to zero
  numsamp_ = 0;

  //----------------------------------------------------------------------
  // define homogeneous direction to compute averages of Smagorinsky constant

  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
  // check if we want to compute averages of Smagorinsky constant
  if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky")
  {
    // store them in parameterlist for access on the element
    modelparams->set<Teuchos::RCP<std::vector<double>>>("dir1coords_", x1coordinates_);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("dir2coords_", x2coordinates_);
  }

  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file

  Teuchos::RCP<std::ofstream> log;

  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);

    if (physicaltype == INPAR::FLUID::loma)
    {
      s.append(".loma_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
      (*log) << "# Statistics for turbulent variable-density flow over a backward-facing step at "
                "low Mach number (first- and second-order moments)\n\n";
    }
    else
    {
      s.append(".flow_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
      (*log) << "# Statistics for turbulent incompressible flow over a backward-facing step "
                "(first- and second-order moments)\n\n";
    }

    log->flush();
  }

  return;
}  // TurbulenceStatisticsBfs::TurbulenceStatisticsBfs

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfs::~TurbulenceStatisticsBfs()
{
  return;
}  // TurbulenceStatisticsBfs::~TurbulenceStatisticsBfs()


//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::DoTimeSample(
    Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> stresses)
{
  // compute squared values of velocity
  squaredvelnp_->Multiply(1.0, *velnp, *velnp, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1line = x1coordinates_->begin();
       x1line != x1coordinates_->end(); ++x1line)
  {
    x1nodnum++;

    for (int x2nodnum = 0; x2nodnum < numx2statlocations_; ++x2nodnum)
    {
      // current x2-coordinate of respective wall
      double x2cwall = x2statlocations_(x2nodnum);

      // current x2-coordinate of supplementary location to respective wall
      double x2csupp = x2supplocations_(x2nodnum);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglep_->PutScalar(0.0);
      // misuse togglev for stresses in u direction
      togglev_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the wall node
        if ((node->X()[0] < (*x1line + 2e-9) and node->X()[0] > (*x1line - 2e-9)) and
            (node->X()[1] < (x2cwall + 2e-5) and node->X()[1] > (x2cwall - 2e-5)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));
          // stresses
          togglev_->ReplaceGlobalValues(1, &one, &(dof[0]));

          countnodes++;
        }
        // this is the supplementary node
        else if ((node->X()[0] < (*x1line + 2e-9) and node->X()[0] > (*x1line - 2e-9)) and
                 (node->X()[1] < (x2csupp + 2e-5) and node->X()[1] > (x2csupp - 2e-5)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
        }
      }

      int countnodesonallprocs = 0;

      discret_->Comm().SumAll(&countnodes, &countnodesonallprocs, 1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative and pressure
        //----------------------------------------------------------------------
        double u;
        velnp->Dot(*toggleu_, &u);
        double p;
        velnp->Dot(*togglep_, &p);
        // stresses
        double tauw;
        stresses->Dot(*togglev_, &tauw);

        //----------------------------------------------------------------------
        // calculate spatial means
        //----------------------------------------------------------------------
        double usm = u / countnodesonallprocs;
        double psm = p / countnodesonallprocs;
        double tauwsm = tauw / countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumu_)(x2nodnum, x1nodnum) += usm;
        (*x1sump_)(x2nodnum, x1nodnum) += psm;
        (*x1sumtauw_)(x2nodnum, x1nodnum) += tauwsm;
      }
    }
  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------
  for (int x1nodnum = 0; x1nodnum < numx1statlocations_; ++x1nodnum)
  {
    // current x1-coordinate
    // caution: if there are supplementary locations in x1-direction, we loop
    //          them first (only DNS geometry)
    double x1c = 1.0e20;
    if (x1nodnum < numx1supplocations_)
    {
      x1c = x1supplocations_(x1nodnum);
    }
    else
    {
      x1c = x1statlocations_(x1nodnum - numx1supplocations_);
    }

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction
    //----------------------------------------------------------------------
    for (std::vector<double>::iterator x2line = x2coordinates_->begin();
         x2line != x2coordinates_->end(); ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the node
        if ((node->X()[0] < (x1c + 2e-5) and node->X()[0] > (x1c - 2e-5)) and
            (node->X()[1] < (*x2line + 2e-9) and node->X()[1] > (*x2line - 2e-9)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
          togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
          togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs = 0;

      discret_->Comm().SumAll(&countnodes, &countnodesonallprocs, 1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity and pressure on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp->Dot(*toggleu_, &u);
        velnp->Dot(*togglev_, &v);
        velnp->Dot(*togglew_, &w);
        velnp->Dot(*togglep_, &p);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->Dot(*toggleu_, &uu);
        squaredvelnp_->Dot(*togglev_, &vv);
        squaredvelnp_->Dot(*togglew_, &ww);
        squaredvelnp_->Dot(*togglep_, &pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr = 1; rr < velnp->MyLength(); ++rr)
        {
          locuv += ((*velnp)[rr - 1] * (*toggleu_)[rr - 1]) * ((*velnp)[rr] * (*togglev_)[rr]);
        }
        discret_->Comm().SumAll(&locuv, &uv, 1);
        for (int rr = 2; rr < velnp->MyLength(); ++rr)
        {
          locuw += ((*velnp)[rr - 2] * (*toggleu_)[rr - 2]) * ((*velnp)[rr] * (*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locuw, &uw, 1);
        for (int rr = 2; rr < velnp->MyLength(); ++rr)
        {
          locvw += ((*velnp)[rr - 1] * (*togglev_)[rr - 1]) * ((*velnp)[rr] * (*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locvw, &vw, 1);

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum, x2nodnum) += u / countnodesonallprocs;
        (*x2sumv_)(x1nodnum, x2nodnum) += v / countnodesonallprocs;
        (*x2sumw_)(x1nodnum, x2nodnum) += w / countnodesonallprocs;
        (*x2sump_)(x1nodnum, x2nodnum) += p / countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum, x2nodnum) += uu / countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum, x2nodnum) += vv / countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum, x2nodnum) += ww / countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum, x2nodnum) += pp / countnodesonallprocs;

        (*x2sumuv_)(x1nodnum, x2nodnum) += uv / countnodesonallprocs;
        (*x2sumuw_)(x1nodnum, x2nodnum) += uw / countnodesonallprocs;
        (*x2sumvw_)(x1nodnum, x2nodnum) += vw / countnodesonallprocs;
      }
    }
  }

  return;
}  // TurbulenceStatisticsBfs::DoTimeSample


//----------------------------------------------------------------------
// sampling of velocity, pressure and temperature values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::DoLomaTimeSample(
    Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> scanp, const double eosfac)
{
  // compute squared values of velocity
  squaredvelnp_->Multiply(1.0, *velnp, *velnp, 0.0);
  squaredscanp_->Multiply(1.0, *scanp, *scanp, 0.0);
  // compute 1/T and (1/T)^2
  for (int rr = 0; rr < squaredscanp_->MyLength(); ++rr)
  {
    if ((*scanp)[rr] > 0)  // temperature in kelvin is always > 0
      (*invscanp_)[rr] = 1 / ((*scanp)[rr]);
  }
  squaredinvscanp_->Multiply(1.0, *invscanp_, *invscanp_, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1line = x1coordinates_->begin();
       x1line != x1coordinates_->end(); ++x1line)
  {
    x1nodnum++;

    for (int x2nodnum = 0; x2nodnum < numx2statlocations_; ++x2nodnum)
    {
      // current x2-coordinate of respective wall
      double x2cwall = x2statlocations_(x2nodnum);

      // current x2-coordinate of supplementary location to respective wall
      double x2csupp = x2supplocations_(x2nodnum);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the wall node
        if ((node->X()[0] < (*x1line + 2e-9) and node->X()[0] > (*x1line - 2e-9)) and
            (node->X()[1] < (x2cwall + 2e-5) and node->X()[1] > (x2cwall - 2e-5)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));

          countnodes++;
        }
        // this is the supplementary node
        else if ((node->X()[0] < (*x1line + 2e-9) and node->X()[0] > (*x1line - 2e-9)) and
                 (node->X()[1] < (x2csupp + 2e-5) and node->X()[1] > (x2csupp - 2e-5)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
        }
      }

      int countnodesonallprocs = 0;

      discret_->Comm().SumAll(&countnodes, &countnodesonallprocs, 1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative, pressure and temperature
        //----------------------------------------------------------------------
        double u;
        velnp->Dot(*toggleu_, &u);
        double p;
        velnp->Dot(*togglep_, &p);
        double T;
        scanp->Dot(*togglep_, &T);

        double rho;
        invscanp_->Dot(*togglep_, &rho);
        // compute density: rho = eosfac/T
        rho *= eosfac;

        //----------------------------------------------------------------------
        // calculate spatial means
        //----------------------------------------------------------------------
        double usm = u / countnodesonallprocs;
        double psm = p / countnodesonallprocs;
        double Tsm = T / countnodesonallprocs;
        double rhosm = rho / countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumu_)(x2nodnum, x1nodnum) += usm;
        (*x1sump_)(x2nodnum, x1nodnum) += psm;
        (*x1sumrho_)(x2nodnum, x1nodnum) += rhosm;
        (*x1sumT_)(x2nodnum, x1nodnum) += Tsm;
      }
    }
  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------
  for (int x1nodnum = 0; x1nodnum < numx1statlocations_; ++x1nodnum)
  {
    // current x1-coordinate
    // caution: if there are supplementary locations in x1-direction, we loop
    //          them first (only DNS geometry)
    double x1c = 1.0e20;
    if (x1nodnum < numx1supplocations_)
    {
      x1c = x1supplocations_(x1nodnum);
    }
    else
    {
      x1c = x1statlocations_(x1nodnum - numx1supplocations_);
    }

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction and calculate pointwise means
    //----------------------------------------------------------------------
    for (std::vector<double>::iterator x2line = x2coordinates_->begin();
         x2line != x2coordinates_->end(); ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the node
        if ((node->X()[0] < (x1c + 2e-5) and node->X()[0] > (x1c - 2e-5)) and
            (node->X()[1] < (*x2line + 2e-9) and node->X()[1] > (*x2line - 2e-9)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
          togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
          togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs = 0;

      discret_->Comm().SumAll(&countnodes, &countnodesonallprocs, 1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity, pressure, density, temperature, and
        // subgrid viscosity on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp->Dot(*toggleu_, &u);
        velnp->Dot(*togglev_, &v);
        velnp->Dot(*togglew_, &w);
        velnp->Dot(*togglep_, &p);

        double T;
        scanp->Dot(*togglep_, &T);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->Dot(*toggleu_, &uu);
        squaredvelnp_->Dot(*togglev_, &vv);
        squaredvelnp_->Dot(*togglew_, &ww);
        squaredvelnp_->Dot(*togglep_, &pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr = 1; rr < velnp->MyLength(); ++rr)
        {
          locuv += ((*velnp)[rr - 1] * (*toggleu_)[rr - 1]) * ((*velnp)[rr] * (*togglev_)[rr]);
        }
        discret_->Comm().SumAll(&locuv, &uv, 1);
        for (int rr = 2; rr < velnp->MyLength(); ++rr)
        {
          locuw += ((*velnp)[rr - 2] * (*toggleu_)[rr - 2]) * ((*velnp)[rr] * (*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locuw, &uw, 1);
        for (int rr = 2; rr < velnp->MyLength(); ++rr)
        {
          locvw += ((*velnp)[rr - 1] * (*togglev_)[rr - 1]) * ((*velnp)[rr] * (*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locvw, &vw, 1);

        double TT;
        squaredscanp_->Dot(*togglep_, &TT);

        double uT;
        double vT;
        double locuT = 0.0;
        double locvT = 0.0;
        for (int rr = 3; rr < velnp->MyLength(); ++rr)
        {
          locuT += ((*velnp)[rr - 3] * (*toggleu_)[rr - 3]) * ((*scanp)[rr] * (*togglep_)[rr]);
        }
        discret_->Comm().SumAll(&locuT, &uT, 1);
        for (int rr = 3; rr < velnp->MyLength(); ++rr)
        {
          locvT += ((*velnp)[rr - 2] * (*togglev_)[rr - 2]) * ((*scanp)[rr] * (*togglep_)[rr]);
        }
        discret_->Comm().SumAll(&locvT, &vT, 1);

        double rho;
        invscanp_->Dot(*togglep_, &rho);
        // compute density: rho = eosfac/T
        rho *= eosfac;
        double rhorho;
        squaredinvscanp_->Dot(*togglep_, &rhorho);
        rhorho *= eosfac * eosfac;

        double rhou;
        double rhov;
        double locrhou = 0.0;
        double locrhov = 0.0;
        for (int rr = 3; rr < velnp->MyLength(); ++rr)
        {
          locrhou += (eosfac * ((*invscanp_)[rr] * (*togglep_)[rr])) *
                     ((*velnp)[rr - 3] * (*toggleu_)[rr - 3]);
        }
        discret_->Comm().SumAll(&locrhou, &rhou, 1);
        for (int rr = 3; rr < velnp->MyLength(); ++rr)
        {
          locrhov += (eosfac * ((*invscanp_)[rr] * (*togglep_)[rr])) *
                     ((*velnp)[rr - 2] * (*togglev_)[rr - 2]);
        }
        discret_->Comm().SumAll(&locrhov, &rhov, 1);


        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum, x2nodnum) += u / countnodesonallprocs;
        (*x2sumv_)(x1nodnum, x2nodnum) += v / countnodesonallprocs;
        (*x2sumw_)(x1nodnum, x2nodnum) += w / countnodesonallprocs;
        (*x2sump_)(x1nodnum, x2nodnum) += p / countnodesonallprocs;

        (*x2sumT_)(x1nodnum, x2nodnum) += T / countnodesonallprocs;
        (*x2sumrho_)(x1nodnum, x2nodnum) += rho / countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum, x2nodnum) += uu / countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum, x2nodnum) += vv / countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum, x2nodnum) += ww / countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum, x2nodnum) += pp / countnodesonallprocs;

        (*x2sumsqT_)(x1nodnum, x2nodnum) += TT / countnodesonallprocs;
        (*x2sumsqrho_)(x1nodnum, x2nodnum) += rhorho / countnodesonallprocs;

        (*x2sumuv_)(x1nodnum, x2nodnum) += uv / countnodesonallprocs;
        (*x2sumuw_)(x1nodnum, x2nodnum) += uw / countnodesonallprocs;
        (*x2sumvw_)(x1nodnum, x2nodnum) += vw / countnodesonallprocs;

        (*x2sumrhou_)(x1nodnum, x2nodnum) += rhou / countnodesonallprocs;
        (*x2sumuT_)(x1nodnum, x2nodnum) += uT / countnodesonallprocs;
        (*x2sumrhov_)(x1nodnum, x2nodnum) += rhov / countnodesonallprocs;
        (*x2sumvT_)(x1nodnum, x2nodnum) += vT / countnodesonallprocs;
      }
    }
  }

  return;
}  // TurbulenceStatisticsBfc::DoLomaTimeSample


//----------------------------------------------------------------------
// sampling of velocity, pressure and scalar values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsBfs::DoScatraTimeSample(
    Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> scanp)
{
  // compute squared values of velocity
  squaredvelnp_->Multiply(1.0, *velnp, *velnp, 0.0);
  squaredscanp_->Multiply(1.0, *scanp, *scanp, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1nodnum = -1;
  //----------------------------------------------------------------------
  // values at lower and upper wall
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1line = x1coordinates_->begin();
       x1line != x1coordinates_->end(); ++x1line)
  {
    x1nodnum++;

    for (int x2nodnum = 0; x2nodnum < numx2statlocations_; ++x2nodnum)
    {
      // current x2-coordinate of respective wall
      double x2cwall = x2statlocations_(x2nodnum);

      // current x2-coordinate of supplementary location to respective wall
      double x2csupp = x2supplocations_(x2nodnum);

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the wall node
        if ((node->X()[0] < (*x1line + 2e-9) and node->X()[0] > (*x1line - 2e-9)) and
            (node->X()[1] < (x2cwall + 2e-5) and node->X()[1] > (x2cwall - 2e-5)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));

          countnodes++;
        }
        // this is the supplementary node
        else if ((node->X()[0] < (*x1line + 2e-9) and node->X()[0] > (*x1line - 2e-9)) and
                 (node->X()[1] < (x2csupp + 2e-5) and node->X()[1] > (x2csupp - 2e-5)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
        }
      }

      int countnodesonallprocs = 0;

      discret_->Comm().SumAll(&countnodes, &countnodesonallprocs, 1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative, pressure and temperature
        //----------------------------------------------------------------------
        double u;
        velnp->Dot(*toggleu_, &u);
        double p;
        velnp->Dot(*togglep_, &p);
        double T;
        scanp->Dot(*togglep_, &T);

        //----------------------------------------------------------------------
        // calculate spatial means
        //----------------------------------------------------------------------
        double usm = u / countnodesonallprocs;
        double psm = p / countnodesonallprocs;
        double Tsm = T / countnodesonallprocs;

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x1sumu_)(x2nodnum, x1nodnum) += usm;
        (*x1sump_)(x2nodnum, x1nodnum) += psm;
        (*x1sumT_)(x2nodnum, x1nodnum) += Tsm;
      }
    }
  }

  //----------------------------------------------------------------------
  // loop locations for statistical evaluation in x1-direction
  //----------------------------------------------------------------------
  for (int x1nodnum = 0; x1nodnum < numx1statlocations_; ++x1nodnum)
  {
    // current x1-coordinate
    // caution: if there are supplementary locations in x1-direction, we loop
    //          them first (only DNS geometry)
    double x1c = 1.0e20;
    if (x1nodnum < numx1supplocations_)
    {
      x1c = x1supplocations_(x1nodnum);
    }
    else
    {
      x1c = x1statlocations_(x1nodnum - numx1supplocations_);
    }

    int x2nodnum = -1;
    //----------------------------------------------------------------------
    // loop nodes in x2-direction and calculate pointwise means
    //----------------------------------------------------------------------
    for (std::vector<double>::iterator x2line = x2coordinates_->begin();
         x2line != x2coordinates_->end(); ++x2line)
    {
      x2nodnum++;

      // toggle vectors are one in the position of a dof of this node,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes = 0;

      for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this is the node
        if ((node->X()[0] < (x1c + 2e-5) and node->X()[0] > (x1c - 2e-5)) and
            (node->X()[1] < (*x2line + 2e-9) and node->X()[1] > (*x2line - 2e-9)))
        {
          std::vector<int> dof = discret_->Dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
          togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
          togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));

          countnodes++;
        }
      }

      int countnodesonallprocs = 0;

      discret_->Comm().SumAll(&countnodes, &countnodesonallprocs, 1);

      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;

      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity, pressure, density, temperature, and
        // subgrid viscosity on this line
        //----------------------------------------------------------------------
        double u;
        double v;
        double w;
        double p;
        velnp->Dot(*toggleu_, &u);
        velnp->Dot(*togglev_, &v);
        velnp->Dot(*togglew_, &w);
        velnp->Dot(*togglep_, &p);

        double T;
        scanp->Dot(*togglep_, &T);

        double uu;
        double vv;
        double ww;
        double pp;
        squaredvelnp_->Dot(*toggleu_, &uu);
        squaredvelnp_->Dot(*togglev_, &vv);
        squaredvelnp_->Dot(*togglew_, &ww);
        squaredvelnp_->Dot(*togglep_, &pp);

        double uv;
        double uw;
        double vw;
        double locuv = 0.0;
        double locuw = 0.0;
        double locvw = 0.0;
        for (int rr = 1; rr < velnp->MyLength(); ++rr)
        {
          locuv += ((*velnp)[rr - 1] * (*toggleu_)[rr - 1]) * ((*velnp)[rr] * (*togglev_)[rr]);
        }
        discret_->Comm().SumAll(&locuv, &uv, 1);
        for (int rr = 2; rr < velnp->MyLength(); ++rr)
        {
          locuw += ((*velnp)[rr - 2] * (*toggleu_)[rr - 2]) * ((*velnp)[rr] * (*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locuw, &uw, 1);
        for (int rr = 2; rr < velnp->MyLength(); ++rr)
        {
          locvw += ((*velnp)[rr - 1] * (*togglev_)[rr - 1]) * ((*velnp)[rr] * (*togglew_)[rr]);
        }
        discret_->Comm().SumAll(&locvw, &vw, 1);

        double TT;
        squaredscanp_->Dot(*togglep_, &TT);

        double uT;
        double vT;
        double locuT = 0.0;
        double locvT = 0.0;
        for (int rr = 3; rr < velnp->MyLength(); ++rr)
        {
          locuT += ((*velnp)[rr - 3] * (*toggleu_)[rr - 3]) * ((*scanp)[rr] * (*togglep_)[rr]);
        }
        discret_->Comm().SumAll(&locuT, &uT, 1);
        for (int rr = 3; rr < velnp->MyLength(); ++rr)
        {
          locvT += ((*velnp)[rr - 2] * (*togglev_)[rr - 2]) * ((*scanp)[rr] * (*togglep_)[rr]);
        }
        discret_->Comm().SumAll(&locvT, &vT, 1);

        //----------------------------------------------------------------------
        // calculate spatial means on this line
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        (*x2sumu_)(x1nodnum, x2nodnum) += u / countnodesonallprocs;
        (*x2sumv_)(x1nodnum, x2nodnum) += v / countnodesonallprocs;
        (*x2sumw_)(x1nodnum, x2nodnum) += w / countnodesonallprocs;
        (*x2sump_)(x1nodnum, x2nodnum) += p / countnodesonallprocs;

        (*x2sumT_)(x1nodnum, x2nodnum) += T / countnodesonallprocs;

        (*x2sumsqu_)(x1nodnum, x2nodnum) += uu / countnodesonallprocs;
        (*x2sumsqv_)(x1nodnum, x2nodnum) += vv / countnodesonallprocs;
        (*x2sumsqw_)(x1nodnum, x2nodnum) += ww / countnodesonallprocs;
        (*x2sumsqp_)(x1nodnum, x2nodnum) += pp / countnodesonallprocs;

        (*x2sumsqT_)(x1nodnum, x2nodnum) += TT / countnodesonallprocs;

        (*x2sumuv_)(x1nodnum, x2nodnum) += uv / countnodesonallprocs;
        (*x2sumuw_)(x1nodnum, x2nodnum) += uw / countnodesonallprocs;
        (*x2sumvw_)(x1nodnum, x2nodnum) += vw / countnodesonallprocs;

        (*x2sumuT_)(x1nodnum, x2nodnum) += uT / countnodesonallprocs;
        (*x2sumvT_)(x1nodnum, x2nodnum) += vT / countnodesonallprocs;
      }
    }
  }

  return;
}  // TurbulenceStatisticsBfc::DoScatraTimeSample


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow over a backward-facing step (first- "
              "and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n\n\n";
    (*log) << std::scientific;

    (*log) << "\n\n\n";
    (*log) << "# lower wall behind step\n";
    (*log) << "#     x1";
    (*log) << "           duxdy         pmean         tauw\n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i = 0; i < x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u = (*x1sumu_)(0, i) / numsamp_;
        double lwx1duxdy = lwx1u / dist;
        double lwx1p = (*x1sump_)(0, i) / numsamp_;
        double lwx1tauw = (*x1sumtauw_)(0, i) / numsamp_;

        (*log) << " " << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1duxdy;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1p;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << lwx1tauw;
        (*log) << "\n";
      }
    }

    if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating ||
        geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
    {
      (*log) << "\n\n\n";
      (*log) << "# upper wall\n";
      (*log) << "#     x1";
      (*log) << "           duxdy         pmean\n";

      // distance from wall to first node off wall
      dist = x2statlocations_(1) - x2supplocations_(1);

      for (unsigned i = 0; i < x1coordinates_->size(); ++i)
      {
        double uwx1u = (*x1sumu_)(1, i) / numsamp_;
        double uwx1duxdy = uwx1u / dist;
        double uwx1p = (*x1sump_)(1, i) / numsamp_;

        (*log) << " " << std::setw(11) << std::setprecision(4) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1duxdy;
        (*log) << "   " << std::setw(11) << std::setprecision(4) << uwx1p;
        (*log) << "\n";
      }
    }

    for (int i = 0; i < numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first
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
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(4) << x1
             << "\n";
      (*log) << "#     x2";
      (*log) << "           umean         vmean         wmean         pmean";
      (*log) << "         urms          vrms          wrms          prms";
      (*log) << "          u'v'          u'w'          v'w'\n";

      for (unsigned j = 0; j < x2coordinates_->size(); ++j)
      {
        double x2u = (*x2sumu_)(i, j) / numsamp_;
        double x2v = (*x2sumv_)(i, j) / numsamp_;
        double x2w = (*x2sumw_)(i, j) / numsamp_;
        double x2p = (*x2sump_)(i, j) / numsamp_;

        double x2urms = 0.0;
        double x2vrms = 0.0;
        double x2wrms = 0.0;
        double x2prms = 0.0;

        if (((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u) > 0.0)
          x2urms = std::sqrt((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u);
        if (((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v) > 0.0)
          x2vrms = std::sqrt((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v);
        if (((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w) > 0.0)
          x2wrms = std::sqrt((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w);
        if (((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p) > 0.0)
          x2prms = std::sqrt((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p);

        double x2uv = (*x2sumuv_)(i, j) / numsamp_ - x2u * x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_ - x2u * x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_ - x2v * x2w;

        (*log) << " " << std::setw(11) << std::setprecision(4) << (*x2coordinates_)[j];
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

}  // TurbulenceStatisticsBfs::DumpStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::DumpLomaStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".loma_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent variable-density flow over a backward-facing step at low "
              "Mach number (first- and second-order moments)";
    (*log) << "\n\n";
#ifdef COMBINE_SAMPLES
    (*log) << "# Statistics are prepared for combinations after restart!!!";
    (*log) << "\n\n";
#endif
    (*log) << "# Caution: The following statistics have to be used carefully:\n";
    (*log) << "#          rhoumean, uTmean, rhovmean, vTmean, rhou'T', rhov'T'\n";
    (*log) << "#          there are not any reference values for rhoumean, uTmean, rhovmean, "
              "vTmean and rhou'T'\n";
    (*log) << "#          it is not clear what is denoted by rhou'T' and rhov'T' in the paper of "
              "Avnacha and Pletcher(2002)\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n\n\n";
    (*log) << std::scientific;

    (*log) << "\n\n\n";
    (*log) << "# lower wall behind step\n";
    (*log) << "#        x1";
    (*log) << "                 duxdy               pmean             rhomean              Tmean\n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i = 0; i < x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u = (*x1sumu_)(0, i) / numsamp_;
        double lwx1duxdy = lwx1u / dist;
        double lwx1p = (*x1sump_)(0, i) / numsamp_;

        double lwx1rho = (*x1sumrho_)(0, i) / numsamp_;
        double lwx1T = (*x1sumT_)(0, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1rho;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1T;
        (*log) << "\n";
      }
    }

    if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating)
    {
      (*log) << "\n\n\n";
      (*log) << "# upper wall\n";
      (*log) << "#        x1";
      (*log)
          << "                 duxdy               pmean             rhomean              Tmean\n";

      // distance from wall to first node off wall
      dist = x2statlocations_(1) - x2supplocations_(1);

      for (unsigned i = 0; i < x1coordinates_->size(); ++i)
      {
        double uwx1u = (*x1sumu_)(1, i) / numsamp_;
        double uwx1duxdy = uwx1u / dist;
        double uwx1p = (*x1sump_)(1, i) / numsamp_;

        double uwx1rho = (*x1sumrho_)(1, i) / numsamp_;
        double uwx1T = (*x1sumT_)(1, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1rho;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1T;
        (*log) << "\n";
      }
    }
    else if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
      dserror("geometry not implemented for loma yet");

    for (int i = 0; i < numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first
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
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(10) << x1
             << "\n";
      (*log) << "#        x2";
      (*log) << "                 umean               vmean               wmean               "
                "pmean             rhomean               Tmean            rhoumean        "
                "rhouTmean            rhovmean        rhovTmean";
      (*log) << "               urms                vrms                wrms                prms   "
                "            rhorms                Trms";
#ifndef COMBINE_SAMPLES
      (*log) << "                u'v'                u'w'                v'w'             rhou'T'  "
                "           rhov'T'\n";
#else
      (*log) << "                u'v'                u'w'                v'w'             rhou'T'  "
                "           rhov'T'";
      (*log) << "                uu                vv                ww             pp         TT  "
                "           rhorho\n";
#endif

      for (unsigned j = 0; j < x2coordinates_->size(); ++j)
      {
        double x2u = (*x2sumu_)(i, j) / numsamp_;
        double x2v = (*x2sumv_)(i, j) / numsamp_;
        double x2w = (*x2sumw_)(i, j) / numsamp_;
        double x2p = (*x2sump_)(i, j) / numsamp_;

        double x2rho = (*x2sumrho_)(i, j) / numsamp_;
        double x2T = (*x2sumT_)(i, j) / numsamp_;

        double x2urms = 0.0;
        double x2vrms = 0.0;
        double x2wrms = 0.0;
        double x2prms = 0.0;

        if (((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u) > 0.0)
          x2urms = std::sqrt((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u);
        if (((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v) > 0.0)
          x2vrms = std::sqrt((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v);
        if (((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w) > 0.0)
          x2wrms = std::sqrt((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w);
        if (((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p) > 0.0)
          x2prms = std::sqrt((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p);

#ifdef COMBINE_SAMPLES
        double x2usq = (*x2sumsqu_)(i, j) / numsamp_;
        double x2vsq = (*x2sumsqv_)(i, j) / numsamp_;
        double x2wsq = (*x2sumsqw_)(i, j) / numsamp_;
        double x2psq = (*x2sumsqp_)(i, j) / numsamp_;
#endif

        // as T and rho are constant in the inflow section
        // <T(rho)^2>-<T(rho)>*<T(rho)> should be zero
        // however, due to small errors, <T(rho)^2>-<T(rho)>*<T(rho)>
        // is only approximately equal zero
        // hence, zero negative values should be excluded
        // as they produce nans
        double x2rhorms = 0.0;
        double x2Trms = 0.0;
        if (std::abs((*x2sumsqrho_)(i, j) / numsamp_ - x2rho * x2rho) > 1e-9)
          x2rhorms = std::sqrt((*x2sumsqrho_)(i, j) / numsamp_ - x2rho * x2rho);
        if (std::abs((*x2sumsqT_)(i, j) / numsamp_ - x2T * x2T) > 1e-9)
          x2Trms = std::sqrt((*x2sumsqT_)(i, j) / numsamp_ - x2T * x2T);

#ifdef COMBINE_SAMPLES
        double x2rhosq = (*x2sumsqrho_)(i, j) / numsamp_;
        double x2Tsq = (*x2sumsqT_)(i, j) / numsamp_;

        double x2uv = (*x2sumuv_)(i, j) / numsamp_;  //-x2u*x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_;  //-x2u*x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_;  //-x2v*x2w;

        double x2rhou = (*x2sumrhou_)(i, j) / numsamp_;  //-x2u*x2rho;
        double x2uT = (*x2sumuT_)(i, j) / numsamp_;      //-x2u*x2T;
        double x2rhov = (*x2sumrhov_)(i, j) / numsamp_;  //-x2v*x2rho;
        double x2vT = (*x2sumvT_)(i, j) / numsamp_;      //-x2v*x2T;
#else
        double x2uv = (*x2sumuv_)(i, j) / numsamp_ - x2u * x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_ - x2u * x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_ - x2v * x2w;

        double x2rhou = (*x2sumrhou_)(i, j) / numsamp_ - x2u * x2rho;
        double x2uT = (*x2sumuT_)(i, j) / numsamp_ - x2u * x2T;
        double x2rhov = (*x2sumrhov_)(i, j) / numsamp_ - x2v * x2rho;
        double x2vT = (*x2sumvT_)(i, j) / numsamp_ - x2v * x2T;
#endif

        double x2rhouppTpp = x2rho * (x2uT - x2u * x2T);
        double x2rhovppTpp = x2rho * (x2vT - x2v * x2T);

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x2coordinates_)[j];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2u;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2v;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2w;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rho;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2T;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhou;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhov;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2urms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2wrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2prms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhorms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2Trms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uv;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uw;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vw;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhouppTpp;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhovppTpp;
#ifdef COMBINE_SAMPLES
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2usq;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vsq;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2wsq;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2psq;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2Tsq;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2rhosq;
#endif
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}  // TurbulenceStatisticsBfs::DumpLomaStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::DumpScatraStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent flow with passive scalar over a backward-facing step "
              "(first- and second-order moments)";
    (*log) << "\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n\n\n";
    (*log) << std::scientific;

    (*log) << "\n\n\n";
    (*log) << "# lower wall behind step\n";
    (*log) << "#        x1";
    (*log) << "                 duxdy               pmean              phimean\n";

    // distance from wall to first node off wall
    double dist = x2supplocations_(0) - x2statlocations_(0);

    for (unsigned i = 0; i < x1coordinates_->size(); ++i)
    {
      if ((*x1coordinates_)[i] > -2e-9)
      {
        double lwx1u = (*x1sumu_)(0, i) / numsamp_;
        double lwx1duxdy = lwx1u / dist;
        double lwx1p = (*x1sump_)(0, i) / numsamp_;

        double lwx1T = (*x1sumT_)(0, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << lwx1T;
        (*log) << "\n";
      }
    }

    if (geotype_ == TurbulenceStatisticsBfs::geometry_LES_flow_with_heating)
    {
      (*log) << "\n\n\n";
      (*log) << "# upper wall\n";
      (*log) << "#        x1";
      (*log) << "                 duxdy               pmean              phimean\n";

      // distance from wall to first node off wall
      dist = x2statlocations_(1) - x2supplocations_(1);

      for (unsigned i = 0; i < x1coordinates_->size(); ++i)
      {
        double uwx1u = (*x1sumu_)(1, i) / numsamp_;
        double uwx1duxdy = uwx1u / dist;
        double uwx1p = (*x1sump_)(1, i) / numsamp_;

        double uwx1T = (*x1sumT_)(1, i) / numsamp_;

        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x1coordinates_)[i];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1duxdy;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << uwx1T;
        (*log) << "\n";
      }
    }
    else if (geotype_ == TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton)
      dserror("geometry not implemented for scatra yet");

    for (int i = 0; i < numx1statlocations_; ++i)
    {
      // current x1-coordinate
      // caution: if there are supplementary locations in x1-direction, we loop
      //          them first
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
      (*log) << "# line in x2-direction at x1 = " << std::setw(11) << std::setprecision(10) << x1
             << "\n";
      (*log) << "#        x2";
      (*log) << "                 umean               vmean               wmean               "
                "pmean               phimean           uphimean           vphimean";
      (*log) << "               urms                vrms                wrms                prms   "
                "             phirms";
      (*log) << "                u'v'                u'w'                v'w'\n";

      for (unsigned j = 0; j < x2coordinates_->size(); ++j)
      {
        double x2u = (*x2sumu_)(i, j) / numsamp_;
        double x2v = (*x2sumv_)(i, j) / numsamp_;
        double x2w = (*x2sumw_)(i, j) / numsamp_;
        double x2p = (*x2sump_)(i, j) / numsamp_;

        double x2T = (*x2sumT_)(i, j) / numsamp_;

        double x2urms = 0.0;
        double x2vrms = 0.0;
        double x2wrms = 0.0;
        double x2prms = 0.0;

        if (((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u) > 0.0)
          x2urms = std::sqrt((*x2sumsqu_)(i, j) / numsamp_ - x2u * x2u);
        if (((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v) > 0.0)
          x2vrms = std::sqrt((*x2sumsqv_)(i, j) / numsamp_ - x2v * x2v);
        if (((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w) > 0.0)
          x2wrms = std::sqrt((*x2sumsqw_)(i, j) / numsamp_ - x2w * x2w);
        if (((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p) > 0.0)
          x2prms = std::sqrt((*x2sumsqp_)(i, j) / numsamp_ - x2p * x2p);

        // as T is constant in the inflow section
        // <T^2>-<T>*<T> should be zero
        // however, due to small errors, <T^2>-<T>*<T>
        // is only approximately equal zero
        // hence, zero negative values should be excluded
        // as they produce nans
        double x2Trms = 0.0;
        if (std::abs((*x2sumsqT_)(i, j) / numsamp_ - x2T * x2T) > 1e-9)
          x2Trms = std::sqrt((*x2sumsqT_)(i, j) / numsamp_ - x2T * x2T);

        double x2uv = (*x2sumuv_)(i, j) / numsamp_ - x2u * x2v;
        double x2uw = (*x2sumuw_)(i, j) / numsamp_ - x2u * x2w;
        double x2vw = (*x2sumvw_)(i, j) / numsamp_ - x2v * x2w;

        double x2uT = (*x2sumuT_)(i, j) / numsamp_ - x2u * x2T;
        double x2vT = (*x2sumvT_)(i, j) / numsamp_ - x2v * x2T;


        (*log) << " " << std::setw(17) << std::setprecision(10) << (*x2coordinates_)[j];
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2u;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2v;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2w;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2p;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2T;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vT;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2urms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2wrms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2prms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2Trms;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uv;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2uw;
        (*log) << "   " << std::setw(17) << std::setprecision(10) << x2vw;
        (*log) << "\n";
      }
    }
    log->flush();
  }

  return;

}  // TurbulenceStatisticsBfs::DumpScatraStatistics


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfs::convertStringToGeoType(const std::string& geotype)
{
  dsassert(geotype != "none", "No geometry supplied");

  geotype_ = TurbulenceStatisticsBfs::none;
  if (geotype == "geometry_DNS_incomp_flow")
    geotype_ = TurbulenceStatisticsBfs::geometry_DNS_incomp_flow;
  else if (geotype == "geometry_LES_flow_with_heating")
    geotype_ = TurbulenceStatisticsBfs::geometry_LES_flow_with_heating;
  else if (geotype == "geometry_EXP_vogel_eaton")
    geotype_ = TurbulenceStatisticsBfs::geometry_EXP_vogel_eaton;
  else
    dserror("(%s) geometry for backward facing step", geotype.c_str());
  return;
}
