/*----------------------------------------------------------------------*/
/*! \file

\brief Write (time and space) averaged values to file for
turbulent flow past a square cylinder

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/
#include <fstream>

#include "turbulence_statistics_sqc.H"

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

    o Create sets for lines in x1- and x2-direction

  o Allocate distributed vector for squares

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsSqc::TurbulenceStatisticsSqc(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::ParameterList& params, const std::string& statistics_outfilename)
    : discret_(actdis), params_(params), statistics_outfilename_(statistics_outfilename)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3)
  {
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");
  }

  // use input parameter HOMDIR to specify sampling
  homdir_ = params_.sublist("TURBULENCE MODEL").get<std::string>("HOMDIR", "not_specified");

  // output to screen
  if (discret_->Comm().MyPID() == 0)
  {
    std::cout
        << "This is the turbulence statistics manager for the flow past a square-section cylinder:"
        << std::endl;

    if (homdir_ == "not_specified")
    {
      std::cout << "Slip-boundary conditions are assumed and" << std::endl;
      std::cout << "sampling is done in time only." << std::endl;
    }
    else if (homdir_ == "z")
    {
      std::cout << "Periodic-boundary conditions are assumed and" << std::endl;
      std::cout << "sampling is done in homogeneous (z-)direction and in time." << std::endl;
    }
    else
      dserror("unknown sampling procedure for square-section cylinder!");
  }

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  squaredvelnp_ = LINALG::CreateVector(*dofrowmap, true);

  toggleu_ = LINALG::CreateVector(*dofrowmap, true);
  togglev_ = LINALG::CreateVector(*dofrowmap, true);
  togglew_ = LINALG::CreateVector(*dofrowmap, true);
  togglep_ = LINALG::CreateVector(*dofrowmap, true);

  // bounds for extension of flow domain in x3-direction
  x3min_ = +10e+19;
  x3max_ = -10e+19;

  //----------------------------------------------------------------------
  // create sets of coordinates for required evaluation lines
  //----------------------------------------------------------------------
  x1ccoordinates_ = Teuchos::rcp(new std::vector<double>);
  x2ccoordinates_ = Teuchos::rcp(new std::vector<double>);
  x2wcoordinates_ = Teuchos::rcp(new std::vector<double>);
  clrcoordinates_ = Teuchos::rcp(new std::vector<double>);
  ctbcoordinates_ = Teuchos::rcp(new std::vector<double>);

  // the criterion allows differences in coordinates by 1e-9
  std::set<double, LineSortCriterion> x1cavcoords;
  std::set<double, LineSortCriterion> x2cavcoords;
  std::set<double, LineSortCriterion> x2wavcoords;
  std::set<double, LineSortCriterion> clravcoords;
  std::set<double, LineSortCriterion> ctbavcoords;

  // only necessary for averaging of Samgorinsky constant
  if (params_.sublist("TURBULENCE MODEL").get<std::string>("PHYSICAL_MODEL", "no_model") ==
      "Dynamic_Smagorinsky")
  {
    if (homdir_ == "z")
    {
      x1coordinates_ = Teuchos::rcp(new std::vector<double>);
      x2coordinates_ = Teuchos::rcp(new std::vector<double>);
    }
  }
  std::set<double, LineSortCriterion> x1avcoords;
  std::set<double, LineSortCriterion> x2avcoords;

  // loop nodes and build sets of lines accessible on this proc
  for (int i = 0; i < discret_->NumMyRowNodes(); ++i)
  {
    DRT::Node* node = discret_->lRowNode(i);

    if (node->X()[1] < (7.0 + 2e-9) && node->X()[1] > (7.0 - 2e-9))
      x1cavcoords.insert(node->X()[0]);
    if (node->X()[0] < (5.0 + 2e-9) && node->X()[0] > (5.0 - 2e-9))
      x2cavcoords.insert(node->X()[1]);
    if (node->X()[0] < (7.5 + 2e-9) && node->X()[0] > (7.5 - 2e-9))
      x2wavcoords.insert(node->X()[1]);
    if ((node->X()[0] < (4.5 + 2e-9) && node->X()[0] > (4.5 - 2e-9)) &&
        (node->X()[1] < (7.5 + 2e-9) && node->X()[1] > (6.5 - 2e-9)))
      clravcoords.insert(node->X()[1]);
    if ((node->X()[0] < (5.5 + 2e-9) && node->X()[0] > (4.5 - 2e-9)) &&
        (node->X()[1] < (7.5 + 2e-9) && node->X()[1] > (7.5 - 2e-9)))
      ctbavcoords.insert(node->X()[0]);

    // only necessary for averaging of Samgorinsky constant
    if (params_.sublist("TURBULENCE MODEL").get<std::string>("PHYSICAL_MODEL", "no_model") ==
        "Dynamic_Smagorinsky")
    {
      if (homdir_ == "z")
      {
        if (node->X()[1] < (0.0 + 2e-9) && node->X()[1] > (0.0 - 2e-9))
          x1avcoords.insert(node->X()[0]);
        if (node->X()[0] < (0.0 + 2e-9) && node->X()[0] > (0.0 - 2e-9))
          x2avcoords.insert(node->X()[1]);
      }
    }


    if (x3min_ > node->X()[2])
    {
      x3min_ = node->X()[2];
    }
    if (x3max_ < node->X()[2])
    {
      x3max_ = node->X()[2];
    }
  }

  // communicate x3mins and x3maxs
  double min;
  discret_->Comm().MinAll(&x3min_, &min, 1);
  x3min_ = min;

  double max;
  discret_->Comm().MaxAll(&x3max_, &max, 1);
  x3max_ = max;

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates on respective lines to
  // all procs
  //--------------------------------------------------------------------
  {
#ifdef PARALLEL
    int myrank = discret_->Comm().MyPID();
#endif
    int numprocs = discret_->Comm().NumProc();

    std::vector<char> sblock;
    std::vector<char> rblock;

#ifdef PARALLEL
    // create an exporter for point to point comunication
    DRT::Exporter exporter(discret_->Comm());
#endif

    // first, communicate coordinates in x1-centerline-direction
    for (int np = 0; np < numprocs; ++np)
    {
      // export set to sendbuffer
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator x1cline = x1cavcoords.begin();
           x1cline != x1cavcoords.end(); ++x1cline)
      {
        DRT::ParObject::AddtoPack(data, *x1cline);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator x1cline = x1cavcoords.begin();
           x1cline != x1cavcoords.end(); ++x1cline)
      {
        DRT::ParObject::AddtoPack(data, *x1cline);
      }
      swap(sblock, data());

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
          x1cavcoords.insert(onecoord);
        }
      }
    }

    // second, communicate coordinates in x2-centerline-direction
    for (int np = 0; np < numprocs; ++np)
    {
      // export set to sendbuffer
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator x2cline = x2cavcoords.begin();
           x2cline != x2cavcoords.end(); ++x2cline)
      {
        DRT::ParObject::AddtoPack(data, *x2cline);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator x2cline = x2cavcoords.begin();
           x2cline != x2cavcoords.end(); ++x2cline)
      {
        DRT::ParObject::AddtoPack(data, *x2cline);
      }
      swap(sblock, data());

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
          x2cavcoords.insert(onecoord);
        }
      }
    }

    // third, communicate coordinates in x2-wakeline-direction
    for (int np = 0; np < numprocs; ++np)
    {
      // export set to sendbuffer
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator x2wline = x2wavcoords.begin();
           x2wline != x2wavcoords.end(); ++x2wline)
      {
        DRT::ParObject::AddtoPack(data, *x2wline);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator x2wline = x2wavcoords.begin();
           x2wline != x2wavcoords.end(); ++x2wline)
      {
        DRT::ParObject::AddtoPack(data, *x2wline);
      }
      swap(sblock, data());

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
          x2wavcoords.insert(onecoord);
        }
      }
    }

    // fourth, communicate coordinates on left/right cylinder boundary
    for (int np = 0; np < numprocs; ++np)
    {
      // export set to sendbuffer
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator clrline = clravcoords.begin();
           clrline != clravcoords.end(); ++clrline)
      {
        DRT::ParObject::AddtoPack(data, *clrline);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator clrline = clravcoords.begin();
           clrline != clravcoords.end(); ++clrline)
      {
        DRT::ParObject::AddtoPack(data, *clrline);
      }
      swap(sblock, data());
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
          clravcoords.insert(onecoord);
        }
      }
    }

    // fifth, communicate coordinates on top/bottom cylinder boundary
    for (int np = 0; np < numprocs; ++np)
    {
      // export set to sendbuffer
      DRT::PackBuffer data;

      for (std::set<double, LineSortCriterion>::iterator ctbline = ctbavcoords.begin();
           ctbline != ctbavcoords.end(); ++ctbline)
      {
        DRT::ParObject::AddtoPack(data, *ctbline);
      }
      data.StartPacking();
      for (std::set<double, LineSortCriterion>::iterator ctbline = ctbavcoords.begin();
           ctbline != ctbavcoords.end(); ++ctbline)
      {
        DRT::ParObject::AddtoPack(data, *ctbline);
      }
      swap(sblock, data());

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
          ctbavcoords.insert(onecoord);
        }
      }
    }

    // only necessary for averaging of Samgorinsky constant
    if (params_.sublist("TURBULENCE MODEL").get<std::string>("PHYSICAL_MODEL", "no_model") ==
        "Dynamic_Smagorinsky")
    {
      if (homdir_ == "z")
      {
        // first, communicate coordinates in x1-direction
        for (int np = 0; np < numprocs; ++np)
        {
          // export set to sendbuffer
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
          swap(sblock, data());

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

        // second, communicate coordinates in x2-centerline-direction
        for (int np = 0; np < numprocs; ++np)
        {
          // export set to sendbuffer
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
          swap(sblock, data());

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
    }
  }

  //----------------------------------------------------------------------
  // push coordinates in vectors
  //----------------------------------------------------------------------
  {
    x1ccoordinates_ = Teuchos::rcp(new std::vector<double>);
    x2ccoordinates_ = Teuchos::rcp(new std::vector<double>);
    x2wcoordinates_ = Teuchos::rcp(new std::vector<double>);
    clrcoordinates_ = Teuchos::rcp(new std::vector<double>);
    ctbcoordinates_ = Teuchos::rcp(new std::vector<double>);

    for (std::set<double, LineSortCriterion>::iterator coord1c = x1cavcoords.begin();
         coord1c != x1cavcoords.end(); ++coord1c)
    {
      x1ccoordinates_->push_back(*coord1c);
    }

    for (std::set<double, LineSortCriterion>::iterator coord2c = x2cavcoords.begin();
         coord2c != x2cavcoords.end(); ++coord2c)
    {
      x2ccoordinates_->push_back(*coord2c);
    }

    for (std::set<double, LineSortCriterion>::iterator coord2w = x2wavcoords.begin();
         coord2w != x2wavcoords.end(); ++coord2w)
    {
      x2wcoordinates_->push_back(*coord2w);
    }

    for (std::set<double, LineSortCriterion>::iterator coordlr = clravcoords.begin();
         coordlr != clravcoords.end(); ++coordlr)
    {
      clrcoordinates_->push_back(*coordlr);
    }

    for (std::set<double, LineSortCriterion>::iterator coordtb = ctbavcoords.begin();
         coordtb != ctbavcoords.end(); ++coordtb)
    {
      ctbcoordinates_->push_back(*coordtb);
    }

    // only necessary for averaging of Samgorinsky constant
    if (params_.sublist("TURBULENCE MODEL").get<std::string>("PHYSICAL_MODEL", "no_model") ==
        "Dynamic_Smagorinsky")
    {
      if (homdir_ == "z")
      {
        x1coordinates_ = Teuchos::rcp(new std::vector<double>);
        x2coordinates_ = Teuchos::rcp(new std::vector<double>);

        for (std::set<double, LineSortCriterion>::iterator coord1 = x1avcoords.begin();
             coord1 != x1avcoords.end(); ++coord1)
        {
          x1coordinates_->push_back(*coord1);
        }

        for (std::set<double, LineSortCriterion>::iterator coord2 = x2avcoords.begin();
             coord2 != x2avcoords.end(); ++coord2)
        {
          x2coordinates_->push_back(*coord2);
        }
      }
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
  x1csumu_ = Teuchos::rcp(new std::vector<double>);
  x1csumu_->resize(size1c, 0.0);
  x2csumu_ = Teuchos::rcp(new std::vector<double>);
  x2csumu_->resize(size2c, 0.0);
  x2w1sumu_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumu_->resize(size2w, 0.0);
  x2w2sumu_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumu_->resize(size2w, 0.0);
  cyllsumu_ = Teuchos::rcp(new std::vector<double>);
  cyllsumu_->resize(sizelr, 0.0);
  cyltsumu_ = Teuchos::rcp(new std::vector<double>);
  cyltsumu_->resize(sizetb, 0.0);
  cylrsumu_ = Teuchos::rcp(new std::vector<double>);
  cylrsumu_->resize(sizelr, 0.0);
  cylbsumu_ = Teuchos::rcp(new std::vector<double>);
  cylbsumu_->resize(sizetb, 0.0);

  x1csumv_ = Teuchos::rcp(new std::vector<double>);
  x1csumv_->resize(size1c, 0.0);
  x2csumv_ = Teuchos::rcp(new std::vector<double>);
  x2csumv_->resize(size2c, 0.0);
  x2w1sumv_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumv_->resize(size2w, 0.0);
  x2w2sumv_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumv_->resize(size2w, 0.0);
  cyllsumv_ = Teuchos::rcp(new std::vector<double>);
  cyllsumv_->resize(sizelr, 0.0);
  cyltsumv_ = Teuchos::rcp(new std::vector<double>);
  cyltsumv_->resize(sizetb, 0.0);
  cylrsumv_ = Teuchos::rcp(new std::vector<double>);
  cylrsumv_->resize(sizelr, 0.0);
  cylbsumv_ = Teuchos::rcp(new std::vector<double>);
  cylbsumv_->resize(sizetb, 0.0);

  x1csumw_ = Teuchos::rcp(new std::vector<double>);
  x1csumw_->resize(size1c, 0.0);
  x2csumw_ = Teuchos::rcp(new std::vector<double>);
  x2csumw_->resize(size2c, 0.0);
  x2w1sumw_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumw_->resize(size2w, 0.0);
  x2w2sumw_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumw_->resize(size2w, 0.0);
  cyllsumw_ = Teuchos::rcp(new std::vector<double>);
  cyllsumw_->resize(sizelr, 0.0);
  cyltsumw_ = Teuchos::rcp(new std::vector<double>);
  cyltsumw_->resize(sizetb, 0.0);
  cylrsumw_ = Teuchos::rcp(new std::vector<double>);
  cylrsumw_->resize(sizelr, 0.0);
  cylbsumw_ = Teuchos::rcp(new std::vector<double>);
  cylbsumw_->resize(sizetb, 0.0);

  x1csump_ = Teuchos::rcp(new std::vector<double>);
  x1csump_->resize(size1c, 0.0);
  x2csump_ = Teuchos::rcp(new std::vector<double>);
  x2csump_->resize(size2c, 0.0);
  x2w1sump_ = Teuchos::rcp(new std::vector<double>);
  x2w1sump_->resize(size2w, 0.0);
  x2w2sump_ = Teuchos::rcp(new std::vector<double>);
  x2w2sump_->resize(size2w, 0.0);
  cyllsump_ = Teuchos::rcp(new std::vector<double>);
  cyllsump_->resize(sizelr, 0.0);
  cyltsump_ = Teuchos::rcp(new std::vector<double>);
  cyltsump_->resize(sizetb, 0.0);
  cylrsump_ = Teuchos::rcp(new std::vector<double>);
  cylrsump_->resize(sizelr, 0.0);
  cylbsump_ = Teuchos::rcp(new std::vector<double>);
  cylbsump_->resize(sizetb, 0.0);

  // second-order moments
  x1csumsqu_ = Teuchos::rcp(new std::vector<double>);
  x1csumsqu_->resize(size1c, 0.0);
  x2csumsqu_ = Teuchos::rcp(new std::vector<double>);
  x2csumsqu_->resize(size2c, 0.0);
  x2w1sumsqu_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumsqu_->resize(size2w, 0.0);
  x2w2sumsqu_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumsqu_->resize(size2w, 0.0);
  cyllsumsqu_ = Teuchos::rcp(new std::vector<double>);
  cyllsumsqu_->resize(sizelr, 0.0);
  cyltsumsqu_ = Teuchos::rcp(new std::vector<double>);
  cyltsumsqu_->resize(sizetb, 0.0);
  cylrsumsqu_ = Teuchos::rcp(new std::vector<double>);
  cylrsumsqu_->resize(sizelr, 0.0);
  cylbsumsqu_ = Teuchos::rcp(new std::vector<double>);
  cylbsumsqu_->resize(sizetb, 0.0);

  x1csumsqv_ = Teuchos::rcp(new std::vector<double>);
  x1csumsqv_->resize(size1c, 0.0);
  x2csumsqv_ = Teuchos::rcp(new std::vector<double>);
  x2csumsqv_->resize(size2c, 0.0);
  x2w1sumsqv_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumsqv_->resize(size2w, 0.0);
  x2w2sumsqv_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumsqv_->resize(size2w, 0.0);
  cyllsumsqv_ = Teuchos::rcp(new std::vector<double>);
  cyllsumsqv_->resize(sizelr, 0.0);
  cyltsumsqv_ = Teuchos::rcp(new std::vector<double>);
  cyltsumsqv_->resize(sizetb, 0.0);
  cylrsumsqv_ = Teuchos::rcp(new std::vector<double>);
  cylrsumsqv_->resize(sizelr, 0.0);
  cylbsumsqv_ = Teuchos::rcp(new std::vector<double>);
  cylbsumsqv_->resize(sizetb, 0.0);

  x1csumsqw_ = Teuchos::rcp(new std::vector<double>);
  x1csumsqw_->resize(size1c, 0.0);
  x2csumsqw_ = Teuchos::rcp(new std::vector<double>);
  x2csumsqw_->resize(size2c, 0.0);
  x2w1sumsqw_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumsqw_->resize(size2w, 0.0);
  x2w2sumsqw_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumsqw_->resize(size2w, 0.0);
  cyllsumsqw_ = Teuchos::rcp(new std::vector<double>);
  cyllsumsqw_->resize(sizelr, 0.0);
  cyltsumsqw_ = Teuchos::rcp(new std::vector<double>);
  cyltsumsqw_->resize(sizetb, 0.0);
  cylrsumsqw_ = Teuchos::rcp(new std::vector<double>);
  cylrsumsqw_->resize(sizelr, 0.0);
  cylbsumsqw_ = Teuchos::rcp(new std::vector<double>);
  cylbsumsqw_->resize(sizetb, 0.0);

  x1csumuv_ = Teuchos::rcp(new std::vector<double>);
  x1csumuv_->resize(size1c, 0.0);
  x2csumuv_ = Teuchos::rcp(new std::vector<double>);
  x2csumuv_->resize(size2c, 0.0);
  x2w1sumuv_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumuv_->resize(size2w, 0.0);
  x2w2sumuv_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumuv_->resize(size2w, 0.0);
  cyllsumuv_ = Teuchos::rcp(new std::vector<double>);
  cyllsumuv_->resize(sizelr, 0.0);
  cyltsumuv_ = Teuchos::rcp(new std::vector<double>);
  cyltsumuv_->resize(sizetb, 0.0);
  cylrsumuv_ = Teuchos::rcp(new std::vector<double>);
  cylrsumuv_->resize(sizelr, 0.0);
  cylbsumuv_ = Teuchos::rcp(new std::vector<double>);
  cylbsumuv_->resize(sizetb, 0.0);

  x1csumuw_ = Teuchos::rcp(new std::vector<double>);
  x1csumuw_->resize(size1c, 0.0);
  x2csumuw_ = Teuchos::rcp(new std::vector<double>);
  x2csumuw_->resize(size2c, 0.0);
  x2w1sumuw_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumuw_->resize(size2w, 0.0);
  x2w2sumuw_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumuw_->resize(size2w, 0.0);
  cyllsumuw_ = Teuchos::rcp(new std::vector<double>);
  cyllsumuw_->resize(sizelr, 0.0);
  cyltsumuw_ = Teuchos::rcp(new std::vector<double>);
  cyltsumuw_->resize(sizetb, 0.0);
  cylrsumuw_ = Teuchos::rcp(new std::vector<double>);
  cylrsumuw_->resize(sizelr, 0.0);
  cylbsumuw_ = Teuchos::rcp(new std::vector<double>);
  cylbsumuw_->resize(sizetb, 0.0);

  x1csumvw_ = Teuchos::rcp(new std::vector<double>);
  x1csumvw_->resize(size1c, 0.0);
  x2csumvw_ = Teuchos::rcp(new std::vector<double>);
  x2csumvw_->resize(size2c, 0.0);
  x2w1sumvw_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumvw_->resize(size2w, 0.0);
  x2w2sumvw_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumvw_->resize(size2w, 0.0);
  cyllsumvw_ = Teuchos::rcp(new std::vector<double>);
  cyllsumvw_->resize(sizelr, 0.0);
  cyltsumvw_ = Teuchos::rcp(new std::vector<double>);
  cyltsumvw_->resize(sizetb, 0.0);
  cylrsumvw_ = Teuchos::rcp(new std::vector<double>);
  cylrsumvw_->resize(sizelr, 0.0);
  cylbsumvw_ = Teuchos::rcp(new std::vector<double>);
  cylbsumvw_->resize(sizetb, 0.0);

  x1csumsqp_ = Teuchos::rcp(new std::vector<double>);
  x1csumsqp_->resize(size1c, 0.0);
  x2csumsqp_ = Teuchos::rcp(new std::vector<double>);
  x2csumsqp_->resize(size2c, 0.0);
  x2w1sumsqp_ = Teuchos::rcp(new std::vector<double>);
  x2w1sumsqp_->resize(size2w, 0.0);
  x2w2sumsqp_ = Teuchos::rcp(new std::vector<double>);
  x2w2sumsqp_->resize(size2w, 0.0);
  cyllsumsqp_ = Teuchos::rcp(new std::vector<double>);
  cyllsumsqp_->resize(sizelr, 0.0);
  cyltsumsqp_ = Teuchos::rcp(new std::vector<double>);
  cyltsumsqp_->resize(sizetb, 0.0);
  cylrsumsqp_ = Teuchos::rcp(new std::vector<double>);
  cylrsumsqp_->resize(sizelr, 0.0);
  cylbsumsqp_ = Teuchos::rcp(new std::vector<double>);
  cylbsumsqp_->resize(sizetb, 0.0);

  //----------------------------------------------------------------------
  // define homogeneous direction to compute averages of Smagorinsky constant

  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
  // check if we want to compute averages of Smagorinsky constant
  if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky")
  {
    if (homdir_ == "z")
    {
      // store them in parameterlist for access on the element
      modelparams->set<Teuchos::RCP<std::vector<double>>>("dir1coords_", x1coordinates_);
      modelparams->set<Teuchos::RCP<std::vector<double>>>("dir2coords_", x2coordinates_);
    }
  }

  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file

  Teuchos::RCP<std::ofstream> log;

  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow past a square-section cylinder "
              "(first- and second-order moments)\n\n";

    log->flush();
  }

  // clear statistics
  this->ClearStatistics();

  return;
}  // TurbulenceStatisticsSqc::TurbulenceStatisticsSqc

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsSqc::~TurbulenceStatisticsSqc()
{
  return;
}  // TurbulenceStatisticsSqc::~TurbulenceStatisticsSqc()


//----------------------------------------------------------------------
// sampling of lift/drag values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsSqc::DoLiftDragTimeSample(double dragforce, double liftforce)
{
  double cdrag = 2.0 * dragforce / (x3max_ - x3min_);
  double clift = 2.0 * liftforce / (x3max_ - x3min_);

  drag_ += cdrag;
  lift_ += clift;
  // TODO: passen rms-Werte?
  dragsq_ += cdrag * cdrag;
  liftsq_ += clift * clift;

  return;
}  // TurbulenceStatisticsSqc::DoLiftDragTimeSample


//----------------------------------------------------------------------
// sampling of velocity/pressure values
//----------------------------------------------------------------------
void FLD::TurbulenceStatisticsSqc::DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp)
{
  // compute squared values of velocity
  squaredvelnp_->Multiply(1.0, *velnp, *velnp, 0.0);

  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  int x1cnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x1-direction and calculate pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x1cline = x1ccoordinates_->begin();
       x1cline != x1ccoordinates_->end(); ++x1cline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if (((node->X()[0] < (*x1cline + 2e-9) && node->X()[0] > (*x1cline - 2e-9)) and
              (node->X()[1] < (7.0 + 2e-9) && node->X()[1] > (7.0 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*x1csumu_)[x1cnodnum] += u / countnodesonallprocs;
      (*x1csumv_)[x1cnodnum] += v / countnodesonallprocs;
      (*x1csumw_)[x1cnodnum] += w / countnodesonallprocs;
      (*x1csump_)[x1cnodnum] += p / countnodesonallprocs;

      (*x1csumsqu_)[x1cnodnum] += uu / countnodesonallprocs;
      (*x1csumsqv_)[x1cnodnum] += vv / countnodesonallprocs;
      (*x1csumsqw_)[x1cnodnum] += ww / countnodesonallprocs;
      (*x1csumsqp_)[x1cnodnum] += pp / countnodesonallprocs;

      (*x1csumuv_)[x1cnodnum] += uv / countnodesonallprocs;
      (*x1csumuw_)[x1cnodnum] += uw / countnodesonallprocs;
      (*x1csumvw_)[x1cnodnum] += vw / countnodesonallprocs;
    }
    x1cnodnum++;
  }

  int x2cnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on centerline in x2-direction (with respect to cylinder
  // center) and calculate pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x2cline = x2ccoordinates_->begin();
       x2cline != x2ccoordinates_->end(); ++x2cline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      if (((node->X()[1] < (*x2cline + 2e-9) && node->X()[1] > (*x2cline - 2e-9)) and
              (node->X()[0] < (5.0 + 2e-9) && node->X()[0] > (5.0 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*x2csumu_)[x2cnodnum] += u / countnodesonallprocs;
      (*x2csumv_)[x2cnodnum] += v / countnodesonallprocs;
      (*x2csumw_)[x2cnodnum] += w / countnodesonallprocs;
      (*x2csump_)[x2cnodnum] += p / countnodesonallprocs;

      (*x2csumsqu_)[x2cnodnum] += uu / countnodesonallprocs;
      (*x2csumsqv_)[x2cnodnum] += vv / countnodesonallprocs;
      (*x2csumsqw_)[x2cnodnum] += ww / countnodesonallprocs;
      (*x2csumsqp_)[x2cnodnum] += pp / countnodesonallprocs;

      (*x2csumuv_)[x2cnodnum] += uv / countnodesonallprocs;
      (*x2csumuw_)[x2cnodnum] += uw / countnodesonallprocs;
      (*x2csumvw_)[x2cnodnum] += vw / countnodesonallprocs;
    }
    x2cnodnum++;
  }

  int x2wnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on first wakeline in x2-direction (x1=7.5) and calculate
  // pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x2wline = x2wcoordinates_->begin();
       x2wline != x2wcoordinates_->end(); ++x2wline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
      if (((node->X()[1] < (*x2wline + 2e-9) && node->X()[1] > (*x2wline - 2e-9)) and
              (node->X()[0] < (7.5 + 2e-9) && node->X()[0] > (7.5 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*x2w1sumu_)[x2wnodnum] += u / countnodesonallprocs;
      (*x2w1sumv_)[x2wnodnum] += v / countnodesonallprocs;
      (*x2w1sumw_)[x2wnodnum] += w / countnodesonallprocs;
      (*x2w1sump_)[x2wnodnum] += p / countnodesonallprocs;

      (*x2w1sumsqu_)[x2wnodnum] += uu / countnodesonallprocs;
      (*x2w1sumsqv_)[x2wnodnum] += vv / countnodesonallprocs;
      (*x2w1sumsqw_)[x2wnodnum] += ww / countnodesonallprocs;
      (*x2w1sumsqp_)[x2wnodnum] += pp / countnodesonallprocs;

      (*x2w1sumuv_)[x2wnodnum] += uv / countnodesonallprocs;
      (*x2w1sumuw_)[x2wnodnum] += uw / countnodesonallprocs;
      (*x2w1sumvw_)[x2wnodnum] += vw / countnodesonallprocs;
    }
    x2wnodnum++;
  }

  x2wnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on second wakeline in x2-direction (x1=11.5) and calculate
  // pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator x2wline = x2wcoordinates_->begin();
       x2wline != x2wcoordinates_->end(); ++x2wline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x2-direction
      if (((node->X()[1] < (*x2wline + 2e-9) && node->X()[1] > (*x2wline - 2e-9)) &&
              (node->X()[0] < (11.5 + 2e-9) && node->X()[0] > (11.5 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*x2w2sumu_)[x2wnodnum] += u / countnodesonallprocs;
      (*x2w2sumv_)[x2wnodnum] += v / countnodesonallprocs;
      (*x2w2sumw_)[x2wnodnum] += w / countnodesonallprocs;
      (*x2w2sump_)[x2wnodnum] += p / countnodesonallprocs;

      (*x2w2sumsqu_)[x2wnodnum] += uu / countnodesonallprocs;
      (*x2w2sumsqv_)[x2wnodnum] += vv / countnodesonallprocs;
      (*x2w2sumsqw_)[x2wnodnum] += ww / countnodesonallprocs;
      (*x2w2sumsqp_)[x2wnodnum] += pp / countnodesonallprocs;

      (*x2w2sumuv_)[x2wnodnum] += uv / countnodesonallprocs;
      (*x2w2sumuw_)[x2wnodnum] += uw / countnodesonallprocs;
      (*x2w2sumvw_)[x2wnodnum] += vw / countnodesonallprocs;
    }
    x2wnodnum++;
  }

  int clrnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on left cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator clrline = clrcoordinates_->begin();
       clrline != clrcoordinates_->end(); ++clrline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if (((node->X()[1] < (*clrline + 2e-9) && node->X()[1] > (*clrline - 2e-9)) &&
              (node->X()[0] < (4.5 + 2e-9) && node->X()[0] > (4.5 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*cyllsumu_)[clrnodnum] += u / countnodesonallprocs;
      (*cyllsumv_)[clrnodnum] += v / countnodesonallprocs;
      (*cyllsumw_)[clrnodnum] += w / countnodesonallprocs;
      (*cyllsump_)[clrnodnum] += p / countnodesonallprocs;

      (*cyllsumsqu_)[clrnodnum] += uu / countnodesonallprocs;
      (*cyllsumsqv_)[clrnodnum] += vv / countnodesonallprocs;
      (*cyllsumsqw_)[clrnodnum] += ww / countnodesonallprocs;
      (*cyllsumsqp_)[clrnodnum] += pp / countnodesonallprocs;

      (*cyllsumuv_)[clrnodnum] += uv / countnodesonallprocs;
      (*cyllsumuw_)[clrnodnum] += uw / countnodesonallprocs;
      (*cyllsumvw_)[clrnodnum] += vw / countnodesonallprocs;
    }
    clrnodnum++;
  }

  int ctbnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on top cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator ctbline = ctbcoordinates_->begin();
       ctbline != ctbcoordinates_->end(); ++ctbline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if (((node->X()[0] < (*ctbline + 2e-9) && node->X()[0] > (*ctbline - 2e-9)) &&
              (node->X()[1] < (7.5 + 2e-9) && node->X()[1] > (7.5 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*cyltsumu_)[ctbnodnum] += u / countnodesonallprocs;
      (*cyltsumv_)[ctbnodnum] += v / countnodesonallprocs;
      (*cyltsumw_)[ctbnodnum] += w / countnodesonallprocs;
      (*cyltsump_)[ctbnodnum] += p / countnodesonallprocs;

      (*cyltsumsqu_)[ctbnodnum] += uu / countnodesonallprocs;
      (*cyltsumsqv_)[ctbnodnum] += vv / countnodesonallprocs;
      (*cyltsumsqw_)[ctbnodnum] += ww / countnodesonallprocs;
      (*cyltsumsqp_)[ctbnodnum] += pp / countnodesonallprocs;

      (*cyltsumuv_)[ctbnodnum] += uv / countnodesonallprocs;
      (*cyltsumuw_)[ctbnodnum] += uw / countnodesonallprocs;
      (*cyltsumvw_)[ctbnodnum] += vw / countnodesonallprocs;
    }
    ctbnodnum++;
  }

  clrnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on right cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator clrline = clrcoordinates_->begin();
       clrline != clrcoordinates_->end(); ++clrline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if (((node->X()[1] < (*clrline + 2e-9) && node->X()[1] > (*clrline - 2e-9)) &&
              (node->X()[0] < (5.5 + 2e-9) && node->X()[0] > (5.5 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*cylrsumu_)[clrnodnum] += u / countnodesonallprocs;
      (*cylrsumv_)[clrnodnum] += v / countnodesonallprocs;
      (*cylrsumw_)[clrnodnum] += w / countnodesonallprocs;
      (*cylrsump_)[clrnodnum] += p / countnodesonallprocs;

      (*cylrsumsqu_)[clrnodnum] += uu / countnodesonallprocs;
      (*cylrsumsqv_)[clrnodnum] += vv / countnodesonallprocs;
      (*cylrsumsqw_)[clrnodnum] += ww / countnodesonallprocs;
      (*cylrsumsqp_)[clrnodnum] += pp / countnodesonallprocs;

      (*cylrsumuv_)[clrnodnum] += uv / countnodesonallprocs;
      (*cylrsumuw_)[clrnodnum] += uw / countnodesonallprocs;
      (*cylrsumvw_)[clrnodnum] += vw / countnodesonallprocs;
    }
    clrnodnum++;
  }

  ctbnodnum = 0;
  //----------------------------------------------------------------------
  // loop nodes on bottom cylinder boundary line and calculate pointwise means
  //----------------------------------------------------------------------
  for (std::vector<double>::iterator ctbline = ctbcoordinates_->begin();
       ctbline != ctbcoordinates_->end(); ++ctbline)
  {
    // toggle vectors are one in the position of a dof for this line,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes contributing to this nodal value on line
    int countnodes = 0;

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the centerline in x1-direction
      if (((node->X()[0] < (*ctbline + 2e-9) && node->X()[0] > (*ctbline - 2e-9)) &&
              (node->X()[1] < (6.5 + 2e-9) && node->X()[1] > (6.5 - 2e-9))) and
          ((homdir_ == "z") or (homdir_ == "not_specified" and
                                   (node->X()[2] < ((x3max_ - x3min_) / 2.0) + 2e-9 &&
                                       node->X()[2] > ((x3max_ - x3min_) / 2.0) - 2e-9))))
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

    if (homdir_ == "z")
    {
      // reduce by 1 due to periodic boundary condition
      countnodesonallprocs -= 1;
    }

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
      (*cylbsumu_)[ctbnodnum] += u / countnodesonallprocs;
      (*cylbsumv_)[ctbnodnum] += v / countnodesonallprocs;
      (*cylbsumw_)[ctbnodnum] += w / countnodesonallprocs;
      (*cylbsump_)[ctbnodnum] += p / countnodesonallprocs;

      (*cylbsumsqu_)[ctbnodnum] += uu / countnodesonallprocs;
      (*cylbsumsqv_)[ctbnodnum] += vv / countnodesonallprocs;
      (*cylbsumsqw_)[ctbnodnum] += ww / countnodesonallprocs;
      (*cylbsumsqp_)[ctbnodnum] += pp / countnodesonallprocs;

      (*cylbsumuv_)[ctbnodnum] += uv / countnodesonallprocs;
      (*cylbsumuw_)[ctbnodnum] += uw / countnodesonallprocs;
      (*cylbsumvw_)[ctbnodnum] += vw / countnodesonallprocs;
    }
    ctbnodnum++;
  }

  return;
}  // TurbulenceStatisticsSqc::DoTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsSqc::DumpStatistics(int step)
{
  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow past a square-section cylinder "
              "(first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";
    (*log) << "\n";

    double liftmean = lift_ / numsamp_;
    double dragmean = drag_ / numsamp_;
    double liftrms = 0.0;
    double dragrms = 0.0;
    if ((liftsq_ / numsamp_ - liftmean * liftmean) > 0.0)
      liftrms = std::sqrt(liftsq_ / numsamp_ - liftmean * liftmean);
    if ((dragsq_ / numsamp_ - dragmean * dragmean) > 0.0)
      dragrms = std::sqrt(dragsq_ / numsamp_ - dragmean * dragmean);

    (*log) << "# lift and drag values\n";
    (*log) << "# mean lift" << std::setw(11) << std::setprecision(4) << liftmean;
    (*log) << "\n";
    (*log) << "# mean drag" << std::setw(11) << std::setprecision(4) << dragmean;
    (*log) << "\n";
    (*log) << "# rms lift " << std::setw(11) << std::setprecision(4) << liftrms;
    (*log) << "\n";
    (*log) << "# rms drag " << std::setw(11) << std::setprecision(4) << dragrms;
    (*log) << "\n";
    (*log) << "\n";
    (*log) << "\n";

    (*log) << "# centerline in x1-direction\n";
    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < x1ccoordinates_->size(); ++i)
    {
      double x1cu = (*x1csumu_)[i] / numsamp_;
      double x1cv = (*x1csumv_)[i] / numsamp_;
      double x1cw = (*x1csumw_)[i] / numsamp_;
      double x1cp = (*x1csump_)[i] / numsamp_;

      double x1curms = 0.0;
      double x1cvrms = 0.0;
      double x1cwrms = 0.0;
      double x1cprms = 0.0;

      if (((*x1csumsqu_)[i] / numsamp_ - x1cu * x1cu) > 0.0)
        x1curms = std::sqrt((*x1csumsqu_)[i] / numsamp_ - x1cu * x1cu);
      if (((*x1csumsqv_)[i] / numsamp_ - x1cv * x1cv) > 0.0)
        x1cvrms = std::sqrt((*x1csumsqv_)[i] / numsamp_ - x1cv * x1cv);
      if (((*x1csumsqw_)[i] / numsamp_ - x1cw * x1cw) > 0.0)
        x1cwrms = std::sqrt((*x1csumsqw_)[i] / numsamp_ - x1cw * x1cw);
      if (((*x1csumsqp_)[i] / numsamp_ - x1cp * x1cp) > 0.0)
        x1cprms = std::sqrt((*x1csumsqp_)[i] / numsamp_ - x1cp * x1cp);

      double x1cuv = (*x1csumuv_)[i] / numsamp_ - x1cu * x1cv;
      double x1cuw = (*x1csumuw_)[i] / numsamp_ - x1cu * x1cw;
      double x1cvw = (*x1csumvw_)[i] / numsamp_ - x1cv * x1cw;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*x1ccoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cu;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cp;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1curms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cvrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cwrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cuv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cuw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cvw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x1cprms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# centerline in x2-direction (with respect to cylinder center)\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < x2ccoordinates_->size(); ++i)
    {
      double x2cu = (*x2csumu_)[i] / numsamp_;
      double x2cv = (*x2csumv_)[i] / numsamp_;
      double x2cw = (*x2csumw_)[i] / numsamp_;
      double x2cp = (*x2csump_)[i] / numsamp_;

      double x2curms = 0.0;
      double x2cvrms = 0.0;
      double x2cwrms = 0.0;
      double x2cprms = 0.0;

      if (((*x2csumsqu_)[i] / numsamp_ - x2cu * x2cu) > 0.0)
        x2curms = std::sqrt((*x2csumsqu_)[i] / numsamp_ - x2cu * x2cu);
      if (((*x2csumsqv_)[i] / numsamp_ - x2cv * x2cv) > 0.0)
        x2cvrms = std::sqrt((*x2csumsqv_)[i] / numsamp_ - x2cv * x2cv);
      if (((*x2csumsqw_)[i] / numsamp_ - x2cw * x2cw) > 0.0)
        x2cwrms = std::sqrt((*x2csumsqw_)[i] / numsamp_ - x2cw * x2cw);
      if (((*x2csumsqp_)[i] / numsamp_ - x2cp * x2cp) > 0.0)
        x2cprms = std::sqrt((*x2csumsqp_)[i] / numsamp_ - x2cp * x2cp);

      double x2cuv = (*x2csumuv_)[i] / numsamp_ - x2cu * x2cv;
      double x2cuw = (*x2csumuw_)[i] / numsamp_ - x2cu * x2cw;
      double x2cvw = (*x2csumvw_)[i] / numsamp_ - x2cv * x2cw;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*x2ccoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cu;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cp;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2curms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cvrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cwrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cuv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cuw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cvw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2cprms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# first wakeline in x2-direction (at x1=7.5)\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < x2wcoordinates_->size(); ++i)
    {
      double x2w1u = (*x2w1sumu_)[i] / numsamp_;
      double x2w1v = (*x2w1sumv_)[i] / numsamp_;
      double x2w1w = (*x2w1sumw_)[i] / numsamp_;
      double x2w1p = (*x2w1sump_)[i] / numsamp_;

      double x2w1urms = 0.0;
      double x2w1vrms = 0.0;
      double x2w1wrms = 0.0;
      double x2w1prms = 0.0;

      if (((*x2w1sumsqu_)[i] / numsamp_ - x2w1u * x2w1u) > 0.0)
        x2w1urms = std::sqrt((*x2w1sumsqu_)[i] / numsamp_ - x2w1u * x2w1u);
      if (((*x2w1sumsqv_)[i] / numsamp_ - x2w1v * x2w1v) > 0.0)
        x2w1vrms = std::sqrt((*x2w1sumsqv_)[i] / numsamp_ - x2w1v * x2w1v);
      if (((*x2w1sumsqw_)[i] / numsamp_ - x2w1w * x2w1w) > 0.0)
        x2w1wrms = std::sqrt((*x2w1sumsqw_)[i] / numsamp_ - x2w1w * x2w1w);
      if (((*x2w1sumsqp_)[i] / numsamp_ - x2w1p * x2w1p) > 0.0)
        x2w1prms = std::sqrt((*x2w1sumsqp_)[i] / numsamp_ - x2w1p * x2w1p);

      double x2w1uv = (*x2w1sumuv_)[i] / numsamp_ - x2w1u * x2w1v;
      double x2w1uw = (*x2w1sumuw_)[i] / numsamp_ - x2w1u * x2w1w;
      double x2w1vw = (*x2w1sumvw_)[i] / numsamp_ - x2w1v * x2w1w;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*x2wcoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w1prms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# second wakeline in x2-direction (at x1=11.5)\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < x2wcoordinates_->size(); ++i)
    {
      double x2w2u = (*x2w2sumu_)[i] / numsamp_;
      double x2w2v = (*x2w2sumv_)[i] / numsamp_;
      double x2w2w = (*x2w2sumw_)[i] / numsamp_;
      double x2w2p = (*x2w2sump_)[i] / numsamp_;

      double x2w2urms = 0.0;
      double x2w2vrms = 0.0;
      double x2w2wrms = 0.0;
      double x2w2prms = 0.0;

      if (((*x2w2sumsqu_)[i] / numsamp_ - x2w2u * x2w2u) > 0.0)
        x2w2urms = std::sqrt((*x2w2sumsqu_)[i] / numsamp_ - x2w2u * x2w2u);
      if (((*x2w2sumsqv_)[i] / numsamp_ - x2w2v * x2w2v) > 0.0)
        x2w2vrms = std::sqrt((*x2w2sumsqv_)[i] / numsamp_ - x2w2v * x2w2v);
      if (((*x2w2sumsqw_)[i] / numsamp_ - x2w2w * x2w2w) > 0.0)
        x2w2wrms = std::sqrt((*x2w2sumsqw_)[i] / numsamp_ - x2w2w * x2w2w);
      if (((*x2w2sumsqp_)[i] / numsamp_ - x2w2p * x2w2p) > 0.0)
        x2w2prms = std::sqrt((*x2w2sumsqp_)[i] / numsamp_ - x2w2p * x2w2p);

      double x2w2uv = (*x2w2sumuv_)[i] / numsamp_ - x2w2u * x2w2v;
      double x2w2uw = (*x2w2sumuw_)[i] / numsamp_ - x2w2u * x2w2w;
      double x2w2vw = (*x2w2sumvw_)[i] / numsamp_ - x2w2v * x2w2w;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*x2wcoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2u;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2v;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2w;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2p;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2urms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2vrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2wrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2uv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2uw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2vw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << x2w2prms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# left cylinder boundary\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < clrcoordinates_->size(); ++i)
    {
      double cyllu = (*cyllsumu_)[i] / numsamp_;
      double cyllv = (*cyllsumv_)[i] / numsamp_;
      double cyllw = (*cyllsumw_)[i] / numsamp_;
      double cyllp = (*cyllsump_)[i] / numsamp_;

      double cyllurms = 0.0;
      double cyllvrms = 0.0;
      double cyllwrms = 0.0;
      double cyllprms = 0.0;

      if (((*cyllsumsqu_)[i] / numsamp_ - cyllu * cyllu) > 0.0)
        cyllurms = std::sqrt((*cyllsumsqu_)[i] / numsamp_ - cyllu * cyllu);
      if (((*cyllsumsqv_)[i] / numsamp_ - cyllv * cyllv) > 0.0)
        cyllvrms = std::sqrt((*cyllsumsqv_)[i] / numsamp_ - cyllv * cyllv);
      if (((*cyllsumsqw_)[i] / numsamp_ - cyllw * cyllw) > 0.0)
        cyllwrms = std::sqrt((*cyllsumsqw_)[i] / numsamp_ - cyllw * cyllw);
      if (((*cyllsumsqp_)[i] / numsamp_ - cyllp * cyllp) > 0.0)
        cyllprms = std::sqrt((*cyllsumsqp_)[i] / numsamp_ - cyllp * cyllp);

      double cylluv = (*cyllsumuv_)[i] / numsamp_ - cyllu * cyllv;
      double cylluw = (*cyllsumuw_)[i] / numsamp_ - cyllu * cyllw;
      double cyllvw = (*cyllsumvw_)[i] / numsamp_ - cyllv * cyllw;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*clrcoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllu;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllp;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllurms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllvrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllwrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylluv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylluw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllvw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyllprms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# top cylinder boundary\n";
    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < ctbcoordinates_->size(); ++i)
    {
      double cyltu = (*cyltsumu_)[i] / numsamp_;
      double cyltv = (*cyltsumv_)[i] / numsamp_;
      double cyltw = (*cyltsumw_)[i] / numsamp_;
      double cyltp = (*cyltsump_)[i] / numsamp_;

      double cylturms = 0.0;
      double cyltvrms = 0.0;
      double cyltwrms = 0.0;
      double cyltprms = 0.0;

      if (((*cyltsumsqu_)[i] / numsamp_ - cyltu * cyltu) > 0.0)
        cylturms = std::sqrt((*cyltsumsqu_)[i] / numsamp_ - cyltu * cyltu);
      if (((*cyltsumsqv_)[i] / numsamp_ - cyltv * cyltv) > 0.0)
        cyltvrms = std::sqrt((*cyltsumsqv_)[i] / numsamp_ - cyltv * cyltv);
      if (((*cyltsumsqw_)[i] / numsamp_ - cyltw * cyltw) > 0.0)
        cyltwrms = std::sqrt((*cyltsumsqw_)[i] / numsamp_ - cyltw * cyltw);
      if (((*cyltsumsqp_)[i] / numsamp_ - cyltp * cyltp) > 0.0)
        cyltprms = std::sqrt((*cyltsumsqp_)[i] / numsamp_ - cyltp * cyltp);

      double cyltuv = (*cyltsumuv_)[i] / numsamp_ - cyltu * cyltv;
      double cyltuw = (*cyltsumuw_)[i] / numsamp_ - cyltu * cyltw;
      double cyltvw = (*cyltsumvw_)[i] / numsamp_ - cyltv * cyltw;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*ctbcoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltu;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltp;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylturms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltvrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltwrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltuv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltuw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltvw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cyltprms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# right cylinder boundary\n";
    (*log) << "#     x2";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < clrcoordinates_->size(); ++i)
    {
      double cylru = (*cylrsumu_)[i] / numsamp_;
      double cylrv = (*cylrsumv_)[i] / numsamp_;
      double cylrw = (*cylrsumw_)[i] / numsamp_;
      double cylrp = (*cylrsump_)[i] / numsamp_;

      double cylrurms = 0.0;
      double cylrvrms = 0.0;
      double cylrwrms = 0.0;
      double cylrprms = 0.0;

      if (((*cylrsumsqu_)[i] / numsamp_ - cylru * cylru) > 0.0)
        cylrurms = std::sqrt((*cylrsumsqu_)[i] / numsamp_ - cylru * cylru);
      if (((*cylrsumsqv_)[i] / numsamp_ - cylrv * cylrv) > 0.0)
        cylrvrms = std::sqrt((*cylrsumsqv_)[i] / numsamp_ - cylrv * cylrv);
      if (((*cylrsumsqw_)[i] / numsamp_ - cylrw * cylrw) > 0.0)
        cylrwrms = std::sqrt((*cylrsumsqw_)[i] / numsamp_ - cylrw * cylrw);
      if (((*cylrsumsqp_)[i] / numsamp_ - cylrp * cylrp) > 0.0)
        cylrprms = std::sqrt((*cylrsumsqp_)[i] / numsamp_ - cylrp * cylrp);

      double cylruv = (*cylrsumuv_)[i] / numsamp_ - cylru * cylrv;
      double cylruw = (*cylrsumuw_)[i] / numsamp_ - cylru * cylrw;
      double cylrvw = (*cylrsumvw_)[i] / numsamp_ - cylrv * cylrw;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*clrcoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylru;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrp;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrurms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrvrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrwrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylruv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylruw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrvw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylrprms;
      (*log) << "   \n";
    }

    (*log) << "\n\n\n";
    (*log) << "# bottom cylinder boundary\n";
    (*log) << "#     x1";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "         urms          vrms          wrms";
    (*log) << "        u'v'          u'w'          v'w'          prms   \n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < ctbcoordinates_->size(); ++i)
    {
      double cylbu = (*cylbsumu_)[i] / numsamp_;
      double cylbv = (*cylbsumv_)[i] / numsamp_;
      double cylbw = (*cylbsumw_)[i] / numsamp_;
      double cylbp = (*cylbsump_)[i] / numsamp_;

      double cylburms = 0.0;
      double cylbvrms = 0.0;
      double cylbwrms = 0.0;
      double cylbprms = 0.0;

      if (((*cylbsumsqu_)[i] / numsamp_ - cylbu * cylbu) > 0.0)
        cylburms = std::sqrt((*cylbsumsqu_)[i] / numsamp_ - cylbu * cylbu);
      if (((*cylbsumsqv_)[i] / numsamp_ - cylbv * cylbv) > 0.0)
        cylbvrms = std::sqrt((*cylbsumsqv_)[i] / numsamp_ - cylbv * cylbv);
      if (((*cylbsumsqw_)[i] / numsamp_ - cylbw * cylbw) > 0.0)
        cylbwrms = std::sqrt((*cylbsumsqw_)[i] / numsamp_ - cylbw * cylbw);
      if (((*cylbsumsqp_)[i] / numsamp_ - cylbp * cylbp) > 0.0)
        cylbprms = std::sqrt((*cylbsumsqp_)[i] / numsamp_ - cylbp * cylbp);

      double cylbuv = (*cylbsumuv_)[i] / numsamp_ - cylbu * cylbv;
      double cylbuw = (*cylbsumuw_)[i] / numsamp_ - cylbu * cylbw;
      double cylbvw = (*cylbsumvw_)[i] / numsamp_ - cylbv * cylbw;

      (*log) << " " << std::setw(11) << std::setprecision(4) << (*ctbcoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbu;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbp;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylburms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbvrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbwrms;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbuv;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbuw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbvw;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << cylbprms;
      (*log) << "   \n";
    }

    log->flush();
  }

  return;

}  // TurbulenceStatisticsSqc::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsSqc::ClearStatistics()
{
  numsamp_ = 0;

  lift_ = 0.0;
  drag_ = 0.0;
  liftsq_ = 0.0;
  dragsq_ = 0.0;

  for (unsigned i = 0; i < x1ccoordinates_->size(); ++i)
  {
    (*x1csumu_)[i] = 0.0;
    (*x1csumv_)[i] = 0.0;
    (*x1csumw_)[i] = 0.0;
    (*x1csump_)[i] = 0.0;

    (*x1csumsqu_)[i] = 0.0;
    (*x1csumsqv_)[i] = 0.0;
    (*x1csumsqw_)[i] = 0.0;
    (*x1csumuv_)[i] = 0.0;
    (*x1csumuw_)[i] = 0.0;
    (*x1csumvw_)[i] = 0.0;
    (*x1csumsqp_)[i] = 0.0;
  }

  for (unsigned i = 0; i < x2ccoordinates_->size(); ++i)
  {
    (*x2csumu_)[i] = 0.0;
    (*x2csumv_)[i] = 0.0;
    (*x2csumw_)[i] = 0.0;
    (*x2csump_)[i] = 0.0;

    (*x2csumsqu_)[i] = 0.0;
    (*x2csumsqv_)[i] = 0.0;
    (*x2csumsqw_)[i] = 0.0;
    (*x2csumuv_)[i] = 0.0;
    (*x2csumuw_)[i] = 0.0;
    (*x2csumvw_)[i] = 0.0;
    (*x2csumsqp_)[i] = 0.0;
  }

  for (unsigned i = 0; i < x2wcoordinates_->size(); ++i)
  {
    (*x2w1sumu_)[i] = 0.0;
    (*x2w1sumv_)[i] = 0.0;
    (*x2w1sumw_)[i] = 0.0;
    (*x2w1sump_)[i] = 0.0;

    (*x2w1sumsqu_)[i] = 0.0;
    (*x2w1sumsqv_)[i] = 0.0;
    (*x2w1sumsqw_)[i] = 0.0;
    (*x2w1sumuv_)[i] = 0.0;
    (*x2w1sumuw_)[i] = 0.0;
    (*x2w1sumvw_)[i] = 0.0;
    (*x2w1sumsqp_)[i] = 0.0;

    (*x2w2sumu_)[i] = 0.0;
    (*x2w2sumv_)[i] = 0.0;
    (*x2w2sumw_)[i] = 0.0;
    (*x2w2sump_)[i] = 0.0;

    (*x2w2sumsqu_)[i] = 0.0;
    (*x2w2sumsqv_)[i] = 0.0;
    (*x2w2sumsqw_)[i] = 0.0;
    (*x2w2sumuv_)[i] = 0.0;
    (*x2w2sumuw_)[i] = 0.0;
    (*x2w2sumvw_)[i] = 0.0;
    (*x2w2sumsqp_)[i] = 0.0;
  }

  for (unsigned i = 0; i < clrcoordinates_->size(); ++i)
  {
    (*cyllsumu_)[i] = 0.0;
    (*cyllsumv_)[i] = 0.0;
    (*cyllsumw_)[i] = 0.0;
    (*cyllsump_)[i] = 0.0;

    (*cyllsumsqu_)[i] = 0.0;
    (*cyllsumsqv_)[i] = 0.0;
    (*cyllsumsqw_)[i] = 0.0;
    (*cyllsumuv_)[i] = 0.0;
    (*cyllsumuw_)[i] = 0.0;
    (*cyllsumvw_)[i] = 0.0;
    (*cyllsumsqp_)[i] = 0.0;

    (*cylrsumu_)[i] = 0.0;
    (*cylrsumv_)[i] = 0.0;
    (*cylrsumw_)[i] = 0.0;
    (*cylrsump_)[i] = 0.0;

    (*cylrsumsqu_)[i] = 0.0;
    (*cylrsumsqv_)[i] = 0.0;
    (*cylrsumsqw_)[i] = 0.0;
    (*cylrsumuv_)[i] = 0.0;
    (*cylrsumuw_)[i] = 0.0;
    (*cylrsumvw_)[i] = 0.0;
    (*cylrsumsqp_)[i] = 0.0;
  }

  for (unsigned i = 0; i < ctbcoordinates_->size(); ++i)
  {
    (*cyltsumu_)[i] = 0.0;
    (*cyltsumv_)[i] = 0.0;
    (*cyltsumw_)[i] = 0.0;
    (*cyltsump_)[i] = 0.0;

    (*cyltsumsqu_)[i] = 0.0;
    (*cyltsumsqv_)[i] = 0.0;
    (*cyltsumsqw_)[i] = 0.0;
    (*cyltsumuv_)[i] = 0.0;
    (*cyltsumuw_)[i] = 0.0;
    (*cyltsumvw_)[i] = 0.0;
    (*cyltsumsqp_)[i] = 0.0;

    (*cylbsumu_)[i] = 0.0;
    (*cylbsumv_)[i] = 0.0;
    (*cylbsumw_)[i] = 0.0;
    (*cylbsump_)[i] = 0.0;

    (*cylbsumsqu_)[i] = 0.0;
    (*cylbsumsqv_)[i] = 0.0;
    (*cylbsumsqw_)[i] = 0.0;
    (*cylbsumuv_)[i] = 0.0;
    (*cylbsumuw_)[i] = 0.0;
    (*cylbsumvw_)[i] = 0.0;
    (*cylbsumsqp_)[i] = 0.0;
  }

  return;
}  // TurbulenceStatisticsSqc::ClearStatistics
