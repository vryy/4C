/*----------------------------------------------------------------------*/
/*! \file

 \brief detailed description in header file levelset_intersection.H

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
 *------------------------------------------------------------------------------------------------*/


#include "levelset_intersection_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_cut/cut_levelsetintersection.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include "../drt_scatra_ele/scatra_ele_parameter_std.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
SCATRA::LEVELSET::Intersection::Intersection()
    : check_lsv_(false), desired_positions_(0), volumeplus_(0.0), volumeminus_(0.0), surface_(0.0)
{
  desired_positions_.reserve(2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::Reset()
{
  volumeplus_ = 0.0;
  volumeminus_ = 0.0;
  surface_ = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::CaptureZeroLevelSet(
    const Teuchos::RCP<const Epetra_Vector>& phi,
    const Teuchos::RCP<const DRT::Discretization>& scatradis, double& volumedomainminus,
    double& volumedomainplus, double& zerosurface,
    std::map<int, GEO::BoundaryIntCells>& elementBoundaryIntCells)
{
  // reset, just to be sure
  Reset();
  volumedomainminus = 0.0;
  volumedomainplus = 0.0;
  zerosurface = 0.0;
  elementBoundaryIntCells.clear();

  // herein the actual capturing happens
  GetZeroLevelSet(*phi, *scatradis, elementBoundaryIntCells);

  // collect contributions from all procs and store in respective variables
  scatradis->Comm().SumAll(&VolumePlus(), &volumedomainplus, 1);
  scatradis->Comm().SumAll(&VolumeMinus(), &volumedomainminus, 1);
  scatradis->Comm().SumAll(&Surface(), &zerosurface, 1);

  // export also interface to all procs
  ExportInterface(elementBoundaryIntCells, scatradis->Comm());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T>
void SCATRA::LEVELSET::Intersection::GetZeroLevelSet(const Epetra_Vector& phi,
    const DRT::Discretization& scatradis, std::map<int, T>& elementBoundaryIntCells,
    bool cut_screenoutput)
{
  // export phi from row to column map
  const Teuchos::RCP<Epetra_Vector> phicol =
      Teuchos::rcp(new Epetra_Vector(*scatradis.DofColMap()));
  LINALG::Export(phi, *phicol);

  // remark: loop over row elements is sufficient
  for (int iele = 0; iele < scatradis.NumMyRowElements(); ++iele)
  {
    // get element from discretization
    const DRT::Element* ele = scatradis.lRowElement(iele);
    const DRT::Element::DiscretizationType distype = ele->Shape();

    // clear vector each loop
    BoundaryIntCellsPerEle<T>().clear();

    // ------------------------------------------------------------------------
    // Prepare cut
    // ------------------------------------------------------------------------
    GEO::CUT::LevelSetIntersection levelset(scatradis.Comm());
    LINALG::SerialDenseMatrix xyze;
    std::vector<double> phi_nodes;
    std::vector<int> nids;
    PrepareCut(ele, scatradis, *phicol, xyze, phi_nodes, nids);

    // check if this element is cut, according to its level-set values
    // -> add it to 'levelset'
    // note: cut is performed in physical space
    if (!levelset.AddElement(1, nids, xyze, ele->Shape(), &phi_nodes[0], false, check_lsv_))
      continue;

    // ------------------------------------------------------------------------
    // call GEO::Cut algorithm and process cut data
    // ------------------------------------------------------------------------
    GEO::CUT::ElementHandle* ehandle = Cut(levelset, xyze, phi_nodes, cut_screenoutput);

    // =========================================================
    // cell is in contact with the interface (cut or touched)
    // =========================================================
    if (ehandle != NULL)
    {
      GEO::CUT::plain_element_set cuteles;

      CollectCutEles(*ehandle, cuteles, distype);

      // ----------------------------------------------------------------------
      // get zero level-set contour
      // ----------------------------------------------------------------------
      GetZeroLevelSetContour(cuteles, xyze, distype);
    }
    // =========================================================
    // element is uncut
    // =========================================================
    else
    {
      double elevol = GEO::ElementVolume(distype, xyze);

      // it is sufficient to check the first node, since the element entirely
      // lies within the plus or minus domain
      if (phi_nodes[0] > 0.0)
        VolumePlus() += elevol;
      else
        VolumeMinus() += elevol;
    }

    // store interface of element
    if (BoundaryIntCellsPerEle<T>().size() > 0)
      elementBoundaryIntCells[ele->Id()] = BoundaryIntCellsPerEle<T>();
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::GetZeroLevelSetContour(
    const GEO::CUT::plain_element_set& cuteles, const LINALG::SerialDenseMatrix& xyze,
    DRT::Element::DiscretizationType distype)
{
  for (GEO::CUT::plain_element_set::const_iterator icutele = cuteles.begin();
       icutele != cuteles.end(); ++icutele)
  {
    // get pointer to cut element
    GEO::CUT::Element* cutele = *icutele;

    GEO::CUT::plain_volumecell_set volcells;
    volcells = cutele->VolumeCells();

    for (GEO::CUT::plain_volumecell_set::const_iterator ivolcell = volcells.begin();
         ivolcell != volcells.end(); ++ivolcell)
    {
      GEO::CUT::VolumeCell* volcell = *ivolcell;
      const GEO::CUT::Point::PointPosition vol_pos = volcell->Position();
      if (IsPointPosition(vol_pos))
      {
        AddToVolume(vol_pos, volcell->Volume());
        // get boundary integration cells for this volume cell
        // we consider only the cells for one position, otherwise we would have the boundary
        // cells twice
        const GEO::CUT::plain_boundarycell_set& bcells = volcell->BoundaryCells();
        for (GEO::CUT::plain_boundarycell_set::const_iterator ibcell = bcells.begin();
             ibcell != bcells.end(); ++ibcell)
        {
          GEO::CUT::BoundaryCell* bcell = *ibcell;

          AddToBoundaryIntCellsPerEle(xyze, *bcell, distype);

          Surface() += bcell->Area();
        }
      }
      else
      {
        AddToVolume(vol_pos, volcell->Volume());
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::AddToBoundaryIntCellsPerEle(
    const LINALG::SerialDenseMatrix& xyze, const GEO::CUT::BoundaryCell& bcell,
    DRT::Element::DiscretizationType distype_ele)
{
  DRT::Element::DiscretizationType distype_bc = bcell.Shape();
  CheckBoundaryCellType(distype_bc);

  const int numnodebc = DRT::UTILS::getNumberOfElementNodes(distype_bc);

  // get physical coordinates of this cell
  LINALG::SerialDenseMatrix coord = bcell.Coordinates();

  // transfer to element coordinates
  LINALG::SerialDenseMatrix localcoord(3, numnodebc, true);

  for (int ivert = 0; ivert < numnodebc; ivert++)
  {
    LINALG::Matrix<3, 1> lcoord;
    LINALG::Matrix<3, 1> pcoord;
    for (int ll = 0; ll < 3; ll++) pcoord(ll, 0) = coord(ll, ivert);

    GEO::currentToVolumeElementCoordinates(distype_ele, xyze, pcoord, lcoord);

    // write as 'physCoord'
    for (int ll = 0; ll < 3; ll++) localcoord(ll, ivert) = lcoord(ll, 0);
  }

  // store boundary element and sum area into surface
  // be careful, we only set physical coordinates
  BoundaryIntCellsPerEle<GEO::BoundaryIntCells>().push_back(
      GEO::BoundaryIntCell(distype_bc, -1, localcoord, Teuchos::null, coord, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::CheckBoundaryCellType(
    DRT::Element::DiscretizationType distype_bc) const
{
  if (distype_bc != DRT::Element::tri3 and distype_bc != DRT::Element::quad4)
  {
    dserror("unexpected type of boundary integration cell: %s",
        DRT::DistypeToString(distype_bc).c_str());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::AddToVolume(GEO::CUT::Point::PointPosition pos, double vol)
{
  switch (pos)
  {
    case GEO::CUT::Point::outside:
      VolumePlus() += vol;
      break;
    case GEO::CUT::Point::inside:
      VolumeMinus() += vol;
      break;
    default:
      /* do nothing for the undecided case */
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::CollectCutEles(GEO::CUT::ElementHandle& ehandle,
    GEO::CUT::plain_element_set& cuteles, DRT::Element::DiscretizationType distype) const
{
  ehandle.CollectElements(cuteles);

  switch (distype)
  {
    case DRT::Element::line2:
    case DRT::Element::hex8:
    {
      if (cuteles.size() != 1) dserror("one cut element expected for linear elements");
      break;
    }
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    {
      if (cuteles.size() != 8) dserror("eight cut elements expected for quadratic elements");
      break;
    }
    default:
    {
      dserror("distype unknown for level set cut algorithm");
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::PrepareCut(const DRT::Element* ele,
    const DRT::Discretization& scatradis, const Epetra_Vector& phicol,
    LINALG::SerialDenseMatrix& xyze, std::vector<double>& phi_nodes,
    std::vector<int>& node_ids) const
{
  const DRT::Element::DiscretizationType distype = ele->Shape();
  unsigned numnode = DRT::UTILS::getNumberOfElementNodes(distype);
  const unsigned probdim = DRT::Problem::Instance()->NDim();

  xyze.Shape(3, numnode);
  switch (distype)
  {
    case DRT::Element::hex8:
      GEO::fillInitialPositionArray<DRT::Element::hex8, 3>(ele, xyze);
      break;
    case DRT::Element::hex20:
      GEO::fillInitialPositionArray<DRT::Element::hex20, 3>(ele, xyze);
      break;
    case DRT::Element::hex27:
      GEO::fillInitialPositionArray<DRT::Element::hex27, 3>(ele, xyze);
      break;
    case DRT::Element::line2:
      switch (probdim)
      {
        case 2:
          GEO::fillInitialPositionArray<DRT::Element::line2, 2>(ele, xyze);
          break;
        case 3:
          GEO::fillInitialPositionArray<DRT::Element::line2, 3>(ele, xyze);
          break;
        default:
          dserror("Unsupported problem dimension! (probdim = %d)", probdim);
          exit(EXIT_FAILURE);
      }
      break;
    default:
      dserror("Unknown elmenet type ( type = %s )", DRT::DistypeToString(distype).c_str());
      break;
  }

  // we assume one dof per node here
  phi_nodes.resize(ele->NumNode(), 0.0);
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // get element location vector
  lm.clear();
  lmowner.clear();
  lmstride.clear();
  ele->LocationVector(scatradis, lm, lmowner, lmstride);
  DRT::UTILS::ExtractMyValues(phicol, phi_nodes, lm);

  // define nodal ID's
  node_ids.resize(numnode, 0.0);
  for (std::size_t i = 0; i < numnode; ++i) node_ids[i] = i;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::ElementHandle* SCATRA::LEVELSET::Intersection::Cut(
    GEO::CUT::LevelSetIntersection& levelset, const LINALG::SerialDenseMatrix& xyze,
    const std::vector<double>& phi_nodes, bool cut_screenoutput) const
{
  try
  {
    levelset.Cut(true, cut_screenoutput);
  }
  catch (std::runtime_error& err)
  {
    std::cerr << "\n--- failed to cut element ---\n"
              << "coordinates:\n";
    xyze.Print(std::cerr);
    std::cerr << "g-function values:\n" << std::setprecision(16);
    std::copy(phi_nodes.begin(), phi_nodes.end(), std::ostream_iterator<double>(std::cerr, ", "));
    std::cerr << "\n";
    run_time_error(err);
  }

  return levelset.GetElement(1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool SCATRA::LEVELSET::Intersection::IsPointPosition(const GEO::CUT::Point::PointPosition& curr_pos,
    const std::vector<GEO::CUT::Point::PointPosition>& desired_pos) const
{
  for (std::vector<GEO::CUT::Point::PointPosition>::const_iterator cit = desired_pos.begin();
       cit != desired_pos.end(); ++cit)
  {
    // OR - combination
    if (curr_pos == *cit) return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<GEO::CUT::Point::PointPosition>&
SCATRA::LEVELSET::Intersection::DesiredPositions()
{
  if (desired_positions_.empty()) desired_positions_.push_back(GEO::CUT::Point::outside);
  return desired_positions_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::SetDesiredPositions(
    const std::vector<GEO::CUT::Point::PointPosition>& desired_pos)
{
  desired_positions_.resize(desired_pos.size(), GEO::CUT::Point::undecided);
  std::copy(desired_pos.begin(), desired_pos.end(), desired_positions_.begin());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::ExportInterface(
    std::map<int, GEO::BoundaryIntCells>& myinterface, const Epetra_Comm& comm)
{
  //-------------------------------
  // prepare parallel communication
  //-------------------------------
  const int myrank = comm.MyPID();
  const int numproc = comm.NumProc();

  int size_one = 1;

  DRT::Exporter exporter(comm);

  // destination proc (the "next" one)
  int dest = myrank + 1;
  if (myrank == (numproc - 1)) dest = 0;

  // source proc (the "last" one)
  int source = myrank - 1;
  if (myrank == 0) source = numproc - 1;

#ifdef DEBUG
  IO::cout << "proc " << myrank << " interface pieces for " << myinterface.size()
           << " elements available before export" << IO::endl;
#endif

  DRT::PackBuffer data;
  packBoundaryIntCells(myinterface, data);
  data.StartPacking();
  packBoundaryIntCells(myinterface, data);

  //-----------------------------------------------------------------
  // pack data (my boundary integration cell groups) for initial send
  //-----------------------------------------------------------------
  std::vector<char> dataSend;
  swap(dataSend, data());

  //-----------------------------------------------
  // send data around in a circle to all processors
  //-----------------------------------------------
  // loop over processors
  for (int num = 0; num < numproc - 1; num++)
  {
    std::vector<int> lengthSend(1, 0);
    lengthSend[0] = dataSend.size();

#ifdef DEBUG
    IO::cout << "--- sending " << lengthSend[0] << " bytes: from proc " << myrank << " to proc "
             << dest << IO::endl;
#endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter.ISend(myrank, dest, &(lengthSend[0]), size_one, length_tag, req_length_data);
    // ... and receive length
    std::vector<int> lengthRecv(1, 0);
    exporter.Receive(source, length_tag, lengthRecv, size_one);
    exporter.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter.ISend(myrank, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
    // ... and receive data
    std::vector<char> dataRecv(lengthRecv[0]);
    exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter.Wait(req_data);

#ifdef DEBUG
    IO::cout << "--- receiving " << lengthRecv[0] << " bytes: to proc " << myrank << " from proc "
             << source << IO::endl;
#endif

    //-----------------------------------------------
    // unpack data (boundary integration cell groups)
    //-----------------------------------------------
    std::map<int, GEO::BoundaryIntCells> interface_recv;

    unpackBoundaryIntCells(dataRecv, interface_recv);

    // add group of cells to my interface map
    /* remark: all groups of boundary integration cells (interface pieces
     * within an element) are collected here */
    for (std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup = interface_recv.begin();
         cellgroup != interface_recv.end(); ++cellgroup)
    {
      myinterface.insert(*cellgroup);
    }

    // make received data the new 'to be sent' data
    dataSend = dataRecv;

    // processors wait for each other
    comm.Barrier();
  }
#ifdef DEBUG
  IO::cout << "proc " << myrank << " interface pieces for " << myinterface.size()
           << " elements available after export" << IO::endl;
#endif
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::packBoundaryIntCells(
    const std::map<int, GEO::BoundaryIntCells>& intcellmap, DRT::PackBuffer& dataSend)
{
  // pack data on all processors
  // loop entries of map (groups of boundary integration cells)
  for (std::map<int, GEO::BoundaryIntCells>::const_iterator cellgroup = intcellmap.begin();
       cellgroup != intcellmap.end(); ++cellgroup)
  {
    // pack data of all boundary integrations cells belonging to an element
    const int elegid = cellgroup->first;
    DRT::ParObject::AddtoPack(dataSend, elegid);

    const int numcells = (cellgroup->second).size();
    DRT::ParObject::AddtoPack(dataSend, numcells);

    for (int icell = 0; icell < numcells; ++icell)
    {
      GEO::BoundaryIntCell cell = cellgroup->second[icell];
      // get all member variables from a single boundary integration cell
      const DRT::Element::DiscretizationType distype = cell.Shape();
      DRT::ParObject::AddtoPack(dataSend, distype);

      // coordinates of cell vertices in (scatra) element parameter space
      //      const Epetra_SerialDenseMatrix& vertices_xi = cell.CellNodalPosXiDomain();
      //      const LINALG::SerialDenseMatrix& vertices_xi = cell.CellNodalPosXiDomain();
      const LINALG::SerialDenseMatrix vertices_xi = cell.CellNodalPosXiDomain();
      DRT::ParObject::AddtoPack(dataSend, vertices_xi);

      // coordinates of cell vertices in physical space
      //      const Epetra_SerialDenseMatrix& vertices_xyz = cell.CellNodalPosXYZ();
      //      const LINALG::SerialDenseMatrix& vertices_xyz = cell.CellNodalPosXYZ();
      const LINALG::SerialDenseMatrix vertices_xyz = cell.CellNodalPosXYZ();
      DRT::ParObject::AddtoPack(dataSend, vertices_xyz);
    }
  }
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void SCATRA::LEVELSET::Intersection::unpackBoundaryIntCells(
    const std::vector<char>& data, std::map<int, GEO::BoundaryIntCells>& intcellmap)
{
  // pointer to current position in a group of cells in local std::string (counts bytes)
  std::vector<char>::size_type posingroup = 0;

  while (posingroup < data.size())
  {
    // extract fluid element gid
    int elegid = -1;
    DRT::ParObject::ExtractfromPack(posingroup, data, elegid);
    if (elegid < 0) dserror("extraction of element gid failed");

    // extract number of boundary integration cells for this element
    int numvecs = -1;
    DRT::ParObject::ExtractfromPack(posingroup, data, numvecs);

    // vector holding group of boundary integration cells belonging to this element
    GEO::BoundaryIntCells intcellvector;

    for (int icell = 0; icell < numvecs; ++icell)
    {
      //--------------------------------------------------------------------
      // extract all member variables for a single boundary integration cell
      //--------------------------------------------------------------------
      // distype of cell
      DRT::Element::DiscretizationType distype;
      int distypeint = -1;
      DRT::ParObject::ExtractfromPack(posingroup, data, distypeint);
      distype = (DRT::Element::DiscretizationType)distypeint;
      if (!(distype == DRT::Element::tri3 || distype == DRT::Element::quad4))
        dserror("unexpected distype %d", distypeint);

      LINALG::SerialDenseMatrix vertices_xi;
      DRT::ParObject::ExtractfromPack(posingroup, data, vertices_xi);

      // coordinates of cell vertices in physical space
      LINALG::SerialDenseMatrix vertices_xyz;
      DRT::ParObject::ExtractfromPack(posingroup, data, vertices_xyz);

      // store boundary integration cells in boundaryintcelllist
      intcellvector.push_back(
          GEO::BoundaryIntCell(distype, -1, vertices_xi, Teuchos::null, vertices_xyz));
    }

    // add group of cells for this element to the map
    intcellmap.insert(std::make_pair(elegid, intcellvector));
  }
  // check correct reading
  if (posingroup != data.size())
    dserror("mismatch in size of data %d <-> %d", (int)data.size(), posingroup);
}


template void SCATRA::LEVELSET::Intersection::GetZeroLevelSet<GEO::BoundaryIntCells>(
    const Epetra_Vector& phi, const DRT::Discretization& scatradis,
    std::map<int, GEO::BoundaryIntCells>& elementBoundaryIntCells, bool cut_screenoutput);
template void SCATRA::LEVELSET::Intersection::GetZeroLevelSet<GEO::BoundaryIntCellPtrs>(
    const Epetra_Vector& phi, const DRT::Discretization& scatradis,
    std::map<int, GEO::BoundaryIntCellPtrs>& elementBoundaryIntCells, bool cut_screenoutput);
