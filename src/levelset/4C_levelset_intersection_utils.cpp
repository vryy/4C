/*----------------------------------------------------------------------*/
/*! \file

 \brief detailed description in header file levelset_intersection.H

\level 2

 *------------------------------------------------------------------------------------------------*/


#include "4C_levelset_intersection_utils.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_cut_integrationcell.hpp"
#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_geometry_element_coordtrafo.hpp"
#include "4C_fem_geometry_element_volume.hpp"
#include "4C_fem_geometry_integrationcell.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_scatra_ele_parameter_std.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ScaTra::LevelSet::Intersection::Intersection()
    : check_lsv_(false), desired_positions_(0), volumeplus_(0.0), volumeminus_(0.0), surface_(0.0)
{
  desired_positions_.reserve(2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::reset()
{
  volumeplus_ = 0.0;
  volumeminus_ = 0.0;
  surface_ = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::CaptureZeroLevelSet(
    const Teuchos::RCP<const Epetra_Vector>& phi,
    const Teuchos::RCP<const Core::FE::Discretization>& scatradis, double& volumedomainminus,
    double& volumedomainplus, double& zerosurface,
    std::map<int, Core::Geo::BoundaryIntCells>& elementBoundaryIntCells)
{
  // reset, just to be sure
  reset();
  volumedomainminus = 0.0;
  volumedomainplus = 0.0;
  zerosurface = 0.0;
  elementBoundaryIntCells.clear();

  // herein the actual capturing happens
  get_zero_level_set(*phi, *scatradis, elementBoundaryIntCells);

  // collect contributions from all procs and store in respective variables
  scatradis->Comm().SumAll(&volume_plus(), &volumedomainplus, 1);
  scatradis->Comm().SumAll(&volume_minus(), &volumedomainminus, 1);
  scatradis->Comm().SumAll(&surface(), &zerosurface, 1);

  // export also interface to all procs
  export_interface(elementBoundaryIntCells, scatradis->Comm());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T>
void ScaTra::LevelSet::Intersection::get_zero_level_set(const Epetra_Vector& phi,
    const Core::FE::Discretization& scatradis, std::map<int, T>& elementBoundaryIntCells,
    bool cut_screenoutput)
{
  // export phi from row to column map
  const Teuchos::RCP<Epetra_Vector> phicol =
      Teuchos::rcp(new Epetra_Vector(*scatradis.DofColMap()));
  Core::LinAlg::Export(phi, *phicol);

  // remark: loop over row elements is sufficient
  for (int iele = 0; iele < scatradis.NumMyRowElements(); ++iele)
  {
    // get element from discretization
    const Core::Elements::Element* ele = scatradis.lRowElement(iele);
    const Core::FE::CellType distype = ele->Shape();

    // clear vector each loop
    boundary_int_cells_per_ele<T>().clear();

    // ------------------------------------------------------------------------
    // Prepare cut
    // ------------------------------------------------------------------------
    Core::Geo::Cut::LevelSetIntersection levelset(scatradis.Comm());
    Core::LinAlg::SerialDenseMatrix xyze;
    std::vector<double> phi_nodes;
    std::vector<int> nids;
    prepare_cut(ele, scatradis, *phicol, xyze, phi_nodes, nids);

    // check if this element is cut, according to its level-set values
    // -> add it to 'levelset'
    // note: cut is performed in physical space
    if (!levelset.add_element(1, nids, xyze, ele->Shape(), phi_nodes.data(), false, check_lsv_))
      continue;

    // ------------------------------------------------------------------------
    // call Core::Geo::Cut algorithm and process cut data
    // ------------------------------------------------------------------------
    Core::Geo::Cut::ElementHandle* ehandle = cut(levelset, xyze, phi_nodes, cut_screenoutput);

    // =========================================================
    // cell is in contact with the interface (cut or touched)
    // =========================================================
    if (ehandle != nullptr)
    {
      Core::Geo::Cut::plain_element_set cuteles;

      collect_cut_eles(*ehandle, cuteles, distype);

      // ----------------------------------------------------------------------
      // get zero level-set contour
      // ----------------------------------------------------------------------
      get_zero_level_set_contour(cuteles, xyze, distype);
    }
    // =========================================================
    // element is uncut
    // =========================================================
    else
    {
      double elevol = Core::Geo::ElementVolume(distype, xyze);

      // it is sufficient to check the first node, since the element entirely
      // lies within the plus or minus domain
      if (phi_nodes[0] > 0.0)
        volume_plus() += elevol;
      else
        volume_minus() += elevol;
    }

    // store interface of element
    if (boundary_int_cells_per_ele<T>().size() > 0)
      elementBoundaryIntCells[ele->Id()] = boundary_int_cells_per_ele<T>();
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::get_zero_level_set_contour(
    const Core::Geo::Cut::plain_element_set& cuteles, const Core::LinAlg::SerialDenseMatrix& xyze,
    Core::FE::CellType distype)
{
  for (Core::Geo::Cut::plain_element_set::const_iterator icutele = cuteles.begin();
       icutele != cuteles.end(); ++icutele)
  {
    // get pointer to cut element
    Core::Geo::Cut::Element* cutele = *icutele;

    Core::Geo::Cut::plain_volumecell_set volcells;
    volcells = cutele->VolumeCells();

    for (Core::Geo::Cut::plain_volumecell_set::const_iterator ivolcell = volcells.begin();
         ivolcell != volcells.end(); ++ivolcell)
    {
      Core::Geo::Cut::VolumeCell* volcell = *ivolcell;
      const Core::Geo::Cut::Point::PointPosition vol_pos = volcell->Position();
      if (is_point_position(vol_pos))
      {
        add_to_volume(vol_pos, volcell->Volume());
        // get boundary integration cells for this volume cell
        // we consider only the cells for one position, otherwise we would have the boundary
        // cells twice
        const Core::Geo::Cut::plain_boundarycell_set& bcells = volcell->BoundaryCells();
        for (Core::Geo::Cut::plain_boundarycell_set::const_iterator ibcell = bcells.begin();
             ibcell != bcells.end(); ++ibcell)
        {
          Core::Geo::Cut::BoundaryCell* bcell = *ibcell;

          add_to_boundary_int_cells_per_ele(xyze, *bcell, distype);

          surface() += bcell->Area();
        }
      }
      else
      {
        add_to_volume(vol_pos, volcell->Volume());
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::add_to_boundary_int_cells_per_ele(
    const Core::LinAlg::SerialDenseMatrix& xyze, const Core::Geo::Cut::BoundaryCell& bcell,
    Core::FE::CellType distype_ele)
{
  Core::FE::CellType distype_bc = bcell.Shape();
  check_boundary_cell_type(distype_bc);

  const int numnodebc = Core::FE::getNumberOfElementNodes(distype_bc);

  // get physical coordinates of this cell
  Core::LinAlg::SerialDenseMatrix coord = bcell.Coordinates();

  // transfer to element coordinates
  Core::LinAlg::SerialDenseMatrix localcoord(3, numnodebc, true);

  for (int ivert = 0; ivert < numnodebc; ivert++)
  {
    Core::LinAlg::Matrix<3, 1> lcoord;
    Core::LinAlg::Matrix<3, 1> pcoord;
    for (int ll = 0; ll < 3; ll++) pcoord(ll, 0) = coord(ll, ivert);

    Core::Geo::currentToVolumeElementCoordinates(distype_ele, xyze, pcoord, lcoord);

    // write as 'physCoord'
    for (int ll = 0; ll < 3; ll++) localcoord(ll, ivert) = lcoord(ll, 0);
  }

  // store boundary element and sum area into surface
  // be careful, we only set physical coordinates
  Core::LinAlg::SerialDenseMatrix dummyMat;
  boundary_int_cells_per_ele<Core::Geo::BoundaryIntCells>().push_back(
      Core::Geo::BoundaryIntCell(distype_bc, -1, localcoord, dummyMat, coord, true));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::check_boundary_cell_type(Core::FE::CellType distype_bc) const
{
  if (distype_bc != Core::FE::CellType::tri3 and distype_bc != Core::FE::CellType::quad4)
  {
    FOUR_C_THROW("unexpected type of boundary integration cell: %s",
        Core::FE::CellTypeToString(distype_bc).c_str());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::add_to_volume(
    Core::Geo::Cut::Point::PointPosition pos, double vol)
{
  switch (pos)
  {
    case Core::Geo::Cut::Point::outside:
      volume_plus() += vol;
      break;
    case Core::Geo::Cut::Point::inside:
      volume_minus() += vol;
      break;
    default:
      /* do nothing for the undecided case */
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::collect_cut_eles(Core::Geo::Cut::ElementHandle& ehandle,
    Core::Geo::Cut::plain_element_set& cuteles, Core::FE::CellType distype) const
{
  ehandle.CollectElements(cuteles);

  switch (distype)
  {
    case Core::FE::CellType::line2:
    case Core::FE::CellType::hex8:
    {
      if (cuteles.size() != 1) FOUR_C_THROW("one cut element expected for linear elements");
      break;
    }
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
    {
      if (cuteles.size() != 8) FOUR_C_THROW("eight cut elements expected for quadratic elements");
      break;
    }
    default:
    {
      FOUR_C_THROW("distype unknown for level set cut algorithm");
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::prepare_cut(const Core::Elements::Element* ele,
    const Core::FE::Discretization& scatradis, const Epetra_Vector& phicol,
    Core::LinAlg::SerialDenseMatrix& xyze, std::vector<double>& phi_nodes,
    std::vector<int>& node_ids) const
{
  const Core::FE::CellType distype = ele->Shape();
  unsigned numnode = Core::FE::getNumberOfElementNodes(distype);
  const unsigned probdim = Global::Problem::Instance()->NDim();

  xyze.shape(3, numnode);
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      Core::Geo::fillInitialPositionArray<Core::FE::CellType::hex8, 3>(ele, xyze);
      break;
    case Core::FE::CellType::hex20:
      Core::Geo::fillInitialPositionArray<Core::FE::CellType::hex20, 3>(ele, xyze);
      break;
    case Core::FE::CellType::hex27:
      Core::Geo::fillInitialPositionArray<Core::FE::CellType::hex27, 3>(ele, xyze);
      break;
    case Core::FE::CellType::line2:
      switch (probdim)
      {
        case 2:
          Core::Geo::fillInitialPositionArray<Core::FE::CellType::line2, 2>(ele, xyze);
          break;
        case 3:
          Core::Geo::fillInitialPositionArray<Core::FE::CellType::line2, 3>(ele, xyze);
          break;
        default:
          FOUR_C_THROW("Unsupported problem dimension! (probdim = %d)", probdim);
          exit(EXIT_FAILURE);
      }
      break;
    default:
      FOUR_C_THROW(
          "Unknown elmenet type ( type = %s )", Core::FE::CellTypeToString(distype).c_str());
      break;
  }

  // we assume one dof per node here
  phi_nodes.resize(ele->num_node(), 0.0);
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // get element location vector
  lm.clear();
  lmowner.clear();
  lmstride.clear();
  ele->LocationVector(scatradis, lm, lmowner, lmstride);
  Core::FE::ExtractMyValues(phicol, phi_nodes, lm);

  // define nodal ID's
  node_ids.resize(numnode, 0.0);
  for (std::size_t i = 0; i < numnode; ++i) node_ids[i] = i;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::ElementHandle* ScaTra::LevelSet::Intersection::cut(
    Core::Geo::Cut::LevelSetIntersection& levelset, const Core::LinAlg::SerialDenseMatrix& xyze,
    const std::vector<double>& phi_nodes, bool cut_screenoutput) const
{
  try
  {
    levelset.Cut(true, cut_screenoutput);
  }
  catch (Core::Exception& err)
  {
    std::cerr << "\n--- failed to cut element ---\n"
              << "coordinates:\n";
    std::cerr << xyze;
    std::cerr << "g-function values:\n" << std::setprecision(16);
    std::copy(phi_nodes.begin(), phi_nodes.end(), std::ostream_iterator<double>(std::cerr, ", "));
    std::cerr << "\n";
    throw;
  }

  return levelset.GetElement(1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool ScaTra::LevelSet::Intersection::is_point_position(
    const Core::Geo::Cut::Point::PointPosition& curr_pos,
    const std::vector<Core::Geo::Cut::Point::PointPosition>& desired_pos) const
{
  for (std::vector<Core::Geo::Cut::Point::PointPosition>::const_iterator cit = desired_pos.begin();
       cit != desired_pos.end(); ++cit)
  {
    // OR - combination
    if (curr_pos == *cit) return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<Core::Geo::Cut::Point::PointPosition>&
ScaTra::LevelSet::Intersection::desired_positions()
{
  if (desired_positions_.empty()) desired_positions_.push_back(Core::Geo::Cut::Point::outside);
  return desired_positions_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::SetDesiredPositions(
    const std::vector<Core::Geo::Cut::Point::PointPosition>& desired_pos)
{
  desired_positions_.resize(desired_pos.size(), Core::Geo::Cut::Point::undecided);
  std::copy(desired_pos.begin(), desired_pos.end(), desired_positions_.begin());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::export_interface(
    std::map<int, Core::Geo::BoundaryIntCells>& myinterface, const Epetra_Comm& comm)
{
  //-------------------------------
  // prepare parallel communication
  //-------------------------------
  const int myrank = comm.MyPID();
  const int numproc = comm.NumProc();

  int size_one = 1;

  Core::Communication::Exporter exporter(comm);

  // destination proc (the "next" one)
  int dest = myrank + 1;
  if (myrank == (numproc - 1)) dest = 0;

  // source proc (the "last" one)
  int source = myrank - 1;
  if (myrank == 0) source = numproc - 1;

#ifdef FOUR_C_ENABLE_ASSERTIONS
  Core::IO::cout << "proc " << myrank << " interface pieces for " << myinterface.size()
                 << " elements available before export" << Core::IO::endl;
#endif

  Core::Communication::PackBuffer data;
  pack_boundary_int_cells(myinterface, data);
  data.StartPacking();
  pack_boundary_int_cells(myinterface, data);

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

#ifdef FOUR_C_ENABLE_ASSERTIONS
    Core::IO::cout << "--- sending " << lengthSend[0] << " bytes: from proc " << myrank
                   << " to proc " << dest << Core::IO::endl;
#endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter.i_send(myrank, dest, lengthSend.data(), size_one, length_tag, req_length_data);
    // ... and receive length
    std::vector<int> lengthRecv(1, 0);
    exporter.Receive(source, length_tag, lengthRecv, size_one);
    exporter.Wait(req_length_data);

    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter.i_send(myrank, dest, dataSend.data(), lengthSend[0], data_tag, req_data);
    // ... and receive data
    std::vector<char> dataRecv(lengthRecv[0]);
    exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter.Wait(req_data);

#ifdef FOUR_C_ENABLE_ASSERTIONS
    Core::IO::cout << "--- receiving " << lengthRecv[0] << " bytes: to proc " << myrank
                   << " from proc " << source << Core::IO::endl;
#endif

    //-----------------------------------------------
    // unpack data (boundary integration cell groups)
    //-----------------------------------------------
    std::map<int, Core::Geo::BoundaryIntCells> interface_recv;

    unpack_boundary_int_cells(dataRecv, interface_recv);

    // add group of cells to my interface map
    /* remark: all groups of boundary integration cells (interface pieces
     * within an element) are collected here */
    for (std::map<int, Core::Geo::BoundaryIntCells>::const_iterator cellgroup =
             interface_recv.begin();
         cellgroup != interface_recv.end(); ++cellgroup)
    {
      myinterface.insert(*cellgroup);
    }

    // make received data the new 'to be sent' data
    dataSend = dataRecv;

    // processors wait for each other
    comm.Barrier();
  }
#ifdef FOUR_C_ENABLE_ASSERTIONS
  Core::IO::cout << "proc " << myrank << " interface pieces for " << myinterface.size()
                 << " elements available after export" << Core::IO::endl;
#endif
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::pack_boundary_int_cells(
    const std::map<int, Core::Geo::BoundaryIntCells>& intcellmap,
    Core::Communication::PackBuffer& dataSend)
{
  // pack data on all processors
  // loop entries of map (groups of boundary integration cells)
  for (std::map<int, Core::Geo::BoundaryIntCells>::const_iterator cellgroup = intcellmap.begin();
       cellgroup != intcellmap.end(); ++cellgroup)
  {
    // pack data of all boundary integrations cells belonging to an element
    const int elegid = cellgroup->first;
    Core::Communication::ParObject::AddtoPack(dataSend, elegid);

    const int numcells = (cellgroup->second).size();
    Core::Communication::ParObject::AddtoPack(dataSend, numcells);

    for (int icell = 0; icell < numcells; ++icell)
    {
      Core::Geo::BoundaryIntCell cell = cellgroup->second[icell];
      // get all member variables from a single boundary integration cell
      const Core::FE::CellType distype = cell.Shape();
      Core::Communication::ParObject::AddtoPack(dataSend, distype);

      // coordinates of cell vertices in (scatra) element parameter space
      const Core::LinAlg::SerialDenseMatrix vertices_xi = cell.cell_nodal_pos_xi_domain();
      Core::Communication::ParObject::AddtoPack(dataSend, vertices_xi);

      // coordinates of cell vertices in physical space
      const Core::LinAlg::SerialDenseMatrix vertices_xyz = cell.CellNodalPosXYZ();
      Core::Communication::ParObject::AddtoPack(dataSend, vertices_xyz);
    }
  }
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void ScaTra::LevelSet::Intersection::unpack_boundary_int_cells(
    const std::vector<char>& data, std::map<int, Core::Geo::BoundaryIntCells>& intcellmap)
{
  // pointer to current position in a group of cells in local std::string (counts bytes)
  std::vector<char>::size_type posingroup = 0;

  while (posingroup < data.size())
  {
    // extract fluid element gid
    int elegid = -1;
    Core::Communication::ParObject::ExtractfromPack(posingroup, data, elegid);
    if (elegid < 0) FOUR_C_THROW("extraction of element gid failed");

    // extract number of boundary integration cells for this element
    int numvecs = -1;
    Core::Communication::ParObject::ExtractfromPack(posingroup, data, numvecs);

    // vector holding group of boundary integration cells belonging to this element
    Core::Geo::BoundaryIntCells intcellvector;

    for (int icell = 0; icell < numvecs; ++icell)
    {
      //--------------------------------------------------------------------
      // extract all member variables for a single boundary integration cell
      //--------------------------------------------------------------------
      // distype of cell
      Core::FE::CellType distype;
      int distypeint = -1;
      Core::Communication::ParObject::ExtractfromPack(posingroup, data, distypeint);
      distype = (Core::FE::CellType)distypeint;
      if (!(distype == Core::FE::CellType::tri3 || distype == Core::FE::CellType::quad4))
        FOUR_C_THROW("unexpected distype %d", distypeint);

      Core::LinAlg::SerialDenseMatrix vertices_xi;
      Core::Communication::ParObject::ExtractfromPack(posingroup, data, vertices_xi);

      // coordinates of cell vertices in physical space
      Core::LinAlg::SerialDenseMatrix vertices_xyz;
      Core::Communication::ParObject::ExtractfromPack(posingroup, data, vertices_xyz);

      // store boundary integration cells in boundaryintcelllist
      Core::LinAlg::SerialDenseMatrix dummyMat;
      intcellvector.push_back(
          Core::Geo::BoundaryIntCell(distype, -1, vertices_xi, dummyMat, vertices_xyz));
    }

    // add group of cells for this element to the map
    intcellmap.insert(std::make_pair(elegid, intcellvector));
  }
  // check correct reading
  if (posingroup != data.size())
    FOUR_C_THROW("mismatch in size of data %d <-> %d", (int)data.size(), posingroup);
}


template void ScaTra::LevelSet::Intersection::get_zero_level_set<Core::Geo::BoundaryIntCells>(
    const Epetra_Vector& phi, const Core::FE::Discretization& scatradis,
    std::map<int, Core::Geo::BoundaryIntCells>& elementBoundaryIntCells, bool cut_screenoutput);
template void ScaTra::LevelSet::Intersection::get_zero_level_set<Core::Geo::BoundaryIntCellPtrs>(
    const Epetra_Vector& phi, const Core::FE::Discretization& scatradis,
    std::map<int, Core::Geo::BoundaryIntCellPtrs>& elementBoundaryIntCells, bool cut_screenoutput);

FOUR_C_NAMESPACE_CLOSE
