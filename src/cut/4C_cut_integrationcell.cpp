/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells

\level 3


*----------------------------------------------------------------------*/
#include "4C_cut_integrationcell.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_facet.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_position.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_geometry_element_volume.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::IntegrationCell::Contains(Core::LinAlg::Matrix<3, 1>& x)
{
  switch (this->Shape())
  {
    case Core::FE::CellType::tet4:
    {
      // find element local position of gauss point
      return Contains<3, Core::FE::CellType::tet4>(x);
    }
    case Core::FE::CellType::hex8:
    {
      return Contains<3, Core::FE::CellType::hex8>(x);
    }
    default:
    {
      FOUR_C_THROW("unknown type of integration cell ");
      break;
    }
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, Core::FE::CellType celltype>
bool Core::Geo::Cut::IntegrationCell::Contains(Core::LinAlg::Matrix<probdim, 1>& x)
{
  const int ncn = Core::FE::num_nodes<celltype>;

  Core::LinAlg::Matrix<probdim, ncn> coords(xyz_);

  Teuchos::RCP<Core::Geo::Cut::Position> pos =
      Core::Geo::Cut::Position::Create(coords, x, celltype);
  pos->Compute();

  return pos->within_limits();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::IntegrationCell::DumpGmsh(std::ofstream& file, int* value)
{
  Output::GmshCellDump(file, Shape(), xyz_, &position_, value);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double Core::Geo::Cut::IntegrationCell::Volume() const
{
  return Core::Geo::ElementVolume(Shape(), xyz_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Line2IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  // not 100% sure what this value really means, but 4 seems more than sufficient.
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Tri3IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Quad4IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Hex8IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  switch (elementshape)
  {
    case Core::FE::CellType::hex8:
      return 6;
    case Core::FE::CellType::hex20:
      return 15;
    case Core::FE::CellType::hex27:
      return 15;
    case Core::FE::CellType::tet4:
      return 6;
    case Core::FE::CellType::tet10:
      return 6;
    case Core::FE::CellType::wedge6:
      return 6;
    case Core::FE::CellType::wedge15:
      return 14;
    case Core::FE::CellType::pyramid5:
      return 6;
    default:
      FOUR_C_THROW("no rule defined for this element type");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Tet4IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  switch (elementshape)
  {
    case Core::FE::CellType::hex8:
      return 6;
    case Core::FE::CellType::hex20:
      return 15;
    case Core::FE::CellType::hex27:
      return 15;
    case Core::FE::CellType::tet4:
      return 6;
    case Core::FE::CellType::tet10:
      return 7;
    case Core::FE::CellType::wedge6:
      return 6;
    case Core::FE::CellType::wedge15:
      return 14;
    case Core::FE::CellType::pyramid5:
      return 6;
    default:
      FOUR_C_THROW("no rule defined for this element type");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Wedge6IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Core::Geo::Cut::Pyramid5IntegrationCell::CubatureDegree(Core::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::IntegrationCell::Print(std::ostream& stream) const
{
  stream << "--- integration cell ( address: " << std::setw(10) << this << " )\n";
  stream << "pos = " << Point::point_position2_string(Position()) << " "
         << "shape = " << Core::FE::CellTypeToString(Shape()) << " "
         << "volume = " << Volume() << "\n";
  for (unsigned i = 0; i < points_.size(); ++i)
  {
    (points_)[i]->Print(stream);
    stream << "\n";
  }
}

FOUR_C_NAMESPACE_CLOSE
