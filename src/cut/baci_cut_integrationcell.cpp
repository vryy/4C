/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells

\level 3


*----------------------------------------------------------------------*/
#include "baci_cut_integrationcell.hpp"

#include "baci_cut_boundarycell.hpp"
#include "baci_cut_facet.hpp"
#include "baci_cut_mesh.hpp"
#include "baci_cut_output.hpp"
#include "baci_cut_position.hpp"
#include "baci_cut_volumecell.hpp"
#include "baci_discretization_geometry_element_volume.hpp"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CORE::GEO::CUT::IntegrationCell::Contains(CORE::LINALG::Matrix<3, 1>& x)
{
  switch (this->Shape())
  {
    case CORE::FE::CellType::tet4:
    {
      // find element local position of gauss point
      return Contains<3, CORE::FE::CellType::tet4>(x);
    }
    case CORE::FE::CellType::hex8:
    {
      return Contains<3, CORE::FE::CellType::hex8>(x);
    }
    default:
    {
      dserror("unknown type of integration cell ");
      break;
    }
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, CORE::FE::CellType celltype>
bool CORE::GEO::CUT::IntegrationCell::Contains(CORE::LINALG::Matrix<probdim, 1>& x)
{
  const int ncn = CORE::FE::num_nodes<celltype>;

  CORE::LINALG::Matrix<probdim, ncn> coords(xyz_);

  Teuchos::RCP<CORE::GEO::CUT::Position> pos =
      CORE::GEO::CUT::Position::Create(coords, x, celltype);
  pos->Compute();

  return pos->WithinLimits();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::IntegrationCell::DumpGmsh(std::ofstream& file, int* value)
{
  OUTPUT::GmshCellDump(file, Shape(), xyz_, &position_, value);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CORE::GEO::CUT::IntegrationCell::Volume() const
{
  return CORE::GEO::ElementVolume(Shape(), xyz_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Line2IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  // not 100% sure what this value really means, but 4 seems more than sufficient.
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Tri3IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Quad4IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Hex8IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  switch (elementshape)
  {
    case CORE::FE::CellType::hex8:
      return 6;
    case CORE::FE::CellType::hex20:
      return 15;
    case CORE::FE::CellType::hex27:
      return 15;
    case CORE::FE::CellType::tet4:
      return 6;
    case CORE::FE::CellType::tet10:
      return 6;
    case CORE::FE::CellType::wedge6:
      return 6;
    case CORE::FE::CellType::wedge15:
      return 14;
    case CORE::FE::CellType::pyramid5:
      return 6;
    default:
      dserror("no rule defined for this element type");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Tet4IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  switch (elementshape)
  {
    case CORE::FE::CellType::hex8:
      return 6;
    case CORE::FE::CellType::hex20:
      return 15;
    case CORE::FE::CellType::hex27:
      return 15;
    case CORE::FE::CellType::tet4:
      return 6;
    case CORE::FE::CellType::tet10:
      return 7;
    case CORE::FE::CellType::wedge6:
      return 6;
    case CORE::FE::CellType::wedge15:
      return 14;
    case CORE::FE::CellType::pyramid5:
      return 6;
    default:
      dserror("no rule defined for this element type");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Wedge6IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CORE::GEO::CUT::Pyramid5IntegrationCell::CubatureDegree(CORE::FE::CellType elementshape) const
{
  return 4;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CORE::GEO::CUT::IntegrationCell::Print(std::ostream& stream) const
{
  stream << "--- integration cell ( address: " << std::setw(10) << this << " )\n";
  stream << "pos = " << Point::PointPosition2String(Position()) << " "
         << "shape = " << CORE::FE::CellTypeToString(Shape()) << " "
         << "volume = " << Volume() << "\n";
  for (unsigned i = 0; i < points_.size(); ++i)
  {
    (points_)[i]->Print(stream);
    stream << "\n";
  }
}

BACI_NAMESPACE_CLOSE
