/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of gauss points access functions

\level 1


*----------------------------------------------------------------------*/

#include "4C_fem_general_utils_gausspoints.hpp"

#include <Intrepid_Cubature.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include <stdexcept>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  namespace
  {
    /// wrapper to intrepid gauss point implementation
    template <class Topology>
    class IntrepidGaussPoints : public GaussPoints
    {
     public:
      explicit IntrepidGaussPoints(int cubDegree)
      {
        // cell type: tetrahedron
        shards::CellTopology cellType = shards::getCellTopologyData<Topology>();

        // retrieve spatial dimension
        int spaceDim = cellType.getDimension();
        // int numNodes = cellType.getNodeCount();

        // create cubature factory
        Intrepid::DefaultCubatureFactory<double> cubFactory;

        // create default cubature
        Teuchos::RCP<Intrepid::Cubature<double>> myCub = cubFactory.create(cellType, cubDegree);

        // retrieve number of cubature points
        int numCubPoints = myCub->getNumPoints();

        cub_points_.resize(numCubPoints, spaceDim);
        cub_weights_.resize(numCubPoints);

        // retrieve cubature points and weights
        myCub->getCubature(cub_points_, cub_weights_);
      }

      int NumPoints() const override { return cub_points_.dimension(0); }

      int NumDimension() const override { return cub_points_.dimension(1); }

      const double* Point(int point) const override { return &cub_points_(point, 0); }

      double Weight(int point) const override { return cub_weights_(point); }

      void print() const override
      {
        // cell type: tetrahedron
        shards::CellTopology cellType = shards::getCellTopologyData<Topology>();

        std::cout << cellType.getName() << " gauss points:\n";
        for (int i = 0; i < NumPoints(); ++i)
        {
          std::cout << "    ";
          for (int j = 0; j < NumDimension(); ++j) std::cout << cub_points_(i, j) << " ";
          std::cout << cub_weights_(i) << "\n";
        }
      }

     private:
      Intrepid::FieldContainer<double> cub_points_;
      Intrepid::FieldContainer<double> cub_weights_;
    };


    /// special case that is not (yet) supported by intrepid
    template <>
    class IntrepidGaussPoints<shards::Pyramid<5>> : public GaussPoints
    {
     public:
      explicit IntrepidGaussPoints<shards::Pyramid<5>>(int cubDegree)
      {
        cub_points_.resize(8, 3);
        cub_weights_.resize(8);

        cub_points_(0, 0) = -0.26318405556971;
        cub_points_(1, 0) = -0.50661630334979;
        cub_points_(2, 0) = -0.26318405556971;
        cub_points_(3, 0) = -0.50661630334979;
        cub_points_(4, 0) = 0.26318405556971;
        cub_points_(5, 0) = 0.50661630334979;
        cub_points_(6, 0) = 0.26318405556971;
        cub_points_(7, 0) = 0.50661630334979;
        cub_points_(0, 1) = -0.26318405556971;
        cub_points_(1, 1) = -0.50661630334979;
        cub_points_(2, 1) = 0.26318405556971;
        cub_points_(3, 1) = 0.50661630334979;
        cub_points_(4, 1) = -0.26318405556971;
        cub_points_(5, 1) = -0.50661630334979;
        cub_points_(6, 1) = 0.26318405556971;
        cub_points_(7, 1) = 0.50661630334979;
        cub_points_(0, 2) = 0.54415184401122;
        cub_points_(1, 2) = 0.12251482265544;
        cub_points_(2, 2) = 0.54415184401122;
        cub_points_(3, 2) = 0.12251482265544;
        cub_points_(4, 2) = 0.54415184401122;
        cub_points_(5, 2) = 0.12251482265544;
        cub_points_(6, 2) = 0.54415184401122;
        cub_points_(7, 2) = 0.12251482265544;

        cub_weights_(0) = 0.10078588207983;
        cub_weights_(1) = 0.23254745125351;
        cub_weights_(2) = 0.10078588207983;
        cub_weights_(3) = 0.23254745125351;
        cub_weights_(4) = 0.10078588207983;
        cub_weights_(5) = 0.23254745125351;
        cub_weights_(6) = 0.10078588207983;
        cub_weights_(7) = 0.23254745125351;
      }

      int NumPoints() const override { return cub_points_.dimension(0); }

      int NumDimension() const override { return cub_points_.dimension(1); }

      const double* Point(int point) const override { return &cub_points_(point, 0); }

      double Weight(int point) const override { return cub_weights_(point); }

      void print() const override
      {
        // cell type: tetrahedron
        shards::CellTopology cellType = shards::getCellTopologyData<shards::Pyramid<5>>();

        std::cout << cellType.getName() << " gauss points:\n";
        for (int i = 0; i < NumPoints(); ++i)
        {
          std::cout << "    ";
          for (int j = 0; j < NumDimension(); ++j) std::cout << cub_points_(i, j) << " ";
          std::cout << cub_weights_(i) << "\n";
        }
      }

     private:
      Intrepid::FieldContainer<double> cub_points_;
      Intrepid::FieldContainer<double> cub_weights_;
    };

  }  // namespace
}  // namespace Core::FE

Core::FE::GaussIntegration::GaussIntegration(Core::FE::CellType distype)
{
  switch (distype)
  {
    case Core::FE::CellType::quad4:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::quad4, 3);
      break;
    case Core::FE::CellType::quad8:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::quad8, 4);
      break;
    case Core::FE::CellType::quad9:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::quad9, 4);
      break;
    case Core::FE::CellType::tri3:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::tri3, 3);
      break;
    case Core::FE::CellType::tri6:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::tri6, 4);
      break;
    case Core::FE::CellType::hex8:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::hex8, 3);
      break;
    case Core::FE::CellType::hex20:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::hex20, 4);
      break;
    case Core::FE::CellType::hex27:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::hex27, 4);
      break;
    case Core::FE::CellType::tet4:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::tet4, 3);
      break;
    case Core::FE::CellType::tet10:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::tet10, 4);
      break;
    case Core::FE::CellType::wedge6:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::wedge6, 3);
      break;
    case Core::FE::CellType::wedge15:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::wedge15, 4);
      break;
    case Core::FE::CellType::pyramid5:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::pyramid5, 3);
      break;
    case Core::FE::CellType::line2:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::line2, 3);
      break;
    case Core::FE::CellType::line3:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::line3, 4);
      break;
    case Core::FE::CellType::nurbs2:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::line2, 3);
      break;
    case Core::FE::CellType::nurbs3:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::line3, 4);
      break;
    case Core::FE::CellType::nurbs4:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::quad4, 3);
      break;
    case Core::FE::CellType::nurbs8:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::hex8, 3);
      break;
    case Core::FE::CellType::nurbs9:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::quad9, 4);
      break;
    case Core::FE::CellType::nurbs27:
      gp_ = GaussPointCache::Instance().Create(Core::FE::CellType::hex27, 4);
      break;
    default:
      FOUR_C_THROW("unsupported element shape");
  }
}

Core::FE::GaussIntegration::GaussIntegration(Core::FE::CellType distype, int degree)
{
  gp_ = GaussPointCache::Instance().Create(distype, degree);
}

Core::FE::GaussPointCache& Core::FE::GaussPointCache::Instance()
{
  static std::unique_ptr<GaussPointCache> instance;
  if (instance == nullptr)
  {
    instance = std::unique_ptr<GaussPointCache>(new GaussPointCache);
  }
  return *instance;
}


Teuchos::RCP<Core::FE::GaussPoints> Core::FE::GaussPointCache::Create(
    Core::FE::CellType distype, int degree)
{
  std::map<std::pair<Core::FE::CellType, int>, Teuchos::RCP<GaussPoints>>::iterator i =
      gp_cache_.find(std::make_pair(distype, degree));
  if (i != gp_cache_.end())
  {
    return i->second;
  }

  // this is expensive and should not be done too often
  Teuchos::RCP<GaussPoints> gp;

  switch (distype)
  {
    case Core::FE::CellType::quad4:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Quadrilateral<4>>(degree));
      break;
    case Core::FE::CellType::quad8:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Quadrilateral<8>>(degree));
      break;
    case Core::FE::CellType::quad9:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Quadrilateral<9>>(degree));
      break;
    case Core::FE::CellType::tri3:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Triangle<3>>(degree));
      break;
    case Core::FE::CellType::tri6:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Triangle<6>>(degree));
      break;
    case Core::FE::CellType::hex8:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Hexahedron<8>>(degree));
      break;
    case Core::FE::CellType::hex20:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Hexahedron<20>>(degree));
      break;
    case Core::FE::CellType::hex27:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Hexahedron<27>>(degree));
      break;
    case Core::FE::CellType::tet4:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Tetrahedron<4>>(degree));
      break;
    case Core::FE::CellType::tet10:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Tetrahedron<10>>(degree));
      break;
    case Core::FE::CellType::wedge6:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Wedge<6>>(degree));
      break;
    case Core::FE::CellType::wedge15:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Wedge<15>>(degree));
      break;
    case Core::FE::CellType::pyramid5:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Pyramid<5>>(degree));
      break;
    case Core::FE::CellType::line2:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Line<2>>(degree));
      break;
    case Core::FE::CellType::line3:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Line<3>>(degree));
      break;
    default:
      FOUR_C_THROW("unsupported element shape");
  }

  gp_cache_[std::make_pair(distype, degree)] = gp;
  return gp;
}

FOUR_C_NAMESPACE_CLOSE
