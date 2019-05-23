/*!----------------------------------------------------------------------
\file drt_utils_gausspoints.cpp

\brief Implementation of gauss points access functions

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/

#include <Intrepid_Cubature.hpp>
#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>
#include <Shards_CellTopology.hpp>

#include <stdexcept>

#include "drt_utils_gausspoints.H"

namespace DRT
{
  namespace UTILS
  {
    namespace
    {
      /// wrapper to intrepid gauss point implementation
      template <class topology>
      class IntrepidGaussPoints : public GaussPoints
      {
       public:
        explicit IntrepidGaussPoints(int cubDegree)
        {
          // cell type: tetrahedron
          shards::CellTopology cellType = shards::getCellTopologyData<topology>();

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

        virtual int NumPoints() const { return cub_points_.dimension(0); }

        virtual int NumDimension() const { return cub_points_.dimension(1); }

        virtual const double* Point(int point) const { return &cub_points_(point, 0); }

        virtual double Weight(int point) const { return cub_weights_(point); }

        virtual void Print() const
        {
          // cell type: tetrahedron
          shards::CellTopology cellType = shards::getCellTopologyData<topology>();

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

        virtual int NumPoints() const { return cub_points_.dimension(0); }

        virtual int NumDimension() const { return cub_points_.dimension(1); }

        virtual const double* Point(int point) const { return &cub_points_(point, 0); }

        virtual double Weight(int point) const { return cub_weights_(point); }

        virtual void Print() const
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
  }    // namespace UTILS
}  // namespace DRT

DRT::UTILS::GaussIntegration::GaussIntegration(DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::quad4:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::quad4, 3);
      break;
    case DRT::Element::quad8:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::quad8, 4);
      break;
    case DRT::Element::quad9:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::quad9, 4);
      break;
    case DRT::Element::tri3:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::tri3, 3);
      break;
    case DRT::Element::tri6:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::tri6, 4);
      break;
    case DRT::Element::hex8:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::hex8, 3);
      break;
    case DRT::Element::hex20:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::hex20, 4);
      break;
    case DRT::Element::hex27:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::hex27, 4);
      break;
    case DRT::Element::tet4:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::tet4, 3);
      break;
    case DRT::Element::tet10:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::tet10, 4);
      break;
    case DRT::Element::wedge6:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::wedge6, 3);
      break;
    case DRT::Element::wedge15:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::wedge15, 4);
      break;
    case DRT::Element::pyramid5:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::pyramid5, 3);
      break;
    case DRT::Element::line2:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::line2, 3);
      break;
    case DRT::Element::line3:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::line3, 4);
      break;
    case DRT::Element::nurbs2:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::line2, 3);
      break;
    case DRT::Element::nurbs3:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::line3, 4);
      break;
    case DRT::Element::nurbs4:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::quad4, 3);
      break;
    case DRT::Element::nurbs8:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::hex8, 3);
      break;
    case DRT::Element::nurbs9:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::quad9, 4);
      break;
    case DRT::Element::nurbs27:
      gp_ = GaussPointCache::Instance().Create(DRT::Element::hex27, 4);
      break;
    default:
      throw std::runtime_error("unsupported element shape");
  }
}

DRT::UTILS::GaussIntegration::GaussIntegration(DRT::Element::DiscretizationType distype, int degree)
{
  gp_ = GaussPointCache::Instance().Create(distype, degree);
}

DRT::UTILS::GaussPointCache* DRT::UTILS::GaussPointCache::instance_;

DRT::UTILS::GaussPointCache& DRT::UTILS::GaussPointCache::Instance()
{
  if (instance_ == NULL)
  {
    instance_ = new GaussPointCache;
  }
  return *instance_;
}

void DRT::UTILS::GaussPointCache::Done()
{
  gp_cache_.clear();
  delete instance_;
  instance_ = NULL;
}

Teuchos::RCP<DRT::UTILS::GaussPoints> DRT::UTILS::GaussPointCache::Create(
    DRT::Element::DiscretizationType distype, int degree)
{
  std::map<std::pair<DRT::Element::DiscretizationType, int>, Teuchos::RCP<GaussPoints>>::iterator
      i = gp_cache_.find(std::make_pair(distype, degree));
  if (i != gp_cache_.end())
  {
    return i->second;
  }

  // this is expensive and should not be done too often
  Teuchos::RCP<GaussPoints> gp;

  switch (distype)
  {
    case DRT::Element::quad4:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Quadrilateral<4>>(degree));
      break;
    case DRT::Element::quad8:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Quadrilateral<8>>(degree));
      break;
    case DRT::Element::quad9:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Quadrilateral<9>>(degree));
      break;
    case DRT::Element::tri3:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Triangle<3>>(degree));
      break;
    case DRT::Element::tri6:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Triangle<6>>(degree));
      break;
    case DRT::Element::hex8:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Hexahedron<8>>(degree));
      break;
    case DRT::Element::hex20:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Hexahedron<20>>(degree));
      break;
    case DRT::Element::hex27:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Hexahedron<27>>(degree));
      break;
    case DRT::Element::tet4:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Tetrahedron<4>>(degree));
      break;
    case DRT::Element::tet10:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Tetrahedron<10>>(degree));
      break;
    case DRT::Element::wedge6:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Wedge<6>>(degree));
      break;
    case DRT::Element::wedge15:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Wedge<15>>(degree));
      break;
    case DRT::Element::pyramid5:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Pyramid<5>>(degree));
      break;
    case DRT::Element::line2:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Line<2>>(degree));
      break;
    case DRT::Element::line3:
      gp = Teuchos::rcp(new IntrepidGaussPoints<shards::Line<3>>(degree));
      break;
    default:
      throw std::runtime_error("unsupported element shape");
  }

  gp_cache_[std::make_pair(distype, degree)] = gp;
  return gp;
}
