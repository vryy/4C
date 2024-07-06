/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_INTEGRATIONCELL_HPP
#define FOUR_C_CUT_INTEGRATIONCELL_HPP

#include "4C_config.hpp"

#include "4C_cut_point.hpp"  // necessary due to enumerator

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class Mesh;
    class VolumeCell;
    class IntegrationCellCreator;

    /*! \class Integration Cell
        \brief Base class for integration cells */
    class IntegrationCell
    {
     public:
      IntegrationCell(Point::PointPosition position, const Core::LinAlg::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : position_(position), xyz_(xyz), points_(points), cell_(cell)
      {
      }

      virtual ~IntegrationCell() = default;
      /// get the dimension of the integration cell
      virtual unsigned n_dim() const = 0;

      /// get the shape of the integration cell
      virtual Core::FE::CellType shape() const = 0;

      virtual int cubature_degree(Core::FE::CellType elementshape) const = 0;

      void dump_gmsh(std::ofstream& file, int* value = nullptr);

      /** \brief calculate the element volume
       *
       *  \note For 1-D elements the element length is returned,
       *        and for 2-D elements the element area.
       *
       *  \author hiermeier \date 11/16 */
      double volume() const;

      const std::vector<Point*>& points() const { return points_; }

      Point::PointPosition position() const { return position_; }

      const Core::LinAlg::SerialDenseMatrix& coordinates() const { return xyz_; }

      bool contains(Core::LinAlg::Matrix<3, 1>& x);

      template <unsigned probdim, Core::FE::CellType celltype>
      bool contains(Core::LinAlg::Matrix<probdim, 1>& x);


      /** Print the integration cells
       *
       *  \author hiermeier \date 02/17 */
      void print(std::ostream& stream) const;
      inline void print() const { print(std::cout); };

     protected:
      Point::PointPosition position_;
      Core::LinAlg::SerialDenseMatrix xyz_;
      std::vector<Point*> points_;
      VolumeCell* cell_;
    };  // class IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// line2 integration cell
    class Line2IntegrationCell : public IntegrationCell
    {
     public:
      Line2IntegrationCell(Point::PointPosition position,
          const Core::LinAlg::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      unsigned n_dim() const override { return Core::FE::dim<Core::FE::CellType::line2>; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::line2; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

    };  // class Line2IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// tri3 integration cell
    class Tri3IntegrationCell : public IntegrationCell
    {
     public:
      Tri3IntegrationCell(Point::PointPosition position, const Core::LinAlg::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      unsigned n_dim() const override { return Core::FE::dim<Core::FE::CellType::tri3>; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::tri3; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

    };  // class Tri3IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// quad4 integration cell
    class Quad4IntegrationCell : public IntegrationCell
    {
     public:
      Quad4IntegrationCell(Point::PointPosition position,
          const Core::LinAlg::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      unsigned n_dim() const override { return Core::FE::dim<Core::FE::CellType::quad4>; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::quad4; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

    };  // class Tri3IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// hex8 integration cell
    class Hex8IntegrationCell : public IntegrationCell
    {
     public:
      Hex8IntegrationCell(Point::PointPosition position, const Core::LinAlg::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned n_dim() const override { return 3; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::hex8; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

      // virtual double Volume() const;

     private:
    };  // class Hex8IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// tet4 integration cell
    class Tet4IntegrationCell : public IntegrationCell
    {
     public:
      Tet4IntegrationCell(Point::PointPosition position, const Core::LinAlg::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned n_dim() const override { return 3; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::tet4; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

      // virtual double Volume() const;
    };

    /*----------------------------------------------------------------------------*/
    /// wedge6 integration cell
    class Wedge6IntegrationCell : public IntegrationCell
    {
     public:
      Wedge6IntegrationCell(Point::PointPosition position,
          const Core::LinAlg::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned n_dim() const override { return 3; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::wedge6; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

      // virtual double Volume() const;
     private:
    };  // class Wedge6IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// pyramid5 integration cell
    class Pyramid5IntegrationCell : public IntegrationCell
    {
     public:
      Pyramid5IntegrationCell(Point::PointPosition position,
          const Core::LinAlg::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned n_dim() const override { return 3; };

      Core::FE::CellType shape() const override { return Core::FE::CellType::pyramid5; }

      int cubature_degree(Core::FE::CellType elementshape) const override;

      // virtual double Volume() const;
     private:
    };  // class Pyramid5IntegrationCell

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
