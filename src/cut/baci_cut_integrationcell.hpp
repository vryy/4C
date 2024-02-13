/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells

\level 3


*----------------------------------------------------------------------*/

#ifndef BACI_CUT_INTEGRATIONCELL_HPP
#define BACI_CUT_INTEGRATIONCELL_HPP

#include "baci_config.hpp"

#include "baci_cut_point.hpp"  // necessary due to enumerator

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class Mesh;
    class VolumeCell;
    class IntegrationCellCreator;

    /*! \class Integration Cell
        \brief Base class for integration cells */
    class IntegrationCell
    {
     public:
      IntegrationCell(Point::PointPosition position, const CORE::LINALG::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : position_(position), xyz_(xyz), points_(points), cell_(cell)
      {
      }

      virtual ~IntegrationCell() = default;
      /// get the dimension of the integration cell
      virtual unsigned Dim() const = 0;

      /// get the shape of the integration cell
      virtual CORE::FE::CellType Shape() const = 0;

      virtual int CubatureDegree(CORE::FE::CellType elementshape) const = 0;

      void DumpGmsh(std::ofstream& file, int* value = nullptr);

      /** \brief calculate the element volume
       *
       *  \note For 1-D elements the element length is returned,
       *        and for 2-D elements the element area.
       *
       *  \author hiermeier \date 11/16 */
      double Volume() const;

      const std::vector<Point*>& Points() const { return points_; }

      Point::PointPosition Position() const { return position_; }

      const CORE::LINALG::SerialDenseMatrix& Coordinates() const { return xyz_; }

      bool Contains(CORE::LINALG::Matrix<3, 1>& x);

      template <unsigned probdim, CORE::FE::CellType celltype>
      bool Contains(CORE::LINALG::Matrix<probdim, 1>& x);


      /** Print the integration cells
       *
       *  \author hiermeier \date 02/17 */
      void Print(std::ostream& stream) const;
      inline void Print() const { Print(std::cout); };

     protected:
      Point::PointPosition position_;
      CORE::LINALG::SerialDenseMatrix xyz_;
      std::vector<Point*> points_;
      VolumeCell* cell_;
    };  // class IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// line2 integration cell
    class Line2IntegrationCell : public IntegrationCell
    {
     public:
      Line2IntegrationCell(Point::PointPosition position,
          const CORE::LINALG::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      unsigned Dim() const override { return CORE::FE::dim<CORE::FE::CellType::line2>; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::line2; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

    };  // class Line2IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// tri3 integration cell
    class Tri3IntegrationCell : public IntegrationCell
    {
     public:
      Tri3IntegrationCell(Point::PointPosition position, const CORE::LINALG::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      unsigned Dim() const override { return CORE::FE::dim<CORE::FE::CellType::tri3>; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::tri3; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

    };  // class Tri3IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// quad4 integration cell
    class Quad4IntegrationCell : public IntegrationCell
    {
     public:
      Quad4IntegrationCell(Point::PointPosition position,
          const CORE::LINALG::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      unsigned Dim() const override { return CORE::FE::dim<CORE::FE::CellType::quad4>; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::quad4; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

    };  // class Tri3IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// hex8 integration cell
    class Hex8IntegrationCell : public IntegrationCell
    {
     public:
      Hex8IntegrationCell(Point::PointPosition position, const CORE::LINALG::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned Dim() const override { return 3; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::hex8; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

      // virtual double Volume() const;

     private:
    };  // class Hex8IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// tet4 integration cell
    class Tet4IntegrationCell : public IntegrationCell
    {
     public:
      Tet4IntegrationCell(Point::PointPosition position, const CORE::LINALG::SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned Dim() const override { return 3; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::tet4; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

      // virtual double Volume() const;
    };

    /*----------------------------------------------------------------------------*/
    /// wedge6 integration cell
    class Wedge6IntegrationCell : public IntegrationCell
    {
     public:
      Wedge6IntegrationCell(Point::PointPosition position,
          const CORE::LINALG::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned Dim() const override { return 3; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::wedge6; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

      // virtual double Volume() const;
     private:
    };  // class Wedge6IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// pyramid5 integration cell
    class Pyramid5IntegrationCell : public IntegrationCell
    {
     public:
      Pyramid5IntegrationCell(Point::PointPosition position,
          const CORE::LINALG::SerialDenseMatrix& xyz, const std::vector<Point*>& points,
          VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      unsigned Dim() const override { return 3; };

      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::pyramid5; }

      int CubatureDegree(CORE::FE::CellType elementshape) const override;

      // virtual double Volume() const;
     private:
    };  // class Pyramid5IntegrationCell

  }  // namespace CUT
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif
