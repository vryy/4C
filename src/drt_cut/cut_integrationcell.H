/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells

\level 3


*----------------------------------------------------------------------*/

#ifndef CUT_INTEGRATIONCELL_H
#define CUT_INTEGRATIONCELL_H

#include "cut_point.H"  // necessary due to enumerator

namespace GEO
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
      IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : position_(position), xyz_(xyz), points_(points), cell_(cell)
      {
      }

      virtual ~IntegrationCell() {}

      /// get the dimension of the integration cell
      virtual unsigned Dim() const = 0;

      /// get the shape of the integration cell
      virtual DRT::Element::DiscretizationType Shape() const = 0;

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const = 0;

      void DumpGmsh(std::ofstream& file, int* value = NULL);

      /** \brief calculate the element volume
       *
       *  \note For 1-D elements the element length is returned,
       *        and for 2-D elements the element area.
       *
       *  \author hiermeier \date 11/16 */
      double Volume() const;

      const std::vector<Point*>& Points() const { return points_; }

      Point::PointPosition Position() const { return position_; }

      const Epetra_SerialDenseMatrix& Coordinates() const { return xyz_; }

      bool Contains(LINALG::Matrix<3, 1>& x);

      template <unsigned probdim, DRT::Element::DiscretizationType celltype>
      bool Contains(LINALG::Matrix<probdim, 1>& x);


      /** Print the integration cells
       *
       *  \author hiermeier \date 02/17 */
      void Print(std::ostream& stream) const;
      inline void Print() const { Print(std::cout); };

     protected:
      Point::PointPosition position_;
      Epetra_SerialDenseMatrix xyz_;
      std::vector<Point*> points_;
      VolumeCell* cell_;
    };  // class IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// line2 integration cell
    class Line2IntegrationCell : public IntegrationCell
    {
     public:
      Line2IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      virtual unsigned Dim() const { return DRT::UTILS::DisTypeToDim<DRT::Element::line2>::dim; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::line2; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

    };  // class Line2IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// tri3 integration cell
    class Tri3IntegrationCell : public IntegrationCell
    {
     public:
      Tri3IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      virtual unsigned Dim() const { return DRT::UTILS::DisTypeToDim<DRT::Element::tri3>::dim; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tri3; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

    };  // class Tri3IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// quad4 integration cell
    class Quad4IntegrationCell : public IntegrationCell
    {
     public:
      Quad4IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell){/* empty construction */};

      virtual unsigned Dim() const { return DRT::UTILS::DisTypeToDim<DRT::Element::quad4>::dim; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::quad4; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

    };  // class Tri3IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// hex8 integration cell
    class Hex8IntegrationCell : public IntegrationCell
    {
     public:
      Hex8IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      virtual unsigned Dim() const { return 3; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::hex8; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

      // virtual double Volume() const;

     private:
    };  // class Hex8IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// tet4 integration cell
    class Tet4IntegrationCell : public IntegrationCell
    {
     public:
      Tet4IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      virtual unsigned Dim() const { return 3; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tet4; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

      // virtual double Volume() const;
    };

    /*----------------------------------------------------------------------------*/
    /// wedge6 integration cell
    class Wedge6IntegrationCell : public IntegrationCell
    {
     public:
      Wedge6IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      virtual unsigned Dim() const { return 3; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::wedge6; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

      // virtual double Volume() const;
     private:
    };  // class Wedge6IntegrationCell

    /*----------------------------------------------------------------------------*/
    /// pyramid5 integration cell
    class Pyramid5IntegrationCell : public IntegrationCell
    {
     public:
      Pyramid5IntegrationCell(Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz,
          const std::vector<Point*>& points, VolumeCell* cell)
          : IntegrationCell(position, xyz, points, cell)
      {
      }

      virtual unsigned Dim() const { return 3; };

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::pyramid5; }

      virtual int CubatureDegree(DRT::Element::DiscretizationType elementshape) const;

      // virtual double Volume() const;
     private:
    };  // class Pyramid5IntegrationCell

  }  // namespace CUT
}  // namespace GEO

#endif
