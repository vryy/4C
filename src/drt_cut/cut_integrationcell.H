/*---------------------------------------------------------------------*/
/*!
\file cut_integrationcell.H

\brief Create and handle integrationcells

\level 3

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_INTEGRATIONCELL_H
#define CUT_INTEGRATIONCELL_H

#include "cut_element.H"
#include "cut_boundingbox.H"

namespace GEO
{
namespace CUT
{
  class Mesh;
  class VolumeCell;
  class IntegrationCellCreator;

/*!
\class Integration Cell
\brief Base class for integrationcells
 */
class IntegrationCell
{
public:

  IntegrationCell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, const std::vector<Point*> & points, VolumeCell * cell )
    : position_( position ),
      xyz_( xyz ),
      points_( points ),
      cell_( cell )
  {}

  virtual ~IntegrationCell() {}

  virtual DRT::Element::DiscretizationType Shape() const = 0;

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const = 0;

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL ) = 0;

  double Volume() const;

  const std::vector<Point*> & Points() const { return points_; }

  Point::PointPosition Position() const { return position_; }

  const Epetra_SerialDenseMatrix & Coordinates() const { return xyz_; }

  bool Contains( LINALG::Matrix<3,1>& x);

  template<DRT::Element::DiscretizationType celltype>
  bool Contains( LINALG::Matrix<3,1>& x);

protected:

  Point::PointPosition position_;
  Epetra_SerialDenseMatrix xyz_;
  std::vector<Point*> points_;
  VolumeCell * cell_;
};

/// hex8 integration cell
class Hex8IntegrationCell : public IntegrationCell
{
public:

  Hex8IntegrationCell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, const std::vector<Point*> & points, VolumeCell * cell )
    : IntegrationCell( position, xyz, points, cell )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::hex8; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const;

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  //virtual double Volume() const;

private:

};

/// tet4 integration cell
class Tet4IntegrationCell : public IntegrationCell
{
public:

  Tet4IntegrationCell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, const std::vector<Point*> & points, VolumeCell * cell )
    : IntegrationCell( position, xyz, points, cell )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tet4; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const;

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  //virtual double Volume() const;

};

/// wedge6 integration cell
class Wedge6IntegrationCell : public IntegrationCell
{
public:

  Wedge6IntegrationCell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, const std::vector<Point*> & points, VolumeCell * cell )
    : IntegrationCell( position, xyz, points, cell )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::wedge6; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const;

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  //virtual double Volume() const;
private:

};

/// pyramid5 integration cell
class Pyramid5IntegrationCell : public IntegrationCell
{
public:

  Pyramid5IntegrationCell( Point::PointPosition position, const Epetra_SerialDenseMatrix & xyz, const std::vector<Point*> & points, VolumeCell * cell )
    : IntegrationCell( position, xyz, points, cell )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::pyramid5; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const;

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  //virtual double Volume() const;
private:

};

}
}

#endif
