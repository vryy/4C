/*-----------------------------------------------------------*/
/*!
\file periodic_boundingbox.H

\brief A class handling a (periodic) bounding box as simulation volume

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef PERIODIC_BOUNDINGBOX_H
#define PERIODIC_BOUNDINGBOX_H

#define BOXOVERLAP  1e-7

#include "../linalg/linalg_fixedsizematrix.H"
#include <Teuchos_RCP.hpp>


// forward declaration
class DiscretizationRuntimeVtuWriter;
class Epetra_Vector;
namespace DRT
{
  class Discretization;
}

namespace GEO
{
namespace MESHFREE
{

/*!
 \brief Construction of periodic boundingbox over entire considered simulation volume
 */
class BoundingBox
{
public:


  BoundingBox();

  virtual ~BoundingBox(){};

  void Init();

  void Setup();

  /// get edge length
  double const& EdgeLength(const int dim) const { return edgelength_[dim]; }

  /// get box
  LINALG::Matrix< 3, 2 > const& Box() const { return box_; }

  /// get flag indicating if periodic boundary conditions are active
  bool HavePBC() const {return havepbc_; }

  /// get const bounding box discretization
  DRT::Discretization const& BoundingBoxDiscret() const { return *boxdiscret_; }

  /// get corner points
  double operator()( int i, int j ) const { return box_(i, j); }

  /// initialize bounding box discretization
  void InitBoundingBoxDiscretization();

  /*!
  \brief shift node (if outside) back in box if periodic boundary conditions
  */
  bool Shift3D( LINALG::Matrix<3,1>& d, LINALG::Matrix<3,1> const X = LINALG::Matrix<3,1>( true ) ) const;

  /*! Check the distance to a reference point position (e.g. node of the
   * same element). If the distance is larger than half of the period
   * length, the point position has been shifted before.
   *
   * Warning: This assumes that the distance between point and reference
   *          point is not larger than half of the period length unless we
   *          shift it. For beam elements, this restricts the element length
   *          to be smaller than this value throughout the entire simulation.
   *          So far, we only check this once in the beginning.
   *
   * Note: this should be equivalent to the previously applied criterion
   *       that the distance between given point and reference point
   *       decreases by either adding or subtracting the period length. */
  bool UnShift3D( LINALG::Matrix<3,1>& d, LINALG::Matrix<3,1> const& ref,
      LINALG::Matrix<3,1> const X = LINALG::Matrix<3,1>( true ) ) const;

  /*!
  \brief get random position inside box
  */
  void RandomPosWithin( LINALG::Matrix< 3, 1 >& pos ) const;

  /*!
   \brief If necessary make the boundingbox larger to include this point as one of the corners of the box
   */
  void AddPoint(const double * x);

  /*!
   \brief Check whether "b" is within this boundingbox
   */
  bool Within( const BoundingBox & b, double norm = 1.0 ) const;

  /*!
   \brief Check the point is within this boundingbox
   */
  bool Within( const double * x, double norm = 1.0 ) const;

  /*!
   \brief Check these points are within this boundingbox
   */
  bool Within( const Epetra_SerialDenseMatrix & xyz, double norm = 1.0 ) const;

  /*!
   \brief Print the corner points of boundingbox on the screen
   */
  void Print();

  double min( int const dim ) const
  {
    return box_( dim, 0 );
  }
  double minx() const
  {
    return box_( 0, 0 );
  }
  double miny() const
  {
    return box_( 1, 0 );
  }
  double minz() const
  {
    return box_( 2, 0 );
  }

  double max( int const dim ) const
  {
    return box_( dim, 1 );
  }
  double maxx() const
  {
    return box_( 0, 1 );
  }
  double maxy() const
  {
    return box_( 1, 1 );
  }
  double maxz() const
  {
    return box_( 2, 1 );
  }

  /*!
   \brief Get the outmost point of the boundingbox
   */
  void UndeformedBoxCornerPointPosition( int i, double * x );
  /*!
   \brief get reference position of corner point i
   */
  LINALG::Matrix< 3, 1 > ReferencePosOfCornerPoint( int i ) const;

  /*!
   \brief get current position of corner point i
   */
  LINALG::Matrix< 3, 1 > CurrentPositionOfCornerPoint( int i ) const;

  /*!
   \brief print box
  */
  void Print( std::ostream& out ) const;

  /*!
   \brief Write output
  */
  void RuntimeOutputStepState( double timen, int stepn ) const;

  /*!
   \brief Apply dirichlet condition accroding to input file
  */
  void ApplyDirichlet( double timen );

  /*!
   \brief init runtime output object for bounding box discretization
  */
  void InitRuntimeOutput();

  //! @name public function dealing with mapping of positions in case of a deforming bounding box
  //! @{

  //! transform from undeformed to global
  void TransformFromUndeformedBoundingBoxSystemToGlobal(
      LINALG::Matrix<3,1> const& xi,
      LINALG::Matrix<3,1>&       x ) const;

  void TransformFromUndeformedBoundingBoxSystemToGlobal(
      double const* xi,
      double *      x ) const;

  //! transform from global to undeformed
  bool TransformFromGlobalToUndeformedBoundingBoxSystem(
      LINALG::Matrix<3,1> const& x,                 ///< input  -> global position
      LINALG::Matrix<3,1>&       xi                 ///< output -> position in undeformed bounding box system
      ) const;
  bool TransformFromGlobalToUndeformedBoundingBoxSystem(
      double const * x,                 ///< input  -> global position
      double *       xi                 ///< output -> position in undeformed bounding box system
      ) const;

  //! @}

protected:
  //! returns init state
  inline const bool& IsInit() const
  { return isinit_; };

  //! returns setup state
  inline const bool& IsSetup() const
  { return issetup_; };

  //! Check the init state
  inline void CheckInit() const
  {
    if (not IsInit())
      dserror("Call Init() first!");
  }

  //! Check the init and setup state
  inline void CheckInitSetup() const
  {
    if (not IsInit() or not IsSetup())
      dserror("Call Init() and Setup() first!");
  }

private:


  /*!
  \brief shift node (if outside) back in box if periodic boundary conditions
  */
  bool Shift1D( int dim, double& d, double const& X = 0.0 ) const;

  /*!
  \brief shift node out of box if it was shifted in previously
  */
  bool UnShift1D( const int dim, double& d, double const& ref, double const& X = 0.0 ) const;


  bool InBetween( double norm, double smin, double smax, double omin, double omax ) const
  {
    double tol = BOXOVERLAP * norm;
    return ( ( omax > smin - tol ) and ( smax > omin - tol ) );
  }

  //! @name private function dealing with mapping of positions in case of a deforming bounding box
  //! @{

  //! evaluate lagrange polynomial that maps from undeformed to global at xi
  void LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobal(
      LINALG::Matrix< 8, 1 >& funct,       ///< to be filled with shape function values
      double                  r,
      double                  s,
      double                  t ) const;

  //! evaluate first derivative of lagrange polynomial that maps from undeformed to global at xi
  void LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobalDeriv1(
      LINALG::Matrix< 3, 8 >& deriv1,      ///< to be filled with shape function derivative values
      double                  r,
      double                  s,
      double                  t ) const;

  //! @}

protected:
  //! @name member variables

  //! indicates if the Init() function has been called
  bool isinit_;

  //! indicates if the Setup() function has been called
  bool issetup_;

private:

  /// discretization with one volume element representing the box ( used e.g. for output)
  Teuchos::RCP<DRT::Discretization> boxdiscret_;
  /// box displacement vector
  Teuchos::RCP<Epetra_Vector> disn_row_;
  Teuchos::RCP<Epetra_Vector> disn_col_;

  bool empty_;
  /// set global flag
  bool havepbc_;
  /// box corners
  LINALG::Matrix<3, 2> box_;
  /// flags for existence of periodic boundary conditions in x, y, z direction
  bool pbconoff_[3];
  ///< box edge lengths in x, y, z direction
  double edgelength_[3];

  //! bounding box discretization runtime vtu writer
  Teuchos::RCP<DiscretizationRuntimeVtuWriter> vtu_writer_ptr_;
};

}
}

#endif
