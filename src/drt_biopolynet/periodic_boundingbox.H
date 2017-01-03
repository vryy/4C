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

  double const& EdgeLength(const int dim) { return edgelength_[dim]; };

  LINALG::Matrix<3,2> const& Box(){return box_;}

  double operator()(int i, int j) const
  {
    return box_(i, j);
  }


  /*!
  \brief shift node (if outside) back in box if periodic boundary conditions
  */
  void Shift1D( const int dim, double& d, double const& X = 0.0 ) const;
  /*!
  \brief shift node (if outside) back in box if periodic boundary conditions
  */
  void Shift3D( LINALG::Matrix<3,1>& d, LINALG::Matrix<3,1> const X = LINALG::Matrix<3,1>( true ) ) const;

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
  void UnShift1D( const int dim, double& d, double const& ref, double const& X = 0.0) const;
  /*!
  \brief shift node out of box if it was shifted in previously
  */
  void UnShift3D( LINALG::Matrix<3,1>& d, LINALG::Matrix<3,1> const& ref,
      LINALG::Matrix<3,1> const X = LINALG::Matrix<3,1>( true ) ) const;

  /*!
   \brief If necessary make the boundingbox larger to include this point as one of the corners of the box
   */
  void AddPoint(const double * x);

  /*!
   \brief Check whether "b" is within this boundingbox
   */
  bool Within(const BoundingBox & b, double norm = 1.0) const;

  /*!
   \brief Check the point is within this boundingbox
   */
  bool Within(const double * x, double norm = 1.0) const;

  /*!
   \brief Check these points are within this boundingbox
   */
  bool Within(const Epetra_SerialDenseMatrix & xyz, double norm = 1.0) const;

  /*!
   \brief Print the corner points of boundingbox on the screen
   */
  void Print();

  double min( int const dim ) const
  {
    return box_(dim, 0);
  }
  double minx() const
  {
    return box_(0, 0);
  }
  double miny() const
  {
    return box_(1, 0);
  }
  double minz() const
  {
    return box_(2, 0);
  }

  double max( int const dim ) const
  {
    return box_(dim, 1);
  }
  double maxx() const
  {
    return box_(0, 1);
  }
  double maxy() const
  {
    return box_(1, 1);
  }
  double maxz() const
  {
    return box_(2, 1);
  }

  /*!
   \brief Get the outmost point of the boundingbox
   */
  void CornerPoint(int i, double * x);

  /*!
   \brief Get the boundingbox
   */
  LINALG::Matrix<3, 2> GetBoundingVolume()
  {
    return box_;
  }

  void Print(std::ostream& out) const;

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

protected:
  //! @name member variables

  //! indicates if the Init() function has been called
  bool isinit_;

  //! indicates if the Setup() function has been called
  bool issetup_;

private:

  bool InBetween(double norm, double smin, double smax, double omin,
      double omax) const
  {
    double tol = BOXOVERLAP * norm;
    return ((omax > smin - tol) and (smax > omin - tol));
  }

  bool empty_;
  // set global flag
  bool havepbc_;
  /// box corners
  LINALG::Matrix<3, 2> box_;
  /// flags for existence of periodic boundary conditions in x, y, z direction
  bool pbconoff_[3];
  ///< box edge lengths in x, y, z direction
  double edgelength_[3];
};

}
}

#endif
