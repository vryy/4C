/*!----------------------------------------------------------------------
\file drt_meshfree_cell.H

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/
#ifndef DRT_MESHFREE_CELL_H
#define DRT_MESHFREE_CELL_H

#include "../drt_meshfree_discret/drt_meshfree_bin.H"
#include "../drt_lib/drt_node.H"

/*==========================================================================*/
// forward declarations
/*==========================================================================*/
namespace DRT
{
//class Node;

namespace MESHFREE
{
class MeshfreeNode;

/*--------------------------------------------------------------------------*/
/*!
 * \brief A virtual class all meshfree cells that are used in DRT::MESHFREE
 *        have to implement
 *
 * This is the pure virtual base class for all meshfree cells to be used
 * with the DRT::MESHFREE::Discretization. Every element (and boundary
 * condition) to be used with the discretization management module has to
 * implement this class. It implements various basic cell methods and
 * stores basic information such as cell to node connectivity.
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2012
 */
/*--------------------------------------------------------------------------*/
template <typename ELEMENT>
class Cell : public DRT::MESHFREE::MeshfreeBin<ELEMENT>
{
public:

  /*========================================================================*/
  //! @name Enums and Friends
  /*========================================================================*/

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  Cell(int id,   //!< (in): A globally unique cell id
       int owner //!< (in): owner processor of the meshfree cell
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a Cell
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  Cell(const DRT::MESHFREE::Cell<ELEMENT>& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual ~Cell();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy the derived class and return pointer to it
   *
   * This method is sort of a copy constructor for a class derived from
   * DRT::Element. It allows to copy construct the derived class without
   * knowing what it actually is using the base class Element.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual DRT::Element* Clone() const = 0;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * Every class implementing ParObject needs a unique id defined at the
   * top of drt_parobject.H
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const = 0;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack this class so it can be communicated
   *
   * \ref Pack and \ref Unpack are used to communicate this cell
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack data from a char vector into this class
   *
   * \ref Pack and \ref Unpack are used to communicate this cell
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Unpack(const std::vector<char>& data);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return the parameters determining the meshfree approximation
   *
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  Teuchos::ParameterList const * const Params() const { return params_.get(); };

  /*========================================================================*/
  //! @name Point methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of points forming the geometry of this cell
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual inline int NumPoint() const { return pointid_.size(); }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return id's of points forming the geometry of this element
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual const int* PointIds() const
  { if (pointid_.size()) return &pointid_[0]; else                return NULL; }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of ptrs to points forming the geometry of this element
   *
   * \warning The pointers to the points are build in
   *          DRT::MESHFREE::MeshfreeDiscretization::FillComplete. A
   *          standalone cell that has not been added to a discretization
   *          (or the discretization has not been called FillComplete)
   *          does not have pointers to points. In this case, the method
   *          returns NULL.
   *
   * \return Ptr to pointers to points of the cell in local nodal ordering.
   *         Returns NULL if pointers to not exist.
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual DRT::Node** Points()
  { if (point_.size()) return (DRT::Node**)(&(point_[0])); else return NULL; }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get const vector of ptrs to points forming the geometry of this element
   *
   * \warning The pointers to the points are build in
   *          DRT::MESHFREE::MeshfreeDiscretization::FillComplete. A
   *          standalone cell that has not been added to a discretization
   *          (or the discretization has not been called FillComplete)
   *          does not have pointers to points. In this case, the method
   *          returns NULL.
   *
   * \return Ptr to pointers to points of the cell in local nodal ordering.
   *         Returns NULL if pointers to not exist.
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual const DRT::Node*const* Points() const
  { if (point_.size()) return (const DRT::Node*const*)(&point_[0]); else return NULL; }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Set a list of point ids this cell is connected to
   *
   * Sets the nodal ids of the points adjacent to this cell and the number
   * of points. This method is used in the construction phase of a
   * discretization. It allows, that cells and points are created separately
   * and be combined later on the way.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  void SetPointIds(
    const int npoint, //!< number of points to be set
    const int* point //!< list of unique global nodal ids
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Set a list of point ids this cell is connected to
   *
   * Here the point ids are directly taken from an input line.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  void SetPointIds(
    const std::string& distype,         //!< discretisation type of dat-file defined cell type
    DRT::INPUT::LineDefinition* linedef //!< line definition of respective cell type
    );

  /*------------------------------------------------------------------------*/
  //! @name Others
  /*------------------------------------------------------------------------*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this cell
   *
   * Prints basic information about this cell to std::ostream. This method would
   * usually be called by the print method of a derived class.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual void Print(std::ostream& os) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Coordinates of cell center computed by point position
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  void CenterAndMaxRadius(LINALG::Matrix<3,1>&  center, double& max_radius);

  /*------------------------------------------------------------------------*/
  //! @name Public methods to be used by DRT::Discretization only
  /*------------------------------------------------------------------------*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointer vector from map of points
   *
   * \warning (public, but to be used by DRT::Discretization ONLY!)
   *
   * The method is used to build the variable point_ in this cell. It is
   * called from DRT::MESHFREE::MeshfreeDiscretization in
   * DRT::MESHFREE::MeshfreeDiscretization::FillComplete() to
   * create the pointers from cells to points (and points to cells)
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool BuildPointPointers(
    std::map<int,Teuchos::RCP<DRT::Node> >& points //!< (in): A map of all points of a discretization
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointer vector from vector of nodal pointers
   *
   * \warning (public, but to be used by DRT::Discretization ONLY!)
   *
   * The method is used to build the variable point_ in this cell. It can be
   * used to explicitly pass the nodal pointers to the cell.
   *
   * The array of pointers is implicitly expected to be of length NumPoint()
   * and contain pointers to points in the correct cell local ordering scheme.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual bool BuildPointPointers(
    DRT::Node** points //!< (in): Pointer to array of pointers to points.
    );

  /*------------------------------------------------------------------------*/
  //@}
  /*------------------------------------------------------------------------*/


private:

  //! \brief List of my point ids, length NumPoint()
  std::vector<int> pointid_;

  //! \brief Pointers to points adjacent to cell in local ordering
  std::vector<DRT::Node*> point_;

  //! meshfree parameter settings
  Teuchos::RCP<Teuchos::ParameterList> params_;

}; // class Cell

} // namespace MESHFREE
} // namespace DRT

//! << operator
template <typename ELEMENT>
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::Cell<ELEMENT>& ele);

#endif  // #ifndef DRT_MESHFREE_CELL_H
