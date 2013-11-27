/*!----------------------------------------------------------------------
\file drt_meshfree_cell.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

------------------------------------------------------------------------------
</pre>

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
class Cell : public DRT::MESHFREE::MeshfreeBin
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
  Cell(const DRT::MESHFREE::Cell& old);

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

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of knots to this cell
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  inline int NumKnot() const
  {
    return knotid_.size();
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return id's of knots adjacent to this cell
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  inline const int* KnotIds() const
  {
    if (knotid_.size()) return &knotid_[0];
    else                return NULL;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of ptrs to knots
   *
   * \warning The pointers to the knots are build in
   *          DRT::Discretization::FillComplete. A standalone
   *          cell that has not been added to a discretization
   *          (or the discretization has not been called FillComplete)
   *          does not have pointers to knots. In this case, the method
   *          returns NULL.
   *
   * \return Ptr to pointers to knots of the cell in local nodal ordering.
   *         Returns NULL if pointers to not exist.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual inline MeshfreeNode** Knots()
  {
    if (knot_.size()) return &knot_[0];
    else              return NULL;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get const vector of ptrs to knots
   *
   * \warning The pointers to the knots are build in
   *          DRT::Discretization::FillComplete. A standalone
   *          cell that has not been added to a discretization
   *          (or the discretization has not been called FillComplete)
   *          does not have pointers to knots. In this case, the method
   *          returns NULL.
   *
   * \return Ptr to pointers to knots of the cell in local nodal ordering.
   *         Returns NULL if pointers to not exist.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  virtual const MeshfreeNode*const* Knots() const
  {
    if (knot_.size()) return (const MeshfreeNode*const*)(&knot_[0]);
    else              return NULL;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Set a list of knot ids this cell is connected to
   *
   * Sets the nodal ids of the knots adjacent to this cell and the number
   * of knots. This method is used in the construction phase of a
   * discretization. It allows, that cells and knots are created separately
   * and be combined later on the way.
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  void SetKnotIds(const int nknot, //!< number of knots to be set
                  const int* knots //!< list of unique global nodal ids
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Set a list of knot ids this cell is connected to
   *
   * Here the knot ids are directly taken from an input line.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  void SetKnotIds(const std::string& distype,         //!< discretisation type of dat-file defined cell type
                  DRT::INPUT::LineDefinition* linedef //!< line definition of respective cell type
    );

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
   * \brief Return the parameters determining the meshfree approximation
   *
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  Teuchos::ParameterList const * const Params() const { return params_.get(); };

  /*------------------------------------------------------------------------*/
  //! @name Public methods to be used by DRT::Discretization only
  /*------------------------------------------------------------------------*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointer vector from map of knots
   *
   * \warning (public, but to be used by DRT::Discretization ONLY!)
   *
   * The method is used to build the variable knot_ in this cell. It is
   * called from DRT::MESHFREE::Discretization in
   * DRT::MESHFREE::Discretization::FillComplete() to
   * create the pointers from cells to knots (and knots to cells)
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/

  virtual bool BuildKnotPointers(std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >& knots //!< (in): A map of all knots of a discretization
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Build pointer vector from vector of nodal pointers
   *
   * \warning (public, but to be used by DRT::Discretization ONLY!)
   *
   * The method is used to build the variable knot_ in this cell. It can be
   * used to explicitly pass the nodal pointers to the cell.
   *
   * The array of pointers is implicitly expected to be of length NumKnot()
   * and contain pointers to knots in the correct cell local ordering scheme.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/

  virtual bool BuildKnotPointers(MeshfreeNode** knots //!< (in): Pointer to array of pointers to knots.
    );

  /*------------------------------------------------------------------------*/
  //@}
  /*------------------------------------------------------------------------*/


private:

  //! \brief List of my nodal ids, length NumKnot()
  std::vector<int> knotid_;

  //! \brief Pointers to adjacent knots in cell local ordering
  std::vector<MeshfreeNode*> knot_;

  //! meshfree parameter settings
  Teuchos::RCP<Teuchos::ParameterList> params_;

}; // class Cell

} // namespace MESHFREE
} // namespace DRT

//! << operator
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::Cell& ele);

#endif  // #ifndef DRT_MESHFREE_CELL_H
