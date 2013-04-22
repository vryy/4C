/*!----------------------------------------------------------------------
\file drt_meshfree_bin.H

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
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*---------------------------------------------------------------------------*/
#ifndef DRT_MESHFREE_MESHFREEBIN_H
#define DRT_MESHFREE_MESHFREEBIN_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"

namespace DRT
{
namespace MESHFREE
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief A meshfree bin adds the possibility to add and delete nodes
 *        from elements dynamically
 *
 * \author Georg Hammerl (hammerl@lnm.mw.tum.de)
 *
 * \date November, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeBin : public DRT::Element
{
public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  MeshfreeBin(int id,   //!< (in): A globally unique bin id
              int owner //!< (in): owner processor of the meshfree bin
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a meshfree bin
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  MeshfreeBin(const DRT::MESHFREE::MeshfreeBin& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeBin();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy the derived class and return pointer to it
   *
   * This method is sort of a copy constructor for a class derived from
   * DRT::Element. It allows to copy construct the derived class without
   * knowing what it actually is using the base class Element.
   *
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual DRT::Element* Clone() const = 0;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * Every class implementing ParObject needs a unique id defined at the
   * top of drt_parobject.H
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const = 0;

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add a node to the meshfree bin
   *
   * \warning It is your own responsibility to make sure that there will not
   *          be any double entries! (This would be disastrous!)
   *
   * Adds entry at the end of nodeid_ and node_ pointers
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void AddNode(
    DRT::Node* nodeptr //!< (in): pointer to node to be added
  )
  {
    nodeid_.push_back(nodeptr->Id());
    node_.push_back(nodeptr);
    return;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add a node to the meshfree bin
   *
   * \warning It is your own responsibility to make sure that there will not
   *          be any double entries! (This would be disastrous!)
   *
   * Adds entry at the end of nodeid_ and node_ pointers
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void AddNode(
    const int gid,     //!< (in): global id of node to be added
    DRT::Node* nodeptr //!< (in): pointer to node to be added
  )
  {
    nodeid_.push_back(gid);
    node_.push_back(nodeptr);
    return;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete a node from the meshfree bin
   *
   * Searches for position of node with specified gid and deletes entry in
   * vectors node_ and nodeid_
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void DeleteNode(int gid //!< (in): global id of node to be deleted
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete all nodes from the meshfree bin
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual inline void DeleteNodes()
  {
    nodeid_.clear();
    node_.clear();
  }

private:

}; // class MeshfreeBin

} // namespace MESHFREE
} // namespace DRT

#endif  // #ifndef DRT_MESHFREE_MESHFREEBIN_H
