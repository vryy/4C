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
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"

namespace DRT
{
namespace MESHFREE
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief The element type class of meshfree bin, creating the same
 *
 * \author Georg Hammerl (hammerl@lnm.mw.tum.de)
 *
 * \date November, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeBinType : public DRT::ElementType
{
public:

  //!< name of specific element type
  std::string Name() const { return "MeshfreeBinType"; }

  //!< returning instance of specific element type
  static MeshfreeBinType & Instance() { return instance_; };

  //!< create parallel object of this element type
  virtual DRT::ParObject* Create( const std::vector<char> & data );

  //!< create element of this element type
  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  //!< create element of this element type
  virtual Teuchos::RCP<DRT::Element> Create( const int id,
                                             const int owner);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ){}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ){}

private:

  //!< static instance of this element type (for self-instantiation)
  static MeshfreeBinType instance_;
};

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
  virtual DRT::Element* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * Every class implementing ParObject needs a unique id defined at the
   * top of drt_parobject.H
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const { return MeshfreeBinType::Instance().UniqueParObjectId(); };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack this class so it can be communicated
   *
   * \ref Pack and \ref Unpack are used to communicate this meshfree bin
   *
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack data from a char vector into this class
   *
   * \ref Pack and \ref Unpack are used to communicate this meshfree bin
   *
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void Unpack(const std::vector<char>& data);

  // return meshfree bin type instance
  virtual DRT::ElementType & ElementType() const { return MeshfreeBinType::Instance(); }

  /*!
  \brief Get number of degrees of freedom of a certain node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

  /*------------------------------------------------------------------------*/
  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)
   *
   * Meshfree bins do not have dofs
   *///
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerElement() const { return 0; }

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

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this meshfree bin
   *
   * Prints basic information about this bin to ostream. This method would
   * usually be called by the print method of a derived class.
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void Print(ostream& os) const;

  /*!
  \brief Evaluate a Neumann boundary condition dummy

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL)
  { return 0; }

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const
  {
      return DRT::Element::dis_none;
  };

private:

}; // class MeshfreeBin

} // namespace MESHFREE
} // namespace DRT

//! << operator
ostream& operator << (ostream& os, const DRT::MESHFREE::MeshfreeBin& ele);

#endif  // #ifndef DRT_MESHFREE_MESHFREEBIN_H
