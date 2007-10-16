/*!----------------------------------------------------------------------
\file drt_cnode.H

\class CONTACT::CNode

\brief A class for a contact node

This class represents a finite element node capable of contact.

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef CNODE_H
#define CNODE_H


#include "../drt_lib/drt_node.H"

using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

/*!
\class CNode

\brief A class for a contact node derived from DRT::Node

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CNode : public DRT::Node
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of CNode
  */
  friend class DRT::Discretization;


  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id     (in): A globally unique node id
  \param coords (in): vector of nodal coordinates, length 3
  \param owner  (in): Owner of this node.
  \param numdof (in): Number of degrees of freedom to this node
  \param dofs   (in): list of global degrees of freedom
  \param isslave(in): flag indicating whether node is slave or master side
  */
  CNode(int id, const double* coords, const int owner, 
        const int numdof, const vector<int>& dofs, const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a CNode

  */
  CNode(const CONTACT::CNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~CNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.
  */
  virtual int UniqueParObjectId() const { return ParObject_CNode;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Print this cnode
  */
  virtual void Print(ostream& os) const ;

  /*!
  \brief Is Node on slave or master side of contact interface
  */
  inline virtual bool IsSlave() const { return isslave_; }

  /*!
  \brief Get no. of degrees of freedom
  */
  inline virtual int NumDof() const { return numdof_; }

  /*!
  \brief Get predefined degrees of freedom
  */
  inline virtual const int* Dofs() const { return &dofs_[0]; }

  //@}

  //! @name Construction


  //@}

protected:

  bool            isslave_;      // true of cnode is on slave side of contact problem
  int             numdof_;       // number of degrees of freedom (length of dofs_)
  vector<int>     dofs_;         // degrees of freedom of this node 
  


}; // class CNode
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CNode& cnode);









#endif  // #ifndef CNODE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
