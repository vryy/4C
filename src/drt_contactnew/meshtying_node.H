/*!----------------------------------------------------------------------
\file meshtying_node.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MESHTYING_NODE_H
#define MESHTYING_NODE_H


#include "../drt_mortar/mortar_node.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

using namespace std;
using namespace Teuchos;


/*!
\brief CONTACT: namespace of the baci contact and meshtying algorithms

*/
namespace CONTACT
{

/*!
\brief A class for a meshtying node derived from DRT::Node

This class represents a finite element node capable of meshtying.

\author popp (popp@lnm.mw.tum.de)
*/
class MtNode : public MORTAR::MortarNode
{
public:

  //! @name Enums and Friends

  /*!
  \brief The Discretization is a friend of MtNode
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
  \param isslave(in): flag indicating whether node is slave or master

  */
  MtNode(int id, const double* coords, const int owner,
         const int numdof, const vector<int>& dofs, const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a MtNode

  */
  MtNode(const CONTACT::MtNode& old);

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::MtNode* Clone() const;

  /*!
  \brief Destructor

  */
  virtual ~MtNode() { return; }

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_lib/drt_parobject.H.

  */
  virtual int UniqueParObjectId() const { return ParObject_MtNode;}

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
  \brief Print this MtNode
  */
  virtual void Print(ostream& os) const ;

  //@}

  //! @name Evaluation methods

  //@}

protected:

  // no new variables as compared to base class

}; // class MtNode
} // namespace CONTACT

// << operator
ostream& operator << (ostream& os, const CONTACT::MtNode& mtnode);

#endif  // #ifndef MESHTYING_NODE_H
#endif  // #ifdef CCADISCRET
