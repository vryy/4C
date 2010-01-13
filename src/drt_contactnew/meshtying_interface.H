/*!----------------------------------------------------------------------
\file mortar_interface.H

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
#ifndef MESHTYING_INTERFACE_H
#define MESHTYING_INTERFACE_H

#include "meshtying_node.H"
#include "meshtying_element.H"
#include "../drt_mortar/mortar_interface.H"

using namespace std;
using namespace Teuchos;

/*!
\brief CONTACT: namespace of the baci contact and meshtying algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A class for a meshtying interface derived from MORTAR::MortarInterface

\author popp (popp@lnm.mw.tum.de)
*/
class MtInterface : public MORTAR::MortarInterface
{
public:

  /*!
  \brief Standard constructor creating empty meshtying interface
 
  This initializes the employed shape function set for lagrangian mutlipliers
  to a specific setting. Throughout the evaluation process, this set will be employed
  for the field of lagrangian multipliers.

  \param id (in): Unique interface id
  \param comm (in): A communicator object
  \param dim (in): Global problem dimension
  \param icontact (in): Global contact parameter list
 
  */
  explicit MtInterface(const int id, const Epetra_Comm& comm, const int dim,
                       const Teuchos::ParameterList& icontact);
    
  /*!
  \brief Destructor

  */
  virtual ~MtInterface() {}

  /*!
  \brief Print this Interface

  */
  virtual void Print(ostream& os) const;

  //! @name Evlauation methods

  /*!
  \brief Add a CONTACT::MtNode to the interface (Filled()==true NOT prerequisite)

  \param mtnode (in): rcp to a meshtying node

  \return Filled()==false

  */
  virtual void AddMtNode(RCP<CONTACT::MtNode> mtnode)
  { idiscret_->AddNode(mtnode); return; }

  /*!
  \brief Add a CONTACT::MtElement to the interface (Filled()==true is prerequisite)

  \param mtele (in): rcp to a meshtying element

  \return Filled()==false

  */
  virtual void AddMtElement(RCP<CONTACT::MtElement> mtele)
  { idiscret_->AddElement(mtele); return; }

  //@}

protected:

  // don't want = operator and cctor
  MtInterface operator = (const MtInterface& old);
  MtInterface(const MtInterface& old);

  // no new variables as compared to base class

}; // class MtInterface
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::MtInterface& interface);


#endif  // #ifndef MESHTYING_INTERFACE_H
#endif  // #ifdef CCADISCRET
