/*!----------------------------------------------------------------------
\file drt_contact_manager_base.H

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
#ifndef DRT_CONTACT_MANAGER_BASE_H
#define DRT_CONTACT_MANAGER_BASE_H

#include "Teuchos_RefCountPtr.hpp"
#include "drt_contact_abstract_strategy.H"
#include "drt_contact_interface.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"

using namespace std;
using namespace Teuchos;

namespace CONTACT
{

// forward declarations

/*!
\brief Main class to control all contact

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class ManagerBase
{
public:

  //! @name Enums and Friends
  //@}

  /*!
  \brief Standard Constructor

  The base class constructor is empty.

  One needs a derived class for a concrete implementation of the contact
  manager into a given FE code environment (see e.g. drt_contact_manager.H
  and drt_contact_manager.cpp for the BACI implementation). This constructor
  then has to be fed with a discretization that is expected to carry at least
  two contact boundary conditions. It extracts all contact boundary conditions
  and constructs one or multiple contact interfaces from them and stores them.

  It also builds the corresponding strategy solver object and stores a reference
  in the strategy_ member variable.
     
  It calls CONTACT::Interface::FillComplete() on all contact interfaces which
  makes the nodes and elements of a contact interfaces redundant on all
  processors that either own a node or an element on the interfaces in the
  input discretization.

  */
  explicit ManagerBase();

  /*!
  \brief Destructor

  */
  virtual ~ManagerBase() {};

  //! @name Access methods

  /*!
  \brief Get Epetra communicator

  */
  inline const Epetra_Comm& Comm() const { return *comm_;}
  
  /*!
  \brief Return the object for the solving strategy.
 
  All necessary steps for the computation algorithm
  have to be specialized in subclasses of AbstractStrategy
 
  */
  inline AbstractStrategy& GetStrategy() { return *strategy_;}


protected:

  // don't want = operator and cctor
  ManagerBase operator = (const ManagerBase& old);
  ManagerBase(const ManagerBase& old);

  
  RCP<Epetra_Comm>               comm_;     // communicator
  RCP<CONTACT::AbstractStrategy> strategy_; // strategy object

}; // class ManagerBase
} // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_MANAGER_BASE_H
#endif  // #ifdef CCADISCRET
