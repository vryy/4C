/*!----------------------------------------------------------------------
\file drt_contact_interface.H

\class CONTACT::Interface

\brief One contact interface

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_CONTACT_INTERFACE_H
#define DRT_CONTACT_INTERFACE_H

#include "Epetra_Comm.h"
#include "../drt_lib/drt_discret.H"
#include "drt_cnode.H"
#include "drt_celement.H"


using namespace std;
using namespace Teuchos;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\class Interface

\brief One contact interface


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class Interface
{
public:

  //! @name Enums and Friends
  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  explicit Interface(const int id, const Epetra_Comm& comm);


  /*!
  \brief Destructor

  */
  virtual ~Interface() {}


  //@}

  /*!
  \brief Print this Interface

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return comm_;}

  /*!
  \brief check whether interface was called FillComplete

  */
  virtual inline bool Filled() const { return idiscret_->Filled();}

  /*!
  \brief Add a CONTACT::CNode to the interface

  */
  inline virtual void AddCNode(RCP<CONTACT::CNode> cnode) 
  { idiscret_->AddNode(cnode); return; }

  /*!
  \brief Add a CONTACT::CElement to the interface

  */
  inline virtual void AddCElement(RCP<CONTACT::CElement> cele) 
  { idiscret_->AddElement(cele); return; }

  /*!
  \brief Finalize construction of contact interface

  */
  virtual void FillComplete();

protected:

  // don't want = operator and cctor
  Interface operator = (const Interface& old);
  Interface(const Interface& old);


  int                           id_;       // unique interface id
  const Epetra_Comm&            comm_;     // communicator
  RCP<DRT::Discretization>      idiscret_; // the discretization of the contact interface


}; // class Interface
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Interface& interface);









#endif  // #ifndef DRT_CONTACT_INTERFACE_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
