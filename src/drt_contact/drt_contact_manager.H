/*!----------------------------------------------------------------------
\file drt_contact_manager.H

\class CONTACT::Manager

\brief Main class to control all contact

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_MANAGER_H
#define DRT_CONTACT_MANAGER_H

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "drt_contact_interface.H"


using namespace std;
using namespace Teuchos;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace CONTACT
{

// forward declarations

/*!
\class Manager

\brief Main class to control all contact


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class Manager
{
public:

  //! @name Enums and Friends
  //@}

  /*!
  \brief Standard Constructor

  The constructor takes a discretization 
  that is expected to have an even number of contact boundary conditions.
  It extracts all contact boundary conditions and constructs one or multiple
  contact interfaces from them and stores them.
  
  It calls CONTACT::Interface::FillComplete() on all contact interfaces which
  makes the nodes and elements of a contact interfaces redundant on all
  processors that either own a node or an element on the interfaces in the
  input discretization.
  
  \param discret (in): A discretization containing contact boundary conditions
  */
  explicit Manager(DRT::Discretization& discret);


  /*!
  \brief Destructor

  */
  virtual ~Manager() {};


  /*!
  \brief Print this Manager

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get Epetra communicator

  */
  inline const Epetra_Comm& Comm() const { return discret_.Comm();};

  /*!
  \brief Get discretization

  */
  inline DRT::Discretization& Discret() const { return discret_;};
  
  /*!
  \brief Initialize contact for next Newton step
  
  */
  void Initialize();
  
  /*!
  \brief Set current deformation state
  
  */
  void SetState(const string& statename, const RCP<Epetra_Vector> vec);
  
  /*!
  \brief Evaluate contact
    
  */
  void Evaluate();
  

protected:

  // don't want = operator and cctor
  Manager operator = (const Manager& old);
  Manager(const Manager& old);


  DRT::Discretization&              discret_;   // the underlying problem discretization
  vector<RCP<CONTACT::Interface> >  interface_; // contact interfaces


}; // class Manager
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Manager& manager);









#endif  // #ifndef DRT_CONTACT_MANAGER_H
#endif  // #ifdef CCADISCRET
