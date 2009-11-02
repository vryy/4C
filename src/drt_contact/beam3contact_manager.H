/*!----------------------------------------------------------------------
\file beam3contact_manager.H

<pre>
Maintainer: Alexander Popp, Christian Cyron
            {popp,cyron}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
// administrative line in order to make a difference between BACI and Ccarat files
#ifdef CCADISCRET

//compile only if beam3 element is complied, too, as beam3 element required for member variables of this class
#ifdef D_BEAM3

//header file only included if not yet included!
#ifndef BEAM3C_MANAGER_H
#define BEAM3C_MANAGER_H

#include "../drt_beam3/beam3.H"
#include "beam3contact.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

using namespace std;
using namespace Teuchos;

namespace CONTACT
{

/*!
\brief contact manager for contact in a beam3 discretization

\author popp (popp@lnm.mw.tum.de), cyron (cyron@lnm.mw.tum.de)
*/
class Beam3cmanager 
{
public:

  //! @name Friends
  
  //no fried classes defined
  
  //@}


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param discret (in): A discretization containing beam elements
  */
  Beam3cmanager(const DRT::Discretization& discret);

  /*!
  \brief Destructor
  */
  virtual ~Beam3cmanager() {};
  
  //@}

  //! @name public evaluation methods

  /*!
  \brief Evaluate beam contact

  This method calls each beam contact pair to compute its contact forces
  and stiffness. Then all entries are assembles into global force resdiual
  and global stiffness matrix.

  */
  virtual int Evaluate();
  
  /*!
  \brief Search for beam contact
  
  The discretization is searched for potential beam element pairs coming
  into contact. For each pair, a temporary Beam3contact class is generated
  which handles penalty force and stiffness computation.
  
  */
  virtual int SearchContact();
  
  //@}

private:
  
  // don't want = operator and cctor
  Beam3cmanager operator = (const Beam3cmanager& old);
  Beam3cmanager(const Beam3cmanager& old);
    
  //! @name member variables

  //! discretizaton
  const DRT::Discretization& discret_;
  
  //! current contact element pairs
  vector<RCP<CONTACT::Beam3contact> > cpairs_;
 
  //@}

}; // class Beam3cmanager
} // namespace CONTACT

#endif  // #ifndef BEAM3C_MANAGER_H
#endif  // #ifdef D_BEAM3
#endif  // #ifdef CCADISCRET 
