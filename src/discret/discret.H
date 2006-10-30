/*!----------------------------------------------------------------------
\file discret.H

\class CCADISCRETIZATION::Discretization

\brief a class to manage one discretization

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DISCRET_H
#define DISCRET_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Comm.h"
using namespace std;
using namespace Teuchos;


/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class Discretization

\brief A class to manage a discretization in parallel

The \ref CCADISCRET::Discretization class supports the ostream& operator <<

\author gee (gee@lnm.mw.tum.de)
*/
class Discretization
{
public:
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param comm : An epetra comm object associated with this discretization
  */
  Discretization(RefCountPtr<Epetra_Comm> comm);

  /*!
  \brief Destructor

  */
  virtual ~Discretization();

  //@}

  // @{ \name Acess methods

  /*!
  \brief Get communicator associated with this class
  */
  inline const Epetra_Comm& Comm() const {return *comm_;}


  //@}

private:

  // don't want = operator
  Discretization operator = (const Discretization& old);


  RefCountPtr<Epetra_Comm> comm_;          // an Epetra_comm

}; // class Discretization
} // namespace CCADISCRET











#endif  // #ifndef DISCRET_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
