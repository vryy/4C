/*!----------------------------------------------------------------------
\file drt_pbcdofset.H

\brief A modified set of degrees of freedom for periodic boundary
       conditions

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef PBCDOFSET_H
#define PBCDOFSET_H

#include <vector>

#include <Epetra_IntVector.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

using namespace std;

#include "drt_dserror.H"
#include "drt_dofset.H"
#include "drt_discret.H"
#include "drt_element.H"
#include "drt_node.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

class Discretization;


/*!
\class DRT::PBCDofSet

\brief A set of degrees of freedom

\note

\author gammi
*/
class PBCDofSet: public DofSet
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor

  */
  PBCDofSet(RefCountPtr<map<int,vector<int> > >  couplednodes);

  /*!
  \brief Destructor

  */
  virtual ~PBCDofSet();

  //@}


  // @{ \name Construction
  
  /// Assign dof numbers to all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const int start);

  //@}

protected:

  // master and slave node connectivity for periodic boundary conditions
  RefCountPtr<map<int,vector<int> > >           perbndcouples_;   
}; // class PBCDofSet
} // namespace DRT



#endif  // #ifndef PBCDOFSET_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
