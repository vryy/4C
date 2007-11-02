/*!----------------------------------------------------------------------
\file drt_pbcdofset.H

\brief A modified set of degrees of freedom for periodic boundary
       conditions. This class is inherited from dofset and replaces the
       method AssignDegreesOfFreedom by a version which uses a map
       of coupled nodes provided by the periodic boundary conditions to
       assign the same degrees of freedom to coupled pairs of master and
       slavenodes. It is absolutely mandatory that for each slave node
       on this proc the master is (ghosted) on this proc, too!

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
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

\author gammi
*/
class PBCDofSet: public DofSet
{
public:



  /*!
  \brief Standard Constructor


  <pre>

  create a dofset that allows coupled nodes for periodic boundary
  conditions                                         gammi 05/07
  
  </pre>

  \param couplednodes (i) list of coupled nodes


  \return void
  
  */
  PBCDofSet(RefCountPtr<map<int,vector<int> > >  couplednodes);

  /*!
  \brief Destructor

  */
  virtual ~PBCDofSet();



  //! @name Construction


  /*!
  \brief Assign dof numbers to all elements and nodes of the discretization.

  <pre>

  Assign dof numbers to all elements and nodes of the discretization.
  This version of AssignDegreesOfFreedom takes care of the numbering
  of coupled degrees of freedom (for periodic boundary conditions,
  several coupled nodes have the same degrees of freedom)   gammi 05/07
  
  </pre>

  \param  dis   (i) discretisation
  \param  start (i) start numbering of dofs with this index


  \return int highest dof number assigned up to now
  
  */
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const int start);

  //@}

protected:

  //!\brief master and slave node connectivity for periodic boundary conditions
  RefCountPtr<map<int,vector<int> > >           perbndcouples_;   
}; // class PBCDofSet
} // namespace DRT



#endif  // #ifndef PBCDOFSET_H
#endif  // #ifdef CCADISCRET
