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

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_utils.H"

#include "../headers/define_sizes.h"

/*!
\brief A set of degrees of freedom

\author gammi
*/
class PBCDofSet: public DRT::DofSet
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


  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new PBCDofSet(*this)); }

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start);

protected:

  /// get number of nodal dofs for this element at this node
  virtual int NumDofPerNode( const DRT::Node & node, unsigned dspos ) const
  {
    if ( slavenodeids_.count( node.Id() )==0 )
    {
      return DRT::DofSet::NumDofPerNode( node, dspos );
    }
    return 0;
  }

private:

  //!\brief master and slave node connectivity for periodic boundary conditions
  RefCountPtr<map<int,vector<int> > >           perbndcouples_;

  std::set<int> slavenodeids_;
}; // class PBCDofSet



#endif  // #ifndef PBCDOFSET_H
#endif  // #ifdef CCADISCRET
