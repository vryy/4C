/*!----------------------------------------------------------------------
\file drt_dofset_pbc.H

\brief A modified set of degrees of freedom for periodic boundary
       conditions. This class is inherited from dofset and replaces the
       method AssignDegreesOfFreedom by a version which uses a map
       of coupled nodes provided by the periodic boundary conditions to
       assign the same degrees of freedom to coupled pairs of master and
       slavenodes. It is absolutely mandatory that for each slave node
       on this proc the master is (ghosted) on this proc, too!

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef PBCDOFSET_H
#define PBCDOFSET_H

#include <vector>
#include <map>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dofset.H"



// forward declarations
namespace DRT
{
class Discretization;


/*!
\brief A set of degrees of freedom

\author gammi
*/
class PBCDofSet: virtual public DRT::DofSet
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
  PBCDofSet(Teuchos::RCP<std::map<int,std::vector<int> > >  couplednodes);

  /*!
  \brief Destructor

  */
  virtual ~PBCDofSet();

  /// Get maximum GID of degree of freedom row map
  virtual int MaxAllGID() const;

  /// Get minimum GID of degree of freedom row map
  virtual int MinAllGID() const;

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new PBCDofSet(*this)); }

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start);

  /// Update the coupled nodes map of dofset
  virtual void SetCoupledNodes(Teuchos::RCP<std::map<int,std::vector<int> > >  couplednodes);

  /// Get coupled nodes map (corresponding col format)
  virtual Teuchos::RCP<std::map<int,std::vector<int> > > GetCoupledNodes(){return perbndcouples_;};

  /// Get connectivity map between slave node and its master node
  virtual Teuchos::RCP<std::map<int,int> > GetSlaveToMasterNodeConnectivity(){return perbnd_slavetomaster_;};

protected:

  /// get number of nodal dofs for this element at this node
  virtual int NumDofPerNode( const DRT::Node & node ) const
  {
    if ( slavenodeids_->count( node.Id() )==0 )
    {
      return DRT::DofSet::NumDofPerNode( node );
    }
    return 0;
  }

  //!\brief master and slave node connectivity for periodic boundary conditions
  Teuchos::RCP<std::map<int,std::vector<int> > >           perbndcouples_;

  //!\brief the largest original GID, to stop the dofset from 'shrinking'
  int myMaxGID_;

  //!\brief the smallest original GID, to stop the dofset from 'shrinking'
  int myMinGID_;

private:

  /// Build the connectivity between slave node and its master node
  void BuildSlaveToMasterNodeConnectivity();

  Teuchos::RCP<std::set<int> > slavenodeids_;

  //!\brief slave node to master node connectivity for periodic boundary conditions (key=slave nid, value=master nid)
  Teuchos::RCP<std::map<int,int> >  perbnd_slavetomaster_;

}; // class PBCDofSet

}

#endif  // #ifndef PBCDOFSET_H
