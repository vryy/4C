/*!----------------------------------------------------------------------
\file drt_dofset_pbc.cpp

\brief A modified set of degrees of freedom for periodic boundary
       conditions

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_dofset_pbc.H"


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
DRT::PBCDofSet::PBCDofSet(Teuchos::RCP<std::map<int,std::vector<int> > >  couplednodes)
  :DofSet(), perbndcouples_(Teuchos::null), myMaxGID_(-1)
{
  SetCoupledNodes(couplednodes);
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
DRT::PBCDofSet::~PBCDofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::PBCDofSet::MaxAllGID() const
{
  return myMaxGID_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::PBCDofSet::MinAllGID() const
{
  return myMinGID_;
}


int DRT::PBCDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // temporarily store the slave node set
  Teuchos::RCP<std::set<int> > tempset = slavenodeids_;
  slavenodeids_ = Teuchos::rcp(new std::set<int>);

  // assign dofs using the empty slave node set. This way the dofrowmap_
  // contains exactly the entries as in a regular dofset
  DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);
  if (pccdofhandling_) dserror("ERROR: Point coupling cinditions not yet implemented for PBCDofSet");

  myMaxGID_ = DRT::DofSet::MaxAllGID();
  myMinGID_ = DRT::DofSet::MinAllGID();

  // restore the slave node set
  slavenodeids_ = tempset;

  // assign dofs for the standard dofset, that is without periodic boundary
  // conditions and with the slave node set back in place
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);


  // loop all master nodes and set the dofs of the slaves to the dofs of the master
  // remark: the previously assigned dofs of slave nodes are overwritten here
  for(std::map<int,std::vector<int> >::iterator master = perbndcouples_->begin();
      master != perbndcouples_->end();
      ++master )
  {
    int master_lid=dis.NodeColMap()->LID(master->first);

    if (master_lid<0)
    {
      dserror("master gid %d not on proc %d, but required by slave %d",
              master->first,
              dis.Comm().MyPID(),
              master->second[0]);
    }

    for (std::vector<int>::iterator slave=master->second.begin();
         slave!=master->second.end();
         ++slave)
    {
      int slave_lid=dis.NodeColMap()->LID(*slave);

      if (slave_lid>-1)
      {
        (*numdfcolnodes_)[slave_lid] = (*numdfcolnodes_)[master_lid];
        (*idxcolnodes_)  [slave_lid] = (*idxcolnodes_)  [master_lid];
      }
      else
      {
#ifdef DEBUG
        if(dis.NodeRowMap()->MyGID(master->first))
        {
          dserror("slave not on proc but master owned by proc\n");
        }
#endif
      }
    }
  }

  return count;
}


/*----------------------------------------------------------------------*
 |  update coupled nodes map                             rasthofer 07/11|
 |                                                       DA wichmann    |
 *----------------------------------------------------------------------*/
void DRT::PBCDofSet::SetCoupledNodes(Teuchos::RCP<std::map<int,std::vector<int> > >  couplednodes)
{
  perbndcouples_=couplednodes;
  slavenodeids_ = Teuchos::rcp(new std::set<int>);

  for( std::map<int,std::vector<int> >::iterator curr = perbndcouples_->begin();
       curr != perbndcouples_->end();
       ++curr )
  {
    std::vector<int> & sids = curr->second;
    std::copy( sids.begin(), sids.end(), std::inserter( *slavenodeids_, slavenodeids_->begin() ) );
  }

  /// Build the connectivity between slave node and its master node
  BuildSlaveToMasterNodeConnectivity();

  return;
}

/*----------------------------------------------------------------------*
 |  Build the connectivity between slave node and its master node       |
 |                                                       schott 05/15   |
 *----------------------------------------------------------------------*/
void DRT::PBCDofSet::BuildSlaveToMasterNodeConnectivity()
{
  perbnd_slavetomaster_ = Teuchos::rcp(new std::map<int, int>);

  for(std::map<int,std::vector<int> >::const_iterator masterslavepair = perbndcouples_->begin();
      masterslavepair != perbndcouples_->end() ; ++masterslavepair)
  {
    // loop slave nodes associated with master
    for(std::vector<int>::const_iterator iter=masterslavepair->second.begin(); iter!=masterslavepair->second.end(); ++iter)
    {
      const int slavegid = *iter;
      (*perbnd_slavetomaster_)[slavegid] = masterslavepair->first;
    }
  }
}
