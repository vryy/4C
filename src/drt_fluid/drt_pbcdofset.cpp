/*!----------------------------------------------------------------------
\file drt_pbcdofset.cpp

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

#include "drt_pbcdofset.H"


#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
PBCDofSet::PBCDofSet(RefCountPtr<map<int,vector<int> > >  couplednodes)
{
  perbndcouples_=couplednodes;

  for( map<int,vector<int> >::iterator curr = perbndcouples_->begin();
       curr != perbndcouples_->end();
       ++curr )
  {
    std::vector<int> & sids = curr->second;
    std::copy( sids.begin(), sids.end(), std::inserter( slavenodeids_, slavenodeids_.begin() ) );
  }
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/07|
 *----------------------------------------------------------------------*/
PBCDofSet::~PBCDofSet()
{
  return;
}


int PBCDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // assign dofs for the standard dofset, that is without periodic boundary conditions
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);

  // loop all master nodes and set the dofs of the slaves to the dofs of the master
  // remark: the previously assigned dofs of slave nodes are overwritten here
  for(map<int,vector<int> >::iterator master = perbndcouples_->begin();
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

    for (vector<int>::iterator slave=master->second.begin();
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
 |  update coupled nodes map                              wichmann 07/11|
 *----------------------------------------------------------------------*/
void PBCDofSet::SetCoupledNodes(RefCountPtr<map<int,vector<int> > >  couplednodes)
{
  perbndcouples_=couplednodes;

  for( map<int,vector<int> >::iterator curr = perbndcouples_->begin();
       curr != perbndcouples_->end();
       ++curr )
  {
    std::vector<int> & sids = curr->second;
    std::copy( sids.begin(), sids.end(), std::inserter( slavenodeids_, slavenodeids_.begin() ) );
  }
}

#endif  // #ifdef CCADISCRET
