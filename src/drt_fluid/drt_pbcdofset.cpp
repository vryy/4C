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
  int count = DRT::DofSet::AssignDegreesOfFreedom(dis,dspos,start);

  // loop all master nodes and set the degrees of freedom of
  // the slaves to the degrees of freedom of the master
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

#endif  // #ifdef CCADISCRET
