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


#endif  // #ifdef CCADISCRET
