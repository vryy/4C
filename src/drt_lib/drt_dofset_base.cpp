/*!----------------------------------------------------------------------
\file drt_dofset_base.cpp
\brief A set of degrees of freedom

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Ulrrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iostream>
#include <algorithm>
#include <numeric>

#include "drt_dofset_base.H"
#include "drt_discret.H"
#include "drt_utils.H"

#include "linalg_utils.H"

#include "../headers/define_sizes.h"

// list of all dof sets
std::list<DRT::DofSetBase*> DRT::DofSetBase::dofsets_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::DofSetBase()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             ukue 04/07|
 *----------------------------------------------------------------------*/
DRT::DofSetBase::~DofSetBase()
{
  dofsets_.remove(this);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetBase::MaxGIDinList()
{
  int count = 0;
  for (std::list<DofSetBase*>::const_iterator i=dofsets_.begin();
       i!=dofsets_.end();
       ++i)
  {
    if (*i==this)
      break;

    // ignore empty (no yet initialized) dof row maps
    if ((*i)->NumGlobalElements()>0)
    {
      count = max((*i)->MaxAllGID() + 1,count);
    }
  }
  return count;
}

#endif  // #ifdef CCADISCRET
