/*----------------------------------------------------------------------*/
/*!
\file post_drt_gid.H

\brief GiD filter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

*/
/*----------------------------------------------------------------------*/
#ifndef POST_DRT_GID_H
#define POST_DRT_GID_H



#include "../post_drt_common/post_drt_common.H"


namespace PostGid
{
  //! Runs the Gid filter
  void runGidFilter(PostProblem &problem);
}

#endif
