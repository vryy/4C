/*!----------------------------------------------------------------------
\file ac_fsi_ls.cpp

\brief cpp-file associated with algorithmic routines for two-way coupled partitioned
       solution approaches to fluid-structure-scalar-scalar interaction
       (FS3I). Specifically related version for multiscale approches. This file thereby holds
       all functions related with the large time scale simulation and
       the small to large time scale 'communication'

\date 2015-07-29

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            089/289-10364

\level 3
----------------------------------------------------------------------*/


#include "ac_fsi.H"


/*----------------------------------------------------------------------*
 | timeloop for small time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
void FS3I::ACFSI::LargeTimeScaleLoop()
{
  while (LargeTimeScaleLoopNotFinished())
  {
    //For now nothing to do!
    dserror("sieht doch gut aus :-)");
  }
}

/*----------------------------------------------------------------------*
 | timeloop for large time scales                            Thon 07/15 |
 *----------------------------------------------------------------------*/
bool FS3I::ACFSI::LargeTimeScaleLoopNotFinished()
{
  return NotFinished() ;
}
