/*----------------------------------------------------------------------*/
/*!
\file inpar_searchtree.H

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_SEARCHTREE_H_
#define INPAR_SEARCHTREE_H_

/*----------------------------------------------------------------------*/
namespace INPAR
{
  /// search tree type
  namespace GEO
  {
    /// specify tree type
    enum TreeType
    {
      Notree,
      Octree3D,
      Quadtree3D,
      Quadtree2D
    };

  }  // namespace GEO

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // /*INPAR_SEARCHTREE_H_*/
#endif  // #ifdef CCADISCRET
