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
#ifndef INPAR_SEARCHTREE_H_
#define INPAR_SEARCHTREE_H_

#include "inpar_parameterlist_utils.H"

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

    /// set the searchtree parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace GEO

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // /*INPAR_SEARCHTREE_H_*/
