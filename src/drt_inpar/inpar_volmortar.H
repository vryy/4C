/*-------------------------------------------------------------------------
\file inpar_volmortar.H

\brief

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_VOLMORTAR_H_
#define INPAR_VOLMORTAR_H_


/*----------------------------------------------------------------------*/
namespace INPAR
{
  /// Input parameters for vol. mortar
  namespace VOLMORTAR
  {
  // Type of integration procedure
  enum IntType
  {
    inttype_segments,       ///< cut procedure of volume meshes
    inttype_elements        ///< fast, elementwise integration
  };

  // Type of cut procedure
  enum CutType
  {
    cuttype_directdivergence, ///< direct divergence for integration
    cuttype_tessellation      ///< tessellation of volume meshes
  };

  } //namespace volmortae
} //namespace inpar

#endif /* INPAR_VOLMORTAR_H_ */
