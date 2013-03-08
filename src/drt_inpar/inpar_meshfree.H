/*!--------------------------------------------------------------------------
\file inpar_meshfree.H

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#ifndef INPAR_MESHFREE_H
#define INPAR_MESHFREE_H

#include "inpar_parameterlist_utils.H"

namespace INPAR
{

/// Input parameters for Meshfree
namespace MESHFREE
{
 /*--------------------------------------------------------------------------*
  *                                                                          *
  * enums to specify meshfree problems                                       *
  *                                                                          *
  * Any changes should be consistent with changes in DRT_validarameters.cpp  *
  *                                                                          *
  *                                                                nis Mar12 *
  *--------------------------------------------------------------------------*/
  //! type of meshfree discretisation
  enum meshfreetype {maxent, particle};
  //! type of prior
  enum priortype {p_gauss};
  //! skewness type of prior
  enum priorskew {p_sym, p_0815};
  //! type of compliance condition
  enum compltype {c_linear, c_stream, c_freesp};
  //! type of gap function/projection
  enum gaptype {g_perp, g_stream};

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Way the neighbourhood search nodes<>elements is managed
   *
   * \param procwise : each proc sends all its proper (row) nodes to all other
   *                   procs where they are checked against (all ghosted)
   *                   knots (broodforce)
   * \param blockwise: nodes and knots are assign to blocks of size of the
   *                   maximum range of influence; nodes are only checked
   *                   against knots in own and neighbour blocks
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  enum NodeAssignType {procwise, blockwise};


}  // namespace MESHFREE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MESHFREE_H
