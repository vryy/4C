/*-----------------------------------------------------------*/
/*!
\file inpar_meshfree.H

\brief inpar meshfree

\maintainer Keijo Nissen

\level 2

*/
/*-----------------------------------------------------------*/

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
   * Any changes should be consistent with changes in drt_validarameters.cpp  *
   *                                                                          *
   *                                                                nis Mar12 *
   *--------------------------------------------------------------------------*/

  //! type of meshfree discretisation
  enum meshfreetype {
    maxent,
    particle,
    geo_decoupled
  };

  //! type of prior
  enum priortype {
    p_gauss
  };

  //! skewness type of prior
  enum priorskew {
    p_sym,
    p_0815
  };

  //! type of compliance condition
  enum compltype {
    c_linear,
    c_stream,
    c_freesp
  };

  //! type of gap function/projection
  enum gaptype {
    g_perp,
    g_stream
  };

  //! type of bin specification
  enum xaabbspectype {
    input,
    dynamic
  };

  //! type of bin specification
  enum binspectype {
    cutoff,
    binsperdir,
    largestele
  };

  //! type of gap function/projection
  enum writebins {
    none,
    rows,
    cols
  };

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
  enum NodeAssignType {
    procwise,
    blockwise
  };


  /// set the meshfree parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace MESHFREE

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MESHFREE_H
