/*-------------------------------------------------------------------------
/*!
\file inpar_volmortar.H

\brief

\level 1

<pre>
\maintainer Philipp Farah
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

  // Type weighting function for quadr. problems
  enum DualQuad
  {
    dualquad_no_mod,   ///< no modification
    dualquad_lin_mod,  ///< linear modification
    dualquad_quad_mod  ///< quadr. modification
  };

  // Type of coupling
  enum CouplingType
  {
    couplingtype_volmortar,  ///< volmortar
    couplingtype_coninter    ///< consist. interpolation
  };

  // Type of coupling
  enum Shapefcn
  {
    shape_dual,  ///< dual shape functions
    shape_std    ///< std. shape functions --> lumped
  };

  /// set the volmortar parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  } //namespace volmortar
} //namespace inpar

#endif /* INPAR_VOLMORTAR_H_ */
