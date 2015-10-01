/*----------------------------------------------------------------------*/
/*!
\file inpar_cut.H

\brief Input parameters for cut library

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_CUT_H
#define INPAR_CUT_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  namespace CUT
  {

  enum VCellGaussPts
  {
    VCellGaussPts_Tessellation,
    VCellGaussPts_MomentFitting,
    VCellGaussPts_DirectDivergence
  };

  enum BCellGaussPts
  {
    BCellGaussPts_Tessellation,
    BCellGaussPts_MomentFitting,
    BCellGaussPts_DirectDivergence
  };

  enum ElementIntegrationType
  {
    EleIntType_Undecided,        ///< for this element, we still have to set the integration type
    EleIntType_StandardUncut,    ///< classical Gauss rule from the underlying element (if element is uncut)
    EleIntType_Tessellation,     ///< element uses Tessellation for integration of the element's volume-cells
    EleIntType_MomentFitting,    ///< element uses MomentFitting for integration of the element's volume-cells
    EleIntType_DirectDivergence  ///< element uses DirectDivergence for integration of the element's volume-cells
  };

  enum NodalDofSetStrategy
  {
    NDS_Strategy_OneDofset_PerNodeAndPosition,  ///< suppress multiple dofset per node and allow just one dofset per node and per dofset position
    NDS_Strategy_full                           ///< full strategy based on found connections of volume-cells in support of the nodal shape function
  };

  /// set the cut parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace CUT

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CUT_H
