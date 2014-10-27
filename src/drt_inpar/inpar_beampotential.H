/*----------------------------------------------------------------------*/
/*!
\file inpar_beampotential.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_BEAMPOTENTIAL_H
#define INPAR_BEAMPOTENTIAL_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for potential-based beam interaction
  namespace BEAMPOTENTIAL
  {
    /// type of potential interaction
    /// (this enum represents the input file parameter BEAMPOTENTIAL_TYPE)
    enum BeamPotentialType
    {
      beampot_surf,         ///< surface potential
      beampot_vol           ///< volume potential
    };

  }  // namespace BEAMPOTENTIAL

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_BEAMPOTENTIAL_H
