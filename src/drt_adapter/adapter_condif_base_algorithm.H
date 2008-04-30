/*----------------------------------------------------------------------*/
/*!
\file adapter_condif_base_algorithm.H

\brief convection-diffusion field base algorithm

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_CONDIF_BASE_ALGORITHM_H
#define ADAPTER_CONDIF_BASE_ALGORITHM_H

#include "adapter_condif.H"


/// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general fluid field interface for FSI
  /*!

  This is the FSI algorithm's view on a fluid algorithm. This pure virtual
  interface contains all the methods any FSI algorithm might want to
  call. The idea is to implement this interface with a concrete adapter class
  for each fluid algorithm we want to use for FSI.

  FSI is quite demanding when it comes to knowledge about the internal details
  of a fluid algorithm. Furthermore there are different coupling versions, all
  of them requiring a slightly different view. Yet these views have a lot in
  common, so a common adapter clas seems appropriate.

  We use this adapter interface instead of deriving from the fluid algorithm
  class. The good thing is that this way we keep control over the exported
  entities. The down side is that the fluid algorithm still has to grant
  access to a lot of internal state variables.

  Currently supported FSI couplings:

  - Dirichlet-Neumann coupling (fixed-point, Newton-Krylov, vector
    extrapolation)

  - Monolithic (Newton with overlapping blocks and block preconditioning)

  \warning Further cleanup is still needed.

  \sa Structure, Ale
  \author u.kue
  \date 11/07
  */

  /// fluid field solver
  class ConDifBaseAlgorithm
  {
  public:

    /// constructor
    explicit ConDifBaseAlgorithm(const Teuchos::ParameterList& prbdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~ConDifBaseAlgorithm();

    /// convection-diffusion field solver
    ConDifTimeIntegration& ConDifField() { return *condif_;}
    //const ConDif& ConDifField() const { return *condif_; }

  private:

    /// convection-diffusion field solver
    Teuchos::RCP<ConDifTimeIntegration> condif_;

  }; // class ConDifBaseAlgorithm

} // namespace ADAPTER

#endif
#endif
