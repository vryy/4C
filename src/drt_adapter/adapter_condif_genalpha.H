/*----------------------------------------------------------------------*/
/*!
\file adapter_condif_genalpha.H

\brief Convection-Diffusion field adapter

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_CONDIF_GENALPHA_H
#define ADAPTER_CONDIF_GENALPHA_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/condif_genalpha_integration.H"

#include "adapter_condif_genalpha.H"
#include "adapter_condif.H"
#include "adapter_utils.H"


/// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {
  /*!
  /// adapter to generalized alpha convection-diffusion algorithm

    \author gjb
    \date 04/08
   */
  class ConDifGenAlpha : public ConDifTimeIntegration
  {
  public:
    ConDifGenAlpha(
      Teuchos::RCP<DRT::Discretization>      dis,
      Teuchos::RCP<LINALG::Solver>           solver,
      Teuchos::RCP<ParameterList>            params,
      Teuchos::RCP<IO::DiscretizationWriter> output
      );

    void PrepareTimeStep();
    void Update();
    void Output();
    void SetVelocityField(int veltype, int velfuncno);
    void SetVelocityField(int veltype, RCP<const Epetra_Vector> extvel);
    void ReadRestart(int step);
    void Integrate();
    void Solve();

  private:

    /// the actual fluid algorithm
    CondifGenAlphaIntegration       condif_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

  };

}

#endif
#endif
