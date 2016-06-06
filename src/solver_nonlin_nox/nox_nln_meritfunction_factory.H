/*-----------------------------------------------------------*/
/*!
\file nox_nln_meritfunction_factory.H

\brief Factory to create a merit function evaluation object.

\maintainer Michael Hiermeier

\date Jun 12, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_MERITFUNCTION_FACTORY_H_
#define NOX_NLN_MERITFUNCTION_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_globaldata.H"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class Utils;

namespace MeritFunction {
  class Generic;
}

namespace NLN {
namespace MeritFunction {
  class Factory
  {
  public:
    //! constructor
    Factory();

    //! destructor
    virtual ~Factory() {};

    //! choose between constrainted and unconstrainted
    Teuchos::RCP<NOX::MeritFunction::Generic> BuildMeritFunction
        (Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData);

    //! unconstrainted factory
    Teuchos::RCP<NOX::MeritFunction::Generic> BuildUnconstrainedMeritFunction
        (Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData);

    //! constrainted factory
    Teuchos::RCP<NOX::MeritFunction::Generic> BuildConstrainedMeritFunction
        (Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData);
  };

  /*! Nonmember function to build a merit function object.

  \relates NOX::NLNSOL::Constraint::MeritFunction::Factory

  */
  Teuchos::RCP<NOX::MeritFunction::Generic> BuildMeritFunction
      (Teuchos::RCP<const NOX::NLN::GlobalData>& noxNlnGlobalData);
} // MeritFunction
} // NLN
} // NOX

#endif /* NOX_NLN_MERITFUNCTION_FACTORY_H_ */
