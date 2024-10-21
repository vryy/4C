// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_TIMADA_CREATE_HPP
#define FOUR_C_STRUCTURE_TIMADA_CREATE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
namespace Solid
{
  // forward declarations
  class TimInt;
  class TimAda;

  /*====================================================================*/
  //! Create auxiliary time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<Solid::TimAda> tim_ada_create(
      const Teuchos::ParameterList& ioflags,     //!< input-output-flags
      const Teuchos::ParameterList& timeparams,  //!< structural dynamic flags
      const Teuchos::ParameterList& sdyn,        //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,     //!< extra flags
      const Teuchos::ParameterList& tap,         //!< adaptive input flags
      Teuchos::RCP<Solid::TimInt> tis            //!< marching time integrator
  );

}  // namespace Solid

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
