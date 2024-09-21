/*----------------------------------------------------------------------*/
/*! \file
\brief Creation of structural time integrators in accordance with user's wishes
\level 1
*/

#ifndef FOUR_C_STRUCTURE_TIMINT_CREATE_HPP
#define FOUR_C_STRUCTURE_TIMINT_CREATE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class Solver;
}

namespace Core::IO
{
  class DiscretizationWriter;
}

/*----------------------------------------------------------------------*/
namespace Solid
{
  // forward declarations
  class TimInt;
  class TimIntImpl;
  class TimIntExpl;

  /*====================================================================*/
  //! Create marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<Solid::TimInt> tim_int_create(
      const Teuchos::ParameterList& timeparams,             //!< time parameters
      const Teuchos::ParameterList& ioflags,                //!< input-output-flags
      const Teuchos::ParameterList& sdyn,                   //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,                //!< extra flags
      Teuchos::RCP<Core::FE::Discretization>& actdis,       //!< discretisation
      Teuchos::RCP<Core::LinAlg::Solver>& solver,           //!< the solver
      Teuchos::RCP<Core::LinAlg::Solver>& contactsolver,    //!< the solver for contact/meshtying
      Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b implicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<Solid::TimIntImpl> tim_int_impl_create(
      const Teuchos::ParameterList& timeparams,             //!< time parameters
      const Teuchos::ParameterList& ioflags,                //!< input-output-flags
      const Teuchos::ParameterList& sdyn,                   //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,                //!< extra flags
      Teuchos::RCP<Core::FE::Discretization>& actdis,       //!< discretisation
      Teuchos::RCP<Core::LinAlg::Solver>& solver,           //!< the solver
      Teuchos::RCP<Core::LinAlg::Solver>& contactsolver,    //!< the contact solver
      Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b explicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<Solid::TimIntExpl> tim_int_expl_create(
      const Teuchos::ParameterList& timeparams,             //!< time parameters
      const Teuchos::ParameterList& ioflags,                //!< input-output-flags
      const Teuchos::ParameterList& sdyn,                   //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,                //!< extra flags
      Teuchos::RCP<Core::FE::Discretization>& actdis,       //!< discretisation
      Teuchos::RCP<Core::LinAlg::Solver>& solver,           //!< the solver
      Teuchos::RCP<Core::LinAlg::Solver>& contactsolver,    //!< the solver for contact/meshtying
      Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< output writer
  );

}  // namespace Solid

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
