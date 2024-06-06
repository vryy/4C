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

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class Discretization;
}

namespace Core::LinAlg
{
  class Solver;
}

namespace Core::IO
{
  class DiscretizationWriter;
}

/*----------------------------------------------------------------------*/
namespace STR
{
  // forward declarations
  class TimInt;
  class TimIntImpl;
  class TimIntExpl;

  /*====================================================================*/
  //! Create marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimInt> TimIntCreate(
      const Teuchos::ParameterList& timeparams,             //!< time parameters
      const Teuchos::ParameterList& ioflags,                //!< input-output-flags
      const Teuchos::ParameterList& sdyn,                   //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,                //!< extra flags
      Teuchos::RCP<Discret::Discretization>& actdis,        //!< discretisation
      Teuchos::RCP<Core::LinAlg::Solver>& solver,           //!< the solver
      Teuchos::RCP<Core::LinAlg::Solver>& contactsolver,    //!< the solver for contact/meshtying
      Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b implicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimIntImpl> TimIntImplCreate(
      const Teuchos::ParameterList& timeparams,             //!< time parameters
      const Teuchos::ParameterList& ioflags,                //!< input-output-flags
      const Teuchos::ParameterList& sdyn,                   //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,                //!< extra flags
      Teuchos::RCP<Discret::Discretization>& actdis,        //!< discretisation
      Teuchos::RCP<Core::LinAlg::Solver>& solver,           //!< the solver
      Teuchos::RCP<Core::LinAlg::Solver>& contactsolver,    //!< the contact solver
      Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b explicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimIntExpl> TimIntExplCreate(
      const Teuchos::ParameterList& timeparams,             //!< time parameters
      const Teuchos::ParameterList& ioflags,                //!< input-output-flags
      const Teuchos::ParameterList& sdyn,                   //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,                //!< extra flags
      Teuchos::RCP<Discret::Discretization>& actdis,        //!< discretisation
      Teuchos::RCP<Core::LinAlg::Solver>& solver,           //!< the solver
      Teuchos::RCP<Core::LinAlg::Solver>& contactsolver,    //!< the solver for contact/meshtying
      Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< output writer
  );

}  // namespace STR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
