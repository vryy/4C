/*----------------------------------------------------------------------*/
/*! \file
\brief Creation of auxiliary time integration scheme for time step size adaptivity
\level 1
*/

/*----------------------------------------------------------------------*/
/* macros */

#ifndef FOUR_C_STRUCTURE_TIMADA_CREATE_HPP
#define FOUR_C_STRUCTURE_TIMADA_CREATE_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_config.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
namespace STR
{
  // forward declarations
  class TimInt;
  class TimAda;

  /*====================================================================*/
  //! Create auxiliary time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimAda> TimAdaCreate(
      const Teuchos::ParameterList& ioflags,     //!< input-output-flags
      const Teuchos::ParameterList& timeparams,  //!< structural dynamic flags
      const Teuchos::ParameterList& sdyn,        //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,     //!< extra flags
      const Teuchos::ParameterList& tap,         //!< adaptive input flags
      Teuchos::RCP<STR::TimInt> tis              //!< marching time integrator
  );

}  // namespace STR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_TIMADA_CREATE_H
