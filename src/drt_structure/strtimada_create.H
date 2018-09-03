/*----------------------------------------------------------------------*/
/*!
\file strtimada_create.H
\brief Creation of auxiliary time integration scheme for time step size adaptivity

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

#ifndef STRTIMADA_CREATE_H
#define STRTIMADA_CREATE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Epetra_MpiComm.h"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"


/*----------------------------------------------------------------------*/
/* Structural Dynamics */
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
#endif /* #ifndef STRTIMADA_CREATE_H */
