/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_factory.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"
#include "nln_operator_factory.H"
#include "nln_operator_fas.H"
#include "nln_operator_ngmres.H"
#include "nln_operator_newton.H"
#include "nln_operator_nonlincg.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Constructor (empty) */
NLNSOL::NlnOperatorFactory::NlnOperatorFactory()
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Create the nonlinear operator */
Teuchos::RCP<NLNSOL::NlnOperatorBase>
NLNSOL::NlnOperatorFactory::Create(const Teuchos::ParameterList& params)
{
  const std::string optype = params.get<std::string>("Nonlinear Operator Type");
  if (optype == "Newton")
  {
    return Teuchos::rcp(new NLNSOL::NlnOperatorNewton());
  }
  else if (optype == "Nonlinear CG")
  {
    return Teuchos::rcp(new NLNSOL::NlnOperatorNonlinCG());
  }
  else if (optype == "AMG with FAS")
  {
    return Teuchos::rcp(new NLNSOL::NlnOperatorFas());
  }
  else if (optype == "NGMRES")
  {
    return Teuchos::rcp(new NLNSOL::NlnOperatorNGmres());
  }
  else
  {
    dserror("Unknown nonlinear operator %s.", optype.c_str());
  }

  return Teuchos::null;
}
