/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer3.cpp

\brief optimizer of the topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double TOPOPT::Optimizer::EvaluateObjective(
    const Teuchos::RCP<const Epetra_Vector> porosity,
    const std::string& action
) const
{
  double value = 0.0;
  // TODO fill!

  if (action=="dissipation")
  {

  }
  else if (action=="pressure inlet")
  {

  }
  else if (action=="pressure drop")
  {

  }
  else
    dserror("evaluation of objective function not defined for %s!",action.c_str());

  return value;
}
