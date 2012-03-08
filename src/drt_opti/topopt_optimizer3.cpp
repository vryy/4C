/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer3.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET


#include "../drt_lib/drt_dserror.H"

#include "topopt_optimizer.H"



// TODO fill!
double TOPOPT::Optimizer::EvaluateObjective(
    Teuchos::RCP<Epetra_Vector> porosity,
    const string& action
)
{
  double value = 0.0;

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



#endif  // #ifdef CCADISCRET
