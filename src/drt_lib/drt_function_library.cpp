/*----------------------------------------------------------------------*/
/*!

\brief Collection of functions that are not problem-specific

\maintainer Sebastian Proell

\level 2

*/
/*----------------------------------------------------------------------*/

#include <Sacado.hpp>
#include "drt_function.H"
#include "drt_function_library.H"


DRT::UTILS::TranslatedFunction::TranslatedFunction(
    Teuchos::RCP<Function> origin, Teuchos::RCP<Function> local)
{
  if (origin->NumberComponents() != nsd)
    dserror("Origin function needs to have exactly %d components but %d were given.", nsd,
        origin->NumberComponents());
  originFunction_ = origin;
  localFunction_ = local;
}

double DRT::UTILS::TranslatedFunction::Evaluate(const int index, const double* x, double t)
{
  if (index < 0 or index >= localFunction_->NumberComponents())
    dserror("Index must be between 0 and %d but is %d.", localFunction_->NumberComponents(), index);

  double new_coord[nsd];
  for (int i = 0; i < nsd; i++)
  {
    new_coord[i] = x[i] - originFunction_->Evaluate(i, x, t);
  }

  return localFunction_->Evaluate(index, new_coord, t);
}
