/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for nonlinear solvers

\maintainer Matthias Mayr

\level 1
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_SOLVER_NONLIN_H
#define INPAR_SOLVER_NONLIN_H

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace NLNSOL
  {
    /// set the nonlinear solver parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace NLNSOL
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_SOLVER_NONLIN_H
