/*-----------------------------------------------------------*/
/*! \file

\brief Central header for forward declaration of Trilinos's classes.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_FORWARD_DECL_H
#define BACI_SOLVER_NONLIN_NOX_FORWARD_DECL_H

#include "baci_config.H"

class Epetra_BlockMap;
class Epetra_Map;
class Epetra_Vector;
class Epetra_Comm;
class Epetra_Operator;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace NOX
{
  class GlobalData;
  class Observer;
  class Utils;
  namespace Abstract
  {
    class Group;
    class Vector;
  }  // namespace Abstract
  namespace Direction
  {
    class Generic;
    class UserDefinedFactory;
  }  // namespace Direction
  namespace Epetra
  {
    class Scaling;
    class Vector;
    class LinearSystem;
    namespace Interface
    {
      class Jacobian;
      class Preconditioner;
      class Required;
    }  // namespace Interface
  }    // namespace Epetra
  namespace LineSearch
  {
    class Generic;
  }  // namespace LineSearch
  namespace MeritFunction
  {
    class Generic;
  }  // namespace MeritFunction
  namespace Solver
  {
    class Generic;
  }  // namespace Solver
  namespace StatusTest
  {
    class Factory;
    class Generic;
  }  // namespace StatusTest
}  // namespace NOX

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

#endif  // BACI_SOLVER_NONLIN_NOX_FORWARD_DECL_H
