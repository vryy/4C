/*-----------------------------------------------------------*/
/*! \file

\brief Central header for forward declaration of Trilinos's classes.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_FORWARD_DECL_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_FORWARD_DECL_HPP

#include "baci_config.hpp"

// Do not try to lint the following names. They are all defined by Trilinos. Including this header
// and not using all identifiers will lead to clang-tidy thinking these are our own names.
// NOLINTBEGIN(readability-identifier-naming)

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

// NOLINTEND(readability-identifier-naming)

#endif
