/*-----------------------------------------------------------*/
/*!
\file str_utils.H

\brief Utility methods for the structural time integration.

\maintainer Michael Hiermeier

\date Dec 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_UTILS_H_
#define SRC_DRT_STRUCTURE_NEW_STR_UTILS_H_

#include <NOX_Abstract_Vector.H>
#include "../solver_nonlin_nox/nox_nln_enum_lists.H"
#include "../drt_inpar/inpar_structure.H"

// forward declaration
class Epetra_Vector;
namespace LINALG {
class Solver;
} // namespace LINALG
namespace NOX {
namespace NLN {
namespace CONSTRAINT {
namespace Interface {
  class Required;
  class Preconditioner;
} // namespace Interface
} // namespace CONSTRAINT
} // namespace NLN
} // namespace NOX

namespace STR {
  class Integrator;

  /*! \brief Assemble a source Epetra_Vector into a target Epetra_Vector
   *
   *  The map of the source vector has to be a sub-map of the target vector.
   *  The entries of the source vector are added to the target vector:
   *
   *      target_vector[gid] = source_scalar * source[gid] + target_scalar * target[gid]
   *
   *  The remaining gids, which are no part of the source map stay untouched! */
  void AssembleVector(const double& scalar_target, Epetra_Vector& target,
      const double& scalar_source, const Epetra_Vector& source);

namespace NLN {
  //! Convert the structural vector types to a corresponding nox norm type
  enum NOX::Abstract::Vector::NormType Convert2NoxNormType(
      const enum INPAR::STR::VectorNorm& normtype);

  /*! Convert the structural model type enums to nox nln solution
   *  type enums
   *
   *  Convert the model type enums to the nox internal solution type
   *  enums. This is necessary, because the nox framework is not
   *  supposed to be restricted to structural problems only.
   */
  void ConvertModelType2SolType(
      std::vector<enum NOX::NLN::SolutionType>& soltypes,
      std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& slinsolvers,
      const std::set<enum INPAR::STR::ModelType>& modeltypes,
      const std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> >& mlinsolvers);

  //! Returns the optimization type of the underlying structural problem
  enum NOX::NLN::OptimizationProblemType
      OptimizationType(const std::vector<enum NOX::NLN::SolutionType>& soltypes);

  //! Set the constraint interfaces
  void CreateConstraintInterfaces(
      std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >& iconstr,
      STR::Integrator& integrator,
      const std::vector<enum NOX::NLN::SolutionType>& soltypes);

  //! Set the constraint preconditioner interfaces
  void CreateConstraintPreconditioner(
      std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Preconditioner> >& iconstr_prec,
      STR::Integrator& integrator,
      const std::vector<enum NOX::NLN::SolutionType>& soltypes);
} // namespace NLN
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_UTILS_H_ */
