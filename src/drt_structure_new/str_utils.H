/*-----------------------------------------------------------*/
/*! \file

\brief Utility methods for the structural time integration.

\maintainer Matthias Mayr

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_UTILS_H_
#define SRC_DRT_STRUCTURE_NEW_STR_UTILS_H_

#include <set>

#include <NOX_Abstract_Vector.H>
#include "../solver_nonlin_nox/nox_nln_enum_lists.H"
#include "../drt_inpar/inpar_structure.H"

// forward declaration
class Epetra_Vector;

namespace LINALG
{
  class Solver;
}  // namespace LINALG

namespace NOX
{
  namespace Epetra
  {
    class Scaling;
  }

  namespace NLN
  {
    namespace CONSTRAINT
    {
      namespace Interface
      {
        class Required;
        class Preconditioner;
      }  // namespace Interface
    }    // namespace CONSTRAINT
  }      // namespace NLN
}  // namespace NOX

namespace STR
{
  class Integrator;

  namespace TIMINT
  {
    class BaseDataSDyn;
    class BaseDataGlobalState;
  }  // namespace TIMINT

  namespace NLN
  {
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
    void ConvertModelType2SolType(std::vector<enum NOX::NLN::SolutionType>& soltypes,
        std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& slinsolvers,
        const std::set<enum INPAR::STR::ModelType>& modeltypes,
        const std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver>>& mlinsolvers);

    /*! \brief Convert the structural model type enumerator to a nox nln solution
     *  type enumerator
     *
     *  \param modeltype (in) : Model type enumerator which has to be converted
     *  \param do_check  (in) : Check if a corresponding solution type exists */
    enum NOX::NLN::SolutionType ConvertModelType2SolType(
        const enum INPAR::STR::ModelType& modeltype, const bool& do_check);

    /*! \brief Convert the structural model type enumerator to a nox nln solution
     *  type enumerator and check if the conversion was successful
     *
     *  \param modeltype (in) : Model type enumerator which has to be converted. */
    inline enum NOX::NLN::SolutionType ConvertModelType2SolType(
        const enum INPAR::STR::ModelType& modeltype)
    {
      return ConvertModelType2SolType(modeltype, true);
    }

    /*! \brief Convert the nox nln solution type enumerator to a structural model
     *  type enumerator
     *
     *  \param soltype (in)   : Solution type enumerator which has to be converted.
     *  \param do_check  (in) : Check if a corresponding model type exists. */
    enum INPAR::STR::ModelType ConvertSolType2ModelType(
        const enum NOX::NLN::SolutionType& soltype, const bool& do_check);

    /*! \brief Convert the structural model type enumerator to a nox nln solution
     *  type enumerator and check if the conversion was successful
     *
     *  \param soltype (in) : Solution type enumerator which has to be converted. */
    inline enum INPAR::STR::ModelType ConvertSolType2ModelType(
        const enum NOX::NLN::SolutionType& soltype)
    {
      return ConvertSolType2ModelType(soltype, true);
    }

    /*! \brief Convert the nox nln statustest quantity type enumerator to a structural model
     *  type enumerator
     *
     *  \param qtype (in)    : Quantity type enumerator which has to be converted.
     *  \param do_check (in) : Check if a corresponding model type exists. */
    enum INPAR::STR::ModelType ConvertQuantityType2ModelType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype, const bool& do_check);

    /*! \brief Convert the nox nln statustest quantity type enumerator to a structural model
     *  type enumerator and check if the conversion was successful
     *
     *  \param qtype (in) : Quantity type enumerator which has to be converted. */
    inline enum INPAR::STR::ModelType ConvertQuantityType2ModelType(
        const enum NOX::NLN::StatusTest::QuantityType& qtype)
    {
      return ConvertQuantityType2ModelType(qtype, true);
    }

    /*! \brief Convert the nox nln statustest quantity type enumerator to a structural model
     *  type enumerator
     *
     *  \param qtype (in)    : Quantity type enumerator which has to be converted. */
    enum INPAR::STR::EleTech ConvertQuantityType2EleTech(
        const enum NOX::NLN::StatusTest::QuantityType& qtype);

    //! Returns the optimization type of the underlying structural problem
    enum NOX::NLN::OptimizationProblemType OptimizationType(
        const std::vector<enum NOX::NLN::SolutionType>& soltypes);

    /// convert structure condition number type to a nox condition number type
    NOX::NLN::LinSystem::ConditionNumber Convert2NoxConditionNumberType(
        const INPAR::STR::ConditionNumber stype);

    //! Set the constraint interfaces
    void CreateConstraintInterfaces(
        std::map<enum NOX::NLN::SolutionType,
            Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>>& iconstr,
        STR::Integrator& integrator, const std::vector<enum NOX::NLN::SolutionType>& soltypes);

    //! Set the constraint preconditioner interfaces
    void CreateConstraintPreconditioner(
        std::map<NOX::NLN::SolutionType,
            Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Preconditioner>>& iconstr_prec,
        STR::Integrator& integrator, const std::vector<enum NOX::NLN::SolutionType>& soltypes);

    //! Create object to scale linear system
    void CreateScaling(Teuchos::RCP<NOX::Epetra::Scaling>& iscale,
        STR::TIMINT::BaseDataSDyn DataSDyn, STR::TIMINT::BaseDataGlobalState GState);
  }  // namespace NLN
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_UTILS_H_ */
