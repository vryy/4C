/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_utils.H

\brief Utility routines for the structural non-linear solver classes.

\maintainer Michael Hiermeier

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_NLN_SOLVER_UTILS_H_
#define SRC_DRT_STRUCTURE_NEW_STR_NLN_SOLVER_UTILS_H_

#include "../drt_inpar/inpar_structure.H"
#include "../solver_nonlin_nox/nox_nln_statustest_factory.H"
#include "../solver_nonlin_nox/nox_nln_enum_lists.H"

// forward declaration
namespace Teuchos
{
  class ParameterList;
}
namespace LINALG
{
  class Solver;
}
namespace STR
{
  namespace TIMINT
  {
    class Base;
    class BaseDataSDyn;
  }  // namespace TIMINT
  namespace NLN
  {
    namespace SOLVER
    {
      //! Check if a xml status test file is specified in the dat file
      bool IsXMLStatusTestFile(const Teuchos::ParameterList& pstatus);

      /*! \brief Create quantity types
       *
       *  This function translates the model type enums and element technology enums
       *  into quantity types.
       *
       *  \author Hiermeier
       */
      void CreateQuantityTypes(std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes,
          const STR::TIMINT::BaseDataSDyn& datasdyn);

      void ConvertModelType2QuantityType(const enum INPAR::STR::ModelType& mt,
          std::vector<enum NOX::NLN::StatusTest::QuantityType>& qt);

      void ConvertEleTech2QuantityType(const enum INPAR::STR::EleTech& et,
          std::vector<enum NOX::NLN::StatusTest::QuantityType>& qt);

      /*! \brief Create a status test parameter list
       *
       *  A status test parameter list for the outer status test is created. The
       *  information comes from the dat file. Actually we consider only convergence
       *  tests which were already in BACI before the NOX framework was introduced.
       *
       *  Feel free to extend the given framework, to generalize it or to use a xml
       *  file, where the path is specified in the "STRUCT NOX/Status Test" sublist.
       *  The last option comes for free and there is no need to modify any code
       *  fragments.
       *
       *  \author Hiermeier
       */
      void SetStatusTestParams(Teuchos::ParameterList& pstatus,
          const STR::TIMINT::BaseDataSDyn& datasdyn,
          const std::set<enum NOX::NLN::StatusTest::QuantityType>& qt);

      //! Split the given tests into and/or combinations
      void SplitAndOrCombo(std::vector<enum NOX::NLN::StatusTest::QuantityType>& combo_or,
          std::vector<enum NOX::NLN::StatusTest::QuantityType>& combo_and,
          const STR::TIMINT::BaseDataSDyn& datasdyn, const std::string& testname,
          const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes);

      /*! \brief Set the combination of different NormF or NormUpdate tests in the status test
       * parameter list
       *
       *  Maybe the current implementation needs a short explanation:
       *  You can specify different combinations in your dat-file. Let's concentrate on
       *  the NormF case and imagine the following status test settings:
       *
       *  NORMCOMBI_RESFINCO           AND        (RESIDUAL and PRESSURE)
       *  NORMCOMBI_EASRES             OR         (RESIDUAL or  EAS)
       *  NORMCOMBI_RESFCONSTR         AND        (RESIDUAL and CONTACT)
       *  NORMCOMBI_RESFPLASTCONSTR    OR         (RESIDUAL or  PLASTICITY)
       *
       *  Following will happen:
       *
       *  (STRUCTURE and PRESSURE and CONTACT) or EAS or PLASTICITY
       *
       *  This means, that the OR-combination is a strong OR and ignores/omits all remaining tests!
       *  In the old implementation the exact behavior was more or less random, because it depended
       * on the order of the implemented tests. In this way the OR-combination is a debugging tool
       * and you can check your results by forcing only one of the corresponding residuals to zero
       *  (i.e. the named partial residual or all remaining parts of the residual and-combinations).
       *
       *  If you want a different behavior, please use a xml file instead.
       *
       *  One last note: In the case you want to commit your dat-file as a test case, your algorithm
       *  should not use any OR-combination, since the goal of your algorithm should be to reduce
       * the whole residual. If you think a part of your residual can not be reduced in a sufficient
       * way, think again and if it stays your opinion, do not check it (by using a xml-file for
       * example). This makes things easier to read and understand. If you use only
       * AND-combinations, you can use the QUANTITY parameter list name option for NormF, NormWRMS
       * and NormUpdate tests. See the NOX::NLN::StatusTest::Factory for more information.
       *
       *  \author Hiermeier
       *  \date Oct 13, 2015
       */
      void SetComboQuantityTestParams(Teuchos::ParameterList& p,
          const STR::TIMINT::BaseDataSDyn& datasdyn, const std::size_t& count,
          const std::string& testname,
          const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes);

      /*! \brief Set the status test corresponding to the given quantity
       *
       *  Create a new sublist corresponding to the count variable and create the specific status
       * test entries in the parameter list for the given quantity
       *
       *  \author Hiermeier
       */
      void SetQuantityTestParams(Teuchos::ParameterList& p,
          const STR::TIMINT::BaseDataSDyn& datasdyn,
          const enum NOX::NLN::StatusTest::QuantityType& qtype, const std::size_t& count,
          const std::string& testname);

      //! Set the status test parameters corresponding to the given quantity
      void SetQuantityTestParams(Teuchos::ParameterList& p,
          const STR::TIMINT::BaseDataSDyn& datasdyn,
          const enum NOX::NLN::StatusTest::QuantityType& qtype, const std::string& testname);

      //! Set the NormUpdate status test parameters
      void SetNormUpdateParams(Teuchos::ParameterList& qlist,
          const enum NOX::NLN::StatusTest::QuantityType& qtype,
          const enum INPAR::STR::ConvNorm& toltype, const double& tol,
          const enum INPAR::STR::VectorNorm& normtype);
      void SetNormUpdateParams(Teuchos::ParameterList& qlist,
          const enum NOX::NLN::StatusTest::QuantityType& qtype, const double& alpha,
          const double& beta, const enum INPAR::STR::ConvNorm& toltype, const double& tol,
          const enum INPAR::STR::VectorNorm& normtype, const bool& isscaled);

      //! Set the NormF status test parameters
      void SetNormFParams(Teuchos::ParameterList& qlist,
          const enum NOX::NLN::StatusTest::QuantityType& qtype,
          const enum INPAR::STR::ConvNorm& toltype, const double& tol,
          const enum INPAR::STR::VectorNorm& normtype);
      void SetNormFParams(Teuchos::ParameterList& qlist,
          const enum NOX::NLN::StatusTest::QuantityType& qtype,
          const enum INPAR::STR::ConvNorm& toltype, const double& tol,
          const enum INPAR::STR::VectorNorm& normtype, const bool& isscaled);

      //! Set the ActiveSet status test parameters
      void SetActiveSetParams(
          Teuchos::ParameterList& qlist, const enum NOX::NLN::StatusTest::QuantityType& qtype);

    }  // namespace SOLVER
  }    // namespace NLN
}  // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_NLN_SOLVER_UTILS_H_ */
