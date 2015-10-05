/*
 * str_nln_solver_FullNewton.H
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#ifndef STR_NLN_SOLVER_FULLNEWTON_H_
#define STR_NLN_SOLVER_FULLNEWTON_H_

#include "str_nln_solver_nox.H"
#include "../solver_nonlin_nox/nox_nln_statustest_factory.H"

namespace Teuchos {
  class ParameterList;
}

namespace STR {
namespace NLN {
namespace SOLVER {
  class FullNewton : public Nox
  {
  public:
    //! constructor
    FullNewton();

    //! destructor
    virtual ~FullNewton() {};

    // derived from the base class
    virtual void Setup();

  protected:
    //! set the full newton parameters in the nox parameter list
    void SetFullNewtonParams();

    //! create quantity types
    //! this function translates the modeltype enums and element technology enums
    //! in quantity types.
    void CreateQuantityTypes(
        std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes) const;

    void SetStatusTestParams(
        const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes);

    void SplitAndOrCombo(
        std::vector<enum NOX::NLN::StatusTest::QuantityType>& combo_or,
        std::vector<enum NOX::NLN::StatusTest::QuantityType>& combo_and,
        const std::string& testname,
        const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes);

    void SetComboQuantityTestParams(
        Teuchos::ParameterList& p,
        const std::size_t& count,
        const std::string& testname,
        const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes);

    void SetQuantityTestParams(
        Teuchos::ParameterList& p,
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const std::size_t& count,
        const std::string& testname
        );

    void SetQuantityTestParams(
        Teuchos::ParameterList& p,
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const std::string& testname
        );

    void SetNormIncrParams(
        Teuchos::ParameterList& qlist,
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const enum INPAR::STR::ConvNorm& toltype,
        const double& tol,
        const enum INPAR::STR::VectorNorm& normtype
        );

    void SetNormIncrParams(
        Teuchos::ParameterList& qlist,
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const double& alpha,
        const double& beta,
        const enum INPAR::STR::ConvNorm& toltype,
        const double& tol,
        const enum INPAR::STR::VectorNorm& normtype,
        const bool& isscaled
        );

    void SetNormFParams(
        Teuchos::ParameterList& qlist,
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const enum INPAR::STR::ConvNorm& toltype,
        const double& tol,
        const enum INPAR::STR::VectorNorm& normtype
        );


    void SetNormFParams(
        Teuchos::ParameterList& qlist,
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const enum INPAR::STR::ConvNorm& toltype,
        const double& tol,
        const enum INPAR::STR::VectorNorm& normtype,
        const bool& isscaled
        );

    void SetActiveSetParams(
        Teuchos::ParameterList& qlist,
        const enum NOX::NLN::StatusTest::QuantityType& qtype
        );

  };  // class FullNewton
} // namespace SOLVER
} // namespace NLN
} // namespace STR


#endif /* STR_NLN_SOLVER_FULLNEWTON_H_ */
