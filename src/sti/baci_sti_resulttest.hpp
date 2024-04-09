/*----------------------------------------------------------------------*/
/*! \file

\brief result testing functionality for scatra-thermo interaction problems

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STI_RESULTTEST_HPP
#define FOUR_C_STI_RESULTTEST_HPP

#include "baci_config.hpp"

#include "baci_lib_resulttest.hpp"

BACI_NAMESPACE_OPEN

namespace STI
{
  // forward declaration
  class Algorithm;
  class Monolithic;

  /*!
    \brief result testing functionality for scatra-thermo interaction problems

    This class provides result testing functionality for quantities associated with
    scatra-thermo interaction in general. Quantities associated with either the scatra
    or the thermo field are not tested by this class, but by field-specific result
    testing classes. Feel free to extend this class if necessary.

    \sa ResultTest
    \author fang
    \date 01/2017
  */
  class STIResultTest : public DRT::ResultTest
  {
   public:
    //! constructor
    STIResultTest(const Teuchos::RCP<STI::Algorithm>&
            sti_algorithm  //!< time integrator for scatra-thermo interaction
    );

    //! test special quantity not associated with a particular element or node
    void TestSpecial(
        INPUT::LineDefinition& res,  //!< input file line containing result test specification
        int& nerr,                   //!< number of failed result tests
        int& test_count              ///< number of result tests
        ) override;

   private:
    //! get special result to be tested
    double ResultSpecial(const std::string& quantity  //! name of quantity to be tested
    ) const;

    //! return time integrator for monolithic scatra-thermo interaction
    const STI::Monolithic& STIMonolithic() const;

    //! time integrator for scatra-thermo interaction
    const Teuchos::RCP<const STI::Algorithm> sti_algorithm_;
  };
}  // namespace STI
BACI_NAMESPACE_CLOSE

#endif
