// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STI_RESULTTEST_HPP
#define FOUR_C_STI_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN

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
  class STIResultTest : public Core::Utils::ResultTest
  {
   public:
    //! constructor
    STIResultTest(const std::shared_ptr<STI::Algorithm>&
            sti_algorithm  //!< time integrator for scatra-thermo interaction
    );

    //! test special quantity not associated with a particular element or node
    void test_special(
        const Core::IO::InputParameterContainer&
            container,   ///< container with expected results as specified in the input file
        int& nerr,       //!< number of failed result tests
        int& test_count  ///< number of result tests
        ) override;

   private:
    //! get special result to be tested
    double result_special(const std::string& quantity  //! name of quantity to be tested
    ) const;

    //! return time integrator for monolithic scatra-thermo interaction
    const STI::Monolithic& sti_monolithic() const;

    //! time integrator for scatra-thermo interaction
    const std::shared_ptr<const STI::Algorithm> sti_algorithm_;
  };
}  // namespace STI
FOUR_C_NAMESPACE_CLOSE

#endif
