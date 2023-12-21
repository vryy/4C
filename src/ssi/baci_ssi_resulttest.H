/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SSI_RESULTTEST_H
#define BACI_SSI_RESULTTEST_H

#include "baci_config.H"

#include "baci_lib_resulttest.H"

BACI_NAMESPACE_OPEN

// forward declaration
namespace DRT
{
  class Node;
}

namespace SSI
{
  // forward declarations
  class SSIBase;
  class SSIMono;

  /*!
    \brief result testing functionality for scalar-structure interaction problems

    This class provides result testing functionality for quantities associated with
    scalar-structure interaction as an overall problem type. Quantities associated
    with either the scalar or the structural field are not tested by this class, but
    by field-specific result testing classes. Feel free to extend this class if necessary.

    \sa ResultTest
    \author fang
    \date 11/2017
  */
  class SSIResultTest : public DRT::ResultTest
  {
   public:
    /*!
     * @brief constructor
     *
     * @param[in] ssi_base  time integrator for scalar-structure interaction
     */
    explicit SSIResultTest(const Teuchos::RCP<const SSI::SSIBase> ssi_base);

    void TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /*!
     * @brief get special result to be tested
     *
     * @param[in] quantity  name of quantity to be tested
     * @return special result
     */
    double ResultSpecial(const std::string& quantity) const;

    //! return time integrator for monolithic scalar-structure interaction
    const SSI::SSIMono& SSIMono() const;

    //! time integrator for scalar-structure interaction
    const Teuchos::RCP<const SSI::SSIBase> ssi_base_;
  };
}  // namespace SSI
BACI_NAMESPACE_CLOSE

#endif
