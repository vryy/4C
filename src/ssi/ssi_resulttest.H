/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef SSI_RESULTTEST_H
#define SSI_RESULTTEST_H

#include "resulttest.H"

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

    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

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
#endif
