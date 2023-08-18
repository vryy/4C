/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure-thermo interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SSTI_RESULTTEST_H
#define BACI_SSTI_RESULTTEST_H

#include "baci_lib_resulttest.H"

namespace SSTI
{
  class SSTIAlgorithm;

  /*!
    @brief result testing functionality for scalar-structure-thermo interaction problems

    This class provides result testing functionality for quantities associated with
    scalar-structure-thermo interaction as an overall problem type. Quantities associated
    with either the scalar or the structural field are not tested by this class, but
    by field-specific result testing classes. Feel free to extend this class if necessary.

  */
  class SSTIResultTest : public DRT::ResultTest
  {
   public:
    //! constructor
    explicit SSTIResultTest(const SSTI::SSTIAlgorithm& ssti_algorithm);

    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /*!
     * @brief get special result to be tested
     *
     * @param[in] quantity  name of quantity to be tested
     * @return special result
     */
    double ResultSpecial(const std::string& quantity) const;

    //! ssti algorithm
    const SSTI::SSTIAlgorithm& ssti_algorithm_;
  };
}  // namespace SSTI
#endif
