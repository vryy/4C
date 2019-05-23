/*----------------------------------------------------------------------*/
/*!
\file sti_resulttest.H

\brief result testing functionality for scatra-thermo interaction problems

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/
#ifndef STI_RESULTTEST_H
#define STI_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"

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
    virtual void TestSpecial(
        DRT::INPUT::LineDefinition& res,  //!< input file line containing result test specification
        int& nerr,                        //!< number of failed result tests
        int& test_count                   ///< number of result tests
    );

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
#endif
