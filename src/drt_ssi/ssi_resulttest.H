/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef SSI_RESULTTEST_H
#define SSI_RESULTTEST_H

#include "../drt_lib/drt_resulttest.H"

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
    //! constructor
    SSIResultTest(const Teuchos::RCP<const SSI::SSIBase>
            ssi_base  //!< time integrator for scalar-structure interaction
    );

    //! test quantity associated with a particular node
    void TestNode(
        DRT::INPUT::LineDefinition& res,  //!< input file line containing result test specification
        int& nerr,                        //!< number of failed result tests
        int& test_count                   //!< number of result tests
        ) override;

    //! test special quantity not associated with a particular element or node
    void TestSpecial(
        DRT::INPUT::LineDefinition& res,  //!< input file line containing result test specification
        int& nerr,                        //!< number of failed result tests
        int& test_count                   //!< number of result tests
        ) override;

   private:
    //! get nodal result to be tested
    double ResultNode(const std::string& quantity,  //!< name of quantity to be tested
        DRT::Node* node                             //!< node carrying the result to be tested
        ) const;

    //! get special result to be tested
    double ResultSpecial(const std::string& quantity  //!< name of quantity to be tested
        ) const;

    //! return time integrator for monolithic scalar-structure interaction
    const SSI::SSIMono& SSIMono() const;

    //! time integrator for scalar-structure interaction
    const Teuchos::RCP<const SSI::SSIBase> ssi_base_;
  };
}  // namespace SSI
#endif
