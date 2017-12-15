/*----------------------------------------------------------------------*/
/*!
\file ssi_monolithic_resulttest.H

\brief result testing functionality for monolithic scalar-structure interaction problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
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
  // forward declaration
  class SSI_Mono;

  /*!
    \brief result testing functionality for monolithic scalar-structure interaction problems

    This class provides result testing functionality for quantities associated with monolithic
    scalar-structure interaction as an overall problem type. Quantities associated with either
    the scalar or the structural field are not tested by this class, but by field-specific result
    testing classes. Feel free to extend this class if necessary.

    \sa ResultTest
    \author fang
    \date 11/2017
  */
  class SSI_Mono_ResultTest : public DRT::ResultTest
  {
    public:

      //! constructor
      SSI_Mono_ResultTest(
          const Teuchos::RCP<const SSI::SSI_Mono>   ssi_mono   //!< time integrator for monolithic scalar-structure interaction
          );

      //! test quantity associated with a particular node
      void TestNode(
          DRT::INPUT::LineDefinition&   res,         //!< input file line containing result test specification
          int&                          nerr,        //!< number of failed result tests
          int&                          test_count   //!< number of result tests
          );

      //! test special quantity not associated with a particular element or node
      void TestSpecial(
          DRT::INPUT::LineDefinition&   res,         //!< input file line containing result test specification
          int&                          nerr,        //!< number of failed result tests
          int&                          test_count   //!< number of result tests
          );

    private:

      //! get nodal result to be tested
      double ResultNode(
          const std::string   quantity,   //!< name of quantity to be tested
          DRT::Node*          node        //!< node carrying the result to be tested
          ) const;

      //! get special result to be tested
      double ResultSpecial(
          const std::string&   quantity   //!< name of quantity to be tested
          ) const;

      //! time integrator for monolithic scalar-structure interaction
      const Teuchos::RCP<const SSI::SSI_Mono> ssi_mono_;
  };
}
#endif
