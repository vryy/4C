/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_factory.H

\maintainer Michael Hiermeier

\date Jul 31, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_STATUSTEST_FACTORY_H_
#define NOX_NLN_STATUSTEST_FACTORY_H_

#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace Teuchos {
  class ParameterList;
}

namespace NOX {
  class Utils;
namespace StatusTest {
  class Generic;
  class Factory;
}
namespace NLN {
namespace StatusTest {
  //! Supported quantity names for distinguished status tests
  //! NOTE: Each name corresponds to only ONE quantity
  enum QuantityType {
    quantity_unknown,
    // check structural quantities
    quantity_structure,
    // check normal contact quantities
    quantity_ncontact,
    // check tangential contact quantities
    quantity_tcontact
  };

  /// Map quantity name to std::string
  static inline const std::string QuantityType2String(const enum QuantityType type)
  {
    switch (type)
    {
    case quantity_structure :
      return "STRUCTURE";
      break;
    case quantity_ncontact :
      return "nCONTACT";
      break;
    case quantity_tcontact :
      return "tCONTACT";
      break;
    case quantity_unknown:
    default :
      return "unknown quantity type";
      break;
    }
    return "";
  };

  inline QuantityType String2QuantityType(const std::string& name)
  {
    QuantityType type = quantity_unknown;
    if (name=="structure" or name=="Structure")
      type = quantity_structure;
    else if (name=="nContact" or name=="Normal Contact")
      type = quantity_ncontact;
    else if (name=="tContact" or name=="Tangential Contact")
      type = quantity_tcontact;

    return type;
  };

  class Factory
  {
  public:
    //! Constructor.
    Factory();

    //! Destructor.
    virtual ~Factory() {};

    //! Returns a outer status test set from a parameter list.
    Teuchos::RCP<NOX::StatusTest::Generic> BuildOuterStatusTests(
        Teuchos::ParameterList& p, const NOX::Utils& utils,
        std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests) const;

  protected:
    /*! \brief Derived buildComboTest.

      Restatement of the base class function, because the internal recursive call
      has to be redefined. Note that the base class input parameter
      <tt>Number of Tests<\tt> is not used in this version.
    */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildComboTest(
        Teuchos::ParameterList& p, const NOX::Utils& u,
        std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
        tagged_tests) const;

    /*! \brief New implementation of the \c BuildNormFTest function.

        The underlying Status Test is capable of a variety of quantities at the same
        time and combines the RelativeNormF and NormF test.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildNormFTest(
        Teuchos::ParameterList& p, const NOX::Utils& u, const bool& relativeNormF=false) const;

    /*! \brief New implementation of the \c BuildNormWRMSTest function.

        The underlying Status Test is capable of a variety of quantities at the same
        time.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildNormWRMSTest(
        Teuchos::ParameterList& p, const NOX::Utils& u) const;

    //! Checks if a tag is present in the param list and adds the test to the tagged_test std::map if true.  Returns true if a tag was present.
    bool CheckAndTagTest(
        const Teuchos::ParameterList& p,
        const Teuchos::RCP<NOX::StatusTest::Generic>& test,
        std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
        tagged_tests) const;
  private:
    //! Throws formated error
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  private:
    /*! \brief Reference to the base class NOX::StatusTest::Factory

        Used for a direct call instead of the member function to prevent too
        many calls of the factory constructor.
     */
    Teuchos::RCP<const NOX::StatusTest::Factory> noxfactory_;

  };  // class Factory
  /*! \brief Nonmember helper function for the NOX::NLN::StatusTest::Factory.

  \relates NOX::NLN::StatusTest::Factory

  */
  Teuchos::RCP<NOX::StatusTest::Generic> BuildOuterStatusTests(Teuchos::ParameterList& p,
      const NOX::Utils& utils, std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
      tagged_tests=0);

} // namespace StatusTest
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_STATUSTEST_FACTORY_H_ */
