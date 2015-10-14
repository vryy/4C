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
#include <boost/algorithm/string/predicate.hpp>

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
  enum QuantityType {
    quantity_unknown,
    // check structural quantities
    quantity_structure,
    // check (semi-smooth) contact quantities
    quantity_contact,
    // check meshtying quantities
    quantity_meshtying,
    // check windkessel quantities
    quantity_windkessel,
    // check Lagrange/penalty enforced constraint quantites
    quantity_lag_pen_constraint,
    // check semi-smooth contact
    quantity_plasticity,
    // check pressure dofs
    quantity_pressure,
    // check eas dofs
    quantity_eas
  };

  /// Map quantity name to std::string
  static inline const std::string QuantityType2String(const enum QuantityType type)
  {
    switch (type)
    {
    case quantity_structure:
      return "STRUCTURE";
      break;
    case quantity_contact:
      return "CONTACT";
      break;
    case quantity_meshtying:
      return "MESHTYING";
      break;
    case quantity_windkessel:
      return "WINDKESSEL";
      break;
    case quantity_lag_pen_constraint:
      return "CONSTRAINT";
      break;
    case quantity_plasticity:
      return "PLASTICITY";
      break;
    case quantity_pressure:
      return "PRESSURE";
      break;
    case quantity_eas:
      return "EAS";
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

    // -- case insensitive comparison --
    if (boost::iequals(name,"structure"))
      type = quantity_structure;
    else if (boost::iequals(name,"contact"))
      type = quantity_contact;
    else if (boost::iequals(name,"meshtying"))
      type = quantity_meshtying;
    else if (boost::iequals(name,"windkessel"))
      type = quantity_windkessel;
    else if (boost::iequals(name,"constraint"))
      type = quantity_lag_pen_constraint;
    else if (boost::iequals(name,"plasticity"))
      type = quantity_plasticity;
    else if (boost::iequals(name,"pressure"))
      type = quantity_pressure;
    else if (boost::iequals(name,"eas"))
      type = quantity_eas;

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
    /*! \brief New implementation of the \c BuildNormFTest function.
     *
     *  The underlying Status Test is capable of a variety of quantities at the same
     *  time and combines the RelativeNormF and NormF test.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildNormFTest(
        Teuchos::ParameterList& p, const NOX::Utils& u, const bool& relativeNormF) const;
    Teuchos::RCP<NOX::StatusTest::Generic> BuildNormFTest(
            Teuchos::ParameterList& p, const NOX::Utils& u) const;

    /*! \brief New implementation of the \c BuildNormUpdateTest function.
     *
     *  The underlying Status Test is capable of a variety of quantities at the same
     *  time.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildNormUpdateTest(
        Teuchos::ParameterList& p, const NOX::Utils& u) const;

    /*! \brief New implementation of the \c BuildNormWRMSTest function.
     *
     *  The underlying Status Test is capable of a variety of quantities at the same
     *  time.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildNormWRMSTest(
        Teuchos::ParameterList& p, const NOX::Utils& u) const;

    /*! \brief Simple active set test
     *
     *  This active set test checks the active set status of the corresponding
     *  quantity. Basically we check if the active set did not change.
     *
     *  OPTIONAL: Based on the implementation of the internal get function, also
     *  a cycling of the active set can be detected and printed to the screen.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildActiveSetTest(
        Teuchos::ParameterList& p, const NOX::Utils& u) const;

    /*! \brief Derived buildComboTest.
     *
     *  Restatement of the base class function, because the internal recursive call
     *  has to be redefined. Note that the base class input parameter
     *  <tt>Number of Tests<\tt> is not used in this version.
     */
    Teuchos::RCP<NOX::StatusTest::Generic> BuildComboTest(
        Teuchos::ParameterList& p, const NOX::Utils& u,
        std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
        tagged_tests) const;

    //! Checks if a tag is present in the parameter list and adds the test to the tagged_test std::map if true.  Returns true if a tag was present.
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
     *
     *  Used for a direct call instead of the member function to prevent too
     *  many calls of the factory constructor.
     */
    Teuchos::RCP<const NOX::StatusTest::Factory> noxfactory_;

  };  // class Factory
  /*! \brief Nonmember helper function for the NOX::NLN::StatusTest::Factory.
   *
   *  \relates NOX::NLN::StatusTest::Factory
   *
   */
  Teuchos::RCP<NOX::StatusTest::Generic> BuildOuterStatusTests(Teuchos::ParameterList& p,
      const NOX::Utils& utils, std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
      tagged_tests=0);

} // namespace StatusTest
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_STATUSTEST_FACTORY_H_ */
