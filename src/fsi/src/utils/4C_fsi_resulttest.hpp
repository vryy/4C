/*----------------------------------------------------------------------*/
/*! \file

\brief testing of fsi specific calculation results

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_RESULTTEST_HPP
#define FOUR_C_FSI_RESULTTEST_HPP

#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace FSI
{
  // forward declarations
  class Monolithic;
  class MonolithicNoNOX;

  /*!
    \brief FSI specific result test class

    Here, additional result tests of quantities that do not belong to a single
    field are tested. Basically, this should work for monolithic and partitioned
    schemes as well.

    Feel free to add further testing functionalities!

    \sa ResultTest
    \author mayr.mt
    \date 11/2012
  */
  class FSIResultTest : public Core::UTILS::ResultTest
  {
   public:
    //! constructor for standard FSI
    FSIResultTest(Teuchos::RCP<FSI::Monolithic>&
                      fsi,  ///< monolithic solver object that was used for the simulation
        const Teuchos::ParameterList& fsidyn  ///< FSI parameter list from input file
    );

    //! constructor for FSI implementation without NOX
    FSIResultTest(
        Teuchos::RCP<FSI::MonolithicNoNOX>
            fsi,  ///< monolithic solver object without NOX that was used for the simulation
        const Teuchos::ParameterList& fsidyn  ///< FSI parameter list from input file
    );

    //! \brief fsi version of nodal value tests
    //!
    //! Possible position flags are "lambdax", "lambday", "lambdaz"
    void test_node(Input::LineDefinition& res,  ///< the line definition from input file
        int& nerr,                              ///< number of tests with errors
        int& test_count                         ///< number of tests performed
        ) override;

    //! \brief fsi version of element value tests
    void test_element(Input::LineDefinition& res,  ///< the line definition from input file
        int& nerr,                                 ///< number of tests with errors
        int& test_count                            ///< number of tests performed
        ) override;


    //! \brief fsi version of special tests
    void test_special(Input::LineDefinition& res,  ///< the line definition from input file
        int& nerr,                                 ///< number of tests with errors
        int& test_count                            ///< number of tests performed
        ) override;

   private:
    //! slave discretisation
    Teuchos::RCP<Core::FE::Discretization> slavedisc_;

    //! Lagrange multiplier living on the slave discretization
    Teuchos::RCP<Epetra_Vector> fsilambda_;

    //! the monolithic solver object itself
    Teuchos::RCP<FSI::Monolithic> fsi_;
  };
}  // namespace FSI
FOUR_C_NAMESPACE_CLOSE

#endif
