/*----------------------------------------------------------------------*/
/*! \file

\brief testing of fsi specific calculation results

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FSI_RESULTTEST_H
#define BACI_FSI_RESULTTEST_H

#include "baci_lib_resulttest.H"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

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
  class FSIResultTest : public DRT::ResultTest
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
    void TestNode(DRT::INPUT::LineDefinition& res,  ///< the line definition from input file
        int& nerr,                                  ///< number of tests with errors
        int& test_count                             ///< number of tests performed
        ) override;

    //! \brief fsi version of element value tests
    void TestElement(DRT::INPUT::LineDefinition& res,  ///< the line definition from input file
        int& nerr,                                     ///< number of tests with errors
        int& test_count                                ///< number of tests performed
        ) override;


    //! \brief fsi version of special tests
    void TestSpecial(DRT::INPUT::LineDefinition& res,  ///< the line definition from input file
        int& nerr,                                     ///< number of tests with errors
        int& test_count                                ///< number of tests performed
        ) override;

   private:
    //! slave discretisation
    Teuchos::RCP<DRT::Discretization> slavedisc_;

    //! Lagrange multiplier living on the slave discretization
    Teuchos::RCP<Epetra_Vector> fsilambda_;

    //! the monolithic solver object itself
    Teuchos::RCP<FSI::Monolithic> fsi_;
  };
}  // namespace FSI
BACI_NAMESPACE_CLOSE

#endif  // FSI_RESULTTEST_H
