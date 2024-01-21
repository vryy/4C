/*----------------------------------------------------------------------*/
/*! \file

\brief testing of thermo calculation results

\level 1

*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 08/09 |
 *----------------------------------------------------------------------*/
#ifndef BACI_THERMO_RESULTTEST_H
#define BACI_THERMO_RESULTTEST_H


/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "baci_config.H"

#include "baci_global_data.H"
#include "baci_io.H"
#include "baci_lib_resulttest.H"
#include "baci_thermo_timint.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                     dano 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  //!
  //! \brief Thermo specific result test class
  //!
  //! \author cd
  class ResultTest : public DRT::ResultTest
  {
   public:
    //! Constructor for time integrators of general kind
    //! \author bborn \date 06/08
    ResultTest(TimInt& tintegrator);

    //! \brief thermo version of nodal value tests
    //!
    //! Possible position flags are "temp",
    //!                             "rate",
    void TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! our discretisation
    Teuchos::RCP<DRT::Discretization> thrdisc_;
    // our solution
    //! global temperature DOFs
    Teuchos::RCP<Epetra_Vector> temp_;
    //! global temperature rate DOFs
    Teuchos::RCP<Epetra_Vector> rate_;
    //! global temperature DOFs
    Teuchos::RCP<Epetra_Vector> flux_;
    //! NOTE: these have to be present explicitly
    //! as they are not part of the problem instance like in fluid3

  };  // DRT::ResultTest

}  // namespace THR


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // THERMO_RESULTTEST_H
