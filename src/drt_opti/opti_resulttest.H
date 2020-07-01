/*---------------------------------------------------------------------*/
/*! \file

\brief Result test for optimization algorithms


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef OPTI_RESULTTEST_H_
#define OPTI_RESULTTEST_H_


#include "../drt_lib/drt_resulttest.H"


class Epetra_SerialDenseVector;
class Epetra_MultiVector;
class Epetra_Vector;

namespace DRT
{
  class Discretization;
}

namespace OPTI
{
  class GCMMA;

  /*!
  \brief Optimization specific result test class

  \author winklmaier
  */
  class OptiResultTest : public DRT::ResultTest
  {
   public:
    /*!
    \brief constructor
    */
    OptiResultTest(const OPTI::GCMMA& optimizer);


    /// our version of nodal value tests
    /*!
    Possible position flags are "velx", "vely", "velz" and
    "pressure". With the obvious meaning.
    */
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

    virtual void TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   private:
    /// pointer to fluid discretization
    const Teuchos::RCP<const DRT::Discretization> optidis_;

    /// optimization variable
    const Teuchos::RCP<const Epetra_Vector> x_;
    /// objective value
    const double obj_;
    /// spatial derivative of objective value
    const Teuchos::RCP<const Epetra_Vector> obj_deriv_;
    /// constraint value(s)
    const Teuchos::RCP<const Epetra_SerialDenseVector> constr_;
    /// spatial derivative of constraint value(s)
    const Teuchos::RCP<const Epetra_MultiVector> constr_deriv_;
  };
}  // namespace OPTI


#endif /* OPTI_RESULTTEST_H_ */
