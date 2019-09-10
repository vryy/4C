/*---------------------------------------------------------------------*/
/*! \file

\brief Result test for the adjoint fluid equations

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_FLUIDADJOINTRESULTTEST_H_
#define TOPOPT_FLUIDADJOINTRESULTTEST_H_


#include "../drt_lib/drt_resulttest.H"


class Epetra_Vector;
namespace DRT
{
  class Discretization;
}

namespace TOPOPT
{
  namespace ADJOINT
  {
    // forward declarations
    class ImplicitTimeInt;


    /*!
    \brief Fluid adjoint specific result test class

    \author winklmaier
    */
    class FluidAdjointResultTest : public DRT::ResultTest
    {
     public:
      /*!
      \brief constructor
      */
      FluidAdjointResultTest(const ImplicitTimeInt& adjointfluid);


      /// our version of nodal value tests
      /*!
      Possible position flags are "velx", "vely", "velz" and
      "pressure". With the obvious meaning.
      */
      virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

     private:
      /// pointer to fluid discretization
      Teuchos::RCP<const DRT::Discretization> fluiddis_;
      /// pointer to unknown vector with values
      const Teuchos::RCP<const Epetra_Vector> mysol_;
    };
  }  // namespace ADJOINT
}  // namespace TOPOPT

#endif /* TOPOPT_FLUIDADJOINTRESULTTEST_H_ */
