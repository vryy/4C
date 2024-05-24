/*----------------------------------------------------------------------*/
/*! \file

\brief xfem based fluid result tests

\level 0

 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_XFLUID_RESULTTEST_HPP
#define FOUR_C_FLUID_XFLUID_RESULTTEST_HPP


#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace FLD
{
  // forward declarations
  class XFluid;
  class XFluidFluid;
  class XFluidFluid;

  /*!
    ResultTest class for XFluid
   */
  class XFluidResultTest : public CORE::UTILS::ResultTest
  {
   public:
    //! ctor for standard XFEM problems
    XFluidResultTest(const FLD::XFluid& xfluid);

    //! ctor for XFF-problems
    XFluidResultTest(const FLD::XFluidFluid& xfluid);

    /// our version of nodal value tests
    /*!
      Possible position flags are "velx", "vely", "velz" and
      "pressure". With the obvious meaning.
     */
    void test_node(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    /// nodal value test (one can specify discretization and corresponding solution here!)
    void test_node(INPUT::LineDefinition& res, int& nerr, int& test_count, int node,
        const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const Epetra_Vector>& velnp);

    /// XFEM discretization
    Teuchos::RCP<const DRT::Discretization> discret_;

    /// solution vector for XFEM discretization
    Teuchos::RCP<const Epetra_Vector> velnp_;

    /// optional additional discretization for the same field (fluid-fluid coupling)
    Teuchos::RCP<const DRT::Discretization> coupl_discret_;

    /// solution vector for additional coupling discretization
    Teuchos::RCP<const Epetra_Vector> coupl_velnp_;

    /// take care of node numbering off-by-one (will be removed soon)
    const bool node_from_zero_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
