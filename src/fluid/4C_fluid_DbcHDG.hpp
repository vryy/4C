/*-----------------------------------------------------------*/
/*! \file

\brief contains utils functions for Dirichlet Boundary Conditions of HDG discretizations

\level 2

*/

#ifndef FOUR_C_FLUID_DBCHDG_HPP
#define FOUR_C_FLUID_DBCHDG_HPP


#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_discretization_utils.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::UTILS
{
  class Dbc;
  // class DbcInfo;
}  // namespace Discret::UTILS

namespace FLD
{
  namespace UTILS
  {
    /** \brief Specialized Dbc evaluation class for HDG discretizations
     *
     *  \author hiermeier \date 10/16 */
    class DbcHdgFluid : public Discret::UTILS::Dbc
    {
     public:
      /// constructor
      DbcHdgFluid(){};

     protected:
      void read_dirichlet_condition(const Teuchos::ParameterList& params,
          const Discret::Discretization& discret, const Core::Conditions::Condition& cond,
          double time, DbcInfo& info, const Teuchos::RCP<std::set<int>>* dbcgids,
          int hierarchical_order) const override;

      void read_dirichlet_condition(const Teuchos::ParameterList& params,
          const Discret::DiscretizationFaces& discret, const Core::Conditions::Condition& cond,
          double time, DbcInfo& info, const Teuchos::RCP<std::set<int>>* dbcgids,
          int hierarchical_order) const;

      void do_dirichlet_condition(const Teuchos::ParameterList& params,
          const Discret::Discretization& discret, const Core::Conditions::Condition& cond,
          double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_IntVector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const override;

      void do_dirichlet_condition(const Teuchos::ParameterList& params,
          const Discret::DiscretizationFaces& discret, const Core::Conditions::Condition& cond,
          double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_IntVector& toggle) const;
    };  // class DbcHDG_Fluid
  }     // namespace UTILS

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
