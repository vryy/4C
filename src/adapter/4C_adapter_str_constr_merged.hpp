// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP
#define FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Solid
{
  namespace Aux
  {
    class MapExtractor;
  }
}  // namespace Solid

namespace Core::LinAlg
{
  class MapExtractor;
}

/*----------------------------------------------------------------------*/
namespace Adapter
{
  /*====================================================================*/
  /*!
   * \brief Adapter to constrained structural time integration.
   * This class wraps one of the standard adapters for structural time
   * integration. The results are modified and/or merged to account for the
   * additional degrees of freedom of the lagrange multipliers.
   *
   * \date 11/08
   */
  class StructureConstrMerged : public FSIStructureWrapper
  {
   public:
    /// Constructor
    StructureConstrMerged(Teuchos::RCP<Structure> stru);

    /// setup this object
    void setup() override;

    /// initial guess of Newton's method
    Teuchos::RCP<const Core::LinAlg::Vector<double>> initial_guess() override;

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Core::LinAlg::Vector<double>> rhs() override;

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Core::LinAlg::Vector<double>> dispnp() const override;

    /// known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Core::LinAlg::Vector<double>> dispn() const override;

    /*! \brief known velocity at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    Teuchos::RCP<const Core::LinAlg::Vector<double>> veln() const override;

    /*! known acceleration at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    Teuchos::RCP<const Core::LinAlg::Vector<double>> accn() const override;

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> dof_row_map() override;

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    void apply_interface_forces_temporary_deprecated(
        Teuchos::RCP<Core::LinAlg::Vector<double>> iforce) override;

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override;

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

    /// update displacement and evaluate elements
    void evaluate(Teuchos::RCP<const Core::LinAlg::Vector<double>>
            dispstepinc  ///< solution increment between time step n and n+1
        ) override;

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() override
    {
      return structure_->get_dbc_map_extractor();
    };

    /// domain map of system matrix
    const Epetra_Map& domain_map() const override;

    /// are there any algebraic constraints?
    bool have_constraint() override { return structure_->have_constraint(); };

    /// Return bool indicating if constraints are defined
    Teuchos::RCP<CONSTRAINTS::ConstrManager> get_constraint_manager() override
    {
      return structure_->get_constraint_manager();
    };

    Inpar::Solid::StcScale get_stc_algo() override { return structure_->get_stc_algo(); };

    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_stc_mat() override
    {
      FOUR_C_THROW("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->get_stc_mat();
    }

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    void update_iter_incr_constr(
        Teuchos::RCP<Core::LinAlg::Vector<double>> lagrincr  ///< Lagrange multiplier increment
        ) override
    {
      structure_->update_iter_incr_constr(lagrincr);
    }

    /// @name Apply interface forces

    //@}

    /// Integrate from t1 to t2
    int integrate() override { return structure_->integrate(); }

   private:
    /// the constraint map setup for full <-> stuct+constr transition
    Teuchos::RCP<Core::LinAlg::MapExtractor> conmerger_;

    /// the complete non-overlapping degree of freedom row map for structure and lagrange
    /// multipliers
    Teuchos::RCP<Epetra_Map> dofrowmap_;

    /// @name local copies of input parameters
    //{@
    Teuchos::RCP<Core::FE::Discretization> discret_;  ///< the discretization
    Teuchos::RCP<Teuchos::ParameterList>
        sdynparams_;  ///< dynamic control flags ... used, but could/should be circumvented
    Teuchos::RCP<Teuchos::ParameterList> xparams_;         ///< eXtra input parameters
    Teuchos::RCP<Core::LinAlg::Solver> solver_;            ///< the linear solver
    Teuchos::RCP<Core::IO::DiscretizationWriter> output_;  ///< the output writer

    //@}

    /// flag indicating if setup() was called
    bool issetup_;

  };  // class StructureConstrained

}  // namespace Adapter

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
