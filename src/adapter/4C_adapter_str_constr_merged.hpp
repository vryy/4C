/*----------------------------------------------------------------------*/
/*! \file

\brief Adapter Layer for Structures with Algebraic Constraints

\level 2

*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP
#define FOUR_C_ADAPTER_STR_CONSTR_MERGED_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"

#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace STR
{
  namespace Aux
  {
    class MapExtractor;
  }
}  // namespace STR

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
    Teuchos::RCP<const Epetra_Vector> initial_guess() override;

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS() override;

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp() const override;

    /// known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn() const override;

    /*! \brief known velocity at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    Teuchos::RCP<const Epetra_Vector> Veln() const override;

    /*! known acceleration at \f$t_{n}\f$
     *
     *  Lagrange multiplier does not have a time derivative. Though we need a map
     *  including the Lagrange multiplier, thus, we include it and set it to zero.
     */
    Teuchos::RCP<const Epetra_Vector> Accn() const override;

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> dof_row_map() override;

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    void apply_interface_forces_temporary_deprecated(Teuchos::RCP<Epetra_Vector> iforce) override;

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override;

    /// direct access to system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

    /// update displacement and evaluate elements
    void evaluate(Teuchos::RCP<const Epetra_Vector>
            dispstepinc  ///< solution increment between time step n and n+1
        ) override;

    //! Return MapExtractor for Dirichlet boundary conditions
    Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() override
    {
      return structure_->GetDBCMapExtractor();
    };

    /// domain map of system matrix
    const Epetra_Map& DomainMap() const override;

    /// are there any algebraic constraints?
    bool HaveConstraint() override { return structure_->HaveConstraint(); };

    /// Return bool indicating if constraints are defined
    Teuchos::RCP<CONSTRAINTS::ConstrManager> get_constraint_manager() override
    {
      return structure_->get_constraint_manager();
    };

    Inpar::STR::StcScale get_stc_algo() override { return structure_->get_stc_algo(); };

    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_stc_mat() override
    {
      FOUR_C_THROW("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->get_stc_mat();
    }

    //! Update iteration
    //! Add residual increment to Lagrange multipliers stored in Constraint manager
    void update_iter_incr_constr(
        Teuchos::RCP<Epetra_Vector> lagrincr  ///< Lagrange multiplier increment
        ) override
    {
      structure_->update_iter_incr_constr(lagrincr);
    }

    /// @name Apply interface forces

    //@}

    /// Integrate from t1 to t2
    int Integrate() override { return structure_->Integrate(); }

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
