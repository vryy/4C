/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_TSI_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_TSI_HPP

#include "4C_config.hpp"

#include "4C_contact_nitsche_strategy.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Coupling;
}

namespace CONTACT
{
  /*!
   \brief Contact solving strategy with Nitsche's method.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   */
  class NitscheStrategyTsi : public NitscheStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategyTsi(const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(
              dof_row_map, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          fix_redistribution_(true)
    {
    }

    //! Shared data constructor
    NitscheStrategyTsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(data_ptr, dof_row_map, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof),
          fix_redistribution_(true)
    {
    }

    void apply_force_stiff_cmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor) override
    {
      FOUR_C_THROW("not implemented");
    }

    void integrate(const CONTACT::ParamsInterface& cparams) override;

    /*! \brief Setup this strategy object (maps, vectors, etc.)

     derived from contact abstract strategy.
     The Nitsche strategy does not have
      */
    void setup(bool redistributed, bool init) override;

    void update_trace_ineq_etimates() override;
    void set_state(const enum Mortar::StateType& statename, const Epetra_Vector& vec) override;

    void set_parent_state(
        const enum Mortar::StateType& statename, const Epetra_Vector& vec) override;

    Teuchos::RCP<const Epetra_Vector> get_rhs_block_ptr(
        const enum CONTACT::VecBlockType& bt) const override;

    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const override;

    //! [derived]
    bool redistribute_contact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) override
    {
      if (fix_redistribution_) return false;
      return CONTACT::AbstractStrategy::redistribute_contact(dis, vel);
    }

    void enable_redistribution() { fix_redistribution_ = false; }

    // don't want = operator and cctor
    NitscheStrategyTsi operator=(const NitscheStrategyTsi& old) = delete;
    NitscheStrategyTsi(const NitscheStrategyTsi& old) = delete;

   protected:
    // create an appropriate vector for the RHS
    Teuchos::RCP<Epetra_FEVector> setup_rhs_block_vec(
        const enum CONTACT::VecBlockType& bt) const override;

    // create an appropriate matrix block
    Teuchos::RCP<Core::LinAlg::SparseMatrix> setup_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt) override;

    // complete matrix block with correct maps
    void complete_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc) override;

    // do not reditribute (during constructor phase)
    bool fix_redistribution_;

    Teuchos::RCP<Epetra_Vector> curr_state_temp_;

    Teuchos::RCP<Epetra_FEVector> ft_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ktt_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ktd_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kdt_;
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
