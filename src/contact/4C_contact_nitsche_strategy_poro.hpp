/*---------------------------------------------------------------------*/
/*! \file
\brief Nitsche poro contact solving strategy

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_PORO_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_PORO_HPP

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
  class NitscheStrategyPoro : public NitscheStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategyPoro(const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::Interface>> interface,
        int dim, Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(
              dof_row_map, NodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof),
          no_penetration_(params.get<bool>("CONTACTNOPEN"))
    {
    }

    //! Shared data constructor
    NitscheStrategyPoro(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof)
        : NitscheStrategy(data_ptr, dof_row_map, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof),
          no_penetration_(params.get<bool>("CONTACTNOPEN"))
    {
    }

    void apply_force_stiff_cmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor) override;

    //  void Integrate(CONTACT::ParamsInterface& cparams);
    void set_state(const enum Mortar::StateType& statename, const Epetra_Vector& vec) override;

    void SetParentState(const enum Mortar::StateType& statename, const Epetra_Vector& vec) override;

    Teuchos::RCP<const Epetra_Vector> get_rhs_block_ptr(
        const enum CONTACT::VecBlockType& bp) const override;

    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> GetMatrixBlockPtr(
        const enum CONTACT::MatBlockType& bp) const;


    // Flag for Poro No Penetration Condition
    bool has_poro_no_penetration() const override { return no_penetration_; }

    // don't want = operator and cctor
    NitscheStrategyPoro operator=(const NitscheStrategyPoro& old) = delete;
    NitscheStrategyPoro(const NitscheStrategyPoro& old) = delete;

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

    bool no_penetration_;

    Teuchos::RCP<Epetra_FEVector> fp_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kpp_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kpd_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kdp_;
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
