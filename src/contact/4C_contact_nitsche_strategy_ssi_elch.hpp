/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy including electrochemistry

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_SSI_ELCH_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_SSI_ELCH_HPP

#include "4C_config.hpp"

#include "4C_contact_nitsche_strategy_ssi.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
   * @brief Contact solving strategy with Nitsche's method.
   *
   * This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   * For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.
   */
  class NitscheStrategySsiElch : public NitscheStrategySsi
  {
   public:
    //! Shared data constructor
    NitscheStrategySsiElch(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof)
        : NitscheStrategySsi(data_ptr, dof_row_map, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof)
    {
    }

    void apply_force_stiff_cmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int step, const int iter, bool predictor) override
    {
      FOUR_C_THROW("not implemented");
    }

    void integrate(const CONTACT::ParamsInterface& cparams) override;

    //! don't want = operator
    NitscheStrategySsiElch operator=(const NitscheStrategySsiElch& old) = delete;

    //! don't want copy constructor
    NitscheStrategySsiElch(const NitscheStrategySsiElch& old) = delete;
  };
}  // namespace CONTACT
FOUR_C_NAMESPACE_CLOSE

#endif
