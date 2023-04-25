/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy including electrochemistry

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_SSI_ELCH_H
#define CONTACT_NITSCHE_STRATEGY_SSI_ELCH_H

#include "contact_nitsche_strategy_ssi.H"

namespace CONTACT
{
  /*!
   * @brief Contact solving strategy with Nitsche's method.
   *
   * This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
   * For a more general documentation of the involved functions refer to CoAbstract Strategy.
   */
  class CoNitscheStrategySsiElch : public CoNitscheStrategySsi
  {
   public:
    //! Shared data constructor
    CoNitscheStrategySsiElch(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof)
        : CoNitscheStrategySsi(data_ptr, DofRowMap, NodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof)
    {
    }

    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor) override
    {
      dserror("not implemented");
    }

    void Integrate(const CONTACT::ParamsInterface& cparams) override;

    //! don't want = operator
    CoNitscheStrategySsiElch operator=(const CoNitscheStrategySsiElch& old) = delete;

    //! don't want copy constructor
    CoNitscheStrategySsiElch(const CoNitscheStrategySsiElch& old) = delete;
  };
}  // namespace CONTACT
#endif  // CONTACT_NITSCHE_STRATEGY_SSI_ELCH_H
