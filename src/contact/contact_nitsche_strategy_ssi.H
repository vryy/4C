/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_STRATEGY_SSI_H
#define CONTACT_NITSCHE_STRATEGY_SSI_H

#include "contact_nitsche_strategy.H"

#include <utility>

namespace CONTACT
{
  /*!
   * @brief Contact solving strategy with Nitsche's method.
   *
   * This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
   * For a more general documentation of the involved functions refer to CoAbstract Strategy.
   */
  class CoNitscheStrategySsi : public CoNitscheStrategy
  {
   public:
    //! Standard constructor
    CoNitscheStrategySsi(const Epetra_Map* dofRowMap, const Epetra_Map* nodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        const Teuchos::RCP<Epetra_Comm>& comm, double alphaf, int maxdof)
        : CoNitscheStrategy(
              dofRowMap, nodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof)
    {
    }

    //! Shared data constructor
    CoNitscheStrategySsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dofRowMap, const Epetra_Map* nodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof)
        : CoNitscheStrategy(data_ptr, dofRowMap, nodeRowMap, params, std::move(interface), dim,
              comm, alphaf, maxdof)
    {
    }

    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor) override
    {
      dserror("not implemented");
    }

    /*!
     * @note currently there is no support for friction; furthermore we can not call
     * 'UpdateTraceIneqEtimates' at the moment since at this stage we do not have the gp
     * concentrations -> so far no adaptive nitsche penalty possible!
     */
    void EvaluateReferenceState() override;

    /*!
     * @brief initialize TraceHE parameter of the mortar elements
     *
     * @note: this approximation is done here since UpdateTraceIneqEtimates() can not be called at
     * the moment
     */
    void InitTraceHE();

    void Integrate(const CONTACT::ParamsInterface& cparams) override;

    void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    void SetParentState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec) override;

    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    /*!
     * @brief get the pointer to the matrix block
     *
     * @param[in] bt   block type of requested matrix block
     * @return pointer to matrix block
     */
    Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) const;

    //! don't want = operator
    CoNitscheStrategySsi operator=(const CoNitscheStrategySsi& old) = delete;
    //! don't want copy constructor
    CoNitscheStrategySsi(const CoNitscheStrategySsi& old) = delete;

   protected:
    Teuchos::RCP<Epetra_FEVector> SetupRhsBlockVec(
        const enum DRT::UTILS::VecBlockType& bt) const override;

    Teuchos::RCP<LINALG::SparseMatrix> SetupMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) override;

    void CompleteMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt, Teuchos::RCP<LINALG::SparseMatrix> kc) override;

    //! current scalar state vector
    Teuchos::RCP<Epetra_Vector> curr_state_scalar_;
    //! ScaTra residual
    Teuchos::RCP<Epetra_FEVector> fs_;
    //! linearization of ScaTra residual w.r.t ScaTra dofs
    Teuchos::RCP<LINALG::SparseMatrix> kss_;
    //! linearization of ScaTra residual w.r.t displacement dofs
    Teuchos::RCP<LINALG::SparseMatrix> ksd_;
    //! linearization of displacement residual w.r.t ScaTra dofs
    Teuchos::RCP<LINALG::SparseMatrix> kds_;
  };
}  // namespace CONTACT
#endif  // CONTACT_NITSCHE_STRATEGY_SSI_H
