/*---------------------------------------------------------------------*/
/*! \file
\brief a derived strategy handling the Lagrange multiplier based TSI contact

\level 3


*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_TSI_LAGRANGE_STRATEGY_H
#define CONTACT_TSI_LAGRANGE_STRATEGY_H

#include "contact_lagrange_strategy.H"
#include "../drt_lib/drt_dserror.H"
#include "contact_defines.H"
#include <Epetra_Operator.h>

// forward declarations
namespace LINALG
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace LINALG

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace CONTACT
{
  // forward declaration
  // class WearInterface;
  /*!
   \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

   This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
   For a more general documentation of the involved functions refer to CoAbstract Strategy.

   */
  class CoTSILagrangeStrategy : public CoLagrangeStrategy
  {
   public:
    /*!
      \brief Standard Constructor

     */
    CoTSILagrangeStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        Teuchos::RCP<const Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
      \brief Destructor

     */
    virtual ~CoTSILagrangeStrategy(){};

    //! @name Access methods

    //@}

    //! @name Evaluation methods

    /*!
      \brief Set current state
      ...Standard Implementation in Abstract Strategy:
      All interfaces are called to set the current deformation state
      (u, xspatial) in their nodes. Additionally, the new contact
      element areas are computed.

      ... + Overloaded Implementation in Poro Lagrange Strategy
      Set structure & fluid velocity and lagrangean multiplier to Contact nodes data container!!!

      \param statetype (in): enumerator defining which quantity to set (see mortar_interface.H for
      an overview) \param vec (in): current global state of the quantity defined by statetype
     */
    virtual void SetState(const enum MORTAR::StateType& statetype, const Epetra_Vector& vec);

    // Overload CONTACT::CoAbstractStrategy::ApplyForceStiffCmt as this is called in the structure
    // --> to early for monolithically coupled algorithms!
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor)
    {
      // structure single-field predictors (e.g.TangDis) may evaluate the structural contact part
      if (predictor) CoAbstractStrategy::ApplyForceStiffCmt(dis, kt, f, step, iter, predictor);
    }

    /*!
      \brief Apply thermo-contact to matrix blocks

      In the TSI case, the contact terms are applied to the global system here.
      The "usual" place, i.e. the
      Evaluate(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
      in the Contact_lagrange_strategy is overloaded to do nothing, since
      in a coupled problem, we need to be very careful, when condensating
      the Lagrange multipliers.

     */
    virtual void Evaluate(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector>& combined_RHS, Teuchos::RCP<ADAPTER::Coupling> coupST,
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> temp);

    /*!
    \brief Overload CONTACT::CoLagrangeStrategy::Recover as this is called in the structure

    --> not enough information available for monolithically coupled algorithms!
    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override { return; };

    virtual void RecoverCoupled(Teuchos::RCP<Epetra_Vector> sinc,  /// displacement  increment
        Teuchos::RCP<Epetra_Vector> tinc,                          /// thermal  increment
        Teuchos::RCP<ADAPTER::Coupling> coupST);

    void StoreNodalQuantities(
        MORTAR::StrategyBase::QuantityType type, Teuchos::RCP<ADAPTER::Coupling> coupST);

    /*!
     \brief Update contact at end of time step

     \param dis (in):  current displacements (-> old displacements)

     */
    virtual void Update(Teuchos::RCP<const Epetra_Vector> dis);

    /*!
     \brief Set time integration parameter from Thermo time integration

     */
    void SetAlphafThermo(const Teuchos::ParameterList& tdyn);


    /*!
    \brief Perform a write restart

    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Hence, a portion of the restart has to be
    performed on the level of the contact algorithm, for short: here's the right place.

    */
    virtual void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const;

    /*!
    \brief Perform a write restart

    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Hence, all the restart action has to be
    performed on the level of the contact algorithm, for short: here's the right place.

    */
    void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis,
        Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    void SetCoupling(Teuchos::RCP<ADAPTER::Coupling> coupST) { coupST_ = coupST; };

    //@}

    // residual and increment norms
    double mech_contact_res_;
    double mech_contact_incr_;
    double thr_contact_incr_;

   protected:
    // don't want = operator and cctor
    CoTSILagrangeStrategy operator=(const CoTSILagrangeStrategy& old) = delete;
    CoTSILagrangeStrategy(const CoTSILagrangeStrategy& old) = delete;

    // time integration
    double tsi_alpha_;

    Teuchos::RCP<Epetra_Vector>
        fscn_;  // structural contact forces of last time step (needed for time integration)
    Teuchos::RCP<Epetra_Vector>
        ftcn_;  // thermal    contact forces of last time step (needed for time integration)
    Teuchos::RCP<Epetra_Vector>
        ftcnp_;  // thermal   contact forces of this time step (needed for time integration)

    Teuchos::RCP<Epetra_Vector> z_thr_;  // current vector of Thermo-Lagrange multipliers at t_n+1
    Teuchos::RCP<Epetra_Map> thr_act_dofs_;  // active thermo dofs
    Teuchos::RCP<Epetra_Map> thr_s_dofs_;    // slave thermo dofs

    Teuchos::RCP<LINALG::SparseMatrix> dinvA_;  // dinv on active displacement dofs (for recovery)
    Teuchos::RCP<LINALG::SparseMatrix> dinvAthr_;  // dinv on active thermal dofs (for recovery)
    // recovery of contact LM
    Teuchos::RCP<LINALG::SparseMatrix>
        kss_a_;  // Part of structure-stiffness (kss) that corresponds to active slave rows
    Teuchos::RCP<LINALG::SparseMatrix>
        kst_a_;  // Part of coupling-stiffness  (kst) that corresponds to active slave rows
    Teuchos::RCP<Epetra_Vector>
        rs_a_;  // Part of structural residual that corresponds to active slave rows

    // recovery of thermal LM
    Teuchos::RCP<LINALG::SparseMatrix>
        ktt_a_;  // Part of structure-stiffness (ktt) that corresponds to active slave rows
    Teuchos::RCP<LINALG::SparseMatrix>
        kts_a_;  // Part of coupling-stiffness  (kts) that corresponds to active slave rows
    Teuchos::RCP<Epetra_Vector>
        rt_a_;  // Part of structural residual that corresponds to active slave rows

    // pointer to TSI coupling object
    Teuchos::RCP<ADAPTER::Coupling> coupST_;
  };  // class CoTSILagrangeStrategy

  namespace UTILS
  {
    //! @name little helpers
    void AddVector(Epetra_Vector& src, Epetra_Vector& dst);
  }  // namespace UTILS
}  // namespace CONTACT


#endif  // #ifndef CONTACT_TSI_LAGRANGE_STRATEGY_H
