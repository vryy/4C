/*!----------------------------------------------------------------------
\file contact_tsi_lagrange_strategy.H
\maintainer Alexander Seitz

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/
#ifndef CONTACT_TSI_LAGRANGE_STRATEGY_H
#define CONTACT_TSI_LAGRANGE_STRATEGY_H

#include "contact_lagrange_strategy.H"
#include "../drt_lib/drt_dserror.H"
#include "contact_defines.H"
#include <Epetra_Operator.h>

namespace LINALG { class SparseMatrix;class BlockSparseMatrixBase;}

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
  }
}

namespace CONTACT
{

// forward declaration
//class WearInterface;
/*!
 \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 \author popp (popp@lnm.mw.tum.de)
 */
class CoTSILagrangeStrategy : public CoLagrangeStrategy
{
public:

  /*!
    \brief Standard Constructor

   */
  CoTSILagrangeStrategy(
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      Teuchos::ParameterList params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
      int dim,
      Teuchos::RCP<Epetra_Comm> comm,
      double alphaf,
      int maxdof);

  /*!
    \brief Destructor

   */
  virtual ~CoTSILagrangeStrategy() {};

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

    \param statename (in): std::string defining which quantity to set (either "displacement", "olddisplacement", "svelocity", "fvelocity", "lm")
    \param vec (in): current global state of the quantity defined by statename
   */
  void SetState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec);

  //Overload CONTACT::CoAbstractStrategy::ApplyForceStiffCmt as this is called in the structure --> to early for monolithically coupled algorithms!
  virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
      Teuchos::RCP<LINALG::SparseOperator>& kt,
      Teuchos::RCP<Epetra_Vector>& f,
      const int step,
      const int iter,
      bool predictor)
  {
    // structure single-field predictors (e.g.TangDis) may evaluate the structural contact part
    if (predictor) CoAbstractStrategy::ApplyForceStiffCmt(dis,kt,f,step,iter,predictor);
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
  virtual void Evaluate(
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
      Teuchos::RCP<Epetra_Vector>& combined_RHS,
      Teuchos::RCP<ADAPTER::Coupling> coupST,
      Teuchos::RCP<Epetra_Vector> dis,
      Teuchos::RCP<Epetra_Vector> temp,
      Teuchos::RCP<const LINALG::MapExtractor> str_dbc,
      Teuchos::RCP<const LINALG::MapExtractor> thr_dbc
      );

  //Overload CONTACT::CoLagrangeStrategy::Recover as this is called in the structure --> no enought information available for monolithically coupled algorithms!
   virtual void Recover(Teuchos::RCP<Epetra_Vector> disi)
   {};

  virtual void RecoverCoupled(
      Teuchos::RCP<Epetra_Vector> sinc,/// displacement  increment
      Teuchos::RCP<Epetra_Vector> tinc,/// thermal  increment
      Teuchos::RCP<ADAPTER::Coupling> coupST);

  //! @name little helpers
  inline void AddVector(Epetra_Vector& src,Epetra_Vector& dst);

  void StoreNodalQuantities(
      MORTAR::StrategyBase::QuantityType type,Teuchos::RCP<ADAPTER::Coupling> coupST);

  /*!
   \brief Update contact at end of time step

   \param dis (in):  current displacements (-> old displacements)

   */
  void Update(Teuchos::RCP<Epetra_Vector> dis,Teuchos::RCP<ADAPTER::Coupling> coupST);

  /*!
   \brief Update contact at end of time step
          don't do it in the standard call from structural time integration

   */
  void Update(Teuchos::RCP<Epetra_Vector> dis){}

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
  void DoWriteRestart(std::map<std::string,Teuchos::RCP<Epetra_Vector> >& restart_vectors,
      bool forcedrestart = false);

  /*!
  \brief Perform a write restart

  A write restart is initiated by the contact manager. However, the manager has no
  direct access to the nodal quantities. Hence, all the restart action has to be
  performed on the level of the contact algorithm, for short: here's the right place.

  */
  void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis);

  void SetCoupling(Teuchos::RCP<ADAPTER::Coupling> coupST){coupST_=coupST;};

  //@}

  // residual and increment norms
  double mech_contact_res_;
  double mech_contact_incr_;
  double thr_contact_incr_;

protected:

  // don't want = operator and cctor
  CoTSILagrangeStrategy operator = (const CoTSILagrangeStrategy& old);
  CoTSILagrangeStrategy(const CoTSILagrangeStrategy& old);

  //time integration
  double tsi_alpha_;

  Teuchos::RCP<Epetra_Vector> fscn_; // structural contact forces of last time step (needed for time integration)
  Teuchos::RCP<Epetra_Vector> ftcn_; // thermal    contact forces of last time step (needed for time integration)
  Teuchos::RCP<Epetra_Vector> ftcnp_; // thermal   contact forces of this time step (needed for time integration)

  Teuchos::RCP<Epetra_Vector> z_thr_;                 // current vector of Thermo-Lagrange multipliers at t_n+1
  Teuchos::RCP<Epetra_Map>    thr_act_dofs_;          // active thermo dofs
  Teuchos::RCP<Epetra_Map>    thr_s_dofs_;            // slave thermo dofs

  Teuchos::RCP<LINALG::SparseMatrix> dinvA_;          // dinv on active displacement dofs (for recovery)
  Teuchos::RCP<LINALG::SparseMatrix> dinvAthr_;       // dinv on active thermal dofs (for recovery)
  // recovery of contact LM
  Teuchos::RCP<LINALG::SparseMatrix> kss_a_;          // Part of structure-stiffness (kss) that corresponds to active slave rows
  Teuchos::RCP<LINALG::SparseMatrix> kst_a_;          // Part of coupling-stiffness  (kst) that corresponds to active slave rows
  Teuchos::RCP<Epetra_Vector>        rs_a_;           // Part of structural residual that corresponds to active slave rows

  // recovery of thermal LM
  Teuchos::RCP<LINALG::SparseMatrix> ktt_a_;          // Part of structure-stiffness (ktt) that corresponds to active slave rows
  Teuchos::RCP<LINALG::SparseMatrix> kts_a_;          // Part of coupling-stiffness  (kts) that corresponds to active slave rows
  Teuchos::RCP<Epetra_Vector>        rt_a_;           // Part of structural residual that corresponds to active slave rows

  // pointer to TSI coupling object
  Teuchos::RCP<ADAPTER::Coupling> coupST_;
}; // class CoTSILagrangeStrategy
} // namespace CONTACT


#endif  // #ifndef CONTACT_TSI_LAGRANGE_STRATEGY_H
