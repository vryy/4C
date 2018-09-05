/*---------------------------------------------------------------------*/
/*!
\file contact_monocoupled_lagrange_strategy.H

\brief This class provides the functionality to use contact with Lagrangian
multipliers for monolithical coupled multifield problems!
Therefore ApplyForceStiffCmt() & Recover() are overloaded by this class and
do nothing, as they are called directly in the structure. To use the contact
the additional methods ApplyForceStiffCmtCoupled() & RecoverCoupled() have
to be called!

\level 3

\maintainer Christoph Ager

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_MONOCOUPLED_LAGRANGE_STRATEGY_H
#define CONTACT_MONOCOUPLED_LAGRANGE_STRATEGY_H

#include "contact_lagrange_strategy.H"
#include "../drt_lib/drt_dserror.H"

namespace ADAPTER
{
  class Coupling;
}

namespace CONTACT
{
  /*!
   \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

   This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
   For a more general documentation of the involved functions refer to CoAbstract Strategy.

   \author popp (popp@lnm.mw.tum.de)
   */
  class MonoCoupledLagrangeStrategy : public CoLagrangeStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    MonoCoupledLagrangeStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MonoCoupledLagrangeStrategy(){};

    // Overload CONTACT::CoAbstractStrategy::ApplyForceStiffCmt as this is called in the structure
    // --> to early for monolithically coupled algorithms!
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor)
    {
      if (has_to_evaluate_ && 0)
        dserror(
            "MonoCoupledLagrangeStrategy::You have to call ApplyForceStiffCmtCoupled() for Contact "
            "Evaluation!");  // what to do in the predictor?
      has_to_evaluate_ = true;
      return;
    };

    // Overload CONTACT::CoLagrangeStrategy::Recover as this is called in the structure --> no
    // enought information available for monolithically coupled algorithms!
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi)
    {
      if (has_to_recover_ && 0)
        dserror(
            "MonoCoupledLagrangeStrategy::You have to call RecoverCoupled() for Contact Recovery!");
      has_to_recover_ = true;
      return;
    };

    //! @name Access methods

    //@}

    //! @name Evaluation methods

    // Alternative Method to CONTACT::CoAbstractStrategy::ApplyForceStiffCmt for monolithically
    // coupled algorithms
    virtual void ApplyForceStiffCmtCoupled(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& k_ss,
        std::map<int, Teuchos::RCP<LINALG::SparseOperator>*> k_sx,
        Teuchos::RCP<Epetra_Vector>& rhs_s, const int step, const int iter, bool predictor);

    // Alternative Method to CONTACT::CoAbstractStrategy::ApplyForceStiffCmt for monolithically
    // coupled algorithms
    virtual void ApplyForceStiffCmtCoupled(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& k_ss, Teuchos::RCP<LINALG::SparseOperator>& k_sx,
        Teuchos::RCP<Epetra_Vector>& rhs_s, const int step, const int iter, bool predictor);

    // Alternative Method to CONTACT::CoLagrangeStrategy::Recover as this is called in the structure
    // --> no enought information available for monolithically coupled algorithms!
    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers here, which had been
    statically condensated during the setup of the global problem!
    Optionally satisfaction or violation of the contact boundary
    conditions can be checked, too.*/
    virtual void RecoverCoupled(
        Teuchos::RCP<Epetra_Vector> disi, std::map<int, Teuchos::RCP<Epetra_Vector>> inc);

    virtual void RecoverCoupled(Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> inc);

    void EvaluateOffDiagContact(Teuchos::RCP<LINALG::SparseOperator>& kteff,
        int Column_Block_Id);  // condensation for all off diagonal matrixes k_s? in monolithically
                               // coupled problems!

   protected:
    // don't want = operator and cctor
    MonoCoupledLagrangeStrategy operator=(const MonoCoupledLagrangeStrategy& old);
    MonoCoupledLagrangeStrategy(const MonoCoupledLagrangeStrategy& old);

    void SaveCouplingMatrices(Teuchos::RCP<LINALG::SparseMatrix> dhat,
        Teuchos::RCP<LINALG::SparseMatrix> mhataam, Teuchos::RCP<LINALG::SparseMatrix> invda);

    std::map<int, Teuchos::RCP<LINALG::SparseOperator>>
        csx_s_;  // offdiagonal coupling stiffness blocks on slave side!

    Teuchos::RCP<LINALG::SparseMatrix> dhat_;
    Teuchos::RCP<LINALG::SparseMatrix> mhataam_;
    Teuchos::RCP<LINALG::SparseMatrix> invda_;

    Teuchos::RCP<Epetra_Vector>
        lambda_;  // current vector of Lagrange multipliers(for poro no pen.) at t_n+1
    Teuchos::RCP<Epetra_Vector>
        lambdaold_;  // old vector of Lagrange multipliers(for poro no pen.) at t_n

    //! pure useage safty flags
    bool has_to_evaluate_;  // checks if ApplyForceStiffCmtCoupled() after every call of
                            // ApplyForceStiffCmt()
    bool has_to_recover_;   // checks if RecoverCoupled() after every call of Recover()

  };  // class MonoCoupledLagrangeStrategy
}  // namespace CONTACT


#endif  // #ifndef CONTACT_MONOCOUPLED_LAGRANGE_STRATEGY_H
