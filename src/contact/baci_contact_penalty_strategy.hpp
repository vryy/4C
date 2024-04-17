/*---------------------------------------------------------------------*/
/*! \file
\brief Penalty contact solving strategy: The contact constrains are enforced
       by a penalty formulation.

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_PENALTY_STRATEGY_HPP
#define FOUR_C_CONTACT_PENALTY_STRATEGY_HPP

#include "baci_config.hpp"

#include "baci_contact_abstract_strategy.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  // forward declaration
  // class WearInterface;
  /*!
   \brief Contact solving strategy with regularization of Lagrangian multipliers,
   also known as Penalty Method or regularization. An Augmented Lagrangian version
   based on the Uzawa algorithm is included, too.

   This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.

   Refer also to the Semesterarbeit of Bernd Budich, 2009

   */
  class PenaltyStrategy : public AbstractStrategy
  {
   public:
    /*!
    \brief Standard constructor

    \param[in] DofRowMap Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    PenaltyStrategy(const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params, std::vector<Teuchos::RCP<CONTACT::Interface>> interface,
        const int spatialDim, const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf,
        const int maxdof);

    /*!
    \brief Shared data constructor

    \param[in] stratData Data container object
    \param[in] DofRowMap Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] interface All contact interface objects
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    PenaltyStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, const int spatialDim,
        const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof);



    //! @name Access methods

    /*!
    \brief Return L2-norm of active constraints

    */
    double ConstraintNorm() const override { return constrnorm_; }

    /*!
    \brief Return L2-norm of slip constraints

    */
    double ConstraintNormTan() { return constrnormtan_; }


    /*!
    \brief Return initial penalty parameter for non-penetration

    */
    double InitialPenalty() override { return initialpenalty_; }

    /*!
    \brief Return initial penalty parameter for tangential direction

    */
    double InitialPenaltyTan() { return initialpenaltytan_; }


    //@}

    //! @name Evaluation methods

    /*!
    \brief Save nodal kappa-coefficients

    Before starting with the time integration, we have to calculate a nodal scaling factor,
    which will compensate the different integration area for computing the nodal weighted
    gap. Omitting this scaling, nodes on edges or boundaries would have a smaller weighted
    gap, even in case of a uniform physical gap. Hence, this scaling is of crucial importance
    for a penalty strategy since the weighted gap determines the lagrangian multipliers.

    */
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override;

    /*!
    \brief Evaluate relative movement of contact bodies in predictor

    This is a tiny control routine for evaluating the relative movement of
    contact bodies in the predictor of an implicit time integration scheme.
    This evaluation (resetting) is ONLY necessary for penalty strategy and
    Uzawa augmented lagrange strategy, thus this tiny routine here.

    */

    void EvaluateRelMovPredict() override;

    /*!
    \brief Initialize general contact variables for next Newton step

    For a penalty strategy this involves the derivative matrix for the regularized lagrange
    multipliers.

    */
    void Initialize() override;

    /*!
    \brief Evaluate contact

    For a penalty strategy this includes the evaluation of regularized forces
    in normal and tangential direction and results in a simple addition of extra
    stiffness contributions to kteff and extra contact forces to feff.

    */
    void EvaluateContact(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override;

    /*!
    \brief Evaluate frictional contact

    This includes the evaluation of of the frictional contact forces.

    */
    void EvaluateFriction(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override;

    /*!
    \brief Reset penalty parameter to intial value

    When applying an Uzawa Augmented Lagrangian version of the penalty approach,
    the penalty parameter is sometimes updated during the Uzawa steps in
    order to accelerate convergence of the constraint norm. This increase
    in penalty stiffness can be dealt with, because at the time it is applied
    the constraint norm is already quite low. Yet, for a new time step, we have
    to come back to the initial penalty parameter. Thus, this method is called
    at the beginning of each time step and resets the penalty parameter to its initial value.

    */
    void ResetPenalty() override;

    void ModifyPenalty() override;

    /*!
    \brief Initialize Uzawa step


    This method is called at the beginning of the second, third, ... Uzawa
    iterarion in order to create an of an out-of-balance force again. First,
    the contact force and stiffness terms are removed from feff and kteff.
    Then the LM and derivatives are updated (Uzawa AugmentedLagrange) and the new
    contact forces and stiffness terms are created by calling Initialize()
    and finally Evaluate().

    */
    void InitializeUzawa(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override;

    /*!
    \brief Compute L2-norm of active constraints

    In a classical penalty approach, the constraint norm is only monitored.
    When applying an Uzawa Augmented Lagrangian version, the constraint norm is the
    relevant stopping criterion of the Uzawa iteration. In order to accelerate
    convergence, a heuristic update formula for the penalty parameter is applied
    in this method, too.

    */
    void UpdateConstraintNorm(int uzawaiter = 0) override;

    /*!
    \brief Store Lagrange multipliers for next Uzawa step

    A method ONLY called for the Uzawa Augmented Lagrangian version of the penalty method.
    At the end of an Uzawa step, the converged Lagrange multiplier value is stored
    in the variable zuzawa_, which is then used in the next Uzawa step.

    */
    void UpdateUzawaAugmentedLagrange() override;

    /*! \brief Compute force terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration*/
    void EvalForce(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Compute force and stiffness terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration*/
    void EvalForceStiff(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Assemble force and stiffness terms to global vector and matrix */
    void Assemble();

    /*! \brief Run at the beginning of the Evaluate() routine
     *         set force evaluation flag
     *
     */
    void PreEvaluate(CONTACT::ParamsInterface& cparams) override;

    /*! \brief Run in the end of the Evaluate() routine to reset
     *         force evaluation flag
     *
     *
     */
    void PostEvaluate(CONTACT::ParamsInterface& cparams) override;


    /*! \brief Return the desired right-hand-side block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint,*/
    Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bt) const override;

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. displ_displ, displ_lm, ...
     *  \param cparams (in): contact parameter interface (read-only) */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const override;

    //@}

    //! @name Empty functions (Lagrange contact)

    // All these functions only have functionality in Lagrange contact simulations,
    // thus they are defined empty here in the case of Penalty contact.
    Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const override { return Teuchos::null; };
    Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const override { return Teuchos::null; };
    bool ActiveSetSemiSmoothConverged() const override { return true; }
    bool ActiveSetConverged() override { return true; }
    int ActiveSetSteps() override { return 0; }
    void ResetActiveSet() override {}
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override { return; };
    void BuildSaddlePointSystem(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override
    {
      dserror(
          "A penalty approach does not have Lagrange multiplier DOFs. So, saddle point system "
          "makes no sense here.");
    }
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override
    {
      dserror(
          "A penalty approach does not have Lagrange multiplier DOFs. So, saddle point system "
          "makes no sense here.");
    }
    void EvalConstrRHS() override {}
    void UpdateActiveSet() override {}
    void UpdateActiveSetSemiSmooth(const bool firstStepPredictor = false) override {}
    bool IsPenalty() const override { return true; };
    void ResetLagrangeMultipliers(
        const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew) override
    {
    }
    bool IsSaddlePointSystem() const override { return false; }
    bool IsCondensedSystem() const override { return false; }
    bool IsNitsche() const override { return false; }

    /*! \brief Recover the current state
     *
     *  The main task of this method is to recover the Lagrange multiplier solution.
     *  The Lagrange multiplier solution will be stored inside the corresponding strategy
     *  and is necessary for different internal evaluation methods. If the Lagrange multiplier
     *  is condensed, this method is the right place to recover it from the displacement solution.
     *  If it is not condensed (saddle-point system) use the ResetLagrangeMultiplier routine
     * instead.
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration \param xold    (in): old solution vector of the NOX solver \param dir     (in):
     * current search direction (in general NOT the actual step, keep in mind that the step length
     * can differ from 1.0) \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 05/2016
     *  \author hiermeier */
    void RunPostComputeX(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
        const Epetra_Vector& dir, const Epetra_Vector& xnew) override
    {
    }
    Teuchos::RCP<const Epetra_Vector> GetLagrMultN(const bool& redist) const override;
    Teuchos::RCP<const Epetra_Vector> GetLagrMultNp(const bool& redist) const override;
    Teuchos::RCP<Epetra_Vector> LagrMultOld() override;
    Teuchos::RCP<const Epetra_Map> LMDoFRowMapPtr(const bool& redist) const override;

   protected:
    //! derived
    std::vector<Teuchos::RCP<CONTACT::Interface>>& Interfaces() override { return interface_; }

    //! derived
    const std::vector<Teuchos::RCP<CONTACT::Interface>>& Interfaces() const override
    {
      return interface_;
    }

    // don't want = operator and cctor
    PenaltyStrategy operator=(const PenaltyStrategy& old) = delete;
    PenaltyStrategy(const PenaltyStrategy& old) = delete;

    std::vector<Teuchos::RCP<Interface>> interface_;  // contact interfaces
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        linzmatrix_;            // global matrix LinZ with derivatives of LM
    double constrnorm_;         // L2-norm of normal contact constraints
    double constrnormtan_;      // L2-norm of tangential contact constraints
    double initialpenalty_;     // initial penalty parameter
    double initialpenaltytan_;  // initial tangential penalty parameter
    bool evalForceCalled_;      //< flag for evaluate force call

    Teuchos::RCP<Epetra_Vector> fc_;               //< contact penalty force
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kc_;  //< contact penalty stiffness

  };  // class PenaltyStrategy
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
