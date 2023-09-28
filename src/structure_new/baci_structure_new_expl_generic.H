/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all explicit time integrators


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_EXPL_GENERIC_H
#define BACI_STRUCTURE_NEW_EXPL_GENERIC_H

#include "baci_solver_nonlin_nox_abstract_prepostoperator.H"
#include "baci_structure_new_integrator.H"

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief A generic fully explicit time integrator for solid dynamics
     *
     * This serves as a base class for all fully explicit time integration schemes.
     */
    class Generic : public STR::Integrator
    {
     public:
      //! constructor
      Generic() = default;

      //! destructor
      ~Generic() override = default;

      //! Setup (has to be implemented by the derived classes)
      void Setup() override;

      //! Apply the right hand side only (derived)
      bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) override;

      //! \brief Apply the stiffness only (derived)
      bool ApplyStiff(const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac) override;

      //! \brief Apply force and stiff at once (derived)
      bool ApplyForceStiff(
          const Epetra_Vector& x, Epetra_Vector& f, CORE::LINALG::SparseOperator& jac) override;

      /*! \brief (derived)
       *
       */
      bool ApplyCorrectionSystem(const enum NOX::NLN::CorrectionType type,
          const std::vector<INPAR::STR::ModelType>& constraint_models, const Epetra_Vector& x,
          Epetra_Vector& f, CORE::LINALG::SparseOperator& jac) override
      {
        return false;
      }

      /*! \brief Calculate characteristic/reference norms for forces (derived)
       */
      [[nodiscard]] double CalcRefNormForce(
          const enum NOX::Abstract::Vector::NormType& type) const override;

      //! return the default step length
      [[nodiscard]] double GetDefaultStepLength() const;

      //! compute the scaling operator for element based scaling using PTC (derived)
      void ComputeJacobianContributionsFromElementLevelForPTC(
          Teuchos::RCP<CORE::LINALG::SparseMatrix>& scalingMatrixOpPtr) override;

      //! Assemble the right hand side
      bool AssembleForce(Epetra_Vector& f,
          const std::vector<INPAR::STR::ModelType>* without_these_models = nullptr) const override;

      //! @name Monolithic update routines
      //! @{
      //! things that should be done before updating (derived)
      void PreUpdate() override{/* do nothing for now */};

      //! (derived)
      void UpdateConstantStateContributions() override;

      virtual void OutputStepstate() { ; };

      //! Update everything on element level after time step and after output (derived)
      void UpdateStepElement() override;

      //! things that should be done after updating (derived)
      void PostUpdate() override;

      /*! \brief Remove contributions from the structural right-hand side stemming
       *  from any condensation operations (typical example is contact) (derived)
       */
      void RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) const override;
      //! @}

      //! @name Attribute access functions
      //@{

      //! Provide Name
      [[nodiscard]] virtual enum INPAR::STR::DynamicType MethodName() const = 0;

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a \f$m\f$-multistep method returns \f$m\f$
      [[nodiscard]] virtual int MethodSteps() const = 0;

      //! Give local order of accuracy of displacement part
      [[nodiscard]] virtual int MethodOrderOfAccuracyDis() const = 0;

      //! Give local order of accuracy of velocity part
      [[nodiscard]] virtual int MethodOrderOfAccuracyVel() const = 0;

      //! Return linear error coefficient of displacements
      [[nodiscard]] virtual double MethodLinErrCoeffDis() const = 0;

      //! Return linear error coefficient of velocities
      [[nodiscard]] virtual double MethodLinErrCoeffVel() const = 0;

      //! @}

     protected:
      //! reset the time step dependent parameters for the element evaluation [derived]
      void ResetEvalParams() override;
    };
  }  // namespace EXPLICIT
}  // namespace STR

namespace NOX
{
  namespace NLN
  {
    namespace PrePostOp
    {
      namespace EXPLICIT
      {
        /*! \brief Explicit time integration helper class
         */
        class Generic : public NOX::NLN::Abstract::PrePostOperator
        {
         public:
          //! constructor
          Generic(const ::STR::EXPLICIT::Generic& expl)
              : expl_(expl), default_step_(expl.GetDefaultStepLength())
          {
          }

          //! destructor
          ~Generic() override = default;

          //! derived
          void runPreIterate(const NOX::Solver::Generic& solver) override;

          //! derived
          void runPreSolve(const NOX::Solver::Generic& nlnSolver) override;

          //! derived
          void runPreComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
              const double& step, const NOX::NLN::Group& curr_grp) override;

          //! derived
          void runPostComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
              const double& step, const NOX::NLN::Group& curr_grp) override;

         private:
          //! reference to the STR::IMPLICIT::Generic object (read-only)
          const ::STR::EXPLICIT::Generic& expl_;

          //! default step length
          const double default_step_;
        };  // class Generic
      }     // namespace EXPLICIT
    }       // namespace PrePostOp
  }         // namespace NLN
}  // namespace NOX

#endif  // BACI_STRUCTURE_NEW_EXPL_GENERIC_H
