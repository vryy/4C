/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all implicit time integrators

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_GENERIC_H_
#define STR_IMPL_GENERIC_H_

#include "str_integrator.H"
#include "../solver_nonlin_nox/nox_nln_abstract_prepostoperator.H"
// required because of used enums
#include <NOX_Abstract_Vector.H>

// forward declaration ...
namespace LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace LINALG

namespace STR
{
  namespace IMPLICIT
  {
    /*! \brief A generic fully implicit time integrator for solid dynamics
     *
     * This serves as a base class for all fully implicit time integration schemes.
     */
    class Generic : public Integrator
    {
     public:
      //! constructor
      Generic();

      //! destructor
      virtual ~Generic(){};

      //! Setup (has to be implemented by the derived classes)
      virtual void Setup();

      /*! \brief Set state variables (derived)
       *
       * \param[in] x State to be set
       */
      virtual void SetState(const Epetra_Vector& x) = 0;

      /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
       *
       * \remark The remaining contributions have been considered in the corresponding model
       *         evaluators. This is due to the fact, that some models use a different
       *         time integration scheme for their terms (e.g. GenAlpha for the structure
       *         and OST for the remaining things).
       *
       * \param[in/out] f Right-hand side vector
       */
      virtual void AddViscoMassContributions(Epetra_Vector& f) const = 0;

      /*! \brief Add the viscous and mass contributions to the Jacobian (TR-rule)
       *
       * \remark The remaining blocks have been considered in the corresponding model
       *         evaluators. This is due to the fact, that some models use a different
       *         time integration scheme for their terms (e.g. GenAlpha for the structure
       *         and OST for the remaining things). Furthermore, constraint/Lagrange
       *         multiplier blocks need no scaling anyway.
       *
       * \param[in/out] jac Jacobian matrix
       */
      virtual void AddViscoMassContributions(LINALG::SparseOperator& jac) const = 0;

      //! derived
      virtual bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) = 0;

      /*! \brief Apply the stiffness only
       *
       * Normally this one is unnecessary, since it makes more sense
       * to evaluate the stiffness and right hand side at once, because of
       * the lower computational overhead. */
      virtual bool ApplyStiff(const Epetra_Vector& x, LINALG::SparseOperator& jac) = 0;

      /*! \brief Apply force and stiff at once
       *
       *  Only one loop over all elements. Especially in the contact case,
       *  the difference between this call and first call ApplyForce and
       *  then ApplyStiff is mentionable because of the projection operations. */
      virtual bool ApplyForceStiff(
          const Epetra_Vector& x, Epetra_Vector& f, LINALG::SparseOperator& jac) = 0;

      bool ApplyCorrectionSystem(const enum NOX::NLN::CorrectionType type,
          const std::vector<INPAR::STR::ModelType>& constraint_models, const Epetra_Vector& x,
          Epetra_Vector& f, LINALG::SparseOperator& jac);

      /*! \brief Remove contributions from the structural right-hand side stemming
       *  from any condensation operations (typical example is contact)
       *
       *  \author hiermeier \date 03/18 */
      void RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) const;

      //! derived
      virtual bool AssembleForce(Epetra_Vector& f,
          const std::vector<INPAR::STR::ModelType>* without_these_models = NULL) const = 0;

      //! derived
      virtual bool AssembleJac(LINALG::SparseOperator& jac,
          const std::vector<INPAR::STR::ModelType>* without_these_models = NULL) const
      {
        return false;
      };

      //! derived
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const = 0;

      //! derived
      virtual void ReadRestart(IO::DiscretizationReader& ioreader) = 0;

      /*! \brief Calculate characteristic/reference norms for forces
       *
       *  The reference norms are used to scale the calculated iterative
       *  displacement norm and/or the residual force norm. For this
       *  purpose we only need the right order of magnitude, so we don't
       *  mind evaluating the corresponding norms at possibly different
       *  points within the time-step (end point, generalized midpoint). */
      virtual double CalcRefNormForce(const enum NOX::Abstract::Vector::NormType& type) = 0;

      //! return time integration factor
      virtual double GetIntParam() const = 0;

      //! return the default step length of the used NOX::LineSearch method
      double GetDefaultStepLength() const;

      //! @name Monolithic update routines
      //! @{
      //! things that should be done before updating (derived)
      virtual void PreUpdate(){/* do nothing for now */};

      //! Update configuration after time step (derived)
      virtual void UpdateStepState() = 0;

      //! Update everything on element level after time step and after output (derived)
      virtual void UpdateStepElement() = 0;

      //! things that should be done after updating (derived)
      virtual void PostUpdate(){/* do nothing for now */};

      //! update constant contributions of the current state for the new time step \f$t_{n+1}\f$
      virtual void UpdateConstantStateContributions(){/* do nothing for some integrators */};
      //! @}

      //! @name Predictor routines (dependent on the implicit integration scheme)
      //!@{

      /*! \brief Predict constant displacements with consistent velocities and accelerations
       *
       * The displacement field is kept constant, however the velocities
       * and accelerations are consistent to the time integration
       * if the constant displacements are taken as correct displacement solution.
       *
       * This method has to be implemented by the individual time
       * integrator since the calculation of consistent velocities and accelerations
       * depends on the actual time integration scheme.
       */
      virtual void PredictConstDisConsistVelAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const = 0;

      /*! \brief Predict displacements based on the assumption of constant velocities.
       *
       * Assuming constant velocities, new displacements are predicted.
       * Calculate consistent accelerations afterwards.
       *
       * This method has to be implemented by the individual time
       * integrator since the calculation of consistent velocities and accelerations
       * depends on the actual time integration scheme.
       *
       * \param[in/out] disnp Displacement vector
       * \param[in/out] velnp Velocity vector
       * \param[in/out] accnp Acceleration vector
       */
      virtual bool PredictConstVelConsistAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const = 0;

      /*! \brief Predict displacements based on the assumption of constant accelerations.
       *
       * Assuming constant accelerations, new velocities and displacements are predicted.
       *
       * This method has to be implemented by the individual time
       * integrator since the calculation of consistent velocities and accelerations
       * depends on the actual time integration scheme.
       *
       * \param[in/out] disnp Displacement vector
       * \param[in/out] velnp Velocity vector
       * \param[in/out] accnp Acceleration vector
       */
      virtual bool PredictConstAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const = 0;
      //!@}

      /*! \brief Set the predictor state flag
       *
       * \param[in] ispredictor_state Predictor state flag
       */
      void SetIsPredictorState(const bool& ispredictor_state);

      //! Get the predictor state flag
      const bool& IsPredictorState() const;

      //! compute the scaling operator for element based scaling using PTC
      void ComputeJacobianContributionsFromElementLevelForPTC(
          Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr);

      /*! \brief Print jacbian into text file for later use in MATLAB
       *
       *  \param[in] NOX group containing the linear system with the Jacobian
       *
       *  \author hiermeier \date 06/17 */
      void PrintJacobianInMatlabFormat(const NOX::NLN::Group& curr_grp) const;

      //! Get the NOX parameter list
      Teuchos::ParameterList& GetNoxParams();

      /// compute the condition number of the structural tangential stiffness matrix if desired
      void ConditionNumber(const NOX::NLN::Group& grp) const;

     protected:
      //! reset the time step dependent parameters for the element evaluation [derived]
      virtual void ResetEvalParams();

     private:
      /*! \brief Flag indicating if the current state is the predictor state.
       *
       *  In the evaluation of the predictor state the SetState() routine is
       *  not allowed to calculate the consistent velocities and accelerations
       *  as usual. This is due to the fact, that the predictor might lead to
       *  velocities and accelerations that are not consistently computed from
       *  the displacements based on the time integration scheme. Instead we
       *  leave the predictor state untouched during the first evaluation. */
      bool ispredictor_state_;

    };  // namespace IMPLICIT
  }     // namespace IMPLICIT
}  // namespace STR

namespace NOX
{
  namespace NLN
  {
    namespace PrePostOp
    {
      namespace IMPLICIT
      {
        /*! \brief Implicit time integration helper class
         *
         *  This class is an implementation of the NOX::NLN::Abstract::PrePostOperator
         *  and is used to modify the computeX() routine of the given NOX::NLN::Group.
         *  It's called by the wrapper class NOX::NLN::GROUP::PrePostOperator.
         *
         *  \author Michael Hiermeier */
        class Generic : public NOX::NLN::Abstract::PrePostOperator
        {
         public:
          //! constructor
          Generic(const ::STR::IMPLICIT::Generic& implicit)
              : impl_(implicit), default_step_(implicit.GetDefaultStepLength()){/* empty */};

          //! destructor
          virtual ~Generic(){};

          /*! \brief Derived function, which is called at the very beginning of a call to
           *  NOX::NLN::Group::computeX()
           *
           *  This method is used to get access to the current direction vector and
           *  to augment/modify the direction vector before the actual solution update is
           *  performed. One possible scenario is the CONTACT::AUG::STEEPESTASCENT::Strategy,
           *  where we calculate the Lagrange multiplier increment in a post-processing
           *  step.
           *
           *  \author hiermeier \date 03/17 */
          void runPreComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
              const double& step, const NOX::NLN::Group& curr_grp) override;

          /*! \brief Derived function, which is called at the end of a call to
           * NOX::NLN::Group::computeX()
           *
           *  This method is used to get access to the current direction vector. The
           *  direction vector is needed for different internal update routines. Two
           *  examples on element level are the EAS recovery, since the EAS DoFs are
           *  condensed and the calculation of the strain increments. The displacement
           *  increment plays also a role in other condensation approaches, like the
           *  mortar dual strategies.
           *
           *  \author hiermeier \date 07/16 */
          void runPostComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
              const double& step, const NOX::NLN::Group& curr_grp) override;

          /*! \brief Derived function, which is called at the very end of a call to
           *  NOX::Solver::Generic::step()
           *
           *  This method gives you the opportunity to do something in the end of
           *  a successful or unsuccessful nonlinear solver step. It is called after
           *  a possible LineSearch or other globalization method has been used.
           *
           *  \author hiermeier \date 03/17 */
          void runPostIterate(const NOX::Solver::Generic& solver) override;

          /*! \brief Derived function, which is called at the very end of a call to
           *  NOX::NLN::Group::applyJacobianInverse()
           *
           *  This method gives you the opportunity to do something in the end of
           *  a successful or unsuccessful linear solver attempt.
           *
           *  \note The result vector is the actual result vector of the internal
           *  linear solver and, accordingly, due to the used sign convention in NOX,
           *  the NEGATIVE direction vector. The sign will be changed again in the
           *  NOX::Epetra::Group::computeNewton method.
           *
           *  \param rhs    : read-only access to the rhs vector
           *  \param result : full access to the result vector
           *  \param xold   : read-only access to the old state vector
           *  \param grp    : read only access to the group object
           *
           *  \author hiermeier \date 12/17 */
          void runPostApplyJacobianInverse(const NOX::Abstract::Vector& rhs,
              NOX::Abstract::Vector& result, const NOX::Abstract::Vector& xold,
              const NOX::NLN::Group& grp) override;

          /// \brief Called at the very beginning of a Newton loop
          /**
           *  \author hiermeier \date 05/18 */
          void runPreSolve(const NOX::Solver::Generic& solver) override;

          /*! \brief Derived function, which is called at the beginning of a call to
           *  NOX::NLN::Group::applyJacobianInverse()
           *
           *  \param rhs    : read-only access to the rhs vector
           *  \param result : full access to the result vector
           *  \param xold   : read-only access to the old state vector
           *  \param grp    : read only access to the group object
           *
           *  \author seitz \date 04/18 */
          virtual void runPreApplyJacobianInverse(const NOX::Abstract::Vector& rhs,
              NOX::Abstract::Vector& result, const NOX::Abstract::Vector& xold,
              const NOX::NLN::Group& grp);

         private:
          /*! \brief Convert NOX::Abstract::Vector to Epetra_Vector
           *
           *  \param[in] vec Vector to be converted
           */
          Epetra_Vector& convert2EpetraVector(NOX::Abstract::Vector& vec) const;
          const Epetra_Vector& convert2EpetraVector(const NOX::Abstract::Vector& vec) const;

          /// get the step length
          bool getStep(double& step, const NOX::Solver::Generic& solver) const;

          /// get the number of necessary system corrections in case of a mod newton direction
          int getNumberOfModifiedNewtonCorrections(const NOX::Solver::Generic& solver) const;

         private:
          //! reference to the STR::IMPLICIT::Generic object (read-only)
          const ::STR::IMPLICIT::Generic& impl_;

          //! default step length
          const double default_step_;

        };  // class Generic
      }     // namespace IMPLICIT
    }       // namespace PrePostOp
  }         // namespace NLN
}  // namespace NOX


#endif /* STR_IMPL_GENERIC_H_ */
