/*-----------------------------------------------------------*/
/*! \file

\brief Generalized Energy Momentum time integrator.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_IMPL_GEMM_H
#define BACI_STRUCTURE_NEW_IMPL_GEMM_H

#include "baci_structure_new_impl_generic.H"

BACI_NAMESPACE_OPEN

namespace STR
{
  namespace IMPLICIT
  {
    class Gemm : public Generic
    {
     public:
      //! constructor
      Gemm();


      //! Setup the class variables
      void Setup() override;

      //! (derived)
      void PostSetup() override;

      //! Set state variables (derived)
      void SetState(const Epetra_Vector& x) override;

      /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
       *
       * \remark Not implemented for the Generalized Energy Momentum time integrator. */
      void AddViscoMassContributions(Epetra_Vector& f) const override;

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       *
       * \remark Not implemented for the Generalized Energy Momentum time integrator. */
      void AddViscoMassContributions(CORE::LINALG::SparseOperator& jac) const override;

      //! Apply the rhs only (derived)
      bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) override;

      //! Apply the stiffness only (derived)
      bool ApplyStiff(const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac) override;

      //! Apply force and stiff at once (derived)
      bool ApplyForceStiff(
          const Epetra_Vector& x, Epetra_Vector& f, CORE::LINALG::SparseOperator& jac) override;

      //! (derived)
      bool AssembleForce(Epetra_Vector& f,
          const std::vector<INPAR::STR::ModelType>* without_these_models = nullptr) const override;

      //! (derived)
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      //! (derived)
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! (derived)
      double CalcRefNormForce(const enum ::NOX::Abstract::Vector::NormType& type) const override;

      //! (derived)
      double GetIntParam() const override;

      //! @name Monolithic update routines
      //! @{
      //! Update configuration after time step (derived)
      void UpdateStepState() override;

      //! Update everything on element level after time step and after output (derived)
      void UpdateStepElement() override;
      //! @}

      //! @name Predictor routines (dependent on the implicit integration scheme)
      //! @{
      /*! Predict constant displacements, consistent velocities and accelerations (derived) */
      void PredictConstDisConsistVelAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const override;

      /*! Predict displacements based on constant velocities and consistent accelerations (derived)
       */
      bool PredictConstVelConsistAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const override;

      /*! Predict displacements based on constant accelerations and consistent velocities (derived)
       */
      bool PredictConstAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const override;
      //! @}

      //! @name access methods
      //@{

      //! Return name
      enum INPAR::STR::DynamicType MethodName() const override { return INPAR::STR::dyna_gemm; }

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a m-multistep method returns m
      int MethodSteps() const override { return 1; }

      //! Give linear order of accuracy of displacement part
      int MethodOrderOfAccuracyDis() const override
      {
        dserror("Not yet available");
        return 0;
      }

      //! Give linear order of accuracy of velocity part
      int MethodOrderOfAccuracyVel() const override
      {
        dserror("Not yet available");
        return 0;
      }

      //! Return linear error coefficient of displacements
      double MethodLinErrCoeffDis() const override
      {
        dserror("Not yet available");
        return 0.0;
      }

      //! 2nd order linear error coefficient of displacements
      double MethodLinErrCoeffDis2() const
      {
        dserror("Not yet available");
        return 0.0;
      }

      //! 3rd order linear error coefficient of displacements
      double MethodLinErrCoeffDis3() const
      {
        dserror("Not yet available");
        return 0.0;
      }

      //! Return linear error coefficient of velocities
      double MethodLinErrCoeffVel() const override
      {
        dserror("Not yet available");
        return 0.0;
      }

      //! 1st order linear error coefficient of velocities
      double MethodLinErrCoeffVel1() const
      {
        dserror("Not yet available");
        return 0.0;
      }

      //! 2nd order linear error coefficient of velocities
      double MethodLinErrCoeffVel2() const
      {
        dserror("Not yet available");
        return 0.0;
      }
      //@}
    };
  }  // namespace IMPLICIT
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_IMPL_GEMM_H
