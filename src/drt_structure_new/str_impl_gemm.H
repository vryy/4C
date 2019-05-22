/*-----------------------------------------------------------*/
/*!

\brief Generalized Energy Momentum time integrator.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_GEMM_H_
#define STR_IMPL_GEMM_H_

#include "str_impl_generic.H"

namespace STR
{
  namespace IMPLICIT
  {
    class Gemm : public Generic
    {
     public:
      //! constructor
      Gemm();

      //! destructor
      virtual ~Gemm(){};

      //! Setup the class variables
      virtual void Setup();

      //! Set state variables (derived)
      virtual void SetState(const Epetra_Vector& x);

      /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
       *
       * \remark Not implemented for the Generalized Energy Momentum time integrator. */
      virtual void AddViscoMassContributions(Epetra_Vector& f) const;

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       *
       * \remark Not implemented for the Generalized Energy Momentum time integrator. */
      virtual void AddViscoMassContributions(LINALG::SparseOperator& jac) const;

      //! Apply the rhs only (derived)
      virtual bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f);

      //! Apply the stiffness only (derived)
      virtual bool ApplyStiff(const Epetra_Vector& x, LINALG::SparseOperator& jac);

      //! Apply force and stiff at once (derived)
      virtual bool ApplyForceStiff(
          const Epetra_Vector& x, Epetra_Vector& f, LINALG::SparseOperator& jac);

      //! (derived)
      bool AssembleForce(Epetra_Vector& f,
          const std::vector<INPAR::STR::ModelType>* without_these_models = NULL) const;

      //! (derived)
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

      //! (derived)
      virtual void ReadRestart(IO::DiscretizationReader& ioreader);

      //! (derived)
      virtual double CalcRefNormForce(const enum NOX::Abstract::Vector::NormType& type);

      //! (derived)
      virtual double GetIntParam() const;

      //! @name Monolithic update routines
      //! @{
      //! Update configuration after time step (derived)
      virtual void UpdateStepState();

      //! Update everything on element level after time step and after output (derived)
      virtual void UpdateStepElement();
      //! @}

      //! @name Predictor routines (dependent on the implicit integration scheme)
      //! @{
      /*! Predict constant displacements, consistent velocities and accelerations (derived) */
      virtual void PredictConstDisConsistVelAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      /*! Predict displacements based on constant velocities and consistent accelerations (derived)
       */
      virtual bool PredictConstVelConsistAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;

      /*! Predict displacements based on constant accelerations and consistent velocities (derived)
       */
      virtual bool PredictConstAcc(
          Epetra_Vector& disnp, Epetra_Vector& velnp, Epetra_Vector& accnp) const;
      //! @}
    };
  }  // namespace IMPLICIT
}  // namespace STR


#endif /* STR_IMPL_GEMM_H_ */
