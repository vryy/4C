/*-----------------------------------------------------------*/
/*! \file

\brief Explicit Euler time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_EXPL_FORWARDEULER_HPP
#define FOUR_C_STRUCTURE_NEW_EXPL_FORWARDEULER_HPP

#include "baci_config.hpp"

#include "baci_structure_new_expl_generic.hpp"

BACI_NAMESPACE_OPEN

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief Explicit Euler time integration for solid dynamics
     *
     */
    class ForwardEuler : public Generic
    {
     public:
      //! constructor
      ForwardEuler();


      //! Setup class variables (derived)
      void Setup() override;

      //! Post setup operation (compute initial equilibrium state) (derived)
      void PostSetup() override;

      //! Set state variables (derived)
      void SetState(const Epetra_Vector& x) override;

      //! return integration factor (derived)
      double GetIntParam() const override { return -1.0; }

      /*! \brief Add the viscous and mass contributions to the right hand side
       */
      void AddViscoMassContributions(Epetra_Vector& f) const override;

      /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
       */
      void AddViscoMassContributions(CORE::LINALG::SparseOperator& jac) const override;

      //! Update configuration after time step (derived)
      void UpdateStepState() override;

      //! (derived)
      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      /*! read restart information of the different time integration schemes
       *  and model evaluators (derived) */
      void ReadRestart(IO::DiscretizationReader& ioreader) override;

      //! @name Attribute access functions
      //@{

      //! Return time integrator name
      [[nodiscard]] enum INPAR::STR::DynamicType MethodName() const override
      {
        return INPAR::STR::dyna_expleuler;
      }

      //! Provide number of steps, e.g. a single-step method returns 1,
      //! a m-multistep method returns m
      [[nodiscard]] int MethodSteps() const override { return 1; }

      //! Give local order of accuracy of displacement part
      [[nodiscard]] int MethodOrderOfAccuracyDis() const override { return 1; }

      //! Give local order of accuracy of velocity part
      [[nodiscard]] int MethodOrderOfAccuracyVel() const override { return 1; }

      /*! \brief Return linear error coefficient of displacements
       *
       *  The local discretization error reads
       *  \f[
       *  e \approx \frac{1}{2}\Delta t_n^2 \ddot{d_n} + HOT(\Delta t_n^3)
       *  \f]
       */
      [[nodiscard]] double MethodLinErrCoeffDis() const override { return 0.5; }

      /*! \brief Return linear error coefficient of velocities
       *
       *  The local discretization error reads
       *  \f[
       *  e \approx \frac{1}{2}\Delta t_n^2 \dddot{d_n} + HOT(\Delta t_n^3)
       *  \f]
       */
      [[nodiscard]] double MethodLinErrCoeffVel() const override { return 0.5; }

      //@}

     private:
      bool modexpleuler_;  //!< modified explicit Euler equation (veln_ instead of vel_ for
                           //!< calc of disn_), default: true

      //! viscous force vector F_viscous F_{viscous;n+1}
      Teuchos::RCP<Epetra_Vector> fvisconp_ptr_;

      //! viscous force vector F_viscous F_{viscous;n}
      Teuchos::RCP<Epetra_Vector> fviscon_ptr_;

      //! pointer to inertial force vector F_{inertial,n+1} at new time
      Teuchos::RCP<Epetra_Vector> finertianp_ptr_;

      //! pointer to inertial force vector F_{inertial,n} at last time
      Teuchos::RCP<Epetra_Vector> finertian_ptr_;
    };
  }  // namespace EXPLICIT
}  // namespace STR

BACI_NAMESPACE_CLOSE

#endif  // BACI_STRUCTURE_NEW_EXPL_FORWARDEULER_H
