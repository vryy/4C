/*----------------------------------------------------------------------*/
/*! \file
 \brief time integration schemes for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_POROMULTI_HPP
#define FOUR_C_SCATRA_TIMINT_POROMULTI_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_bdf2.hpp"
#include "4C_scatra_timint_genalpha.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_ost.hpp"
#include "4C_scatra_timint_stat.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace SCATRA
{
  class ScaTraTimIntPoroMulti : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    ScaTraTimIntPoroMulti(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    // -----------------------------------------------------------------
    // general methods
    // -----------------------------------------------------------------

    /// initialize algorithm
    void Init() override;

    //! update the solution after convergence of the nonlinear iteration.
    void Update() override{};

    //! set the nodal L2-flux
    virtual void set_l2_flux_of_multi_fluid(Teuchos::RCP<const Epetra_MultiVector> multiflux);

    //! set solution field of the multiphase fluid
    virtual void set_solution_field_of_multi_fluid(Teuchos::RCP<const Epetra_Vector> phinp_fluid,
        Teuchos::RCP<const Epetra_Vector> phin_fluid);

    //! set the velocity field (zero or field by function)
    virtual void SetVelocityField(const int nds)
    {
      FOUR_C_THROW(
          "SetVelocityField(...) cannot be used for transport within a multiphase porous medium!"
          " Use SetSolutionFields(...) instead!");
    };

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void SetVelocityField(
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity/press. vector
        Teuchos::RCP<const Epetra_Vector> acc,      //!< acceleration vector
        Teuchos::RCP<const Epetra_Vector> vel,      //!< velocity vector
        Teuchos::RCP<const Epetra_Vector> fsvel,    //!< fine-scale velocity vector
        const int nds,  //!< number of the dofset the velocity/pressure state belongs to
        const bool setpressure =
            false  //!< flag whether the fluid pressure needs to be known for the scatra
    )
    {
      FOUR_C_THROW(
          "SetVelocityField(...) cannot be used for transport within a multiphase porous medium!"
          " Use SetSolutionFields(...) instead!");
    };

    //! write state vectors (phinp and convective velocity) to BINIO
    void OutputState() override;

    //! problem specific output
    void output_problem_specific() override;

    //! add parameters depending on the problem
    void add_problem_specific_parameters_and_vectors(
        Teuchos::ParameterList& params  //!< parameter list
        ) override;

   protected:
    //! output of oxygen partial pressure
    void output_oxygen_partial_pressure();
    //! do we employ L2-projection for reconstruction of velocity field
    bool L2_projection_;
  };


  class ScaTraTimIntPoroMultiOST : public ScaTraTimIntPoroMulti, public TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiOST(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void Init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

  };  // class TimIntPoroMultiOST


  class ScaTraTimIntPoroMultiBDF2 : public ScaTraTimIntPoroMulti, public TimIntBDF2
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiBDF2(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void Init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

  };  // class TimIntPoroMultiBDF2


  class ScaTraTimIntPoroMultiGenAlpha : public ScaTraTimIntPoroMulti, public TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void Init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

  };  // class TimIntPoroMultiGenAlpha


  class ScaTraTimIntPoroMultiStationary : public ScaTraTimIntPoroMulti, public TimIntStationary
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! initialize time integration scheme
    void Init() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

  };  // class TimIntPoroMultiStationary
}  // namespace SCATRA



FOUR_C_NAMESPACE_CLOSE

#endif
