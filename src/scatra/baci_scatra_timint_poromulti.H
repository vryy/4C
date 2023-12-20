/*----------------------------------------------------------------------*/
/*! \file
 \brief time integration schemes for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_POROMULTI_H
#define BACI_SCATRA_TIMINT_POROMULTI_H

#include "baci_config.H"

#include "baci_scatra_timint_bdf2.H"
#include "baci_scatra_timint_genalpha.H"
#include "baci_scatra_timint_implicit.H"
#include "baci_scatra_timint_ost.H"
#include "baci_scatra_timint_stat.H"

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


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
    virtual void SetL2FluxOfMultiFluid(Teuchos::RCP<const Epetra_MultiVector> multiflux);

    //! set solution field of the multiphase fluid
    virtual void SetSolutionFieldOfMultiFluid(Teuchos::RCP<const Epetra_Vector> phinp_fluid,
        Teuchos::RCP<const Epetra_Vector> phin_fluid);

    //! set the velocity field (zero or field by function)
    virtual void SetVelocityField(const int nds)
    {
      dserror(
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
      dserror(
          "SetVelocityField(...) cannot be used for transport within a multiphase porous medium!"
          " Use SetSolutionFields(...) instead!");
    };

    //! write state vectors (phinp and convective velocity) to BINIO
    void OutputState() override;

    //! problem specific output
    void OutputProblemSpecific() override;

    //! add parameters depending on the problem
    void AddProblemSpecificParametersAndVectors(Teuchos::ParameterList& params  //!< parameter list
        ) override;

   protected:
    //! output of oxygen partial pressure
    void OutputOxygenPartialPressure();
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



BACI_NAMESPACE_CLOSE

#endif  // SCATRA_TIMINT_POROMULTI_H
