/*----------------------------------------------------------------------*/
/*!
 \file scatra_timint_poromulti.H

 \brief time integration schemes for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_SCATRA_TIMINT_POROMULTI_H_
#define SRC_DRT_SCATRA_SCATRA_TIMINT_POROMULTI_H_

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"
#include "../drt_scatra/scatra_timint_genalpha.H"
#include "../drt_scatra/scatra_timint_stat.H"


namespace SCATRA
{
  class ScaTraTimIntPoroMulti : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    ScaTraTimIntPoroMulti(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// virtual destructor
    virtual ~ScaTraTimIntPoroMulti(){};

    // -----------------------------------------------------------------
    // general methods
    // -----------------------------------------------------------------

    /// initialize algorithm
    virtual void Init();

    //! update the solution after convergence of the nonlinear iteration.
    virtual void Update(const int num = 0  //!< field number
    )
    {
      return;
    };

    //! set the nodal L2-flux
    virtual void SetL2FluxOfMultiFluid(
        Teuchos::RCP<const Epetra_MultiVector> multiflux, const int nds_flux);

    //! set solution field of the multiphase fluid
    virtual void SetSolutionFieldOfMultiFluid(Teuchos::RCP<const Epetra_Vector> phinp_fluid,
        Teuchos::RCP<const Epetra_Vector> phin_fluid, const int nds_phi_fluid);

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
    virtual void OutputState();

    //! problem specific output
    void OutputProblemSpecific();

    //! add parameters depending on the problem
    virtual void AddProblemSpecificParametersAndVectors(
        Teuchos::ParameterList& params  //!< parameter list
    );

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
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~ScaTraTimIntPoroMultiOST();

    //! initialize time integration scheme
    void Init();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

  };  // class TimIntPoroMultiOST


  class ScaTraTimIntPoroMultiBDF2 : public ScaTraTimIntPoroMulti, public TimIntBDF2
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiBDF2(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~ScaTraTimIntPoroMultiBDF2();

    //! initialize time integration scheme
    void Init();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

  };  // class TimIntPoroMultiBDF2


  class ScaTraTimIntPoroMultiGenAlpha : public ScaTraTimIntPoroMulti, public TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~ScaTraTimIntPoroMultiGenAlpha();

    //! initialize time integration scheme
    void Init();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

  };  // class TimIntPoroMultiGenAlpha


  class ScaTraTimIntPoroMultiStationary : public ScaTraTimIntPoroMulti, public TimIntStationary
  {
   public:
    //! Standard Constructor
    ScaTraTimIntPoroMultiStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~ScaTraTimIntPoroMultiStationary();

    //! initialize time integration scheme
    void Init();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

  };  // class TimIntPoroMultiStationary
}  // namespace SCATRA



#endif /* SRC_DRT_SCATRA_SCATRA_TIMINT_POROMULTI_H_ */
